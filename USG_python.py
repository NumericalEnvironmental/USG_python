##########################################################################
#
# USG_python.py - a 2D, unstructured grid confined groundwater flow model
# for a 1-layer aquifer with variable geometry
#
##########################################################################

from numpy import *
from pandas import *
from scipy.spatial import Delaunay
from scipy.spatial import KDTree
from scipy.interpolate import griddata
from scipy.spatial import Voronoi,voronoi_plot_2d                
import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull
from scipy.sparse import *
from scipy.sparse.linalg import *


### classes ###


class Connections:
    def __init__(self, vor, cells_df):
        self.nodes = vor.ridge_points                   # nodes[0] and nodes[1] are the connected nodes
        node1_loc = empty((len(self.nodes), 2), float)
        node2_loc = empty((len(self.nodes), 2), float)        
        vert1_loc = empty((len(self.nodes), 2), float)
        vert2_loc = empty((len(self.nodes), 2), float)
        for i, node in enumerate(self.nodes):
            # fetch Voronoi regions associated with nodes[0] and nodes[1]
            r1 = cells_df['vor_region'][node[0]]     
            r2 = cells_df['vor_region'][node[1]]
            # find the Voronoi vertices index numbers that are associated with the connections between nodes[0] and nodes[1]
            s1 = vor.regions[r1]
            s2 = vor.regions[r2]
            vert_indx = array(list(set(s1) & set(s2)))
            # note coordinate sets for nodes and Voronoi vertices
            node1_loc[i] = vor.points[node[0]]
            node2_loc[i] = vor.points[node[1]]            
            vert1_loc[i] = vor.vertices[vert_indx[0]]
            vert2_loc[i] = vor.vertices[vert_indx[1]]
        # interfacial areas of connections
        width = Distance(vert1_loc, vert2_loc)
        self.area = width * 0.5 * (array(cells_df['b'][self.nodes.T[0]]) + array(cells_df['b'][self.nodes.T[1]]))               
        # distance from node_1 to node_2, used to approximate gradient
        self.d = Distance(node1_loc, node2_loc)
        self.conduct = 0.5 * (array(cells_df['K'][self.nodes.T[0]]) + array(cells_df['K'][self.nodes.T[1]])) * self.area/self.d
        self.WriteConnections()             # summarize
    def WriteConnections(self):
        # write connection attributes to output file
        output_file = open('connections.csv','w')
        output_file.writelines(['node_1', ',', 'node_2', ',', 'area', ',', 'd', ',', 'conductance', '\n'])
        for i, node in enumerate(self.nodes):
            output_line = [str(node[0]), ',', str(node[1]), ',', str(self.area[i]), ',', str(self.d[i]), ',', str(self.conduct[i]), '\n']
            output_file.writelines(output_line)
        output_file.close()


class LHS_Matrix:
    def __init__(self, cells_df, connections, params, dt):
        # compute all of the terms on the mass balance equation set left-hand-side; diagonal terms will be subject to change as time step changes
        row_index = []
        col_index = []
        data = []
        # diagonal terms
        for i in xrange(len(cells_df)):         # diagonal term setup isn't vectorized because of sparse matrix structure
            row_index.append(i)
            col_index.append(i)
            data.append(cells_df['vol'][i]*cells_df['Ss'][i]/dt + params.gamma*cells_df['sum_conduct'][i])
        # off-diagonal terms
        for i in xrange(len(connections.nodes)):
            node_1 = connections.nodes[i, 0]
            node_2 = connections.nodes[i, 1]
            row_index.append(node_1)
            col_index.append(node_2)
            data.append(-params.gamma * connections.conduct[i])
            row_index.append(node_2)
            col_index.append(node_1)
            data.append(-params.gamma * connections.conduct[i])
        # summarize
        self.row_index = array(row_index)
        self.col_index = array(col_index)
        self.data = array(data)
    def Diagonal(self, cells_df, params, dt):
        # update the diagonal terms in the left-hand-side matrix of equations
        self.data[:len(cells_df)] = cells_df['vol']*cells_df['Ss']/dt + params.gamma*cells_df['sum_conduct']


class Params:
    def __init__(self):
        # miscellaneous grid and model run parameters
        line_input = []        
        input_file = open('parameters.txt','r')
        for line in input_file: line_input.append(line.split())
        input_file.close()
        self.fixed_head_area = float(line_input[0][1])              # default planar area of fixed-head cells (very large)
        self.hull_area = float(line_input[1][1])                    # default area of cells on convex hull
        self.min_del = float(line_input[2][1])                      # minimum node-to-node distance (during grid generation)
        self.num_refine = int(line_input[3][1])                     # levels of grid refinement (via Deleanay triangulation)
        self.num_quad_level = int(line_input[4][1])                 # levels of quad delineation of nodes (for better node numbering)
        self.interp_q = int(line_input[5][1])                       # interpolate recharge to newly-generated nodes
        self.show_map = int(line_input[6][1])                       # show node map after grid generation, but before running flow model
        self.gamma = float(line_input[7][1])                        # central-difference time stepping weighting factor
        self.dt_init = float(line_input[8][1])                      # initial time step
        self.dt_min = float(line_input[9][1])                       # minimum time step (if dt is less than this value, simulation will abort)
        self.dt_max = float(line_input[10][1])                      # maximum time step
        self.t_end = float(line_input[11][1])                       # end-time of simulation
        self.dh_max = float(line_input[12][1])                      # maximum head chnage (anywhere in grid) per time step
        self.dt_decrease = float(line_input[13][1])                 # factor by which to decrease dt if max_dh is violated
        self.dt_increase = float(line_input[14][1])                 # factor by which to increase time step unitl dt_max is reached


### support functions ###


def Distance(u, v):
    # return distance between [u1, u2] and [v1, v2]; can be vectorized
    # credit goes to: http://stackoverflow.com/questions/17936587/in-numpy-find-euclidean-distance-between-each-pair-from-two-arrays, user2357112
    d = (u-v)**2
    d = d.sum(axis=-1)
    return sqrt(d)
    
    
def Mesher(tri, params):
    # return midpoints along connecting triangles to use as new node points
    nodes_tree = KDTree(tri.points)                                         # set up KDTree for input points
    elements = tri.simplices                                                # indices of points forming triangles
    doublets_0 = array([elements[:, 0], elements[:, 1]]).T                  # convert to connecting vertex pairs
    doublets_1 = array([elements[:, 1], elements[:, 2]]).T
    doublets_2 = array([elements[:, 2], elements[:, 0]]).T
    doublets = concatenate((doublets_0, doublets_1, doublets_2), axis=0)
    middles = tri.points[doublets].mean(axis=1)                             # midpoints between vertices
    middles_tuple = tuple([tuple(row) for row in middles])                  # remove duplicates by converting to hashable tuples and then convert back to arrays
    unique_middles_tuple = set(middles_tuple)
    unique_middles = array([array(row) for row in unique_middles_tuple])
    d_near,i_near = nodes_tree.query(unique_middles)                        # remove those new points that are too close to prior existing points
    new_pts = unique_middles[array(d_near) >= params.min_del]
    return new_pts


def Splitter(segment):
    # return midpoints along a segmented boundary
    nodes_tree = KDTree(segment)                                            # set up KDTree for input points    
    d_near,i_near = nodes_tree.query(segment, k=2, eps=0, p=2, distance_upper_bound=inf)
    nearest_pt = i_near[:, 1]
    x_mids = 0.5 * (segment[:, 0] + segment[nearest_pt, 0])
    y_mids = 0.5 * (segment[:, 1] + segment[nearest_pt, 1])
    middles = array([x_mids, y_mids]).T
    middles_tuple = tuple([tuple(row) for row in middles])                  # remove duplicates by converting to hashable tuples and then convert back to arrays
    unique_middles_tuple = set(middles_tuple)
    unique_middles = array([array(row) for row in unique_middles_tuple])
    return unique_middles


def PurgeNeighbors(points, params):
    # randomly remove points in a set that are too close to other points
    purged = False
    while not purged:
        nodes_tree = KDTree(points)
        d_near,i_near = nodes_tree.query(points, k=2, eps=0, p=2, distance_upper_bound=inf)
        nearest_pt = i_near[:, 1]
        nearest_d = d_near[:, 1]
        close_mark = (nearest_d < params.min_del) * 1
        if sum(close_mark) == 0:
            purged = True
        else:
            pt_delete = (where(close_mark == 1)[0])
            points = delete(points, pt_delete[0], 0)
    return points


def Gridder(params):

    # read in point data and use data frame(s) to interpolate/fill in grid points
    cells_df = read_csv('wells.txt', sep='\t')

    for i in xrange(params.num_refine):

        print 'Grid refinement step no.',i+1

        # create new internal points (those not part of a segmented boundary)
        pts_internal_df = cells_df[(cells_df['fixed'] < 2)]
        pts_internal_locs = array([pts_internal_df['x'], pts_internal_df['y']])
        tri = Delaunay(pts_internal_locs.T)
        new_pts_internal = Mesher(tri, params)

        # remove new_pts_internal that are too close to other points
        new_pts_internal = PurgeNeighbors(new_pts_internal, params)

        # interpolate fixed-head flag (it's boolean, so use nearest-neighbor method); initialize data frame to hold results
        fixed = griddata(pts_internal_locs.T, pts_internal_df['fixed'], new_pts_internal, method='nearest')
        new_cells = {'x':new_pts_internal.T[0], 'y':new_pts_internal.T[1], 'fixed':array(fixed)}
        new_cells_df = DataFrame(new_cells)
        new_cells_df = new_cells_df[['x', 'y', 'fixed']]

        # create new segmented boundary points
        pts_segment_df = cells_df[(cells_df['fixed'] == 2)]
        pts_segment_locs = array([pts_segment_df['x'], pts_segment_df['y']])
        segment = pts_segment_locs.T    
        new_pts_segment = Splitter(segment)

        # add new segmented boundary points to data frame
        segment_matrix = array([new_pts_segment.T[0], new_pts_segment.T[1], zeros(len(new_pts_segment), int) + 2])
        new_cells_df = new_cells_df.append(DataFrame(segment_matrix.T, columns=['x', 'y', 'fixed']))
        new_cells_df.reset_index(drop = True, inplace = True)

        # interpolate to estimate K, Ss, b, and h at new node locations
        pts_all = array([cells_df['x'], cells_df['y']])
        pts_new = array([new_cells_df['x'], new_cells_df['y']])
        new_cells_df['K'] = griddata(pts_all.T, cells_df['K'], pts_new.T, method='linear')
        new_cells_df['Ss'] = griddata(pts_all.T, cells_df['Ss'], pts_new.T, method='linear')
        new_cells_df['b'] = griddata(pts_all.T, cells_df['b'], pts_new.T, method='linear')
        new_cells_df['h'] = griddata(pts_all.T, cells_df['h'], pts_new.T, method='linear')    

        # distributed recharge terms and wells
        if params.interp_q: new_cells_df['q'] = griddata(pts_all.T, cells_df['q'], pts_new.T, method='linear')
        else: new_cells_df['q'] = 0.
        new_cells_df['Q'] = 0.              # this is always node-specific for initial data points

        # tracking flag (time-series tracking is only applied to marked initial data points)
        new_cells_df['track'] = 0

        # merge with earlier cells set iteration
        new_cells_df = new_cells_df[['x', 'y', 'b', 'h', 'K', 'Ss', 'Q', 'q', 'fixed', 'track']]
        cells_df = cells_df.append(new_cells_df)
        cells_df.reset_index(drop = True, inplace = True)

    return cells_df


def Quadrants(nodes_df):
    # divide set of nodes into quadrants
    Xm = nodes_df['x'].median()
    Ym = nodes_df['y'].median()    
    quad = zeros(len(nodes_df), float)
    quad = quad + ((nodes_df['x']<=Xm)&(nodes_df['y']>Ym))*0 + ((nodes_df['x']>Xm)&(nodes_df['y']>Ym))*1 \
        + ((nodes_df['x']<=Xm)&(nodes_df['y']<=Ym))*2 + ((nodes_df['x']>Xm)&(nodes_df['y']<=Ym))*3
    return quad


def DecToN(num, n):
    # convert integer num from decimal to base-n numbering system (with n <= 9)
    if num:
        answer = []
        while num > 0:
            num, r = divmod(num, n)
            answer.insert(0,str(r))
        return int(''.join(answer))
    else: return 0


def SortNodes(cells_df, params):
    # heirarchically split the grid into quadrants to group nearby cells together - this will reduce sparse matrix bandwidth
    for i in xrange(params.num_quad_level):
        if not i: cells_df['quad'] = Quadrants(cells_df)      # initialize quadrant assignments (top level)
        else:
            for j in xrange(4**i):
                # step through quadrants defined at previous level (i-1)
                quad_index = DecToN(j, 4) * 10.**(-(i-1))
                # the next statement selectively only operates on those df rows with 'quad' = current definition of quad_index
                cells_df.ix[cells_df['quad']==quad_index, 'quad'] += 10.**(-i) * Quadrants(cells_df.ix[cells_df['quad']==quad_index])
    # sort by 'quad' and re-index; new indices (row numbers) = node numbers that will be used to define connections and populate sparse matrix
    cells_df = cells_df.sort('quad')
    cells_df.reset_index(drop = True, inplace = True)    
    return cells_df


def PolygonArea(p):
    # classic cross-product planimeter approach for calculating area of polygon by vertices (p = array of points)
    pc = append(p,[p[0]],axis=0)
    s = (pc[1:-1,0]*(pc[2:,1]-pc[0:-2,1])).sum()
    s += pc[-1,0]*(pc[1,1]-pc[-2,1])
    return abs(s/2)


def ProcessGeometry(vor, cells_df, node_locs, params):
    # compute planar areas and volumes of cells
    area = empty(len(cells_df), float)
    for i, p_region in enumerate(vor.point_region):
        vert_index_set = vor.regions[p_region]
        vert_coord_set = vor.vertices[vert_index_set]
        area[i] = PolygonArea(vert_coord_set) 
    cells_df['plan_area'] = area
    cells_df['Q'] += cells_df['q'] * area                                                       # lump areal recharge with pumping wells into single term
    hull = ConvexHull(node_locs)
    for hull_cell in hull.vertices: cells_df.loc[hull_cell, 'plan_area'] = params.hull_area     # assign default hull cell areas
    cells_df['plan_area'] += (cells_df['fixed'] > 0) * params.fixed_head_area                   # mark constant head cells
    cells_df['vol'] = cells_df['plan_area'] * cells_df['b']                                     # calculate cell volumes
    return cells_df


def SumConductances(cells_df, connections):
    # return the sum of the conductance terms associated with each cell
    sum_conduct = zeros(len(cells_df), float)
    for i in xrange(len(connections.nodes)):
        node_1 = connections.nodes[i, 0]
        node_2 = connections.nodes[i, 1]
        sum_conduct[node_1] += connections.conduct[i]
        sum_conduct[node_2] += connections.conduct[i]
    return sum_conduct


def RHS_Vector(cells_df, connections):
    # construct explicit matrix (run for each time step)
    b = cells_df['Q'] - cells_df['sum_conduct'] * cells_df['h']
    for i in xrange(len(connections.nodes)):
        node_1 = connections.nodes[i, 0]
        node_2 = connections.nodes[i, 1]
        b[node_1] += connections.conduct[i] * cells_df['h'][node_2]
        b[node_2] += connections.conduct[i] * cells_df['h'][node_1]
    return b


def BuildGrid(params):
    
    # read in point data and populate unstructured grid
    cells_df = Gridder(params)

    # optimize cell/node order
    cells_df = SortNodes(cells_df, params)
    node_locs = array([cells_df['x'], cells_df['y']]).T
    print 'Processed node numbering scheme.'

    # calculate Voronoi polygons
    vor = Voronoi(node_locs)
    if params.show_map:
        voronoi_plot_2d(vor)
        plt.show()
    cells_df['vor_region'] = vor.point_region           # note index number of Voronoi region associated with each node
    print 'Created Voronoi diagram; computed interfacial areas and connections.'

    # populate Connections object, which holds arrays for node-to-node connections,including interfacial areas and distances
    connections = Connections(vor, cells_df)
    cells_df['sum_conduct'] = SumConductances(cells_df, connections)
    print 'Mapped cell connections.'

    # process cell geometry
    cells_df = ProcessGeometry(vor, cells_df, node_locs, params)
    print 'Processed internal cell geometry.'
    cells_df.to_csv('init_cells.csv')               # write grid setup and initial conditions to output file

    return cells_df, connections


def FlowModel(params, cells_df, connections):

    # initialize ...
    print 'Running flow model ...'
    t = 0.
    dt = params.dt_init
    track_cells = array(where(cells_df['track']==1))[0]             # note the cell numbers for nodes that serve as monitor wells
    monitor_well = []
    t_monitor = []

    # set up left-hand-side of equation set with LHS_Matrix class; diagonal terms will change as the time step changes
    implicit = LHS_Matrix(cells_df, connections, params, dt)
    A = csr_matrix( (array(implicit.data),(array(implicit.row_index), array(implicit.col_index))), shape=(len(cells_df), len(cells_df)) )

    while t < params.t_end:

        converged = False

        while converged == False:

            # modify diganonal according to time step size
            implicit.Diagonal(cells_df, params, dt)
            
            # update sparse equation matrix
            A = csr_matrix( (array(implicit.data),(array(implicit.row_index), array(implicit.col_index))), shape=(len(cells_df), len(cells_df)) )

            # construct explicit vector
            b = RHS_Vector(cells_df, connections)

            # solve equations
            dh = spsolve(A,b)

            # determine convergence (all dH must be <= params.dh_max
            converged = 1 - sign(sum((abs(dh) > params.dh_max) * 1))
            if not converged: dt *= params.dt_decrease
            assert(dt > params.dt_min)

        # update values
        t += dt
        
        t_monitor.append(t)
        cells_df['h'] += dh

        # append to time series data frame
        row = []
        for monitor in track_cells:
            row.append(cells_df['h'][monitor])
        monitor_well.append(row)
        
        # update time step
        dt *= params.dt_increase
        dt = min(dt, params.dt_max, params.t_end - t)

    # write final head distribution to file
    cells_df.to_csv('final_cells.csv')              

    # write time series data frame to file
    t_monitor = array(t_monitor)
    monitor_well = array(monitor_well)
    monitor_df = DataFrame(monitor_well)
    col_label = []
    for i in track_cells: col_label.append('cell_' + str(i))
    monitor_df.columns = col_label
    monitor_df.insert(0, 'time', t_monitor)
    monitor_df.to_csv('monitor_df.csv')
   
    return


### script ###


def USG():

    # read model parameters
    params = Params()

    # generate model grid
    cells_df, connections = BuildGrid(params)

    # run groundwater flow model
    FlowModel(params, cells_df, connections)

    print 'Completed.'


### run script ###
USG()
























