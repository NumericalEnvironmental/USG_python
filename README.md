# USG_python
A python script that generates an (irregular) numerical grid from (irreguar) scattered data points with assigned properties, then subsequently runs a single=layer, confined groundwater flow model, per instructions. Input files include:

(1) wells.txt - initial scattered properties file, including:

x --> x-location

y --> y-location
b --> aquifer thickness
h --> initial head
K --> hydraulic conductivity
Ss --> specific storage
Q --> extraction(-) or injection(+) rate, volume/time
q --> recharge, length/time
fixed --> 0=variable-head cell, 1=regular fixed-head cell, 2=fixed-head cell that is part of a linear segment (special interpolation)
track --> mark this cell as a monitor well (for time series output) with a "1"

(2) parameters.txt - miscellaneous gridding and model run parameters, as indicated (see internal code comments)

More info can be found here: (link will be posted shortly)
