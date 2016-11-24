# USG_python
A Python 2.7 script that generates an (irregular) numerical grid from (irreguar) scattered data points with assigned properties, then subsequently runs a single=layer, confined groundwater flow model, per instructions. Input files include:

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

More info can be found here: https://numericalenvironmental.wordpress.com/2016/11/13/an-unstructured-grid-finite-difference-groundwater-flow-model-in-python/

THIS CODE/SOFTWARE IS PROVIDED IN SOURCE OR BINARY FORM "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
