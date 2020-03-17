from Core_calculation.creation_network import Network
from Core_calculation.force_balance import *
from Core_calculation.tensile_test import *
from scipy.spatial import Voronoi, voronoi_plot_2d
import matplotlib.pyplot as plt
#from Several_networks_test import *
from Plotting.network_plotting import *
import os
from datetime import date

## PARAMETERS
dimension=2 #dimension of the problem
complexity_network=50 #number of random seed points
length_domain=1.0
min_distance = 0.0001*length_domain
space_discretization = 0.001*length_domain
Ef=1.0
A=1.
disturbance=0.02
traction_distance = 0.1*length_domain
#iteration = 15



data_path = '../Data/reg_Voronoi/'

today = date.today()

new_dir = data_path+today.strftime("%b-%d-%Y")+'_'+'%04d' % len(os.listdir(data_path))
os.mkdir(new_dir)
path = new_dir



## EXPERIMENT
creation="reg_Voronoi"
constitutive = 'linear2'
scheme='nonlinear'
side = 'right'
plot = True
video = True
phase = 'only_one'
stress_rep = True
details = True


network = Network(dimension, complexity_network, length_domain, min_distance, Ef, A, disturbance, creation, path)

network = network.set_fibers(creation, path)
plot_geometry(network)
plt.show()
print len(network.ridge_vertices)-2*len(network.vertices), len(network.ridge_vertices), len(network.vertices)
#suppriemr tous les points au desus de la ligne sans raccamoder puis supprimer toutes les alones fibres. Dans le cut updown. le cut side works.

test_1 = Tensile_test(constitutive, scheme, side, space_discretization, traction_distance, plot, video, phase, path)

network = test_1.full_test(network, path,details)


