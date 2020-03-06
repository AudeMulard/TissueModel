from creation_network import Network
from force_balance import *
from tensile_test import *
from scipy.spatial import Voronoi, voronoi_plot_2d
import matplotlib.pyplot as plt
from Several_networks_test import *
from network_plotting import *
import os
from datetime import date

## PARAMETERS
dimension=2 #dimension of the problem
complexity_network=50 #number of random seed points
length_domain=1.0
min_distance = 0.0001*length_domain
space_discretization = 0.01*length_domain
Ef=1.0
A=1.
disturbance=0.02
traction_distance = 0.1*length_domain
#iteration = 15



data_path = '../Data/influence_points/'

today = date.today()

new_dir = data_path+today.strftime("%b-%d-%Y")+'_'+'%04d' % len(os.listdir(data_path))
os.mkdir(new_dir)
path = new_dir



## EXPERIMENT
creation="disturbed_grid"
constitutive = 'linear2'
scheme='nonlinear'
side = 'right'
plot = False
video = False
phase = 'only_one'
stress_rep = True
details = True


network = Network(dimension, complexity_network, length_domain, min_distance, Ef, A, disturbance, creation, path)


network = network.set_fibers(creation, path)

print len(network.vertices), len(network.ridge_vertices),len(network.ridge_vertices)-2*len(network.interior_nodes)
#print network.interior_nodes
plot_network(network)
plt.show()

test_1 = Tensile_test(constitutive, scheme, side, space_discretization, traction_distance, plot, video, 0, path)

network = test_1.full_test(network, path,details)

"""
number_networks = 10

test_1 = Tensile_test(constitutive, scheme, side, space_discretization, traction_distance, plot, video, 0, path)

plot_impact_number_points(number_networks,dimension, length_domain, min_distance, Ef, A, B, creation, test_1, path,False)

"""

