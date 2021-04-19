from Network_generation.creation_network import Network
from Network_generation.add_ons_network import Cell
from Core_calculation.tensile_test import *
from Plotting.information_network import *
from scipy.spatial import Voronoi, voronoi_plot_2d
import matplotlib.pyplot as plt
#from Several_networks_test import *
from Plotting.network_plotting import *
import os
from datetime import date
import matplotlib.patches as patches
import cProfile
import numpy 
## PARAMETERS
dimension=2 #dimension of the problem
complexity_network=200 #number of random seed points
length_domain=(1.0,1.0,1.0)
min_distance = 0.0001*length_domain[0]
space_discretization = 0.1*length_domain[0]
beam_Young = 33000
beam_poisson = 0.3
beam_profile = 0.01
disturbance=0.02
traction_distance = 1.0*length_domain[0]
hyperstatic_param = 2
element_size = 0.0005
connector_coeff = 0.01

### SET DATA FILE

data_path = '../Data_1/default/'

today = date.today()

new_dir = data_path+today.strftime("%b-%d-%Y")+'_'+'%04d' % len(os.listdir(data_path))
os.mkdir(new_dir)
path = new_dir



### CREATE NETWORK
creation="growth_network"
generation = 'grid'
constitutive = 'spring'
side = 'right'
plot = True
video = True
phase = 'only_one'
stress_rep = True
details = True





network = Network(dimension, complexity_network, length_domain, min_distance, beam_Young, beam_poisson, beam_profile, connector_coeff, disturbance, hyperstatic_param, creation, generation, path)
network = network.set_fibers(path)
plot_geometry(network)
plt.show()
test_1 = Tensile_test(constitutive, side, space_discretization, traction_distance,element_size, plot, video, path,details)
test_1.save_parameters(network,path)
os.system("abaqus cae script=new_solver_4.py")

