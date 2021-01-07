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
complexity_network=100 #number of random seed points
length_domain=(1.0,1.0,1.0)
min_distance = 0.0001*length_domain[0]
space_discretization = 0.1*length_domain[0]
beam_Young = 1000
beam_poisson = 0.3
beam_profile = 0.1
connector_coeff = 1.0
disturbance=0.02
traction_distance = 1.0*length_domain[0]
hyperstatic_param = dimension
element_size = 0.001


### SET DATA FILE

path = '../Data/Study_networks/understanding/'
#data_path = '../Data_1/Study_networks/Nov-17-2020_0073/'
today = date.today()



current_dir = os.getcwd()

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

test_1 = Tensile_test(constitutive, side, space_discretization, traction_distance,element_size, plot, video, path,details)
os.chdir(path)
#filenames =  fnmatch.filter(sorted_ls('.'), 'network_vertices_initial_*.csv')
network = load_network_info(912)
os.chdir(current_dir)

for connector_coeff in [0.006]:#,0.005,0.004,0.003]:
	#network = load_network_info(int(filename[-9:-4]))
	network.connector_coeff = connector_coeff
	network.save_network('initial',path)
	test_1.save_parameters(network,path)
	os.system("abaqus cae script=new_solver_2.py")

