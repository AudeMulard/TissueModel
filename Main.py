from Network_generation.creation_network import Network
from Core_calculation.tensile_test import *
from Plotting.information_network import *
import matplotlib.pyplot as plt
from Plotting.network_plotting import *
import os
from datetime import date


sys.path.append('C:/Users/am2548/TissueModel/Documents/PhD/Code/TissueModel/') # solving package
sys.path.append('D:/twoD/plotting_Codes/') # plotting package

## PARAMETERS
generation = 'Bridson_sampling' # method to generate initial points
creation='Voronoi' # geometry of the network
dimension=2 #dimension of the problem
complexity_network=1900 #number of random seed points
length_domain=(1.0,1.0,1.0) # this change won't work for growth network, will need to be changed
min_distance = 7.5e-3*length_domain[0] #minimum distance between points

beam_Young = 33000 # Young modulus of the fibers
beam_poisson = 0.3 # Poisson modulus of the fibers
beam_profile = 0.01 # diameter of the fibers
connector_coeff = 0.0001 # connector coefficient between fibers

traction_distance = length_domain[0] # traction distance of the test
hyperstatic_param = 0 # ratio of number of fibers compared to number of nodes
element_size = 0.0025 # approximate length of the element in the finite element method
side = 'right'
space_discretization = traction_distance/100. # intervals to record data during the simulation
material_density=5 # material density for the explicit simulation

### SET DATA FILE
data_path = '../Data/default/'
today = date.today()
new_dir = data_path+today.strftime("%b-%d-%Y")+'_'+'%04d' % len(os.listdir(data_path))
os.mkdir(new_dir)
path = new_dir
current_dir = os.getcwd()

import time
### CREATE NETWORK
"""
network = Network(dimension, complexity_network, length_domain, min_distance, beam_Young, beam_poisson, beam_profile, connector_coeff, hyperstatic_param, creation, generation, path)
network = network.set_fibers(path)
plot_geometry(network);plt.show() #plots the geometry of the network

### CREATE TEST
test_1 = Tensile_test(side, space_discretization, traction_distance, element_size, '.')
test_1.save_parameters(network,path)

### LAUNCH TEST
os.system("abaqus cae script=SOLVER.py") #noGUI: no graphical interface of abaqus. if script instead, the grpahical interface of abaqus opens
"""

complexity = [50,100,200,400,500,600,800,900,1300,1500,2000,2500]
#hyperstatic_param = [2.8,3.]
for i in range(len(complexity)):
	print(complexity[i])
	network = Network(dimension, complexity[i], length_domain, min_distance, beam_Young, beam_poisson, beam_profile,connector_coeff, hyperstatic_param, creation, generation, path)
	network = network.set_fibers(path)
	length = []
	network = network.create_ridge_node_list()
	for i in range(len(network.list_nodes_ridges)):
		length.append(len(network.list_nodes_ridges[i]))
	print(np.mean(length))

"""
for k in range(4):
	generation = 'Bridson_sampling' # method to generate initial points
	creation='growth_network'
	dimension = 2
	for complexity_network in [4000]:
		network = Network(dimension, complexity_network, length_domain, min_distance, beam_Young, beam_poisson, beam_profile,connector_coeff, hyperstatic_param, creation, generation, path)
		network = network.set_fibers(path)
		print(len(network.vertices))
		element_size = max(network.min_distance,0.005)
		print(element_size)
		test_1 = Tensile_test(side, space_discretization, traction_distance, element_size, path)
		test_1.save_parameters(network,path)
		os.chdir('../Data/')# insert where the network geometry was recorded
		os.chdir(sorted_ls('.')[-1])
		os.chdir(sorted_ls('.')[-1])
		filenames=sorted(fnmatch.filter(sorted_ls('.'), 'network_vertices_initial_*.csv'))
		print(int(filenames[-1][-9:-4]))
		if int(filenames[-1][-9:-4])!=int(len(network.vertices)):
			time.sleep(30)
		filenames=sorted(fnmatch.filter(sorted_ls('.'), 'network_vertices_initial_*.csv'))
		print(int(filenames[-1][-9:-4]))
		os.chdir(current_dir)
		os.system("abaqus cae noGUI=SOLVER.py")


for k in range(5):
	generation = 'Bridson_sampling' # method to generate initial points
	creation='growth_network' # geometry of the network
	dimension=2 #dimension of the problem
	complexity_network=1900 #number of random seed points
	#complexity_network=50
	network = Network(dimension, complexity_network, length_domain, min_distance, beam_Young, beam_poisson, beam_profile,connector_coeff, hyperstatic_param, creation, generation, path)
	network = network.set_fibers(path)
	print(len(network.vertices))
	#element_size = network.min_distance/3.
	element_size = 0.01
	test_1 = Tensile_test(side, space_discretization, traction_distance, element_size, path)
	test_1.save_parameters(network,path)
	os.system("abaqus cae noGUI=SOLVER.py")

for k in range(5):
	generation = 'Bridson_sampling' # method to generate initial points
	creation='Voronoi' # geometry of the network
	dimension=2 #dimension of the problem
	complexity_network=1100 #number of random seed points
	#complexity_network=50
	network = Network(dimension, complexity_network, length_domain, min_distance, beam_Young, beam_poisson, beam_profile,connector_coeff, hyperstatic_param, creation, generation, path)
	network = network.set_fibers(path)
	print(len(network.vertices))
	element_size = 0.005
	#element_size = network.min_distance/3.
	test_1 = Tensile_test(side, space_discretization, traction_distance, element_size, path)
	test_1.save_parameters(network,path)
	os.system("abaqus cae noGUI=SOLVER.py")

for k in range(10):
	generation = 'Bridson_sampling' # method to generate initial points
	creation='Voronoi' # geometry of the network
	dimension=3 #dimension of the problem
	complexity_network=1900 #number of random seed points
	#complexity_network=50
	network = Network(dimension, complexity_network, length_domain, min_distance, beam_Young, beam_poisson, beam_profile,connector_coeff, hyperstatic_param, creation, generation, path)
	network = network.set_fibers(path)
	print(len(network.vertices))
	#element_size = network.min_distance/3.
	element_size = 0.005
	test_1 = Tensile_test(side, space_discretization, traction_distance, element_size, path)
	test_1.save_parameters(network,path)
	os.system("abaqus cae noGUI=SOLVER.py")




for k in range(10):
	generation = 'Bridson_sampling' # method to generate initial points
	creation='growth_network' # geometry of the network
	dimension=3 #dimension of the problem
	complexity_network=8000 #number of random seed points
	#complexity_network=50
	network = Network(dimension, complexity_network, length_domain, min_distance, beam_Young, beam_poisson, beam_profile,connector_coeff, hyperstatic_param, creation, generation, path)
	network = network.set_fibers(path)
	print(len(network.vertices))
	#element_size = network.min_distance/3.
	element_size = 0.005
	test_1 = Tensile_test(side, space_discretization, traction_distance, element_size, path)
	test_1.save_parameters(network,path)
	os.system("abaqus cae noGUI=SOLVER.py")


generation = 'Bridson_sampling' # method to generate initial points
creation='Voronoi' # geometry of the network
dimension=2 #dimension of the problem
complexity_network=300 #number of random seed points
#complexity_network=50
network = Network(dimension, complexity_network, length_domain, min_distance, beam_Young, beam_poisson, beam_profile,connector_coeff, hyperstatic_param, creation, generation, path)
network = network.set_fibers(path)
print(len(network.vertices))
element_size = 0.01
test_1 = Tensile_test(side, space_discretization, traction_distance, element_size, path)
test_1.save_parameters(network,path)
print(0.1**2/network.beam_profile**2 * np.sqrt(2*material_density/(network.beam_young*np.pi)) *10)
os.system("abaqus cae script=SOLVER.py")

generation = 'Bridson_sampling' # method to generate initial points
creation='growth_network' # geometry of the network
dimension=2 #dimension of the problem
complexity_network=300 #number of random seed points
#complexity_network=50
network = Network(dimension, complexity_network, length_domain, min_distance, beam_Young, beam_poisson, beam_profile,connector_coeff, hyperstatic_param, creation, generation, path)
network = network.set_fibers(path)
print(len(network.vertices))
element_size = network.min_distance/3.
test_1 = Tensile_test(side, space_discretization, traction_distance, element_size, path)
test_1.save_parameters(network,path)
print(tensile_test.element_size**2/network.beam_profile**2 * np.sqrt(2*material_density/(network.beam_young*np.pi)) *10)
os.system("abaqus cae noGUI=SOLVER.py")
"""