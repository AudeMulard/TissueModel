from Network_generation.creation_network import Network
from Network_generation.add_ons_network import Cell
from Core_calculation.force_balance import *
from Core_calculation.tensile_test import *
from scipy.spatial import Voronoi, voronoi_plot_2d
import matplotlib.pyplot as plt
#from Several_networks_test import *
from Plotting.network_plotting import *
import os
from datetime import date
import matplotlib.patches as patches
import fnmatch,time

## PARAMETERS
dimension=2 #dimension of the problem
complexity_network=100 #number of random seed points
length_domain=(1.0,1.0,1.0)
min_distance = 0.0001*length_domain[0]
space_discretization = 0.01*length_domain[0]
k_tension=1.0
k_compression = 1.0
A=1.
disturbance=0.02
traction_distance = 0.5*length_domain[0]
#iteration = 15



## EXPERIMENT
creation="Voronoi"
generation = 'random'
constitutive = 'spring'
side = 'right'
plot = True
video = True
phase = 'only_one'
stress_rep = True
details = True

data_path = '../Data/Study_networks/' 

today = date.today()

new_dir = data_path+today.strftime("%b-%d-%Y")+'_'+'%04d' % len(os.listdir(data_path))
os.mkdir(new_dir)
path = new_dir

list_modes = [['Voronoi','random'],['Voronoi','grid'],['Voronoi','regular'],['growth_network','grid']]

test_1 = Tensile_test(constitutive, side, space_discretization, traction_distance, plot, video, path,details)

for mode in list_modes:
	print(mode)
	
	##################### Parameter study (comparison on same network) #########################
	creation = mode[0]
	generation = mode[1]
	"""
	network = Network(dimension, complexity_network, length_domain, min_distance, k_tension, k_compression, A, disturbance, creation,generation, path,name='%s_%s' % (creation, generation))
	network = network.set_fibers(path)

	
	
	#### Changing geometry
	# Changing complexity
	for complexity in [10,50,100,200,500,1000]:
		network = Network(dimension, complexity, length_domain, min_distance, k_tension, k_compression, A, disturbance, creation,generation, path)
		network = network.set_fibers(path)
		network = test_1.full_test(network, path,test_1.details,name='%s_%s_complexity_%04d' % (creation, generation,complexity))

	if mode[1] == 'grid':
		for disturbance in [0.,0.0001,0.001,0.01,0.1,1.]:
			network = Network(dimension, complexity, length_domain, min_distance, k_tension, k_compression, A, disturbance, creation,generation, path)
			network = network.set_fibers(path)
			network = test_1.full_test(network, path,test_1.details,name='%s_%s_disturbance_%04d' % (creation,generation,disturbance))
	"""
	#### Time calculation
	with open(os.path.join(path,'time.txt'), 'a') as writeFile:
		writeFile.write('%s \n' % mode)
		writeFile.write('%d \n' % time.time())
		start = time.time()
		network = Network(dimension, complexity_network, length_domain, min_distance, k_tension, k_compression, A, disturbance, creation,generation, path)
		network = network.set_fibers(path)
		writeFile.write('Generation time: %d \n' % (time.time() - start))
		test_1 = Tensile_test(constitutive, side, space_discretization, traction_distance, plot, video, path,details)
		network = test_1.full_test(network, path,test_1.details,name='%s_%s_time' % (creation,generation))
		writeFile.write('Test time: %d \n' % (time.time() - start))
		writeFile.close()
	"""
	#### Without changing network
	# Change of k_tension, k_compression
	network = Network(dimension, complexity_network, length_domain, min_distance, k_tension, k_compression, A, disturbance, creation,generation, path)
	network = network.set_fibers(network.creation, path)

	for k_compression in [100.,10.,1.,0.1,0.01,0.001]:
		network.k_compression = k_compression
		network.k_tension = k_compression
		network.vertices = np.array(network.vertices_ini)
		network = test_1.full_test(network, path,test_1.details,name='%s_%s_k_%03d' % (creation, generation,k_compression))
	# Change of ratio k_compression on k_tension
	network = Network(dimension, complexity_network, length_domain, min_distance, k_tension, k_compression, A, disturbance, creation,generation, path)
	network = network.set_fibers(network.creation, path)

	for k_compression in [100.,10.,1.,0.1,0.01,0.001]:
		network.k_compression = k_compression
		network.vertices = np.array(network.vertices_ini)
		network = test_1.full_test(network, path,test_1.details,name='%s_%s_ratio_%03d' % (creation, generation,k_compression))
	"""



