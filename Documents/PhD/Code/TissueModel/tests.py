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
import fnmatch

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
_plot = True
video = True
phase = 'only_one'
stress_rep = True
details = True

list_modes = [['Voronoi','random'],['Voronoi','regular'],['growth_network','grid']]#,['Voronoi','grid']

for mode in list_modes:
	creation = mode[0]
	generation = mode[1]
	# Make results directory

	data_path = '../Data/Testing_%s_%s/' % (creation,generation)

	try:
		os.makedirs(data_path)
	except OSError:
		if not os.path.isdir(data_path):
			raise

	today = date.today()

	new_dir = data_path+today.strftime("%b-%d-%Y")+'_'+'%04d' % len(os.listdir(data_path))
	os.mkdir(new_dir)
	path = new_dir

	with open(os.path.join(path,'testing.txt'), 'a') as writeFile:
		writeFile.write('%s \n' % mode)
	############################### REPEATABILITY ######################################
	
	# Repeatability for same network
	with open(os.path.join(path,'testing.txt'), 'a') as writeFile:
		writeFile.write('Repeatability for the same network')
	from Testing.Identical_network import *
	network = Network(dimension, complexity_network, length_domain, min_distance, k_tension, k_compression, A, disturbance, creation,generation, path)
	test_1 = Tensile_test(constitutive, side, space_discretization, traction_distance, plot, video, path,details)
	network_def_tests(network,path,test_1)
	identical_check(path)
	
	# Repeatability for different networks, same parameters
	with open(os.path.join(path,'testing.txt'), 'a') as writeFile:
		writeFile.write('Repeatability for different networks, same parameters')
	from Testing.Identical_parameters import *
	network = Network(dimension, complexity_network, length_domain, min_distance, k_tension, k_compression, A, disturbance, creation,generation, path)
	test_1 = Tensile_test(constitutive, side, space_discretization, traction_distance, _plot, video, path,details)
	network_def_idparams(network,path,test_1)
	idparams_check(path)
	
	# Elasticity test
	with open(os.path.join(path,'testing.txt'), 'a') as writeFile:
		writeFile.write('Elasticity test')
	from Testing.Elasticity import *
	network = Network(dimension, complexity_network, length_domain, min_distance, k_tension, k_compression, A, disturbance, creation,generation, path)
	test_1 = Tensile_test(constitutive, side, space_discretization, traction_distance, _plot, video, path,details)
	network_def_elasticity(network,path,test_1)
	elasticity_check(path)

	# Add the change of step size
	with open(os.path.join(path,'testing.txt'), 'a') as writeFile:
		writeFile.write('Step size change test')
	from Testing.Step_size import *
	network = Network(dimension, complexity_network, length_domain, min_distance, k_tension, k_compression, A, disturbance, creation,generation, path)
	test_1 = Tensile_test(constitutive, side, space_discretization, traction_distance, _plot, video, path,details)
	network_def_tests(network,path,test_1)
	identical_check(path)
	

	############################### BIOLOGICAL MEANING ######################################
	with open(os.path.join(path,'testing.txt'), 'a') as writeFile:
		writeFile.write('Biological meaning')
	for i in range(10):
		network = Network(dimension, complexity_network, length_domain, min_distance, k_tension, k_compression, A, disturbance, creation,generation, path)
		network = network.set_fibers(path)
		test_1 = Tensile_test(constitutive, side, space_discretization, traction_distance, _plot, video, path,details)
		network = test_1.full_test(network, path,test_1.details,name='%02d' % i)

	# Comparison of stress-strain curve to experiments + Evolution of geometry: look at orientations of fibers

	from Testing.Biological_meaning import *

	check_bio(path)

