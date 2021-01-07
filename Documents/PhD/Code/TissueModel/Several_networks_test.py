
from scipy.spatial import Voronoi, voronoi_plot_2d
import matplotlib.pyplot as plt
import numpy


def standard_dev_complexity(complexity_network, end_strain, number_points,number_networks,dimension, complexity, length_domain, min_distance, beam_Young, beam_poisson, beam_profile, truss_young, truss_poisson, truss_area,disturbance, hyperstatic_param, creation, generation, path):
	try:
		network = one_network_test(complexity_network,dimension, complexity, length_domain, min_distance, beam_Young, beam_poisson, beam_profile, truss_young, truss_poisson, truss_area,disturbance, hyperstatic_param, creation, generation, path)
		"""current_path=os.getcwd()
		os.chdir('../Data/')
		os.chdir(sorted_ls('.')[-1])
		os.chdir(sorted_ls('.')[-1])
		test_number = fnmatch.filter(sorted_ls('.'), 'network_vertices_01_00_*.csv')[-1][-13:-4]
		strain,stress=stress_strain_curve(test_number,network)
		end_strain.append(stress[-1])
		number_points.append(len(network.vertices))
		os.chdir(current_path)"""
	except (numpy.linalg.linalg.LinAlgError, ValueError):
		pass
	return end_strain, number_points
	

def one_network_test(complexity_network,dimension, complexity, length_domain, min_distance, beam_Young, beam_poisson, beam_profile, truss_young, truss_poisson, truss_area,disturbance, hyperstatic_param, creation, generation, path):
	network = Network(dimension,complexity_network, length_domain, min_distance, beam_Young, beam_poisson, beam_profile, truss_young, truss_poisson, truss_area,disturbance, hyperstatic_param, creation, generation, path)
	network = network.set_fibers(path)
	test_1 = Tensile_test(constitutive, side, space_discretization, traction_distance,element_size, plot, video, path,details)
	test_1.save_parameters(network,path)
	os.system("abaqus cae noGUI=new_solver_1.py")
	return network

def plot_impact_number_points(number_networks,dimension, complexity, length_domain, min_distance, beam_Young, beam_poisson, beam_profile, truss_young, truss_poisson, truss_area,disturbance, hyperstatic_param, creation, generation, path):
	end_strain = []
	number_points = []
	for complexity in range(50,100):
		complexity_network = 10*complexity
		print('Complexity:', complexity_network)
		end_strain, number_points = standard_dev_complexity(complexity_network,end_strain, number_points,number_networks,dimension, complexity, length_domain, min_distance, beam_Young, beam_poisson, beam_profile, truss_young, truss_poisson, truss_area,disturbance, hyperstatic_param, creation, generation, path)
		with open('strain_points.csv', 'w') as writeFile:
			writer = csv.writer(writeFile)
			writer.writerows([number_points])
			writer.writerows([end_strain])
		#os.chdir(current_path)

from Network_generation.creation_network import Network
from Network_generation.add_ons_network import Cell
from Core_calculation.force_balance import *
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
space_discretization = 0.01*length_domain[0]
beam_Young = 33000
beam_poisson = 0.3
beam_profile = 0.01
truss_young = 1000.
truss_poisson = 0.3
truss_area = 0.1
disturbance=0.02
traction_distance = 0.3*length_domain[0]
hyperstatic_param = dimension
element_size = 0.01


data_path = '../Data/validity/'

today = date.today()

new_dir = data_path+today.strftime("%b-%d-%Y")+'_'+'%04d' % len(os.listdir(data_path))
os.mkdir(new_dir)
path = new_dir



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

number_networks = 450
test_1 = Tensile_test(constitutive, side, space_discretization, traction_distance,element_size, plot, video, path,details)
plot_impact_number_points(number_networks,dimension, complexity_network, length_domain, min_distance, beam_Young, beam_poisson, beam_profile, truss_young, truss_poisson, truss_area,disturbance, hyperstatic_param, creation, generation, path)