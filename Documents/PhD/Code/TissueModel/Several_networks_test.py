
from scipy.spatial import Voronoi, voronoi_plot_2d
import matplotlib.pyplot as plt
import numpy


def standard_dev_complexity(complexity_network, number_networks,dimension, length_domain, min_distance, Ef, A, B, creation, tensile_test_1, path, end_strain, number_points,details):
	try:
		network = one_network_test(complexity_network,dimension, length_domain, min_distance, Ef, A, B, creation,tensile_test_1, path,details)
		end_strain.append(network.stress[-1])
		number_points.append(len(network.vertices))
	except (numpy.linalg.linalg.LinAlgError, ValueError):
		pass
	return end_strain, number_points
	

def one_network_test(complexity_network,dimension, length_domain, min_distance, Ef, A, B, creation, tensile_test_1, path,details):
	network = Network(dimension, complexity_network, length_domain, min_distance, Ef, A, B, creation, path)
	network = network.set_fibers(creation, path)
	network = tensile_test_1.full_test(network,path,details)
	return network

def plot_impact_number_points(number_networks,dimension, length_domain, min_distance, Ef, A, B, creation, tensile_test_1, path,details):
	end_strain = []
	number_points = []
	for complexity_network in range(200,500):
		print('Complexity:', complexity_network)
		end_strain, number_points = standard_dev_complexity(complexity_network,number_networks,dimension, length_domain, min_distance, Ef, A, B, creation, tensile_test_1, path,end_strain, number_points,details)
		with open(os.path.join(path,'strain_points.csv'), 'w') as writeFile:
			writer = csv.writer(writeFile)
			writer.writerows([number_points])
			writer.writerows([end_strain])

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
complexity_network=150 #number of random seed points
length_domain=1.0
min_distance = 0.0001*length_domain
space_discretization = 0.01*length_domain
Ef=1.0
A=1.
disturbance=0.02
traction_distance = 0.3*length_domain
#iteration = 15



data_path = '../Data/influence_points/'

today = date.today()

new_dir = data_path+today.strftime("%b-%d-%Y")+'_'+'%04d' % len(os.listdir(data_path))
os.mkdir(new_dir)
path = new_dir



## EXPERIMENT
creation="Voronoi"
constitutive = 'linear2'
scheme='nonlinear'
side = 'right'
plot = False
video = True
phase = 'only_one'
stress_rep = True
details = True

number_networks = 450
test_1 = Tensile_test(constitutive, scheme, side, space_discretization, traction_distance, plot, video, phase, path)
plot_impact_number_points(number_networks,dimension, length_domain, min_distance, Ef, A, disturbance, creation, test_1, path,details)
