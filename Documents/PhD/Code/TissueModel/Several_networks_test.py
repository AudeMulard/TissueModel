from creation_network import Network
from force_balance import *
from tensile_test import *
from scipy.spatial import Voronoi, voronoi_plot_2d
import matplotlib.pyplot as plt
import numpy


def standard_dev_complexity(complexity_network, number_networks,dimension, length_domain, min_distance, Ef, A, B, creation, tensile_test_1, path):
	end_strain = []
	for i in range(number_networks):
		print 'Network', i
		try:
			network = one_network_test(complexity_network,dimension, length_domain, min_distance, Ef, A, B, creation,tensile_test_1, path)
			if str(network.stress[tensile_test_1.iterations]) != 'nan':
				print network.stress[tensile_test_1.iterations]
				end_strain.append(network.stress[tensile_test_1.iterations])
				network.save_network('final_%s_%s' % (complexity_network,i), path)
			plt.close()
		except (numpy.linalg.linalg.LinAlgError, ValueError):
			print('Singular Matrix')
			pass
	return numpy.mean(end_strain), numpy.std(end_strain)
	

def one_network_test(complexity_network,dimension, length_domain, min_distance, Ef, A, B, creation, tensile_test_1, path):
	network = Network(dimension, complexity_network, length_domain, min_distance, Ef, A, B, creation, path)
	network = network.create_network(creation)
	network = network.set_fibers(creation, path)
	print len(network.vertices)
	network = tensile_test_1.full_test(network,path)
	return network

def plot_impact_number_points(number_networks,dimension, length_domain, min_distance, Ef, A, B, creation, tensile_test_1, path):
	mean = []
	standard_dev = []
	for complexity_network in [30,40,50,60, 70]:
		length_domain = complexity_network/10.
		tensile_test_1.traction_distance = 0.3*length_domain
		print 'Complexity:', complexity_network
		mean1, standard_dev1=standard_dev_complexity(complexity_network,number_networks,dimension, length_domain, min_distance, Ef, A, B, creation, tensile_test_1, path)
		mean.append(mean1)
		standard_dev.append(standard_dev1)
		print mean, standard_dev
	return mean, standard_dev


