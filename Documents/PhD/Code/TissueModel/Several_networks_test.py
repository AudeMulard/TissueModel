from creation_network import Network
from force_balance import *
from tensile_test import *
from scipy.spatial import Voronoi, voronoi_plot_2d
import matplotlib.pyplot as plt
import numpy


def standard_dev_complexity(complexity_network, number_networks,dimension, length_domain, min_distance, Ef, A, B, creation, tensile_test_1, path, end_strain, number_points):
	for i in range(number_networks):
		print('Network', i)
		try:
			network = one_network_test(complexity_network,dimension, length_domain, min_distance, Ef, A, B, creation,tensile_test_1, path)
			#if str(network.stress[tensile_test_1.iterations]) != 'nan':
			print(network.stress[-1])
			end_strain.append(network.stress[-1])
			number_points.append(len(network.vertices))
			#network.save_network('final_%s_%s' % (len(network.vertices),i), path)
			#plt.close()
		except (numpy.linalg.linalg.LinAlgError, ValueError):
			print('Singular Matrix')
			pass
	return end_strain, number_points
	

def one_network_test(complexity_network,dimension, length_domain, min_distance, Ef, A, B, creation, tensile_test_1, path):
	network = Network(dimension, complexity_network, length_domain, min_distance, Ef, A, B, creation, path)
	network = network.create_network(creation)
	network = network.set_fibers(creation, path)
	network = tensile_test_1.full_test(network,path)
	return network

def plot_impact_number_points(number_networks,dimension, length_domain, min_distance, Ef, A, B, creation, tensile_test_1, path):
	end_strain = []
	number_points = []
	for complexity_network in [10,30,50,75,90,100, 110,120,130,140,150]:
		#length_domain = complexity_network/10.
		#tensile_test_1.traction_distance = 0.3*length_domain
		print('Complexity:', complexity_network)
		end_strain, number_points = standard_dev_complexity(complexity_network,number_networks,dimension, length_domain, min_distance, Ef, A, B, creation, tensile_test_1, path,end_strain, number_points)
		with open(os.path.join(path,'strain_points.csv'), 'w') as writeFile:
				writer = csv.writer(writeFile)
				writer.writerows([number_points])
				writer.writerows([end_strain])
	return end_strain, number_points


