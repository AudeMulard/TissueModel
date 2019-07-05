from creation_network import Network
from force_balance import *
from tensile_test import *
from scipy.spatial import Voronoi, voronoi_plot_2d
import matplotlib.pyplot as plt
import numpy


def standard_dev_complexity(complexity_network, number_networks,dimension, length_domain, min_distance, Ef, A, B, creation, tensile_test_1):
	end_strain = []
	for i in range(number_networks):
		print 'Network', i
		try:
			network = one_network_test(complexity_network,dimension, length_domain, min_distance, Ef, A, B, creation,tensile_test_1)
			end_strain.append(network.stress[tensile_test_1.iterations])
			network.plot_network_extension()
			plt.savefig("%d%d.png" % (i,complexity_network))
			plt.close()
		except (np.linalg.linalg.LinAlgError):
			print('Singular Matrix')
			pass
	return numpy.mean(end_strain), numpy.std(end_strain)
	

def one_network_test(complexity_network,dimension, length_domain, min_distance, Ef, A, B, creation, tensile_test_1):
	network = Network(dimension, complexity_network, length_domain, min_distance, Ef, A, B, creation)
	network = network.create_network(creation)
	network = network.set_fibers(creation)
	network = tensile_test_1.full_test(network)
	return network

def plot_impact_number_points(number_networks,dimension, length_domain, min_distance, Ef, A, B, creation, tensile_test_1):
	mean = []
	standard_dev = []
	for complexity_network in [10,20,30,40,50,60, 70, 80, 90, 100, 125,150,300,500]:
		print 'Complexity:', complexity_network
		mean1, standard_dev1=standard_dev_complexity(complexity_network,number_networks,dimension, length_domain, min_distance, Ef, A, B, creation, tensile_test_1)
		mean.append(mean1)
		standard_dev.append(standard_dev1)
		print mean, standard_dev
	return mean, standard_dev

"""
import numpy
xdata = numpy.array([0, 0.1,0.2,0.3,0.4])
ydata = numpy.array([0,0.146*10E-3,1.135*10E-3,2.987*10E-3,4.722*10E-3])
ax.scatter(xdata,ydata, color = 'red')
"""


