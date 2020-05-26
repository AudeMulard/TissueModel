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
complexity_network=50 #number of random seed points
length_domain=(1.0,1.0,1./3.)
min_distance = 0.0001*length_domain[0]
space_discretization = 0.01*length_domain[0]
k_tension=1.0
k_compression = 1.0
A=1.
disturbance=0.02
traction_distance = 0.1*length_domain[0]
#iteration = 15


data_path = '../Data/default/'

today = date.today()

new_dir = data_path+today.strftime("%b-%d-%Y")+'_'+'%04d' % len(os.listdir(data_path))
os.mkdir(new_dir)
path = new_dir



## EXPERIMENT
creation="Voronoi"
generation = 'regular'
constitutive = 'spring'
side = 'right'
plot = True
video = True
phase = 'only_one'
stress_rep = True
details = True

network = Network(dimension, complexity_network, length_domain, min_distance, k_tension, k_compression, A, disturbance, creation, generation, path)
network = network.set_fibers(path)
print len(network.ridge_vertices)-2*len(network.interior_nodes)
#plot_geometry(network)
#plt.show()
test_1 = Tensile_test(constitutive, side, space_discretization, traction_distance, plot, video, path,details)
network = test_1.full_test(network, path,details)
"""



#print Network.__mro__()
for i in range(20):
	network = Network(dimension, complexity_network, length_domain, min_distance, k_tension, k_compression, A, disturbance, creation, generation, path)
	network = network.set_fibers(creation, path)
	print(len(network.ridge_vertices)-2*len(network.interior_nodes)), len(network.vertices), len(network.boundary_nodes_right), len(network.boundary_nodes_left)
	#plot_geometry(network)
	#plt.show()
	lengths, cos_theta_square,list_angle = [],[],[]
	for ridge in network.ridge_vertices:
		# Minimum length
		length = numpy.sqrt(length_square(network,network.vertices[ridge[0]]-network.vertices[ridge[1]]))
		lengths.append(length)
		# Minimum & max angle
		cos_theta_square.append(numpy.dot(network.vertices[ridge[0]]-network.vertices[ridge[1]], [1.,0.])**2/length**2)
		list_angle.append(numpy.arccos(numpy.sqrt(cos_theta_square[-1])))

	list_angle = sorted(list_angle)
	lengths = sorted(lengths)

	diff_angle, diff_lengths = [],[]
	for i in range(len(list_angle)-1):
		diff_angle.append(list_angle[i+1]-list_angle[i])
		diff_lengths.append(lengths[i+1]-lengths[i])

	print min(diff_angle), min(diff_lengths)

	test_1 = Tensile_test(constitutive, side, space_discretization, traction_distance, plot, video, path,details)
	#network = test_1.full_test(network, path,details)

	F = write_vector_F(network, constitutive)
	#print F
	J = write_matrix_J(network, constitutive)
	#if np.linalg.cond(J) > 10e10:
	#	print('STOP SIMULATION: ERROR')
	print np.linalg.cond(J)




for k_compression in [1.0,0.1, 0.05, 0.01, 0.005, 0.001, 0.0005, 0.0001]:
	print 'ratio fibers: ', k_compression
	network.k_compression = k_compression
	network = test_1.full_test(network, path,details,name=k_compression)
	network.vertices = np.array(network.vertices_ini)

os.chdir(path)

filenames=fnmatch.filter(os.listdir('.'), 'stress_strain_*.csv')


print filenames
fig, ax1 = plt.subplots()

ax1.set_xlabel('strain')
ax1.set_ylabel('stress')

for filename in filenames:
	with open(filename, 'r') as readFile:
		reader = csv.reader(readFile)
		curve = np.array(list(reader))
		stress = [float(i) for i in curve[1]]
		strain=[float(i) for i in curve[0]]
	ax1.plot(strain, stress, marker='o',linestyle='dashed', markersize=5., label=filename[-7:-4]) 
	
ax1.legend()
#ax1.legend(loc='upper left')
plt.savefig('stress_strain.pdf')
plt.show()

os.chdir('../../../TissueModel/')
"""
