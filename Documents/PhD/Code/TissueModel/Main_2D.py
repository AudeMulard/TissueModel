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
dimension=3 #dimension of the problem
complexity_network=200 #number of random seed points
length_domain=(1.0,1.0,1.0)
min_distance = 0.0001*length_domain[0]
space_discretization = 0.1*length_domain[0]
beam_Young = 33000
beam_poisson = 0.3
beam_profile = 0.01
disturbance=0.02
traction_distance = 1.0*length_domain[0]
hyperstatic_param = 3
element_size = 0.01
connector_coeff = 1.

### SET DATA FILE

data_path = '../Data_1/default/'

today = date.today()

new_dir = data_path+today.strftime("%b-%d-%Y")+'_'+'%04d' % len(os.listdir(data_path))
os.mkdir(new_dir)
path = new_dir



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





network = Network(dimension, complexity_network, length_domain, min_distance, beam_Young, beam_poisson, beam_profile, connector_coeff, disturbance, hyperstatic_param, creation, generation, path)
network = network.set_fibers(path)
test_1 = Tensile_test(constitutive, side, space_discretization, traction_distance,element_size, plot, video, path,details)
test_1.save_parameters(network,path)
#os.system("abaqus cae script=new_solver_2.py")

"""
### TENSILE TEST WITH DIFFERENT BC
test_1 = Tensile_test(constitutive, side, space_discretization, traction_distance, element_size, plot, video, path,details)
test_1.save_parameters(network,path)

os.system("abaqus cae script=new_solver_1.py")
#os.chdir(path)

test_number_1 = fnmatch.filter(sorted_ls('.'), 'network_vertices_01_00_*.csv')[0][-13:-4]
strain_1, stress_1=stress_strain_curve(test_number_1,network)
plt.plot(strain_1,stress_1,label='u2=0')
test_number_2 = fnmatch.filter(sorted_ls('.'), 'network_vertices_01_00_*.csv')[1][-13:-4]
strain_2, stress_2=stress_strain_curve(test_number_2,network)
plt.plot(strain_2,stress_2,label='free u2')

plt.show()

### RESULTS OF TEST

from Plotting.information_network import load_network_info, sorted_ls, stress_strain_curve
#os.chdir(path)
filenames = fnmatch.filter(os.listdir('.'), 'stress_data_*.csv')

test_number = int(filenames[0][-13:-4])
strain,stress = stress_strain_curve(test_number,network)
#strain = [i/test_1.iterations * test_1.traction_distance for i in range(test_1.iterations+1)]
print(stress,strain)
fig_stress_strain,axss = plt.subplots()
axssd =axss.twinx()
axss.plot(strain, stress)
plot_second_der(axssd,strain,stress,color='c')
plt.show()





cell1 = Cell((0.5,0.5),0.1)
cell2 = Cell((1.5,0.5),0.1)
network = cell1.add_cell(network)
network = cell2.add_cell(network)
plot_geometry(network)
plt.show()

#cell = Cell(0.5,0.5,0.1)


test_1 = Tensile_test(constitutive, side, space_discretization, traction_distance, plot, video, path,details)
network = test_1.full_test(network, path,details)
"""

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
