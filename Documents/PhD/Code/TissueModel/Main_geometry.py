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
complexity_network=100 #number of random seed points
length_domain=1.0
min_distance = 0.0001*length_domain
space_discretization = 0.01*length_domain
Ef=1.0
A=1.
disturbance=0.02
traction_distance = 0.1*length_domain
#iteration = 15



data_path = '../Data/disturbed_grid/'

today = date.today()

new_dir = data_path+today.strftime("%b-%d-%Y")+'_'+'%04d' % len(os.listdir(data_path))
os.mkdir(new_dir)
path = new_dir


## EXPERIMENT
creation="disturbed_grid"
constitutive = 'linear2'
scheme='nonlinear'
side = 'right'
plot = True
video = True
phase = 'only_one'
stress_rep = True
details = True


network = Network(dimension, complexity_network, length_domain, min_distance, Ef, A, disturbance, creation, path)

network = network.set_fibers(creation, path)
"""
os.chdir('../Data/disturbed_grid/')
if len(sys.argv) != 1:
	os.chdir(sys.argv[1])
else:
	os.chdir(sorted_ls('.')[-1])
	print(os.getcwd())
print sorted_ls('.')
network=load_network_info('.','initial')
"""
plot_geometry(network)
plt.show()

for node in network.interior_nodes:
	cond = True
	list_ridge=network.list_nodes_ridges[node]
	if network.vertices[node][0]<=float(network.length)/10.:
		for ridge_number in list_ridge:
			ridge = network.ridge_vertices[ridge_number]
			if network.vertices[ridge[0]][0]< network.vertices[node][0] or network.vertices[ridge[1]][0]<network.vertices[node][0]:
				cond = False
				print cond
		if cond == True:
			print True
			network.vertices=np.append(network.vertices,[[0.,network.vertices[node][1]]],axis=0)
			network.ridge_vertices= np.append(network.ridge_vertices,[[len(network.vertices)-1,node]],axis=0)

plot_geometry(network)
plt.show()

test_1 = Tensile_test(constitutive, scheme, side, space_discretization, traction_distance, plot, video, 0, path)


network = test_1.full_test(network, path,details)


