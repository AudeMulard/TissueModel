from creation_network import Network
from force_balance import *
from tensile_test import *
from scipy.spatial import Voronoi, voronoi_plot_2d
import matplotlib.pyplot as plt

## PARAMETERS
dimension=3 #dimension of the problem
complexity_network=20 #number of random seed points
length_domain=1.0
min_distance = length_domain * 0.025
defo = 0.1*length_domain
Ef=1.0
A=1.4E-8
B=3.8
iteration = 15


#np.set_printoptions(precision=2)

## EXPERIMENT
creation="Voronoi"
constitutive = 'linear2'
scheme='nonlinear'
side = 'right'
plot = True
video = False

x = Network(dimension, complexity_network, length_domain, min_distance, Ef, A, B, creation)

x = x.create_network(creation)
if creation == "Voronoi":
	last_network = 'network.txt'
	myfile = open(last_network, 'w')
	myfile.write(str(x.vertices.tolist()))
	myfile.write('\n')
	myfile.write(str(x.ridge_vertices))
	myfile.close()

x = x.set_fibers(creation)
x = new_bc(x, defo, side)
x.plot_network()
x = solve_force_balance(x, defo, constitutive, scheme, side)

#x = full_test(x, defo, constitutive, scheme, side, iteration, plot,video)


#x.plot_network_extension()
x.plot_network()
plt.show()

