from creation_network import Network
from force_balance import *
from tensile_test import *
from scipy.spatial import Voronoi, voronoi_plot_2d
import matplotlib.pyplot as plt
from Several_networks_test import *

## PARAMETERS
dimension=2 #dimension of the problem
complexity_network=300 #number of random seed points
length_domain=1.0
min_distance = 0.005*length_domain
space_discretization = 0.005#*length_domain
Ef=1.4E-8
A=1.4E-8
B=3.8
traction_distance = 0.1*length_domain
#iteration = 15


#np.set_printoptions(precision=2)

## EXPERIMENT
creation="Voronoi"
constitutive = 'linear2'
scheme='nonlinear'
side = 'right'
plot = True
video = True
phase = 'traction'


## One tensile_test

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

x.plot_network()
print 'plotted'
plt.savefig("initial.png")


fig = plt.figure()
ax = fig.gca()
test_1 = Tensile_test(constitutive, scheme, side, space_discretization, traction_distance, plot, video, phase, ax)

x = test_1.full_test(x)

x.plot_network_extension()
#x.plot_network()


plt.savefig("Final.png")
#plt.show()

# Comparison to collagen

"""
import numpy
xdata = numpy.array([0, 0.1,0.2,0.3,0.4])
ydata = numpy.array([0,0.146*10E-3,1.135*10E-3,2.987*10E-3,4.722*10E-3])
ax.scatter(xdata,ydata, color = 'red')
"""

## Impact of number of points on standard deviation

"""
fig = plt.figure()
ax = fig.gca()
tensile_test_1 = Tensile_test(constitutive, scheme, side, space_discretization, traction_distance, plot, video, phase, ax)
print plot_impact_number_points(100,dimension, length_domain, min_distance, Ef, A, B, creation, tensile_test_1)
"""
