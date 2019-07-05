from creation_network import Network
from force_balance import *
from tensile_test import *
from scipy.spatial import Voronoi, voronoi_plot_2d
import matplotlib.pyplot as plt


## PARAMETERS
dimension=2 #dimension of the problem
complexity_network=50 #number of random seed points
length_domain=1.0
min_distance = length_domain * 0.025
defo = 0.025 #*length_domain
Ef=1.0
A=1.4E-8
B=3.8
iteration = 5


## EXPERIMENT
creation="Voronoi"
constitutive = 'linear2'
scheme='nonlinear'
side = 'right'
plot = True
video = False

fig = plt.figure()
ax = fig.gca()

x = Network(dimension, complexity_network, length_domain, min_distance, Ef, A, B, creation)
x = x.create_network(creation)
x = x.set_fibers(creation)



for value in [0.5,1.0,2.0,10.0]:
	print 'Ef',value
	x.Ef = value
	x = full_test(x, defo, constitutive, scheme, side, iteration, plot,video,'traction',ax)

import numpy
xdata = numpy.array([0, 0.1,0.2,0.3,0.4])
ydata = numpy.array([0,0.146*10E-3,1.135*10E-3,2.987*10E-3,4.722*10E-3])
ax.scatter(xdata,ydata, color = 'red')


plt.savefig('50points_3d.png')
plt.show()
