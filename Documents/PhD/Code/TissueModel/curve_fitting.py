import numpy
import scipy.optimize as optimization
from creation_network import Network
from force_balance import *
from tensile_test import *

xdata = numpy.array([0, 0.1,0.2,0.3,0.4])
ydata = numpy.array([0,0.146,1.135,2.987,4.722])
sigma = numpy.array([1.0,1.0,1.0,1.0,1.0,1.0])


## Collagen parameters
Ef=1.0
A=1.4E-8
B=3.8

x0 = numpy.array([1.0, 1.4E-8,3.8])

## PARAMETERS
dimension=3 #dimension of the problem
complexity_network=30 #number of random seed points
length_domain=1.0
min_distance = length_domain * 0.025
defo = 0.025*length_domain
iteration = 5

## EXPERIMENT
creation="Voronoi"
constitutive = 'linear2'
scheme='nonlinear'
side = 'right'
plot = True
video = False

network = Network(dimension, complexity_network, length_domain, min_distance, Ef, A, B, creation)

network = network.create_network(creation)
network = network.set_fibers(creation)
print 'Network created'

fig = plt.figure()
ax = fig.gca()
#network = full_test(network, defo, constitutive, scheme, side, iteration, plot,video,'traction',ax)
def find_values(network,Ef,A,B):
	network.Ef = Ef
	network.A = A
	network.B = B
	network = full_test(network, defo, constitutive, scheme, side, iteration, plot,video,'traction',ax)
	return network.stress, network.strain

def func(x,Ef,A,B):
	network.stress, network.strain = find_values(network, Ef,A,B)
	if x==0.0:
		return network.stress[0]
	else:
		return network.stress[4]

"""
	print network.strain
	print network.strain.index(x)
	print network.stress[network.strain.index(x)]
	return network.stress[network.strain.index(x)]

for value in xdata:
	print network.strain, network.stress, func(value,1.0, 1.4E-8,3.8)
"""
print optimization.curve_fit(func, xdata, ydata, x0, sigma)

