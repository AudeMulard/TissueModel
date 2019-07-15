from creation_network import Network
from force_balance import *
from tensile_test import *
from scipy.spatial import Voronoi, voronoi_plot_2d
import matplotlib.pyplot as plt
from Several_networks_test import *

## PARAMETERS
dimension=2 #dimension of the problem
complexity_network=400 #number of random seed points
length_domain=10.0
min_distance = 0.1#*length_domain
space_discretization = 0.05*length_domain
Ef=1.0
A=1.4E-8
B=3.8
traction_distance = 0.3*length_domain
#iteration = 15

path = '../Data/poster/network_size'
#np.set_printoptions(precision=2)

## EXPERIMENT
creation="Voronoi"
constitutive = 'linear2'
scheme='nonlinear'
side = 'right'
plot = True
video = False
phase = 'traction'


## The code stops when there is no numerical value
import sys
"""
if not sys.warnoptions:
    import os, warnings
    warnings.simplefilter("error") # Change the filter in this process
    os.environ["PYTHONWARNINGS"] = "error"
"""
## One tensile_test

"""

x = Network(dimension, complexity_network, length_domain, min_distance, Ef, A, B, creation, path)

x = x.set_fibers(creation, path)
print x.mean_length



x.plot_network()

print 'plotted'
plt.savefig("initial.pdf")

test_1 = Tensile_test(constitutive, scheme, side, space_discretization, traction_distance, plot, video, phase, path)

x = test_1.full_test(x, path)


"""

# Comparison to collagen

"""
import numpy
xdata = numpy.array([0, 0.1,0.2,0.3,0.4])
ydata = numpy.array([0,0.146*10E-3,1.135*10E-3,2.987*10E-3,4.722*10E-3])
ax.scatter(xdata,ydata, color = 'red')
"""

## Impact of number of points on standard deviation

fig = plt.figure()
ax = fig.gca()
test_1 = Tensile_test(constitutive, scheme, side, space_discretization, traction_distance, plot, video, phase, ax)
print plot_impact_number_points(100,dimension, length_domain, min_distance, Ef, A, B, creation, test_1, path)


## Creep test

"""
x = Network(dimension, complexity_network, length_domain, min_distance, Ef, A, B, creation)

x = x.set_fibers(creation)


x.plot_network()
print 'Network created'
plt.savefig("initial.pdf")


fig = plt.figure()
ax = fig.gca()
test_traction_1 = Tensile_test(constitutive, scheme, side, space_discretization, 0.3*length_domain, plot, video, 'traction', ax)

x = test_traction_1.full_test(x)
print 'First traction', x.strain, x.stress

x.plot_network_extension()
plt.savefig("First_traction.pdf")

test_compression_1 = Tensile_test(constitutive, scheme, side, space_discretization, -0.2*length_domain, plot, video, 'compression', ax)

x = test_compression_1.full_test(x)
print 'First compression', x.strain, x.stress

x.plot_network_extension()
plt.savefig("First_compression.pdf")

test_traction_2 = Tensile_test(constitutive, scheme, side, space_discretization, 0.7*length_domain, plot, video, 'traction', ax)

x = test_traction_2.full_test(x)
print 'Second traction', x.strain, x.stress

x.plot_network_extension()
plt.savefig("Second_traction.pdf")

test_compression_2 = Tensile_test(constitutive, scheme, side, space_discretization, -0.8*length_domain, plot, video, 'compression', ax)

x = test_compression_2.full_test(x)
print 'Second compression', x.strain, x.stress

x.plot_network_extension()
plt.savefig("Second_compression.pdf")


fig = plt.figure()
ax = fig.gca()
ax.scatter(x.strain[0:11], x.stress[0:11], color = 'red', label='First traction')
ax.scatter(x.strain[12:22], x.stress[12:22], color='green', label='First compression')
ax.scatter(x.strain[23:34], x.stress[23:34], color='blue', label='Second traction')
ax.scatter(x.strain[35:46], x.stress[35:46], color='black', label = 'Second compression')
plt.legend()
plt.savefig('stress_strain.pdf')

fig = plt.figure()
ax = fig.gca()
plt.ylabel('Extension', fontsize=40)
#plt.xlabel('Time', fontsize=40)
ax.plot([0.,1.,2.,3.,4.],[0.0, 0.3*length_domain, 0.3*length_domain-0.2*length_domain,0.3*length_domain-0.2*length_domain+0.7*length_domain,0.3*length_domain-0.2*length_domain+0.7*length_domain-0.8*length_domain], linestyle='solid')
plt.savefig('cycle.png')
"""

