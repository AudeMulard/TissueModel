from creation_network import Network
from force_balance import *
from tensile_test import *
from scipy.spatial import Voronoi, voronoi_plot_2d
import matplotlib.pyplot as plt

## PARAMETERS
dim=2 #dimension of the problem
complexity_network=50 #number of random seed points
length_domain=1.0
min_distance = length_domain * 0.025
defo = 0.01*length_domain
Ef=1.
A=1.4E-8
B=3.8
iteration = 1
plot = False

np.set_printoptions(precision=2)

## EXPERIMENT
creation="Voronoi"
constitutive = 'exponential'
scheme='nonlinear'
side = 'right'


x = Network(dim, complexity_network, length_domain, min_distance, Ef, A, B, creation)

x= x.set_fibers(creation)

x = full_test(x, defo, constitutive, scheme, side, iteration, plot)

x.plot_network_extension()
plt.show()

