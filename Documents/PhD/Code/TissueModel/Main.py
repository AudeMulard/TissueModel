from creation_network import Network
from force_balance import *
from macro_curve import *
from scipy.spatial import Voronoi, voronoi_plot_2d
import matplotlib.pyplot as plt

## PARAMETERS
dim=2 #dimension of the problem
complexity_network=50 #number of random seed points
length_domain=1.0
min_distance = length_domain * 0.025
defo = 0.5*length_domain
Ef=1.
A=1.4E-8
B=3.8

## EXPERIMENT
creation="1 straight line"
constitutive = 'exponential'
scheme='nonlinear'
side = 'right'
iteration = 10

x = Network(dim, complexity_network, length_domain, min_distance, Ef, A, B, creation)

x= x.set_fibers(creation)
#print x.vertices
#x.plot_network()
#plt.savefig('initial_network.png')

np.set_printoptions(precision=2)
#print x.vertices



#x = solve_force_balance(x,defo, constitutive = 'constant', scheme='nonlinear')
#x.plot_network()
#plt.savefig('final_network.png')

#plt.savefig('extension_network.png')
#plt.plot([0,0],[0,1]); plt.plot([1,1],[0,1]); plt.plot([0,1],[0,0]); plt.plot([0,1],[1,1])
#plt.plot(x.vertices[:,0],x.vertices[:,1])
#plt.legend()

x = plot_strain_stress_curve(x, defo, constitutive, scheme, side, iteration)

#x.plot_network_extension()
plt.show()

# singular matrix problem: alone point in interior_nodes
