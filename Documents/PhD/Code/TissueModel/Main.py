from creation_network import Network
from force_balance import *
from scipy.spatial import Voronoi, voronoi_plot_2d
import matplotlib.pyplot as plt


dim=2 #dimension of the problem
complexity_network=50 #number of random seed points
length_domain=1.0
min_distance = length_domain * 0.025
defo = 0.3*length_domain
Ef=1.
A=1.4E-8
B=3.8

x = Network(dim, complexity_network, length_domain, min_distance, Ef, A, B, 'small grid')

x= x.set_fibers()
#print x.vertices
x.plot_network()
plt.savefig('initial_network.png')

np.set_printoptions(precision=2)
#print x.vertices



x = solve_force_balance(x,defo, constitutive = 'exponential', scheme='nonlinear')
x.plot_network()
plt.savefig('final_network.png')

#x.plot_network_extension()
#plt.savefig('extension_network.png')
#plt.plot([0,0],[0,1]); plt.plot([1,1],[0,1]); plt.plot([0,1],[0,0]); plt.plot([0,1],[1,1])
#plt.plot(x.vertices[:,0],x.vertices[:,1])
#plt.legend()



#plt.show()

# singular matrix problem: alone point in interior_nodes
