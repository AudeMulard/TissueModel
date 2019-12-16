import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import dual_annealing, minimize, basinhopping, differential_evolution
#np.set_printoptions(precision=0)
import creation_network



### Set boundary conditions on exterior nodes

def new_bc(network, defo, side):
	new_positions = network.vertices
	if side == 'left':
		for k in network.boundary_nodes_left:
		        new_positions[k,0] = new_positions[k,0] - defo
	if side == 'right':
		for k in network.boundary_nodes_right:
		        new_positions[k,0] = new_positions[k,0] + defo
	if side == 'both':
		for k in network.boundary_nodes_left:
		        new_positions[k,0] = new_positions[k,0] - defo
		for k in network.boundary_nodes_right:
		        new_positions[k,0] = new_positions[k,0] + defo
#	for k in network.interior_nodes:
#		new_positions[k] = new_positions[k] * 1.05 # essential movement to not have a singular matrix at first step
	network.vertices=new_positions
	return network
"""
# Def of norm square
def length(network,x):	
	if network.dimension == 2:
		return np.sqrt(x[0]**2+x[1]**2)
	elif network.dimension == 3:
		return x[0]**2+x[1]**2+x[2]**2

### Calculate elastic energy 

# Calculate energy of one ridge
def elas_energy_total(vertices_interior):
	energy =0
	for ridge in network.ridge_vertices:
		if ridge[0] in network.interior_nodes:
			i = network.interior_nodes.index(ridge[0])
			if ridge[1] in network.interior_nodes:
				j = network.interior_nodes.index(ridge[1])
				energy += network.Ef*(length(network,network.vertices_interior[i]-network.vertices_interior[j])-length(network,network.vertices_ini[ridge[0]]-network.vertices_ini[ridge[1]]))**2
			else:
				energy += network.Ef*(length(network,network.vertices_interior[i]-network.vertices[ridge[1]])-length(network,network.vertices_ini[ridge[0]]-network.vertices_ini[ridge[1]]))**2
		else:
			if ridge[1] in network.interior_nodes:
				j = network.interior_nodes.index(ridge[1])
				energy += network.Ef*(length(network,network.vertices[ridge[0]]-network.vertices_interior[j])-length(network,network.vertices_ini[ridge[0]]-network.vertices_ini[ridge[1]]))**2
			else:
				energy += network.Ef*(length(network,network.vertices[ridge[0]]-network.vertices[ridge[1]])-length(network,network.vertices_ini[ridge[0]]-network.vertices_ini[ridge[1]]))**2
	return energy

"""
### Find the minimum of energy

def energy_minimum(network):
	vertices_interior = []
	for node in network.interior_nodes:
		vertices_interior.append(network.vertices[node][0])
		vertices_interior.append(network.vertices[node][1])
	downing = [-0.02] * (2 * len(network.interior_nodes))
	#lw = [min(min(network.vertices[:,0]),min(network.vertices[:,1]))] * 2 * len(network.interior_nodes)
	lw = np.array(vertices_interior)+np.array(downing)
	adding = [0.02] * (2 * len(network.interior_nodes))
	up = np.array(vertices_interior)+np.array(adding)
	#ret = basinhopping(network.elas_energy_total, vertices_interior)
	#ret = dual_annealing(network.elas_energy_total,bounds = list(zip(lw,up)))
	#ret = minimize(network.elas_energy_total, vertices_interior, method = 'L-BFGS-B')
	ret = differential_evolution(network.elas_energy_total,bounds = list(zip(lw,up)))
	#ret = shgo(network.elas_energy_total, bounds = list(zip(lw,up)))
	print ret.x
	i=0
	for node in network.interior_nodes:
		network.vertices[node][0] = ret.x[2*i]
		network.vertices[node][1] = ret.x[2*i+1]
		i+=1
	return network


def iterative_newton(network, constitutive):
	max_iter = 100
	epsilon = 1e-8
	for k in range(max_iter):
		F = write_vector_F(network, constitutive)
		J = write_matrix_J(network, constitutive)
		diff = np.linalg.solve(J,-F)
		for i in range(len(network.interior_nodes)):
			j=network.interior_nodes[i]
			network.vertices[j,0] = network.vertices[j,0] + diff[2*i]
			network.vertices[j,1] = network.vertices[j,1] + diff[2*i+1]
		if np.linalg.norm(diff) < epsilon:
			print('convergence!, nre iter:', k)
			break
		elif k== max_iter-1:
			raise ValueError 
	return network



