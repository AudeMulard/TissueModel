import numpy as np
from force_balance import *
from creation_network import Network


# Calulate the force on boundary point

def calulate_macro_stress(network, constitutive):
	stress = 0
	for node in network.boundary_nodes_right:
		for j in network.list_nodes_ridges[node]:
			if node==network.ridge_vertices[j][0]:
				stress += write_force(network, node, network.ridge_vertices[j][1], constitutive)
			if node==network.ridge_vertices[j][1]:
				stress += write_force(network, node, network.ridge_vertices[j][0], constitutive)
	for node in network.boundary_nodes_left:
		for j in network.list_nodes_ridges[node]:
			if node==network.ridge_vertices[j][0]:
				stress -= write_force(network, node, network.ridge_vertices[j][1], constitutive)
			if node==network.ridge_vertices[j][1]:
				stress -= write_force(network, node, network.ridge_vertices[j][0], constitutive)
	stress = 1/network.length**2*np.linalg.norm(stress)
	return stress

# sum them up and divide by volume to give constraint.

def full_test(network, defo, constitutive, scheme, side, iteration, plot, video):
	stress = [0.0]
	strain = [0.0]
	for i in range(iteration):
		print 'Step', i
		network = new_bc(network, defo, side)
		if i==0 and constitutive != 'constant':
			network = linear_scheme(network)
		network=solve_force_balance(network, defo, constitutive, scheme, side)
		if plot == True:
			stress.append(calulate_macro_stress(network, constitutive))
			strain.append(defo*(i+1))
		if video == True:
			network.plot_network_extension()
			plt.savefig("step%d.png" % i)
	if plot == True:
		fig = plt.figure()
		ax = fig.gca()
		ax.scatter(strain, stress)
		plt.savefig('stress_strain.png')
	return network
