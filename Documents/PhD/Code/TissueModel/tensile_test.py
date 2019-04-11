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
		print stress
	for node in network.boundary_nodes_left:
		for j in network.list_nodes_ridges[node]:
			if node==network.ridge_vertices[j][0]:
				stress -= write_force(network, node, network.ridge_vertices[j][1], constitutive)
			if node==network.ridge_vertices[j][1]:
				stress -= write_force(network, node, network.ridge_vertices[j][0], constitutive)
		print stress
	stress = 1/network.length**2*np.linalg.norm(stress)
	return stress

# sum them up and divide by volume to give constraint.

def full_test(network, defo, constitutive, scheme, side, iteration, plot):
	stress = [0.0]
	strain = [0.0]
	for i in range(iteration):
		network=solve_force_balance(network, defo, constitutive, scheme, side)
		if plot == True:
			stress.append(calulate_macro_stress(network, constitutive))
			strain.append(defo*(i+1))
	if plot == True:
		plt.scatter(strain, stress)
	return network
