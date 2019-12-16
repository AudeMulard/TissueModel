import numpy as np
from force_balance import *
from creation_network import Network
import csv, os
from network_plotting import *


class Tensile_test:
	def __init__(self, constitutive, scheme, side, space_discretization, traction_distance, plot, video, phase, path):
		self.constitutive = constitutive
		self.scheme = scheme
		self.side = side
		self.space_discretization = space_discretization
		self.traction_distance = traction_distance
		self.iterations = abs(int(traction_distance / space_discretization))
		self.plot = plot
		self.video = video
		self.phase = phase

# Calulate the force on boundary point

	def calculate_macro_stress(self,network):
		stress = 0
		for node in network.boundary_nodes_right:
			for j in network.list_nodes_ridges[node]:
				if node==network.ridge_vertices[j][0]:
					i = node
					k = network.ridge_vertices[j][1]
					stress += network.Ef*(np.sqrt(length_square(network,network.vertices[i]-network.vertices[k]))-np.sqrt(length_square(network,network.vertices_ini[i]-network.vertices_ini[k]))) *network.vertices[node][0]#write_force(network, node, network.ridge_vertices[j][1], self.constitutive
				if node==network.ridge_vertices[j][1]:
					i = node
					k = network.ridge_vertices[j][0]
					stress += network.Ef*(np.sqrt(length_square(network,network.vertices[i]-network.vertices[k]))-np.sqrt(length_square(network,network.vertices_ini[i]-network.vertices_ini[k]))) *network.vertices[node][0]#write_force(network, node, network.ridge_vertices[j][0], self.constitutive)[0]*network.vertices[node][0]
		for node in network.boundary_nodes_left:
			for j in network.list_nodes_ridges[node]:
				if node==network.ridge_vertices[j][0]:
					i = node
					k = network.ridge_vertices[j][1]
					stress -= network.Ef*(np.sqrt(length_square(network,network.vertices[i]-network.vertices[k]))-np.sqrt(length_square(network,network.vertices_ini[i]-network.vertices_ini[k]))) *network.vertices[node][0]#write_force(network, node, network.ridge_vertices[j][1], self.constitutive)[0]*network.vertices[node][0]
				if node==network.ridge_vertices[j][1]:
					i = node
					k = network.ridge_vertices[j][0]
					stress -= network.Ef*(np.sqrt(length_square(network,network.vertices[i]-network.vertices[k]))-np.sqrt(length_square(network,network.vertices_ini[i]-network.vertices_ini[k]))) *network.vertices[node][0]#write_force(network, node, network.ridge_vertices[j][0], self.constitutive)[0]*network.vertices[node][0]
		stress = stress/(network.length)
		return stress

# sum them up and divide by volume to give constraint.

	def full_test(self, network, path):
		for i in range(self.iterations):
	#		print 'Step', i
			network = new_bc(network, self.space_discretization*self.traction_distance/abs(self.traction_distance), self.side)
			network=solve_force_balance(network, self.space_discretization, self.constitutive, self.scheme, self.side, i)
			network.stress.append(self.calculate_macro_stress(network))
			#network.strain.append(network.strain[-1]+self.space_discretization*self.traction_distance/abs(self.traction_distance))
			network.strain.append((max(network.vertices[:,0])-network.length)/network.length)
			last_network = 'stress_strain_%d.csv' % len(network.vertices)
			with open(os.path.join(path,last_network), 'w') as writeFile:
				writer = csv.writer(writeFile)
				writer.writerows([network.strain])
				writer.writerows([network.stress])
			writeFile.close()
			if self.plot == True:
				network.save_network(i, path)
		return network
