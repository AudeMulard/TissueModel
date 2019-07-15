import numpy as np
from force_balance import *
from creation_network import Network
import csv, os


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
					stress += write_force(network, node, network.ridge_vertices[j][1], self.constitutive)
				if node==network.ridge_vertices[j][1]:
					stress += write_force(network, node, network.ridge_vertices[j][0], self.constitutive)
		for node in network.boundary_nodes_left:
			for j in network.list_nodes_ridges[node]:
				if node==network.ridge_vertices[j][0]:
					stress -= write_force(network, node, network.ridge_vertices[j][1], self.constitutive)
				if node==network.ridge_vertices[j][1]:
					stress -= write_force(network, node, network.ridge_vertices[j][0], self.constitutive)
		stress = np.linalg.norm(stress)/(network.length/network.mean_length)
		return stress

# sum them up and divide by volume to give constraint.

	def full_test(self, network, path):
		for i in range(self.iterations):
	#		print 'Step', i
			network = new_bc(network, self.space_discretization*self.traction_distance/abs(self.traction_distance), self.side)
			if i==0 and self.constitutive != 'constant':
				network = linear_scheme(network)
			network=solve_force_balance(network, self.space_discretization, self.constitutive, self.scheme, self.side, i)
			if self.plot == True:
				network.stress.append(self.calculate_macro_stress(network))
				network.strain.append(network.strain[-1]+self.space_discretization*self.traction_distance/abs(self.traction_distance))
				last_network = 'stress_strain_%d.csv' % network.complexity
				with open(os.path.join(path,last_network), 'w') as writeFile:
					writer = csv.writer(writeFile)
					writer.writerows([network.strain])
					writer.writerows([network.stress])
				writeFile.close()	
			if self.video == True:
				network.save_network(i, path)
		return network
