import numpy as np
from force_balance import *
from creation_network import Network


class Tensile_test:
	def __init__(self, constitutive, scheme, side, space_discretization, traction_distance, plot, video, phase, ax):
		self.constitutive = constitutive
		self.scheme = scheme
		self.side = side
		self.space_discretization = space_discretization
		self.traction_distance = traction_distance
		self.iterations = int(traction_distance / space_discretization)
		self.plot = plot
		self.video = video
		self.phase = phase
		self.ax = ax

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
		stress = 1/network.length**2*np.linalg.norm(stress)
		return stress

# sum them up and divide by volume to give constraint.

	def full_test(self, network):
		network.stress = [0.0]
		network.strain = [0.0]
		for i in range(self.iterations):
	#		print 'Step', i
			network = new_bc(network, self.space_discretization, self.side)
			if i==0 and self.constitutive != 'constant':
				network = linear_scheme(network)
			network=solve_force_balance(network, self.space_discretization, self.constitutive, self.scheme, self.side, i)
			if self.plot == True:
				network.stress.append(self.calculate_macro_stress(network))
				network.strain.append(self.space_discretization*(i+1))
			if self.video == True:
				network.plot_network_extension()
				plt.savefig("%sstep%d.png" % (self.phase,i))
		if self.plot == True:
			fig = plt.figure()
			ax = fig.gca()
			self.ax.scatter(network.strain, network.stress)
			plt.savefig('stress_strain.png')
		return network
