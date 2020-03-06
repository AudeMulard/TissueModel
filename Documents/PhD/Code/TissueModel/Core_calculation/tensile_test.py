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
					#stress += network.Ef*(np.sqrt(length_square(network,network.vertices[i]-network.vertices[k]))-np.sqrt(length_square(network,network.vertices_ini[i]-network.vertices_ini[k]))) *network.vertices[node][0]
					stress += write_force(network, node, network.ridge_vertices[j][1], self.constitutive)[0]*network.vertices[node][0]
					#stress += network.Ef*(np.sqrt(length_square(network,network.vertices[i]-network.vertices[k]))-np.sqrt(length_square(network,network.vertices_ini[i]-network.vertices_ini[k])))* (network.vertices[i][0]-network.vertices[k][0]) /np.sqrt(length_square(network,network.vertices[i]-network.vertices[k]))*network.vertices[node][0]
				if node==network.ridge_vertices[j][1]:
					i = node
					k = network.ridge_vertices[j][0]
					#stress += network.Ef*(np.sqrt(length_square(network,network.vertices[i]-network.vertices[k]))-np.sqrt(length_square(network,network.vertices_ini[i]-network.vertices_ini[k]))) *network.vertices[node][0]#
					stress += write_force(network, node, network.ridge_vertices[j][0], self.constitutive)[0]*network.vertices[node][0]
					#stress += network.Ef*(np.sqrt(length_square(network,network.vertices[i]-network.vertices[k]))-np.sqrt(length_square(network,network.vertices_ini[i]-network.vertices_ini[k])))* (network.vertices[i][0]-network.vertices[k][0]) /np.sqrt(length_square(network,network.vertices[i]-network.vertices[k]))*network.vertices[node][0]
		for node in network.boundary_nodes_left:
			for j in network.list_nodes_ridges[node]:
				if node==network.ridge_vertices[j][0]:
					i = node
					k = network.ridge_vertices[j][1]
					#stress -= network.Ef*(np.sqrt(length_square(network,network.vertices[i]-network.vertices[k]))-np.sqrt(length_square(network,network.vertices_ini[i]-network.vertices_ini[k]))) *network.vertices[node][0]
					stress+=write_force(network, node, network.ridge_vertices[j][1], self.constitutive)[0]*network.vertices[node][0]
					#stress += network.Ef*(np.sqrt(length_square(network,network.vertices[i]-network.vertices[k]))-np.sqrt(length_square(network,network.vertices_ini[i]-network.vertices_ini[k])))* (network.vertices[i][0]-network.vertices[k][0]) /np.sqrt(length_square(network,network.vertices[i]-network.vertices[k]))*network.vertices[node][0]
				if node==network.ridge_vertices[j][1]:
					i = node
					k = network.ridge_vertices[j][0]
					#stress -= network.Ef*(np.sqrt(length_square(network,network.vertices[i]-network.vertices[k]))-np.sqrt(length_square(network,network.vertices_ini[i]-network.vertices_ini[k]))) *network.vertices[node][0]#
					stress+=write_force(network, node, network.ridge_vertices[j][0], self.constitutive)[0]*network.vertices[node][0]
					#stress += network.Ef*(np.sqrt(length_square(network,network.vertices[i]-network.vertices[k]))-np.sqrt(length_square(network,network.vertices_ini[i]-network.vertices_ini[k])))* (network.vertices[i][0]-network.vertices[k][0]) /np.sqrt(length_square(network,network.vertices[i]-network.vertices[k]))*network.vertices[node][0]
		stress = stress/(network.length)
		return stress

# sum them up and divide by volume to give constraint.

	def full_test(self, network, path,details):
		network.save_network('temp',path)
		i = 0
		while (max(network.vertices[:,0])-network.length) <= self.traction_distance:
			current_disp = (max(network.vertices[:,0])-network.length)
			print 'Displacement: ', current_disp, i
		#for ite in range(1,self.iterations):
	#		print 'Step', i"""
			result = False
			tries = 0
			#network = new_bc(network, self.space_discretization*self.traction_distance/abs(self.traction_distance), self.side)
			#network=solve_force_balance(network, self.space_discretization, self.constitutive, self.scheme, self.side, ite,details)
			space_discretization = self.space_discretization
			while result == False and tries <=10:
				try:
					network = new_bc(network, space_discretization*self.traction_distance/abs(self.traction_distance), self.side)
					network=solve_force_balance(network, space_discretization, self.constitutive, self.scheme, self.side, current_disp,details)
					result = True
					network.save_network('temp',path)
				except (ValueError,RuntimeWarning):
					print 'New discretization: ', space_discretization/2, ' on step ', current_disp
					space_discretization = space_discretization/2
					with open('network_vertices.csv','r') as readFile:
						reader = csv.reader(readFile)
						list_vertices = np.array(list(reader))
						network.vertices=list_vertices.astype(float)
					tries +=1
					continue
			space_discretization = self.space_discretization
			network.stress.append(self.calculate_macro_stress(network))
			network.strain.append((max(network.vertices[:,0])-network.length)/network.length)
			#last_network = 'stress_strain_%03d.csv' % network.complexity
			last_network = 'stress_strain_%s.csv' % network.creation
			with open(os.path.join(path,last_network), 'w') as writeFile:
				writer = csv.writer(writeFile)
				writer.writerows([network.strain])
				writer.writerows([network.stress])
			writeFile.close()
			if self.plot == True:
				if self.phase!='only_one':
					if self.traction_distance>0:
						j = (2*self.phase)*self.iterations+current_disp
					elif self.traction_distance<0:
						j = (2*self.phase+1)*self.iterations+current_disp
					network.save_network(j, path)
				else:
					network.save_network(i, path)
					i+=1
		return network
