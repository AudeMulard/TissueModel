import numpy as np
from Core_calculation.force_balance import *
from Network_generation.creation_network import Network
import csv, os
from Plotting.network_plotting import *


class Tensile_test:
	def __init__(self, constitutive, side, space_discretization, traction_distance, _plot, video, path,details):
		self.constitutive = constitutive
		self.side = side
		self.space_discretization = space_discretization
		self.traction_distance = traction_distance
		self.iterations = abs(int(traction_distance / space_discretization))
		self._plot = _plot
		self.video = video
		self.details = details

	# Adds the parameters of the test to the parameters file

	def save_parameters(self,network,path):
		filename = 'parameters_%03d.csv' % len(network.vertices)
		with open(os.path.join(path,filename), 'a') as writeFile:
			writer = csv.writer(writeFile)
			writer.writerow(["constitutive",self.constitutive])
			writer.writerow(["side",self.side])
			writer.writerow(["space discretization",self.space_discretization])
			writer.writerow(["traction_distance",self.traction_distance])
		writeFile.close()

	# Calulate the macroscopic stress of the network by adding the stress of all the boundary fibers.

	def calculate_macro_stress(self,network):
		stress = 0
		for node in network.boundary_nodes_right:
			for j in network.list_nodes_ridges[node]:
				if node==network.ridge_vertices[j][0]:
					i = node
					k = network.ridge_vertices[j][1]
					stress += write_force(network, node, network.ridge_vertices[j][1], self.constitutive)[0]*network.vertices[node][0]
				if node==network.ridge_vertices[j][1]:
					i = node
					k = network.ridge_vertices[j][0]
					stress += write_force(network, node, network.ridge_vertices[j][0], self.constitutive)[0]*network.vertices[node][0]
		for node in network.boundary_nodes_left:
			for j in network.list_nodes_ridges[node]:
				if node==network.ridge_vertices[j][0]:
					i = node
					k = network.ridge_vertices[j][1]
					stress+=write_force(network, node, network.ridge_vertices[j][1], self.constitutive)[0]*network.vertices[node][0]
				if node==network.ridge_vertices[j][1]:
					i = node
					k = network.ridge_vertices[j][0]
					stress+=write_force(network, node, network.ridge_vertices[j][0], self.constitutive)[0]*network.vertices[node][0]
		stress = stress/(network.length[0]) # divide the sum by length to have stress and not force
		return stress

	def one_step(self, network, path, details,i,length_ini):
		current_disp = (max(network.vertices[:,0])-network.length[0])
		print('Displacement: ', current_disp, i)
		result = False
		tries = 0
		space_discretization = self.space_discretization
		if network.dimension == 2: #separates fibers depending on if they are in tension or compression state to be able to assign the right constant k
			index=0
			for ridge in network.ridge_vertices:
				length_new = length_square(network,network.vertices[ridge[0]]-network.vertices[ridge[1]])
				if length_ini[network.ridge_vertices.index(ridge)] <= length_new:
					network.state_ridge[index]='tension'
				else:
					network.state_ridge[index]='compression'
				index+=1
		while result == False and tries <=10: # solve the balance, if it does not work with initial discretization, reduce the discretization until it works. (more than 10 divides probably means there is an issue with the network)
			#try:
			network = new_bc(network, space_discretization*self.traction_distance/abs(self.traction_distance), self.side) #change boundary positions
			network= solve_force_balance(network, self.constitutive, details) #solve the equilibrium of new configuration
			result = True
			network.save_network('temp',path) # save a configuration in case the next discretization fails, so that the code can begin from the beginning of the step
			"""except (ValueError):#,RuntimeWarning):
				print('New discretization: ', space_discretization/2, ' on step ', current_disp)
				space_discretization = space_discretization/2
				with open(os.path.join(path,'network_vertices.csv'),'r') as readFile:
					reader = csv.reader(readFile)
					list_vertices = np.array(list(reader))
					network.vertices=list_vertices.astype(float)
				tries +=1
				continue"""
		space_discretization = self.space_discretization
		network.stress.append(self.calculate_macro_stress(network))
		network.strain.append((max(network.vertices[:,0])-network.length[0])/network.length[0])
		return network

	def full_test(self, network, path,details,**kw):
		self.save_parameters(network,path)
		network.save_network('temp',path) # save a first temporary version of the network
		i = 0
		network.strain = []
		network.stress = []
		length_ini = []
		for ridge in network.ridge_vertices:
			length_ini.append(length_square(network,network.vertices_ini[ridge[0]]-network.vertices_ini[ridge[1]]))
		if kw.get('name')!=None and kw['name'][-5:] == 'compr':
			while (max(network.vertices[:,0])-network.length[0]) >= 0.0:
				network = self.one_step(network, path, details,i,length_ini)
				if kw.get('name'):
					last_network = 'stress_strain_%03d_%s.csv' % (len(network.vertices),kw['name'])
				else:
					last_network = 'stress_strain_%03d.csv' % len(network.vertices)
				with open(os.path.join(path,last_network), 'w') as writeFile:
					writer = csv.writer(writeFile)
					writer.writerows([network.strain])
					writer.writerows([network.stress])
				writeFile.close()
				if self._plot == True:
					if kw.get('name')!=None:
						network.save_network(i, path, name = kw['name'])
					else:
						network.save_network(i, path)
					i+=1
		else:
			while (max(network.vertices[:,0])-network.length[0]) <= self.traction_distance:
				network = self.one_step(network, path, details,i,length_ini)
				if kw.get('name')!=None:
					last_network = 'stress_strain_%03d_%s.csv' % (len(network.vertices),kw['name'])
				else:
					last_network = 'stress_strain_%03d.csv' % len(network.vertices)
				with open(os.path.join(path,last_network), 'w') as writeFile:
					writer = csv.writer(writeFile)
					writer.writerows([network.strain])
					writer.writerows([network.stress])
				writeFile.close()
				if self._plot == True:
					if kw.get('name')!=None:
						network.save_network(i, path, name = kw['name'])
					else:
						network.save_network(i, path)
					i+=1
		return network
