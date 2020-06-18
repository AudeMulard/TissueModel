import numpy as np
from scipy.spatial import Voronoi, voronoi_plot_2d
import Network_generation.network_types
#import seaborn as sns
import matplotlib.pyplot as plt
import csv
import os
import Core_calculation.force_balance
from Plotting.network_plotting import *
# Class of the network, with different initial positions and all the corrections needed
import fnmatch
import numba


class Network:
	def __init__(self, dimension, complexity_network, length_domain, min_distance, k_tension, k_compression, A, B, creation, generation, path, **kw):
		self.complexity = complexity_network
		self.length=length_domain
		self.dimension = dimension
		self.k_tension = k_tension
		self.k_compression = k_compression
		self.A = A
		self.disturbance = B
		self.min_distance = min_distance
		self.creation = creation
		self.generation = generation
		self.lengths_ini = []
		if kw.get('name')!=None: self.kw = kw
		else: self.kw = kw

		self.strain = [0.]
		self.stress = [0.]
		
	def create_network(self):
		self.vertices, self.ridge_vertices = Network_generation.network_types.select_network(self)
		self.interior_nodes = []
		self.boundary_nodes_left=[]
		self.boundary_nodes_right=[]
		return self

# Deletion of the first point because it creates intersections in the network
	def delete_first_point(self):
		#self = self.create_ridge_node_list()
		ridges_to_delete = []
		for i in range(len(self.ridge_vertices)):
			if self.ridge_vertices[i][0]==-1 or self.ridge_vertices[i][1]==-1:
				ridges_to_delete.append(i)
		#ridges_to_delete = self.list_nodes_ridges[-1]
		ridges_to_delete=np.array(sorted(ridges_to_delete, reverse=True))
		for k in ridges_to_delete:
			self.ridge_vertices=np.delete(self.ridge_vertices, int(k), axis=0)
		self = self.delete_points_with_two_ridges()
		return self

	def delete_points_with_two_ridges(self):
		self = self.create_ridge_node_list()
		for i in range(len(self.list_nodes_ridges)):
			if len(self.list_nodes_ridges[i])==2:
				self.delete_two_ridge_points(i)
		self = self.delete_doubles()
		self = self.create_ridge_node_list()
		self = self.delete_alone_points()
		return self


	def delete_point(self, i):
		self.vertices = np.delete(self.vertices, i, axis=0)
		for ridge in self.ridge_vertices:
			for i in range(2):
				if i < ridge[i]:
					ridge[i] -= 1
		return self

	def delete_two_ridge_points(self,i):
		list_ridges = self.list_nodes_ridges[i]
		ridges = self.ridge_vertices
		if ridges[list_ridges[0]][0] == i and ridges[list_ridges[1]][0] == i:
			ridges[list_ridges[0]][0] = ridges[list_ridges[1]][1]
			ridges[list_ridges[1]][0] = ridges[list_ridges[0]][1]
		elif ridges[list_ridges[0]][1] == i and ridges[list_ridges[1]][0] == i:
			ridges[list_ridges[0]][1] = ridges[list_ridges[1]][1]
			ridges[list_ridges[1]][0] = ridges[list_ridges[0]][0]
		elif ridges[list_ridges[0]][0] == i and ridges[list_ridges[1]][1] == i:
			ridges[list_ridges[0]][0] = ridges[list_ridges[1]][0]
			ridges[list_ridges[1]][1] = ridges[list_ridges[0]][1]
		elif ridges[list_ridges[0]][1] == i and ridges[list_ridges[1]][1] == i:
			ridges[list_ridges[0]][1] = ridges[list_ridges[1]][0]
			ridges[list_ridges[1]][1] = ridges[list_ridges[0]][0]
		self.ridge_vertices == ridges
		return self

# Sort the nodes in three parts: two boundary lists and the interior nodes list

	def sort_nodes(self):
		nodes = self.vertices
		boundary_nodes_left =[]
		boundary_nodes_right=[]
		interior_nodes=[]
		vertices_interior = []
		for i in range(len(nodes)):
        		if nodes[i,0] <= 0.0:
        		        boundary_nodes_left.append(i)
        		elif nodes[i,0] >= self.length[0]:
        		        boundary_nodes_right.append(i)
        		else:
        		        interior_nodes.append(i)
		self.interior_nodes = interior_nodes
		self.boundary_nodes_left = boundary_nodes_left
		self.boundary_nodes_right = boundary_nodes_right
		return self

# Create lists where at each point of the network, you have all the ridges associated

	def create_ridge_node_list(self):
		nodes=self.vertices
		ridges = self.ridge_vertices
		list_nodes_ridges=[[] for i in range(len(nodes))]
		for i in range(len(ridges)):
			list_nodes_ridges[ridges[i][0]].append(i)
			list_nodes_ridges[ridges[i][1]].append(i)
		self.list_nodes_ridges=list_nodes_ridges
		return self

# Merge nodes close to each other

	def merge_nodes(self):
		nodes = self.vertices
		nodes_to_delete = []
		for i in range(len(nodes)):
			if i in self.interior_nodes:
				for j in range(i):
					distance = np.sqrt((nodes[i,0]-nodes[j,0])**2+(nodes[i,1]-nodes[j,1])**2)
					if distance < self.min_distance:
					    for k in range(len(self.ridge_vertices)):
					        if self.ridge_vertices[k][0]==j:
					            self.ridge_vertices[k][0]=i
					        if self.ridge_vertices[k][1]==j:
					            self.ridge_vertices[k][1]=i
					    if j not in nodes_to_delete:
					        nodes_to_delete = np.append(nodes_to_delete,j)
		nodes_to_delete=np.array(sorted(nodes_to_delete, reverse=True))
		for point in nodes_to_delete:
			nodes=np.delete(nodes,int(point),0)
		self.vertices = nodes
		# Renumber points after deleting some
		for ridge in self.ridge_vertices:
			for i in range(2):
				r=0
				for node in nodes_to_delete:
					if node < ridge[i]:
						r+=1
				ridge[i] = ridge[i] - r
		ridges_to_delete = []
		for i in range(len(self.ridge_vertices)):
			if self.ridge_vertices[i][0]==self.ridge_vertices[i][1]:
				ridges_to_delete = np.append(ridges_to_delete,i)
		ridges_to_delete=np.array(sorted(ridges_to_delete, reverse=True))
		for ridge in ridges_to_delete:
			self.ridge_vertices=np.delete(self.ridge_vertices,int(ridge),0)
		return self

# Cut network to make it belong to the domain area, first top and bottom then on the sides
	def cut_network_updown(self):
		nodes = self.vertices
		nodes_to_delete = []
		ridges_to_delete = []
		if self.dimension == 2:
			interval = [1]
		elif self.dimension == 3:
			interval = [1,2]
		# Select in one list all the points outside of the domain
		for coord in interval:
			for i in range(len(nodes)):
        			if nodes[i,coord] < 0.0 or nodes[i,coord] > self.length[coord]:
        			    nodes_to_delete = np.append(nodes_to_delete,i)
		nodes_to_delete = set(nodes_to_delete)
		# delete ridges that are completely out of the domain
		for i in range(len(self.ridge_vertices)):
			if self.ridge_vertices[i][0] in nodes_to_delete and self.ridge_vertices[i][1] in nodes_to_delete:
				ridges_to_delete = np.append(ridges_to_delete, i)
		ridges_to_delete=np.array(sorted(ridges_to_delete, reverse=True))
		for ridge in ridges_to_delete:
			self.ridge_vertices=np.delete(self.ridge_vertices,ridge, axis=0)
		if self.creation == 'growth_network':
			self = self.create_ridge_node_list()
			for node in nodes_to_delete:
				node = int(node)
				if len(self.list_nodes_ridges[node])==2:
					list_ridges = self.list_nodes_ridges[node]
					ridges = self.ridge_vertices
					if ridges[list_ridges[0]][0] == node and ridges[list_ridges[1]][0] == node:
						ridges[list_ridges[0]][0] = ridges[list_ridges[1]][1]
						ridges[list_ridges[1]][0] = ridges[list_ridges[0]][1]
					elif ridges[list_ridges[0]][1] == node and ridges[list_ridges[1]][0] == node:
						ridges[list_ridges[0]][1] = ridges[list_ridges[1]][1]
						ridges[list_ridges[1]][0] = ridges[list_ridges[0]][0]
					elif ridges[list_ridges[0]][0] == node and ridges[list_ridges[1]][1] == node:
						ridges[list_ridges[0]][0] = ridges[list_ridges[1]][0]
						ridges[list_ridges[1]][1] = ridges[list_ridges[0]][1]
					elif ridges[list_ridges[0]][1] == node and ridges[list_ridges[1]][1] == node:
						ridges[list_ridges[0]][1] = ridges[list_ridges[1]][0]
						ridges[list_ridges[1]][1] = ridges[list_ridges[0]][0]
		# delete ridges attached to the points to delete
		ridges_to_delete=[]
		for node in nodes_to_delete:
			node= int(node)
			for k in range(len(self.ridge_vertices)):
			    if self.ridge_vertices[k][0]==node:
			        ridges_to_delete = np.append(ridges_to_delete,k)
			    elif self.ridge_vertices[k][1]==node:
			        ridges_to_delete = np.append(ridges_to_delete,k)
		ridges_to_delete=np.array(sorted(ridges_to_delete, reverse=True))
		for ridge in ridges_to_delete:
			self.ridge_vertices=np.delete(self.ridge_vertices,ridge, axis=0)
		# delete the exterior points
		nodes_to_delete=np.array(sorted(nodes_to_delete, reverse=True))
		for point in nodes_to_delete:
			nodes=np.delete(nodes,int(point),0)
		# Renumber points in the ridge_vertices after deleting some
		for ridge in self.ridge_vertices:
			for i in range(2):
				r=0
				for node in nodes_to_delete:
					if node < ridge[i]:
						r+=1
				ridge[i] = ridge[i] - r
		self.vertices=nodes
		if self.creation == 'Voronoi' and self.generation != 'regular':
			for k in range(3):
				self = self.create_ridge_node_list()
				for i in range(len(self.list_nodes_ridges)):
					if len(self.list_nodes_ridges[i])==2:
						self.delete_two_ridge_points(i)
				self = self.delete_doubles()
				self = self.create_ridge_node_list()
				self = self.delete_alone_points()
				## Supress alone fibres and points
				while ridges_to_delete.tolist() !=[]:
					ridges_to_delete=[]
					nodes_to_delete=[]
					self = self.create_ridge_node_list()
					self = self.sort_nodes()
					for i in range(len(self.interior_nodes)):
						j = self.interior_nodes[i]
						if len(self.list_nodes_ridges[j])<=1:
							nodes_to_delete = np.append(nodes_to_delete,j)
							ridges_to_delete= np.append(ridges_to_delete,self.list_nodes_ridges[j])
					nodes = self.vertices
					ridges_to_delete=np.array(sorted(ridges_to_delete, reverse=True))
					for ridge in ridges_to_delete:
						self.ridge_vertices=np.delete(self.ridge_vertices,ridge, axis=0)
					nodes_to_delete=np.array(sorted(nodes_to_delete, reverse=True))
					for point in nodes_to_delete:
						nodes=np.delete(nodes,int(point),0)
					# Renumber points after deleting some
					for ridge in self.ridge_vertices:
						for i in range(2):
							r=0
							for node in nodes_to_delete:
								if node < ridge[i]:
									r+=1
							ridge[i] = ridge[i] - r
					self.vertices=nodes
		return self

	def cut_network_side(self):
		nodes = self.vertices
		nodes_to_delete_left = []
		nodes_to_delete_right = []
		ridges_to_delete = []
		for i in range(len(nodes)):
        		if nodes[i,0] < 0.0:
        		    nodes_to_delete_left = np.append(nodes_to_delete_left,i)
        		if nodes[i,0] > self.length[0]:
        		    nodes_to_delete_right = np.append(nodes_to_delete_right,i)
		nodes_to_delete = np.append(nodes_to_delete_left, nodes_to_delete_right)
		for i in range(len(self.ridge_vertices)):
			if self.ridge_vertices[i][0] in nodes_to_delete and self.ridge_vertices[i][1] in nodes_to_delete:
				ridges_to_delete = np.append(ridges_to_delete, i)
		ridges_to_delete=np.array(sorted(ridges_to_delete, reverse=True))
		for ridge in ridges_to_delete:
			self.ridge_vertices=np.delete(self.ridge_vertices,int(ridge), axis=0)
		for node in nodes_to_delete_left:
			node=int(node)
			for k in range(len(self.ridge_vertices)):
			    if self.ridge_vertices[k][0]==node:
			        if self.dimension == 2:
			            y = nodes[node,1]+(nodes[self.ridge_vertices[k][1],1]-nodes[node,1])/(nodes[self.ridge_vertices[k][1],0]-nodes[node,0])*(0.0-nodes[node,0])
			            nodes=np.append(nodes,[[0.0,y]], axis=0)
			        elif self.dimension == 3:
			            y = nodes[node,1]+(nodes[self.ridge_vertices[k][1],1]-nodes[node,1])/(nodes[self.ridge_vertices[k][1],0]-nodes[node,0])*(0.0-nodes[node,0])
			            z = nodes[node,2]+(nodes[self.ridge_vertices[k][1],2]-nodes[node,2])/(nodes[self.ridge_vertices[k][1],0]-nodes[node,0])*(0.0-nodes[node,0])
			            nodes=np.append(nodes,[[0.0,y,z]], axis=0)
			        self.ridge_vertices[k][0]=len(nodes)-1
			    elif self.ridge_vertices[k][1]==node:
			        if self.dimension == 2:
			            y = nodes[node,1]+(nodes[self.ridge_vertices[k][0],1]-nodes[node,1])/(nodes[self.ridge_vertices[k][0],0]-nodes[node,0])*(0.0-nodes[node,0])
			            nodes=np.append(nodes,[[0.0,y]], axis=0)
			        if self.dimension == 3:
			            y = nodes[node,1]+(nodes[self.ridge_vertices[k][0],1]-nodes[node,1])/(nodes[self.ridge_vertices[k][0],0]-nodes[node,0])*(0.0-nodes[node,0])
			            z = nodes[node,2]+(nodes[self.ridge_vertices[k][0],2]-nodes[node,2])/(nodes[self.ridge_vertices[k][0],0]-nodes[node,0])*(0.0-nodes[node,0])
			            nodes=np.append(nodes,[[0.0,y,z]], axis=0)
			        self.ridge_vertices[k][1]=len(nodes)-1
		for node in nodes_to_delete_right:
			node= int(node)
			for k in range(len(self.ridge_vertices)):
			    if self.ridge_vertices[k][0]==node:
			        if self.dimension==2:
			            y = nodes[node,1]+(nodes[self.ridge_vertices[k][1],1]-nodes[node,1])/(nodes[self.ridge_vertices[k][1],0]-nodes[node,0])*(self.length[0]-nodes[node,0])
			            nodes=np.append(nodes,[[self.length[0],y]], axis=0)
			        if self.dimension == 3:
			            y = nodes[node,1]+(nodes[self.ridge_vertices[k][1],1]-nodes[node,1])/(nodes[self.ridge_vertices[k][1],0]-nodes[node,0])*(self.length[0]-nodes[node,0])
			            z = nodes[node,2]+(nodes[self.ridge_vertices[k][1],2]-nodes[node,2])/(nodes[self.ridge_vertices[k][1],0]-nodes[node,0])*(self.length[0]-nodes[node,0])
			            nodes=np.append(nodes,[[self.length[0],y,z]], axis=0)
			        self.ridge_vertices[k][0]=len(nodes)-1
			    elif self.ridge_vertices[k][1]==node:
			        if self.dimension == 2:
			            y = nodes[node,1]+(nodes[self.ridge_vertices[k][0],1]-nodes[node,1])/(nodes[self.ridge_vertices[k][0],0]-nodes[node,0])*(self.length[0]-nodes[node,0])
			            nodes=np.append(nodes,[[self.length[0],y]], axis=0)
			        if self.dimension == 3:
			            y = nodes[node,1]+(nodes[self.ridge_vertices[k][0],1]-nodes[node,1])/(nodes[self.ridge_vertices[k][0],0]-nodes[node,0])*(self.length[0]-nodes[node,0])
			            z = nodes[node,2]+(nodes[self.ridge_vertices[k][0],2]-nodes[node,2])/(nodes[self.ridge_vertices[k][0],0]-nodes[node,0])*(self.length[0]-nodes[node,0])
			            nodes=np.append(nodes,[[self.length[0],y,z]], axis=0)
			        self.ridge_vertices[k][1]=len(nodes)-1
		nodes_to_delete=np.array(sorted(nodes_to_delete, reverse=True))
		for point in nodes_to_delete:
			nodes=np.delete(nodes,int(point),0)
		# Renumber points after deleting some
		for ridge in self.ridge_vertices:
			for i in range(2):
				r=0
				for node in nodes_to_delete:
					if node < ridge[i]:
						r+=1
				ridge[i] = ridge[i] - r 
		self.vertices=nodes
		return self
	
	def cut_network(self):
		self=self.cut_network_side()
		self=self.cut_network_updown()
		return self

	def delete_doubles(self):
		ridges = {tuple(np.sort(node)) for node in self.ridge_vertices}
		ridges = set(ridges)
		self.ridge_vertices = [list(l) for l in ridges]
		return self

	def delete_alone_points(self):
		self = self.create_ridge_node_list()
		nodes_to_delete = []
		for i in range(len(self.vertices)):
			if len(self.list_nodes_ridges[i])==0:
				nodes_to_delete = np.append(nodes_to_delete,i)
		nodes_to_delete=np.array(sorted(nodes_to_delete, reverse=True))
		for point in nodes_to_delete:
			self.vertices=np.delete(self.vertices,int(point),0)
		# Renumber points after deleting some
		for ridge in self.ridge_vertices:
			for i in range(2):
				r=0
				for node in nodes_to_delete:
					if node < ridge[i]:
						r+=1
				ridge[i] = ridge[i] - r
		return self
					
	def delete_single_ridge_points(self):
		self = self.delete_doubles()
		self = self.sort_nodes()
		self =self.create_ridge_node_list()
		ridges_to_delete = []
		for i in range(len(self.vertices)):
			if len(self.list_nodes_ridges[i])==1 and i in self.interior_nodes:
				ridges_to_delete.append(self.list_nodes_ridges[i][0])
		ridges_to_delete=np.array(sorted(ridges_to_delete, reverse=True))
		for ridge in ridges_to_delete:
			self.ridge_vertices=np.delete(self.ridge_vertices,ridge, axis=0)
		self = self.delete_alone_points()
		return self

	def distribution_length_fiber(self, path):
		lengths = []
		for ridge in self.ridge_vertices:
			lengths.append(np.linalg.norm(self.vertices[ridge[0]]-self.vertices[ridge[1]]))
		self.lengths_ini = lengths
		self.mean_length = np.mean(lengths)
		return self

	def save_network(self, step, path,**kw):
		if step == 'initial':
			if kw.get('name')!=None:
				filename = 'parameters_%03d_%s.csv' % (len(self.vertices),kw['name'])
			else:
				filename = 'parameters_%03d.csv' % len(self.vertices)
			filename = 'parameters_%03d.csv' % len(self.vertices)
			with open(os.path.join(path,filename), 'w') as writeFile:
				writer = csv.writer(writeFile)
				writer.writerow(["dimension",self.dimension])
				writer.writerow(["complexity",self.complexity])
				writer.writerow(["length",self.length])
				writer.writerow(["k_tension",self.k_tension])
				writer.writerow(["k_compression",self.k_compression])
				writer.writerow(["merge_distance",self.min_distance])
				writer.writerow(["number of nodes",len(self.vertices)])
				writer.writerow(["number of springs",len(self.ridge_vertices)])
				writer.writerow(["disturbance",self.disturbance])
				writer.writerow(["A",self.A])
				writer.writerow(["creation", self.creation])
				writer.writerow(["generation", self.generation])
			writeFile.close()
			last_network = 'network_ridge_vertices_%03d.csv' % len(self.vertices)
			with open(os.path.join(path,last_network), 'w') as writeFile:
				writer = csv.writer(writeFile)
				writer.writerows(self.ridge_vertices)
			writeFile.close()
			if self.kw.get('name')!=None:
				last_network = 'network_vertices_initial_%03d_%s.csv' % (len(self.vertices),self.kw['name'])
			else:
				last_network = 'network_vertices_initial_%03d.csv' % len(self.vertices)
			with open(os.path.join(path,last_network), 'w') as writeFile:
				writer = csv.writer(writeFile)
				writer.writerows(self.vertices)
			writeFile.close()
		elif step=='temp':
			last_network = 'network_vertices.csv'
			with open(os.path.join(path,last_network), 'w') as writeFile:
				writer = csv.writer(writeFile)
				writer.writerows(self.vertices)
			writeFile.close()
		else:
			if kw.get('name')!=None:
				last_network = 'network_vertices_%03d_%03d_%s.csv' % (int(step),len(self.vertices),kw['name'])
			else:
				last_network = 'network_vertices_%03d_%03d.csv' % (int(step),len(self.vertices))
			with open(os.path.join(path,last_network), 'w') as writeFile:
				writer = csv.writer(writeFile)
				writer.writerows(self.vertices)
			writeFile.close()
		
	def add_boundary_nodes(self):
		self = self.create_ridge_node_list()
		for node in self.interior_nodes:
			cond = True
			list_ridge=self.list_nodes_ridges[node]
			if self.vertices[node][0]<=float(self.length)/10.:
				for ridge_number in list_ridge:
					ridge = self.ridge_vertices[ridge_number]
					if self.vertices[ridge[0]][0]< self.vertices[node][0] or self.vertices[ridge[1]][0]<self.vertices[node][0]:
						cond = False
				if cond == True:
					self.vertices=np.append(self.vertices,[[0.,self.vertices[node][1]]],axis=0)
					self.ridge_vertices= np.append(self.ridge_vertices,[[len(self.vertices)-1,node]],axis=0)
		return self

# Global function that applies all the corrections to the network

	def set_fibers(self, path):
		if self.creation == 'old':
			filename = fnmatch.filter(os.listdir('.'), 'network_vertices_initial_087.csv')
			with open(filename[0],'r') as readFile:
				reader = csv.reader(readFile)
				list_vertices = np.array(list(reader))
				self.vertices=list_vertices.astype(float)
			filename = fnmatch.filter(os.listdir('.'), 'network_ridge_vertices_087.csv')
			with open(filename[0], 'r') as readFile:
				reader = csv.reader(readFile)
				list_ridge_vertices=np.array(list(reader))
				self.ridge_vertices=list_ridge_vertices.astype(int)
				self.ridge_vertices = [list(l) for l in self.ridge_vertices]
		elif self.creation == 'Voronoi' and self.generation == 'regular':
			self = self.create_network()
			self = self.delete_first_point()
			self = self.cut_network()
			self = self.delete_alone_points()
			self = self.delete_single_ridge_points()
			self = self.delete_doubles()
		elif self.creation == 'growth_network':
			self = self.create_network()
			self = self.delete_doubles()
			self = self.cut_network_updown()
			if self.dimension != 3:
				self = self.delete_doubles()
				self = self.cut_network_updown()
				self = self.delete_alone_points()
				self = self.delete_single_ridge_points()
				self.min_distance = 0.01
				self = self.merge_nodes()
				#self=self.delete_points_with_two_ridges()
				self = self.sort_nodes()
				self = self.create_ridge_node_list()
		elif self.creation != 'Voronoi' and self.creation != "growth_network":# and creation != 'reg_Voronoi':
			self = self.create_network()
			#self = self.delete_first_point()
			#self = self.cut_network()
			plot_geometry(self)
		elif self.creation == 'Voronoi' and self.generation != 'regular':
			self = self.create_network()
			while len(self.boundary_nodes_right) <= min(self.complexity/10,8) or len(self.boundary_nodes_left) <= min(self.complexity/10,8):
				self = self.create_network()
				if self.creation == 'Voronoi': self = self.cut_network();self = self.delete_first_point()
				elif self.creation == 'growth_network': self = self.cut_network_updown()
				self = self.delete_alone_points()
				self = self.delete_doubles()
				self = self.sort_nodes()
			print('number of boundary_nodes: ', len(self.boundary_nodes_right),len(self.boundary_nodes_left))
		if self.generation!= 'regular' and self.creation != 'old' and self.dimension!=3 and self.creation != '1 straight line':
			while len(self.ridge_vertices)-self.dimension*len(self.interior_nodes)<=1:
				self.min_distance +=0.001
				self = self.delete_doubles()
				self = self.delete_points_with_two_ridges()
				self = self.sort_nodes()
				self = self.create_ridge_node_list()
				self = self.merge_nodes()
				self = self.sort_nodes()
			print('hyperstatic number: ',len(self.ridge_vertices)-self.dimension*len(self.interior_nodes))
			print('number of boundary_nodes: ', len(self.boundary_nodes_right),len(self.boundary_nodes_left))
			print('min_distance: ', self.min_distance)
			self = self.delete_points_with_two_ridges()
		if self.creation!='1 straight line':
			self = self.delete_alone_points()
			self = self.delete_single_ridge_points()
			self = self.delete_alone_points()
		#self = self.create_ridge_node_list()
		self.vertices_ini = np.array(self.vertices.tolist())
		self = self.sort_nodes()
		#network.vertices_ini=np.array(self.vertices.tolist())
		self = self.distribution_length_fiber(path)
		self.save_network('initial', path)
		self = self.create_ridge_node_list()
		self.state_ridge = ['tension']*len(self.ridge_vertices)
		self.ridge_vertices = [list(l) for l in self.ridge_vertices]
		return self

	def length_ridge(self,x):	
		if self.dimension == 2:
			return np.sqrt(x[0]**2+x[1]**2)
		elif self.dimension == 3:
			return x[0]**2+x[1]**2+x[2]**2

	def elas_energy_total(self,vertices_interior):
		energy =0
		for ridge in self.ridge_vertices:
			if ridge[0] in self.interior_nodes:
				i = self.interior_nodes.index(ridge[0])
				if ridge[1] in self.interior_nodes:
					j = self.interior_nodes.index(ridge[1])
					energy += self.Ef*(self.length_ridge(vertices_interior[2*i:2*i+2]-vertices_interior[2*j:2*j+2])-self.length_ridge(self.vertices_ini[ridge[0]]-self.vertices_ini[ridge[1]]))**2
				else:
					energy += self.Ef*(self.length_ridge(vertices_interior[2*i:2*i+2]-self.vertices[ridge[1]])-self.length_ridge(self.vertices_ini[ridge[0]]-self.vertices_ini[ridge[1]]))**2
			else:
				if ridge[1] in self.interior_nodes:
					j = self.interior_nodes.index(ridge[1])
					energy += self.Ef*(self.length_ridge(self.vertices[ridge[0]]-vertices_interior[2*j:2*j+2])-self.length_ridge(self.vertices_ini[ridge[0]]-self.vertices_ini[ridge[1]]))**2
				else:
					energy += self.Ef*(self.length_ridge(self.vertices[ridge[0]]-self.vertices[ridge[1]])-self.length_ridge(self.vertices_ini[ridge[0]]-self.vertices_ini[ridge[1]]))**2
		#print energy
		return energy

	def set_fibers_explanation(self, path):
		import matplotlib.pyplot as plt
		from matplotlib.collections import LineCollection
		fig = plt.figure()
		ax1 = fig.add_subplot(221)
		ax2 = fig.add_subplot(222, sharex=ax1, sharey=ax1)
		ax3 = fig.add_subplot(223, sharex=ax1, sharey=ax1)
		ax4 = fig.add_subplot(224, sharex=ax1, sharey=ax1)
		Seeds=np.random.rand(self.complexity,self.dimension)*self.length[0]
		ax1.scatter(Seeds[:,0],Seeds[:,1])
		ax1.set_ylim([-0.1*self.length[1],self.length[1]*1.1])
		ax1.set_xlim([-0.1*self.length[0],self.length[0]*1.1])
		ax1.set_title('Seeds')
		####
		voronoi = Voronoi(Seeds)
		self.vertices = Voronoi(Seeds).vertices
		self.ridge_vertices = Voronoi(Seeds).ridge_vertices
		self = self.delete_first_point()
		ax2.scatter(self.vertices[:,0],self.vertices[:,1], s =2.)
		line_segments = []
		for simplex in self.ridge_vertices:
		        simplex = np.asarray(simplex)
		        line_segments.append([(x, y) for x, y in self.vertices[simplex]])
		lc = LineCollection(line_segments,linestyle='solid')
		ax2.set_title('Voronoi tesselation')
		ax2.add_collection(lc)
		####
		self = self.cut_network()
		ax3.scatter(self.vertices[:,0],self.vertices[:,1], s =2.)
		line_segments = []
		for simplex in self.ridge_vertices:
		        simplex = np.asarray(simplex)
		        line_segments.append([(x, y) for x, y in self.vertices[simplex]])
		lc = LineCollection(line_segments,linestyle='solid')
		ax3.set_title('Network cut')
		ax3.add_collection(lc)
		####
		self.min_distance = 0.035
		self = self.merge_nodes()
		ax4.scatter(self.vertices[:,0],self.vertices[:,1], s =2.)
		line_segments = []
		for simplex in self.ridge_vertices:
		        simplex = np.asarray(simplex)
		        line_segments.append([(x, y) for x, y in self.vertices[simplex]])
		lc = LineCollection(line_segments,linestyle='solid')
		ax4.set_title('Close points merged')
		ax4.add_collection(lc)
		plt.subplots_adjust(hspace=0.5)
		plt.savefig('../Data/default/creation_network.pdf')


