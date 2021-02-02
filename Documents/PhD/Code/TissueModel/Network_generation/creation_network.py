import numpy as np
import csv, os, fnmatch

# Class of the network, with different initial positions and all the corrections needed

class Network:
	def __init__(self, dimension, complexity_network, length_domain, min_distance, beam_Young, beam_poisson, beam_profile, connector_coeff, disturbance, hyperstatic_param, creation, generation, path, **kw):
		self.complexity = complexity_network
		self.length=length_domain
		self.dimension = dimension
		self.beam_young = beam_Young
		self.beam_poisson = beam_poisson
		self.beam_profile = beam_profile
		self.connector_coeff = connector_coeff
		self.disturbance = disturbance
		self.min_distance = min_distance
		self.creation = creation
		self.generation = generation
		self.lengths_ini = []
		self. hyperstatic_param =  hyperstatic_param
		if kw.get('name')!=None: self.kw = kw
		else: self.kw = kw

		self.strain = [0.]
		self.stress = [0.]
		
	def create_network(self):
		from scipy.spatial import Voronoi, voronoi_plot_2d
		import Network_generation.network_types
		self.vertices, self.ridge_vertices = Network_generation.network_types.select_network(self)
		self.interior_nodes = []
		self.boundary_nodes_left=[]
		self.boundary_nodes_right=[]
		return self

# Deletion of the first point because it creates intersections in the network
	def delete_first_point(self):
		ridges_to_delete = []
		for i in range(len(self.ridge_vertices)):
			if self.ridge_vertices[i][0]==-1 or self.ridge_vertices[i][1]==-1:
				ridges_to_delete.append(i)
		ridges_to_delete=np.array(sorted(ridges_to_delete, reverse=True))
		for k in ridges_to_delete:
			self.ridge_vertices=np.delete(self.ridge_vertices, int(k), axis=0)
		self = self.delete_points_with_two_ridges()
		return self

	def delete_points_with_two_ridges(self):
		self = self.sort_nodes()
		self = self.create_ridge_node_list()
		length_lists=[]
		for node in self.interior_nodes:
			length_lists.append(len(self.list_nodes_ridges[node]))
		while 2 in length_lists:
			self = self.create_ridge_node_list()
			for i in range(len(self.list_nodes_ridges)):
				if len(self.list_nodes_ridges[i])==2:
					self.delete_two_ridge_points(i)
			self = self.delete_doubles()
			self = self.create_ridge_node_list()
			self = self.delete_alone_points()
			self = self.sort_nodes()
			length_lists=[]
			for node in self.interior_nodes:
				length_lists.append(len(self.list_nodes_ridges[node]))
		return self


	def delete_point(self, i):
		self.vertices = np.delete(self.vertices, i, axis=0)
		for ridge in self.ridge_vertices:
			for i in range(2):
				if i < ridge[i]:
					ridge[i] -= 1
		return self

	def delete_two_ridge_points(self,i):
		self = self.create_ridge_node_list()
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
		self = self.create_ridge_node_list()
		nodes = np.array(self.vertices)
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
		self = self.create_ridge_node_list()
		nodes = self.vertices
		nodes_to_delete = []
		for i in range(len(nodes)):
			for j in range(i):
				if self.dimension == 2:
					distance = np.sqrt((nodes[i,0]-nodes[j,0])**2+(nodes[i,1]-nodes[j,1])**2)
				elif self.dimension == 3:
					distance = np.sqrt((nodes[i, 0] - nodes[j, 0]) ** 2 + (nodes[i, 1] - nodes[j, 1]) ** 2 +(nodes[i, 2] - nodes[j, 2]) ** 2  )
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
		self = self.create_ridge_node_list()
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
			self.ridge_vertices=np.delete(self.ridge_vertices,int(ridge), axis=0)
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
			self.ridge_vertices=np.delete(self.ridge_vertices,int(ridge), axis=0)
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
						self.ridge_vertices=np.delete(self.ridge_vertices,int(ridge), axis=0)
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
		self = self.create_ridge_node_list()
		nodes = self.vertices
		self = self.delete_doubles()
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
					elif self.dimension == 3:
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
					elif self.dimension == 3:
						y = nodes[node,1]+(nodes[self.ridge_vertices[k][1],1]-nodes[node,1])/(nodes[self.ridge_vertices[k][1],0]-nodes[node,0])*(self.length[0]-nodes[node,0])
						z = nodes[node,2]+(nodes[self.ridge_vertices[k][1],2]-nodes[node,2])/(nodes[self.ridge_vertices[k][1],0]-nodes[node,0])*(self.length[0]-nodes[node,0])
						nodes=np.append(nodes,[[self.length[0],y,z]], axis=0)
					self.ridge_vertices[k][0]=len(nodes)-1
				elif self.ridge_vertices[k][1]==node:
					if self.dimension == 2:
						y = nodes[node,1]+(nodes[self.ridge_vertices[k][0],1]-nodes[node,1])/(nodes[self.ridge_vertices[k][0],0]-nodes[node,0])*(self.length[0]-nodes[node,0])
						nodes=np.append(nodes,[[self.length[0],y]], axis=0)
					elif self.dimension == 3:
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

	def delete_boundary_to_boundary_fibre(self):
		self = self.create_ridge_node_list()
		ridges_to_delete=[]
		for i in range(len(self.ridge_vertices)):
			ridge = self.ridge_vertices
			if ridge[i][0] in self.boundary_nodes_left and ridge[i][1] in self.boundary_nodes_right:
				ridges_to_delete = np.append(ridges_to_delete, i)
			elif ridge[i][1] in self.boundary_nodes_left and ridge[i][0] in self.boundary_nodes_right:
				ridges_to_delete = np.append(ridges_to_delete, i)
			elif ridge[i][0] in self.boundary_nodes_left and ridge[i][1] in self.boundary_nodes_left:
				ridges_to_delete = np.append(ridges_to_delete, i)
			elif ridge[i][0] in self.boundary_nodes_right and ridge[i][1] in self.boundary_nodes_right:
				ridges_to_delete = np.append(ridges_to_delete, i)
		ridges_to_delete=np.array(sorted(ridges_to_delete, reverse=True))
		for ridge in ridges_to_delete:
			self.ridge_vertices=np.delete(self.ridge_vertices,int(ridge), axis=0)
		return self
	
	def delete_two_fibres_boundary_points(self):
		ridges_to_delete=[]
		self = self.create_ridge_node_list()
		for node in self.boundary_nodes_left:
			if len(self.list_nodes_ridges[node])>1:
				ridges_to_delete.append(self.list_nodes_ridges[node][1:len(self.list_nodes_ridges[node])])
		for node in self.boundary_nodes_right:
			if len(self.list_nodes_ridges[node])>1:
				ridges_to_delete.append(self.list_nodes_ridges[node][1:len(self.list_nodes_ridges[node])])
		from functools import reduce
		if len(ridges_to_delete)>0: ridges_to_delete=reduce(lambda x,y: x+y, ridges_to_delete)
		ridges_to_delete=np.array(sorted(ridges_to_delete, reverse=True))
		for ridge in ridges_to_delete:
			self.ridge_vertices=np.delete(self.ridge_vertices,int(ridge), axis=0)
		return self

	def create_attached_nodes_list(self):
		self = self.create_ridge_node_list()
		nodes=self.vertices
		ridges = self.ridge_vertices
		list_attached_nodes=[[] for i in range(len(nodes))]
		for i in range(len(ridges)):
			list_attached_nodes[ridges[i][0]].append(ridges[i][1])
			list_attached_nodes[ridges[i][1]].append(ridges[i][0])
		self.list_attached_nodes=list_attached_nodes
		return self 

	def delete_doubles(self):
		self = self.create_ridge_node_list()
		ridges = {tuple(np.sort(node)) for node in self.ridge_vertices}
		ridges = set(ridges)
		self.ridge_vertices = [list(l) for l in ridges]
		return self

	def delete_doubles_growth_network(self):
		self = self.create_ridge_node_list()
		self.ridge_vertices = np.array(self.ridge_vertices)
		if self.creation == 'growth_network':
			self = self.create_attached_nodes_list()
			ridges_to_delete=[]
			for k in range(len(self.vertices)):
				attached_nodes=self.list_attached_nodes[k]
				for i in range(len(attached_nodes)):
					node1=attached_nodes[i]
					for j in range(i, len(attached_nodes)):
						node2=attached_nodes[j]
						if [node1,node2] in self.ridge_vertices:
							if np.linalg.det(np.array([self.vertices[node1]-self.vertices[k], self.vertices[node1]-self.vertices[node2]]))<10e-10:
								if self.vertices[node2][0]<self.vertices[k][0]<self.vertices[node1][0] or self.vertices[node2][0]>self.vertices[k][0]>self.vertices[node1][0]:
									try:
										ridges_to_delete.append(self.ridge_vertices.tolist().index([node1,node2]))
									except ValueError:
										continue
						if [node2,node1] in self.ridge_vertices:
							if np.linalg.det(np.array([self.vertices[node1]-self.vertices[k], self.vertices[node1]-self.vertices[node2]]))<10e-10:
								if self.vertices[node2][0]<self.vertices[k][0]<self.vertices[node1][0] or self.vertices[node2][0]>self.vertices[k][0]>self.vertices[node1][0]:
									try:
										ridges_to_delete.append(self.ridge_vertices.tolist().index([node2,node1]))
									except ValueError:
										continue
			ridges_to_delete=np.array(sorted(set(ridges_to_delete), reverse=True))
			for ridge in ridges_to_delete:
				self.ridge_vertices=np.delete(self.ridge_vertices,int(ridge), axis=0)
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
		length_lists=[]
		for node in self.interior_nodes:
			length_lists.append(len(self.list_nodes_ridges[node]))
		while 1 in length_lists:
			ridges_to_delete = []
			for i in range(len(self.vertices)):
				if len(self.list_nodes_ridges[i])==1 and i in self.interior_nodes:
					ridges_to_delete.append(self.list_nodes_ridges[i][0])
			ridges_to_delete=np.array(sorted(ridges_to_delete, reverse=True))
			for ridge in ridges_to_delete:
				self.ridge_vertices=np.delete(self.ridge_vertices,ridge, axis=0)
			self =self.create_ridge_node_list()
			length_lists=[]
			for node in self.interior_nodes:
				length_lists.append(len(self.list_nodes_ridges[node]))
		self = self.delete_alone_points()
		return self

	def distribution_length_fiber(self, path):
		lengths = []
		for ridge in self.ridge_vertices:
			lengths.append(np.linalg.norm(self.vertices[ridge[0]]-self.vertices[ridge[1]]))
		self.lengths_ini = lengths
		self.mean_length = np.mean(lengths)
		return self

	def delete_crossings(self):
		ridges_to_delete=[]
		ridge_vertices=[[0,0]]
		buffer=10e-6
		for i in range(len(self.ridge_vertices)):
			ridge=self.ridge_vertices[i]
			for ridge1 in self.ridge_vertices:
				if ridge[0]!=ridge1[0] or ridge[1]!=ridge1[1]:
					a1 = (self.vertices[ridge[0]][1]-self.vertices[ridge[1]][1])/(self.vertices[ridge[0]][0]-self.vertices[ridge[1]][0])
					a2 = (self.vertices[ridge1[0]][1]-self.vertices[ridge1[1]][1])/(self.vertices[ridge1[0]][0]-self.vertices[ridge1[1]][0])
					b1 = self.vertices[ridge[0]][1]-a1*self.vertices[ridge[0]][0]
					b2 = self.vertices[ridge1[0]][1]-a2*self.vertices[ridge1[0]][0]
					if a1!=a2:
						x_inter=(b2-b1)/(a1-a2)
						y_inter = a1*x_inter+b1
						if self.vertices[ridge[0]][0]+buffer<x_inter<self.vertices[ridge[1]][0]-buffer or self.vertices[ridge[1]][0]+buffer<x_inter<self.vertices[ridge[0]][0]-buffer:
							if self.vertices[ridge[0]][1]<=y_inter<=self.vertices[ridge[1]][1] or self.vertices[ridge[1]][1]<=y_inter<=self.vertices[ridge[0]][1]:
								if self.vertices[ridge1[0]][0]+buffer<x_inter<self.vertices[ridge1[1]][0]-buffer or self.vertices[ridge1[1]][0]+buffer<x_inter<self.vertices[ridge1[0]][0]-buffer:
									if self.vertices[ridge1[0]][1]<=y_inter<=self.vertices[ridge1[1]][1] or self.vertices[ridge1[1]][1]<=y_inter<=self.vertices[ridge1[0]][1]:
										ridges_to_delete.append(i)
										if [x_inter,y_inter] not in self.vertices.tolist():
											self.vertices=np.append(self.vertices,[[x_inter,y_inter]],axis=0)
											ridge_vertices=np.append(ridge_vertices,[[len(self.vertices)-1,ridge[0]]],axis=0)
											ridge_vertices=np.append(ridge_vertices,[[len(self.vertices)-1,ridge[1]]],axis=0)
											ridge_vertices=np.append(ridge_vertices,[[len(self.vertices)-1,ridge1[0]]],axis=0)
											ridge_vertices=np.append(ridge_vertices,[[len(self.vertices)-1,ridge1[1]]],axis=0)
		ridges_to_delete = list(set([x for x in ridges_to_delete if ridges_to_delete.count(x) > 1]))
		ridges_to_delete=np.array(sorted(set(ridges_to_delete), reverse=True))
		for ridge in ridges_to_delete:
			self.ridge_vertices=np.delete(self.ridge_vertices,ridge,axis=0)
		if ridge_vertices[1:]!=[]:
			self.ridge_vertices=np.append(self.ridge_vertices,ridge_vertices[1:],axis=0)
		return self

	def save_network(self, step, path,**kw):
		filenames=sorted(fnmatch.filter(os.listdir(path), 'network_vertices_initial_*.csv'))
		number_network = len(filenames)
		if step == 'initial':
			if kw.get('name')!=None:
				filename = 'parameters_%05d_%s.csv' % (len(self.vertices),kw['name'])
			else:
				filename = 'parameters_%02d_%05d.csv' % (number_network,len(self.vertices))
			with open(os.path.join(path,filename), 'w') as writeFile:
				writer = csv.writer(writeFile)
				writer.writerow(["dimension",self.dimension])
				writer.writerow(["complexity",self.complexity])
				writer.writerow(["length",self.length])
				writer.writerow(["merge_distance",self.min_distance])
				writer.writerow(["number of nodes",len(self.vertices)])
				writer.writerow(["number of springs",len(self.ridge_vertices)])
				writer.writerow(["disturbance",self.disturbance])
				writer.writerow(["hyperstatic_param",self.hyperstatic_param])
				writer.writerow(["creation", self.creation])
				writer.writerow(["generation", self.generation])
				writer.writerow(["beam Young", self.beam_young])
				writer.writerow(["beam Poisson", self.beam_poisson])
				writer.writerow(["beam Profile", self.beam_profile])
				writer.writerow(["connector coeff", self.connector_coeff])
			writeFile.close()
			last_network = 'network_ridge_vertices_%02d_%05d.csv' % (number_network,len(self.vertices))
			with open(os.path.join(path,last_network), 'w') as writeFile:
				writer = csv.writer(writeFile)
				writer.writerows(self.ridge_vertices)
			writeFile.close()
			if self.kw.get('name')!=None:
				last_network = 'network_vertices_initial_%05d_%s.csv' % (len(self.vertices),self.kw['name'])
			else:
				last_network = 'network_vertices_initial_%02d_%05d.csv' % (number_network,len(self.vertices))
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
				last_network = 'network_vertices_%03d_%05d_%s.csv' % (int(step),len(self.vertices),kw['name'])
			else:
				last_network = 'network_vertices_%03d_%05d.csv' % (int(step),len(self.vertices))
			with open(os.path.join(path,last_network), 'w') as writeFile:
				writer = csv.writer(writeFile)
				writer.writerows(self.vertices)
			writeFile.close()
		
	def add_boundary_nodes(self):
		self = self.create_ridge_node_list()
		for node in self.interior_nodes:
			cond = True
			list_ridge=self.list_nodes_ridges[node]
			if self.vertices[node][0]<=float(self.length[1])/5.:
				for ridge_number in list_ridge:
					ridge = self.ridge_vertices[ridge_number]
					if self.vertices[ridge[0]][0]< self.vertices[node][0] or self.vertices[ridge[1]][0]<self.vertices[node][0]:
						cond = False
				if cond == True:
					if self.dimension==2: self.vertices=np.append(self.vertices,[[0.,self.vertices[node][1]]],axis=0)
					elif self.dimension==3: self.vertices=np.append(self.vertices,[[0.,self.vertices[node][1],self.vertices[node][2]]],axis=0)
					self.ridge_vertices= np.append(self.ridge_vertices,[[len(self.vertices)-1,node]],axis=0)
			if self.vertices[node][0]>=float(self.length[1])*4/5.:
				for ridge_number in list_ridge:
					ridge = self.ridge_vertices[ridge_number]
					if self.vertices[ridge[0]][0]< self.vertices[node][0] or self.vertices[ridge[1]][0]<self.vertices[node][0]:
						cond = False
				if cond == True:
					if self.dimension==2: self.vertices=np.append(self.vertices,[[1.,self.vertices[node][1]]],axis=0)
					elif self.dimension==3: self.vertices=np.append(self.vertices,[[1.,self.vertices[node][1],self.vertices[node][2]]],axis=0)
					self.ridge_vertices= np.append(self.ridge_vertices,[[len(self.vertices)-1,node]],axis=0)
		return self

	def add_boundary_nodes_adaptedtoGN(self):	#also checks if some nodes not linked to boundary
		self = self.create_ridge_node_list()
		for node in self.interior_nodes:
			self = self.create_attached_nodes_list()
			cond = False
			list_connected_points = np.array([self.vertices[i] for i in self.list_attached_nodes[node]])
			if self.vertices[node][0]>=max(list_connected_points[:,0]) and self.vertices[node][0]>=0.9:	#bc right
				if self.dimension==2: self.vertices=np.append(self.vertices,[[1.,self.vertices[node][1]]],axis=0)
				elif self.dimension==3: self.vertices=np.append(self.vertices,[[1.,self.vertices[node][1],self.vertices[node][2]]],axis=0)
				self.ridge_vertices= np.append(self.ridge_vertices,[[len(self.vertices)-1,node]],axis=0)
			if self.vertices[node][0]<=min(list_connected_points[:,0]) and self.vertices[node][0]<=0.1:
				if self.dimension==2: self.vertices=np.append(self.vertices,[[0.,self.vertices[node][1]]],axis=0)
				elif self.dimension==3: self.vertices=np.append(self.vertices,[[0.,self.vertices[node][1],self.vertices[node][2]]],axis=0)
				self.ridge_vertices= np.append(self.ridge_vertices,[[len(self.vertices)-1,node]],axis=0)
		return self


# Global function that applies all the corrections to the network

	def generate_network(self, path):
		self.min_distance = 0.0001*self.length[0]
		if self.creation == 'old':
			filename = fnmatch.filter(os.listdir('.'), 'network_vertices_initial_00_00400.csv')
			with open(filename[0],'r') as readFile:
				reader = csv.reader(readFile)
				list_vertices = list(reader)
				self.vertices=[np.array(list_vertices[i]).astype(float) for i in range(len(list_vertices))]
				self.vertices = np.array([list(l) for l in self.vertices if len(l)>0])
			filename = fnmatch.filter(os.listdir('.'), 'network_ridge_vertices_00400.csv')
			with open(filename[0], 'r') as readFile:
				reader = csv.reader(readFile)
				list_ridge_vertices=np.array(list(reader),dtype=object)
				self.ridge_vertices=[np.array(list_ridge_vertices[i]).astype(int) for i in range(len(list_ridge_vertices))]
				self.ridge_vertices = [list(l) for l in self.ridge_vertices if len(l)>0]
		elif self.creation == 'Voronoi' and self.generation == 'regular':
			self = self.create_network()
			self = self.delete_first_point()
			self = self.cut_network()
			self = self.delete_alone_points()
			self = self.delete_single_ridge_points()
			self = self.delete_doubles()
		elif self.creation == 'growth_network':
			from Plotting.network_plotting import plot_geometry
			self = self.create_network()
			self = self.delete_doubles()
			self = self.cut_network()
			self = self.delete_alone_points()
			self = self.delete_single_ridge_points()
			self.minimum_distance = 0.01
			self = self.merge_nodes()
			self = self.delete_single_ridge_points()
			self=self.delete_points_with_two_ridges()
			self = self.delete_single_ridge_points()
			self = self.sort_nodes()
			self = self.delete_boundary_to_boundary_fibre()
			self = self.delete_single_ridge_points()
			self = self.sort_nodes()
			self = self.create_ridge_node_list()
		elif self.creation != 'Voronoi' and self.creation != "growth_network":# and creation != 'reg_Voronoi':
			self = self.create_network()
			plot_geometry(self)
		elif self.creation == 'Voronoi' and self.generation == "grid":
			self = self.create_network()
			self = self.delete_first_point()
			self = self.cut_network()
			self = self.sort_nodes()
			while len(self.boundary_nodes_left)+len(self.boundary_nodes_right)<np.sqrt(self.complexity):
				print(len(self.boundary_nodes_left)+len(self.boundary_nodes_right),2*np.sqrt(self.complexity))
				self = self.add_boundary_nodes()
				self = self.sort_nodes()
			print(len(self.boundary_nodes_left),len(self.boundary_nodes_right))
		elif self.creation == 'Voronoi' and self.generation != 'regular' and self.generation != 'grid':
			self = self.create_network()
			while len(self.boundary_nodes_right) <= min(self.complexity/10,2) or len(self.boundary_nodes_left) <= min(self.complexity/10,2):
				self = self.create_network()
				if self.creation == 'Voronoi': self = self.cut_network();self = self.delete_first_point()
				elif self.creation == 'growth_network': self = self.cut_network_updown()
				self = self.delete_alone_points()
				self = self.merge_nodes()
				self = self.delete_doubles()
				self = self.sort_nodes()
				print('number of boundary_nodes: ', len(self.boundary_nodes_right),len(self.boundary_nodes_left))
		#if len(self.boundary_nodes_left)>0 and len(self.boundary_nodes_right)>0:
		if self.generation!= 'regular' and self.creation != 'old' and self.creation != '1 straight line':# and self.creation != 'growth_network':
			while len(self.ridge_vertices)-float(self.hyperstatic_param)*len(self.interior_nodes)<=0:
				#if len(self.boundary_nodes_left)>0 and len(self.boundary_nodes_right)>0:
					self.min_distance +=0.001
					print('minimum distance: ', self.min_distance)
					print('hyperstatic number: ',len(self.ridge_vertices)-float(self.hyperstatic_param)*len(self.interior_nodes))
					self = self.delete_doubles()
					self = self.delete_points_with_two_ridges()
					self = self.sort_nodes()
					self = self.merge_nodes()
					self = self.sort_nodes()
					#self = self.add_boundary_nodes()
					self = self.delete_points_with_two_ridges()
					self = self.sort_nodes()
					self = self.delete_boundary_to_boundary_fibre()
					self = self.sort_nodes()
					self = self.delete_two_fibres_boundary_points()
				#else:
				#	print('No boundary nodes')
				#	break
			print('hyperstatic number: ',len(self.ridge_vertices)-float(self.hyperstatic_param)*len(self.interior_nodes))
			print('number of boundary_nodes: ', len(self.boundary_nodes_right),len(self.boundary_nodes_left))
			print('min_distance: ', self.min_distance)
		if self.creation!='1 straight line':
			self = self.delete_alone_points()
			self = self.delete_single_ridge_points()
			self = self.delete_alone_points()
		self = self.merge_nodes()
		self.vertices_ini = np.array(self.vertices.tolist())
		self = self.sort_nodes()
		self = self.add_boundary_nodes()
		self.min_distance = 10e-3
		for k in range(4):
			print(k)
			if self.creation == 'growth_network': self = self.add_boundary_nodes_adaptedtoGN()
			self = self.add_boundary_nodes()
			if self.dimension == 2: self = self.delete_crossings()
			self = self.delete_single_ridge_points()
			self = self.sort_nodes()
			self = self.delete_points_with_two_ridges()
			self = self.delete_doubles()
			self = self.merge_nodes()
			self = self.sort_nodes()
			self = self.delete_boundary_to_boundary_fibre()
			if self.dimension == 2: self = self.delete_doubles_growth_network()
			from Plotting.network_plotting import plot_geometry
			import matplotlib.pyplot as plt
		#else:
		#	print('No boundary nodes')
		self.state_ridge = ['tension']*len(self.ridge_vertices)
		self.ridge_vertices = [list(l) for l in self.ridge_vertices]
		self.save_network('initial', path)
		return self

	def set_fibers(self,path):
		self = self.generate_network(path)
		while len(self.boundary_nodes_left)==0 or len(self.boundary_nodes_right)==0: 
			print('No boundary nodes,starting again')
			self = self.generate_network(path)
		return self
