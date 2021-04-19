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
		if self.creation == 'growth_network' and self.dimension ==2:
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
		print('3.5 ',len(self.ridge_vertices), len(self.vertices))
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
			from numpy import array
			#self.vertices = array([array([ 1.        , -0.14131224,  2.01291918]), array([1.        , 4.81691333, 5.00903468]), array([1.        , 0.56002956, 1.02526791]), array([0.34225941, 0.38934292, 0.31279518]), array([ 1.        ,  0.76776869, -2.30605185]), array([ 1.        ,  1.37457455, -1.03588721]), array([1.        , 0.46913882, 0.13638542]), array([1.        , 1.11399981, 1.13588783]), array([1.        , 0.62986422, 0.71640301]), array([1.        , 1.54005956, 0.55651846]), array([ 1.        , -0.30954281,  1.93686526]), array([0.28867379, 0.3716109 , 0.3564105 ]), array([1.        , 0.62341191, 0.73436899]), array([ 1.        ,  0.57479703, -0.06480281]), array([0.4299194 , 0.33099411, 0.71401939]), array([0.73596098, 0.6068841 , 0.53255922]), array([0.16629807, 0.28897343, 0.44760019]), array([0.79570057, 0.76937491, 0.7439834 ]), array([ 1.        ,  1.43413902, -5.90286115]), array([1.        , 2.11403839, 7.62853439]), array([  1.        , -11.24875549,  -2.69990818]), array([ 1.        , -0.19432905,  0.38055463]), array([0.28122633, 0.36307508, 0.34744822]), array([ 1.        , -3.72099735,  4.85104445]), array([0.43752615, 0.46964893, 0.42022496]), array([1.        , 0.4821586 , 0.79336295]), array([0.27501147, 0.3526628 , 0.38382613]), array([0.42209649, 0.56557854, 0.20815474]), array([0.34313806, 0.49565864, 0.18825196]), array([0.78345139, 0.67172222, 0.8116053 ]), array([ 1.        , -0.43991924, -0.09460948]), array([1.        , 0.10945877, 0.8514108 ]), array([0.11749761, 0.37124424, 0.67003405]), array([ 1.        ,  0.89162275, -0.17222027]), array([0.70767294, 0.3860642 , 0.9569931 ]), array([0.1478201 , 0.42121404, 0.30809927]), array([1.        , 4.30122148, 1.54976541]), array([0.71863526, 0.57154476, 0.59235986]), array([1.        , 3.86552585, 6.85365374]), array([0.25647157, 0.33470272, 0.31765833]), array([ 1.        ,  0.56047809, -1.02360412]), array([1.        , 7.62751385, 3.21300663]), array([ 1.        ,  6.26368662, -0.26985564]), array([0.10490397, 0.16532974, 0.37240364]), array([ 1.        ,  1.85369528, -1.48971603]), array([ 1.        , -2.19454015,  2.66256886]), array([ 1.        ,  4.48937555, -7.31433446]), array([0.24565866, 0.5319537 , 0.25967429]), array([0.65792266, 0.79794489, 0.5076317 ]), array([0.61477463, 0.73502498, 0.5413749 ]), array([ 1.        ,  0.51121696, -1.24103822]), array([1.        , 1.21345788, 0.54094173]), array([0.26417911, 0.28359229, 0.71024428]), array([-0.0351766 ,  0.62203527,  0.47049268]), array([ 1.        ,  2.75720863, -1.45343809]), array([ 1.        ,  4.57950825, -3.64024766]), array([0.49436576, 0.60736224, 0.60394032]), array([0.98763744, 0.7821412 , 0.70951529]), array([ 1.        , -0.10912358,  0.89839669]), array([0.8157497 , 0.57989641, 0.8805405 ]), array([1.        , 1.14780834, 0.27946668]), array([0.27186074, 0.35803104, 0.707144  ]), array([0.27342002, 0.72512309, 0.28358295]), array([0.41210924, 0.53189966, 0.31462061]), array([ 1.        , -0.10609421,  1.50889914]), array([1.        , 1.12789225, 0.1833722 ]), array([0.30950993, 0.52386193, 0.28607199]), array([0.65839562, 0.7658841 , 0.64823161]), array([1.        , 0.8646841 , 0.88252156]), array([-0.19942984,  0.50312156,  0.31207341]), array([0.78925128, 0.68068633, 0.83122386]), array([0.8208108 , 0.62082084, 0.76592232]), array([ 1.00000000e+00,  8.10346024e+04, -5.96972676e+04]), array([ 1.        , -0.42237035,  0.77353909]), array([0.49702654, 0.65428991, 0.33198483]), array([0.41374514, 0.37165719, 0.47890751]), array([0.22049866, 0.59224012, 0.24579343]), array([0.55701689, 0.64135835, 0.48780837]), array([0.65865606, 0.46193078, 0.60913528]), array([ 1.        , -4.1253005 , -8.46096968]), array([ 1.        , -1.26245636, -0.26703234]), array([ 1.        ,  0.21718146, -0.05085768]), array([0.56086145, 0.65567863, 0.51398948]), array([1.        , 1.61035866, 1.26869511]), array([1.        , 0.72757748, 0.56803931]), array([ 1.        , -4.89336502,  0.7403044 ]), array([0.28892971, 0.36211023, 0.41526705]), array([1.        , 4.14226232, 0.09357856]), array([0.30634534, 0.20932914, 0.65414439]), array([1.        , 0.16538294, 0.37399579]), array([ 1.        ,  0.52137967, -6.65393521]), array([0.82809521, 0.60695094, 0.50572011]), array([ 1.        ,  0.37787139, -0.27509968]), array([ 1.        , -2.11133236, -2.43216018]), array([  1.        , -17.20787977,   2.46774629]), array([0.73685566, 0.49670495, 0.29596188]), array([ 1.        ,  0.38762394, -0.00230266]), array([0.48420687, 0.34664356, 0.65904063]), array([ 1.        , -3.06515838,  6.01201853]), array([1.        , 1.24832502, 0.17057429]), array([0.52825162, 0.35787097, 0.38959727]), array([0.76398365, 0.78547523, 0.6238551 ]), array([0.57145538, 0.56291895, 0.52859752]), array([ 1.        , -1.16435476, -0.75366395]), array([0.49927159, 0.41077818, 0.71618993]), array([1.        , 1.35224528, 3.18514668]), array([0.64692935, 0.29809416, 0.90719757]), array([1.        , 0.19834651, 0.95327616]), array([  1.        , -12.97511998,  13.54254324]), array([1.        , 0.1411284 , 1.58723109]), array([0.32566903, 0.41401254, 0.40093059]), array([1.        , 0.28413894, 0.14350783]), array([ 1.        ,  1.16697175, -0.34445575]), array([0.74187377, 0.5739529 , 0.38786967]), array([ 1.        , -0.74196991,  2.28906434]), array([ 1.        ,  5.26614997, -3.84031148]), array([ 1.        , -2.19844689,  3.39771273]), array([0.65406314, 0.60310696, 0.37687178]), array([1.        , 0.82164508, 1.17142882]), array([1.        , 3.76521984, 8.40692236]), array([0.37517607, 0.62929416, 0.32851358]), array([0.49217101, 0.54880999, 0.39332046]), array([0.27490363, 0.19949484, 0.53634058]), array([0.60303651, 0.11582554, 0.66327453]), array([  1.        , -23.19154494,  -4.72905559]), array([ 1.        ,  0.7445371 , -1.46274852]), array([ 1.        ,  2.78615652, -1.96720951]), array([ 1.        , -0.17804728, -0.12172181]), array([ 1.        , -0.4729228 , -0.22845005]), array([0.62011301, 0.43024069, 0.44157247]), array([0.44306311, 0.40148928, 0.51172545]), array([ 1.        , -1.69953396,  2.20315126]), array([ 1.        , -0.7867652 ,  1.43525356]), array([ 1.        ,  2.07168065, -0.40598591]), array([1.        , 0.63476625, 6.63590547]), array([1.        , 0.77240289, 0.28873466]), array([1.        , 1.09081765, 1.29843049]), array([0.21457228, 0.04803875, 0.51205137]), array([ 1.        ,  1.31238188, -3.67778618]), array([0.82378296, 0.77045051, 0.32498663]), array([ 1.        ,  0.84896732, -0.91332248]), array([1.        , 0.14960969, 0.29966057]), array([0.20310044, 0.46177802, 0.43470993]), array([0.19031819, 0.45239858, 0.40543426]), array([1.        , 0.86064291, 0.27196012]), array([1.        , 0.31807614, 0.80042194]), array([0.53392029, 0.88175637, 0.48591099]), array([1.        , 1.61779948, 0.72757132]), array([1.        , 0.13000401, 0.34469542]), array([ 1.        , -2.31069462,  0.40916031]), array([1.        , 1.3714676 , 1.46890177]), array([0.43669546, 0.65032245, 0.5669187 ]), array([ 1.        , 12.07438156, 15.24196871]), array([ 1.        ,  1.97146631, -3.71332763]), array([1.        , 0.81594575, 0.56111075]), array([0.45060716, 0.51149992, 0.6142986 ]), array([1.        , 3.73527676, 6.15956119]), array([  1.        , -15.97375338,  -8.5094825 ]), array([ 1.        , -2.36924056, -0.83740676]), array([1.        , 2.02674065, 1.17971612]), array([ 1.        ,  1.00282388, -0.22670725]), array([0.16305831, 0.51883025, 0.27880075]), array([  1.        , -18.70655032, -28.07345816]), array([0.49115601, 0.31464235, 0.33915302]), array([  1.        , -11.68271183,  -9.07873205]), array([1.        , 0.33199924, 0.95089209]), array([ 1.        , -1.55777873,  0.59294611]), array([1.        , 0.17735697, 0.54518977]), array([1.        , 0.01446016, 0.17543484]), array([0.60758069, 0.4209716 , 0.55272726]), array([1.        , 0.84236914, 0.83553995]), array([1.        , 3.58000803, 3.80561009]), array([ 1.        ,  0.84334386, -1.18283412]), array([1.        , 0.23806584, 2.22155478]), array([1.        , 2.07085799, 2.05401509]), array([ 1.        , -0.03065749,  0.23060631]), array([ 1.        , -0.41294113,  1.76083517]), array([1.        , 0.98682558, 0.73186138]), array([1.        , 0.36117434, 0.30339991]), array([ 1.        , -5.05118487,  7.81666645]), array([ 1.        , -2.10471204, -3.56767104]), array([ 1.        , -4.16895966, -2.72685728]), array([0.55228144, 0.62371967, 0.45556039]), array([0.6878947 , 0.36558807, 0.49313861]), array([0.79052924, 0.67574405, 0.46797547]), array([ 1.       ,  0.8675691, -0.1341263]), array([1.        , 0.42734971, 2.00506244]), array([ 1.        , -1.08426906,  0.77786484]), array([ 1.        ,  0.0312219 , -0.02749693]), array([0.75941413, 0.32899267, 0.65103231]), array([1.        , 2.8420736 , 5.57622615]), array([ 1.        ,  5.76610313, -1.1636538 ]), array([ 1.        ,  1.89725272, -2.02854622]), array([1.        , 0.50638453, 0.10369249]), array([ 1.        ,  0.32930387, -0.3786238 ]), array([ 1.        , -0.47690617, -0.5123451 ]), array([0.72880193, 0.587258  , 0.62674939]), array([ 1.        ,  2.653393  , -2.94583177]), array([ 1.        ,  1.67426079, -2.07135249]), array([0.65028112, 0.52040844, 0.57292843]), array([ 1.        ,  0.11188942, 17.13286482]), array([1.        , 0.01880301, 0.0324151 ]), array([0.52465078, 0.43654193, 0.60233175]), array([0.43762597, 0.54233062, 0.53565962]), array([0.38948932, 0.39659423, 0.57324004]), array([ 1.        , -1.5065243 ,  0.84911107]), array([1.        , 2.7907563 , 4.96840416]), array([ 1.        ,  1.31960953, -4.67156172]), array([   1.        ,  -60.2308953 , -141.09473237]), array([1.        , 0.04373941, 0.03563928]), array([1.        , 0.21013333, 2.28585403]), array([ 1.        ,  1.09991143, -0.30388425]), array([0.85501035, 0.71562757, 0.58777274]), array([1.        , 8.34526   , 4.03254712]), array([ 1.        , -1.03011388, -1.21115429]), array([0.68177109, 0.4043535 , 0.79834342]), array([0.32428214, 0.557952  , 0.02065461]), array([0.79352948, 0.87850355, 0.67406597]), array([0.5002353 , 0.61543314, 0.45220888]), array([1.        , 1.66207836, 0.99751545]), array([ 1.        , 20.17285559, 12.66785013]), array([0.6618719 , 0.32382581, 0.6331845 ]), array([0.3920383 , 0.5849604 , 0.55762882]), array([1.        , 0.05095901, 2.93961112]), array([0.56007266, 0.6666609 , 0.43134199]), array([1.        , 5.13382324, 0.72858836]), array([1.        , 1.13144748, 0.05842851]), array([1.        , 0.56680732, 0.01288796]), array([ 1.        ,  1.61874858, -0.24382543]), array([1.        , 2.59972392, 2.73421155]), array([ 1.        ,  7.98726757, -1.73116728]), array([ 1.        ,  0.76744747, -1.23309956]), array([0.78663644, 0.89136634, 0.7046433 ]), array([1.        , 6.19882076, 4.2267132 ]), array([  1.        ,  72.62521397, -38.23977476]), array([ 1.        , -0.74377854,  1.52438165]), array([ 1.        , -4.8664516 ,  6.28400958]), array([0.50591626, 0.55814295, 0.50149552]), array([0.75335727, 0.26340939, 0.3851052 ]), array([1.        , 0.67488445, 0.57624304]), array([0.52254073, 0.87394187, 0.44800087]), array([1.        , 1.73136587, 2.59467366]), array([1.        , 3.52010592, 2.14691331]), array([ 1.        , -1.91176312, -0.66646613]), array([1.        , 3.04647852, 2.15433209]), array([ 1.        , -0.04495035,  0.50754347]), array([ 1.        , 56.08484085, -1.17768867]), array([0.67708654, 0.62060172, 0.32969649]), array([ 1.        , -1.61845942,  0.53831414]), array([ 1.        , -2.23744007,  4.92580436]), array([1.        , 0.05106827, 1.90074811]), array([1.        , 0.25604753, 2.04176815]), array([ 1.        , -0.02600798, -0.76385248]), array([0.67819498, 0.53512248, 0.50903649]), array([1.        , 8.39039674, 9.38120764]), array([1.02222556, 1.01283365, 0.67017903]), array([1.        , 0.09713548, 0.76930098]), array([1.        , 0.70046799, 0.2379698 ]), array([ 1.        , -0.22765303,  1.76095828]), array([ 1.        , -1.96307169,  1.61914266]), array([1.        , 2.52687344, 0.8707168 ]), array([0.56110728, 0.51869728, 0.47364333]), array([ 1.        ,  0.77987712, -0.10697404]), array([0.51139569, 0.11432973, 0.72143262]), array([0.56386959, 0.32173515, 0.68896038]), array([ 1.        , -0.07126444,  4.0206719 ]), array([ 1.        , -0.59099314,  1.20119951]), array([1.        , 0.8114914 , 1.04493551]), array([1.        , 0.81751684, 1.01859461]), array([ 1.        , -4.50687157,  4.39761736]), array([ 1.        ,  1.04627559, -1.56800801]), array([ 1.        ,  1.0242117 , -1.09303896]), array([1.        , 0.35702925, 2.82086575]), array([ 1.        , -0.11320196, -0.61693295]), array([ 1.        , -0.11589805,  0.60847908]), array([1.        , 0.72200293, 0.00639289]), array([0.92969974, 0.51232137, 0.3125155 ]), array([0.3528075 , 0.40546924, 0.55956521]), array([ 1.        , -0.85035667,  1.66828768]), array([0.55313064, 0.62688279, 0.46134338]), array([1.01622455, 0.62298533, 0.6989393 ]), array([1.        , 2.49271814, 0.3232829 ]), array([ 1.        ,  0.07620506, -0.51971455]), array([0.58639382, 0.62675701, 0.60792897]), array([1.        , 0.60717842, 0.66940215]), array([0.64477293, 0.50758148, 0.57032256]), array([0.87689594, 0.51609378, 0.92927565]), array([1.        , 0.11833208, 0.00172752]), array([ 1.        ,  1.8045481 , -3.10990513]), array([ 1.        , -0.33383606,  0.73399217]), array([ 1.        ,  2.15420086, -4.78585712]), array([  1.        , -43.28067395, -19.7345204 ]), array([1.        , 1.23920668, 0.48147165]), array([0.73843523, 0.60214687, 0.65933474]), array([1.        , 0.37930594, 1.44879388]), array([1.        , 1.6639609 , 1.20042534]), array([1.        , 1.04653246, 2.2598793 ]), array([0.65450603, 0.94535796, 0.53700563]), array([1.        , 0.70660392, 1.15860071]), array([ 1.        ,  0.78773119, -0.64143185]), array([0.63492617, 0.60618059, 0.71648155]), array([1.        , 0.36626559, 0.73325144]), array([0.60748959, 0.73002992, 0.52841999]), array([ 1.        , -1.99094716, -1.73927139]), array([0.78632428, 0.37184968, 0.71476619]), array([1.25250652, 0.82215357, 0.71589333]), array([0.6670228 , 0.95929117, 0.67636659]), array([0.49633931, 0.85594899, 0.36071282]), array([ 1.        , -5.79244759,  4.15875683]), array([0.77445846, 0.6547088 , 0.54823567]), array([1.        , 3.23553271, 3.06397813]), array([0.89939663, 0.53924896, 0.66240987]), array([ 1.        ,  0.95927064, -2.33730915]), array([ 1.        , -0.16322033,  0.96848305]), array([1.        , 1.17089207, 1.14915932]), array([ 1.        , -0.79821544,  1.9574441 ]), array([1.        , 3.32599613, 1.54079287]), array([ 1.        , -1.01556723,  1.79426955]), array([0.69156298, 0.61654185, 0.59245849]), array([ 1.        , -2.46934561,  0.62057891]), array([1.        , 0.65586231, 1.58888057]), array([1.        , 1.01391466, 1.21780552]), array([1.        , 1.94527298, 0.72885228]), array([ 1.        , -5.71115588, -0.67407967]), array([1.        , 0.42094208, 0.65762358]), array([ 1.        ,  2.27965496, -0.12874935]), array([1.        , 0.06108848, 0.44194409]), array([1.        , 1.89577444, 2.80964822]), array([ 1.        , -1.56734862, -0.31798952]), array([1.        , 0.15370704, 0.91057925]), array([1.        , 0.22666518, 1.22821767]), array([1.        , 0.80941788, 1.19663286]), array([ 1.        , -2.29039083,  1.25803891]), array([1.        , 4.02917905, 1.28523699]), array([ 1.        ,  0.17937879, -0.17547594]), array([0.81828932, 0.42275709, 0.79047203]), array([1.        , 0.38821396, 0.22937448]), array([0.77548393, 0.83923364, 0.66486935]), array([ 1.        , -3.8072131 ,  5.03737269]), array([ 1.        ,  0.08813177, -0.38346178]), array([1.        , 0.15523773, 0.96209975]), array([1.        , 1.86831235, 0.31478996]), array([ 1.        , -0.41232839,  3.36479128]), array([1.        , 2.39426497, 2.80440797]), array([0.90562981, 0.84940836, 0.83914652]), array([1.        , 0.90875276, 0.59275192]), array([1.        , 0.08298973, 0.83914282]), array([ 1.        , -0.10318189,  0.8702681 ]), array([ 1.        , -0.52539404,  0.8513201 ]), array([ 1.        , -0.10062047, -0.15813409]), array([0.7875641 , 0.74678094, 0.70095783]), array([ 1.        , -1.21360179,  3.62085575]), array([ 1.        , -5.10423332,  3.63115382]), array([1.        , 0.06470307, 0.64468861]), array([ 1.        , -0.14331575, -2.03523857]), array([1.        , 0.22545522, 0.50423102]), array([1.        , 0.77353936, 0.13137833]), array([ 1.        , -0.00876118,  1.01316893]), array([1.        , 0.39490154, 0.36402905]), array([1.        , 0.54202135, 1.94573461]), array([ 1.        , -0.0294532 ,  0.10613414]), array([ 1.        , 16.35491582,  2.15745943]), array([ 1.        ,  0.94151213, -0.7910782 ]), array([1.        , 0.69800352, 0.09526804]), array([1.        , 0.12319316, 2.20586702]), array([ 1.        ,  0.48042587, -0.32558581]), array([1.        , 0.71865669, 0.81850113]), array([1.        , 0.05371406, 0.52113673]), array([1.        , 0.84477097, 2.02155531]), array([ 1.        , -0.42722305,  1.22866023]), array([ 1.        , -6.21813319,  9.16980202]), array([ 1.        , -1.68908169, -0.8161756 ]), array([ 1.        ,  0.16725455, -0.81216095]), array([1.        , 0.71186438, 2.78010942]), array([0.94222791, 1.09072047, 0.73680264]), array([1.        , 0.27394323, 2.51779089]), array([ 1.        , -1.59319886,  0.89077893]), array([1.        , 0.32761263, 0.13035366]), array([ 1.        ,  1.04772108, -1.17638881]), array([ 1.        , -0.38243467, -0.07965983]), array([1.        , 1.03988761, 0.23611084]), array([1.        , 1.74763933, 0.99265721]), array([1.        , 0.57674707, 0.8059142 ]), array([ 1.        ,  0.77344363, -0.18167568]), array([  1.        ,  14.41857303, -21.95053003]), array([0.92657573, 0.77490041, 0.8174383 ]), array([ 1.        , -1.03797385, -0.3937006 ]), array([1.        , 2.737113  , 0.05810584]), array([ 1.        , -0.64966268,  0.17985058]), array([1.        , 1.72501718, 0.11396755]), array([1.        , 0.33766934, 0.33475558]), array([ 1.        ,  0.47779648, -0.12559595]), array([1.        , 1.07646258, 0.10398689]), array([0.78320891, 0.84928103, 0.74638959]), array([1.        , 1.03036935, 0.59312324]), array([0.82575168, 0.43464161, 0.80814584]), array([ 1.        ,  1.05985282, -1.47309642]), array([ 1.        ,  1.53090886, -0.17522798]), array([ 1.        , -0.74653243,  0.87723077]), array([0.81089087, 0.4758125 , 0.63473613]), array([ 1.        , -1.52076817, -1.69718277]), array([1.        , 0.90756976, 0.21114983]), array([1.        , 0.0349472 , 0.10279024]), array([1.        , 1.31862289, 1.99649614]), array([1.        , 1.47576653, 1.29252204]), array([1.        , 1.88874017, 1.46557626]), array([ 1.        , -1.00308143, -0.79790986]), array([1.        , 0.71214888, 1.22083472]), array([ 1.        , -0.51568374,  0.28309243]), array([ 1.        , -1.00235887,  0.52393501]), array([ 1.        , -0.35994325,  1.26341936]), array([0.93925852, 0.82875324, 0.89831091]), array([1.        , 0.83186629, 0.1414568 ]), array([ 1.        , -0.24327174,  1.40579002]), array([ 1.        ,  0.00962897, -0.40525752]), array([ 1.        , -0.3769887 , -2.27791441]), array([1.        , 7.52704041, 5.2793143 ]), array([ 1.        , -0.30280684,  1.78785157]), array([1.        , 0.1109287 , 0.89767666]), array([1.        , 0.3324481 , 0.50549245]), array([ 1.        , -2.31394608,  1.35756848]), array([ 1.        , -0.39368658,  1.17568341]), array([1.        , 0.13765535, 1.70805747]), array([ 1.        ,  0.48817981, -0.72128734]), array([1.        , 2.48612106, 0.74874211]), array([ 1.        , -0.02547329,  0.81972986]), array([1.        , 0.59800856, 0.78549863]), array([1.        , 0.60936021, 0.28939042]), array([1.        , 5.9195122 , 2.46271126]), array([ 1.        , -0.71351021,  2.00134376]), array([ 1.        , -0.08148031,  0.81627148]), array([ 1.        , -0.14778538,  1.55527023]), array([ 1.        , -0.15715993, -0.41731276]), array([ 1.        ,  0.03756828, -0.18630163]), array([1.        , 1.28595756, 0.11800314]), array([1.        , 0.15400746, 0.10614905]), array([1.        , 0.37662845, 0.62762752]), array([1.        , 0.39153302, 0.85825651]), array([1.        , 0.22961855, 0.97345687]), array([0.919194  , 0.51147062, 0.31161369]), array([1.        , 0.16485434, 1.612375  ]), array([ 1.        , -1.82300731, 10.40295209]), array([1.        , 0.79982124, 0.13813832]), array([ 1.        , -1.12016666,  1.68834769]), array([ 1.        ,  2.34236382, -1.50970187]), array([1.        , 9.10340275, 3.71322731]), array([1.        , 2.54908323, 0.22384165]), array([1.        , 0.49306911, 0.73197816]), array([ 1.        , -1.82045545,  4.76760331]), array([1.        , 0.31049558, 0.61819234]), array([0.94683902, 0.86094109, 0.94664849]), array([ 1.        , -3.27345449,  2.74169895]), array([1.        , 0.6902818 , 0.48361414]), array([0.96996897, 0.55330264, 0.93622195]), array([ 1.        , -0.53897415,  0.21931723]), array([0.88245487, 0.77304405, 0.60258983]), array([1.        , 0.28022379, 1.24357567]), array([1.        , 0.63866034, 0.88029031]), array([ 1.        , -3.84333837,  1.82005827]), array([ 1.        ,  0.79980339, -0.53478872]), array([ 1.        ,  0.6526323 , -1.04319818]), array([1.        , 1.07551494, 0.21342046]), array([1.        , 0.98392978, 0.59268443]), array([1.        , 1.06892063, 1.04082519]), array([1.        , 0.81343528, 0.52855969]), array([1.        , 0.93031168, 0.95481775]), array([1.        , 1.26816887, 2.36662513]), array([0.83400058, 0.44777878, 0.82768249]), array([ 1.        ,  3.22857543, -0.33124528]), array([ 1.        ,  1.84641352, -0.83754711]), array([ 1.        ,  0.77613277, -2.10859522]), array([1.01997402, 0.85912428, 0.94258184]), array([ 1.        , -0.37357541,  0.66933295]), array([0.75260074, 0.83673931, 0.68805963]), array([1.        , 0.66521817, 0.77102881]), array([1.        , 2.29138775, 3.50448125]), array([ 1.        ,  0.94132608, -1.71975888]), array([1.        , 1.32994921, 4.38371575]), array([ 1.        , -0.70328679,  1.4826183 ]), array([1.        , 0.96594289, 0.24863242]), array([0.82444318, 0.87025878, 0.66350225]), array([1.        , 1.08488315, 0.59182652]), array([1.        , 1.45641741, 0.72403879]), array([0.6770886 , 0.77365766, 0.59345796]), array([ 1.        , -0.47103988,  1.06035932]), array([ 1.        , -1.20690118, -4.12905294]), array([ 1.        ,  0.2840387 , -1.53815773]), array([ 1.        , -0.79358435, -0.07637077]), array([1.        , 0.01672126, 0.25033246]), array([0.91217923, 0.4137023 , 0.97765307]), array([ 1.        , -0.23821692,  0.53660958]), array([1.        , 0.02772281, 0.49641483]), array([ 1.        ,  8.35729864, -3.48247831]), array([ 1.        , -0.53556384,  6.27892036]), array([ 1.       , -0.093962 ,  3.0981085]), array([ 1.        , -0.71269785, -0.13478057]), array([ 1.        ,  0.12040887, -0.85350748]), array([ 1.        , -0.80053258,  1.16586056]), array([1.        , 0.88229512, 0.28816744]), array([1.        , 0.14236638, 0.31404425]), array([1.        , 0.1640387 , 0.72974499]), array([ 1.        , -0.95949293,  1.76803777]), array([1.        , 0.20304265, 1.06397229]), array([ 1.        , -0.36708757,  1.10706536]), array([ 1.        ,  0.19800382, -0.38229245]), array([1.        , 0.24491226, 0.29143052]), array([ 1.        ,  0.47902794, -0.13305322]), array([ 1.        ,  1.94639763, -0.68206447]), array([ 1.        , -0.3175356 , -0.29195597]), array([1.        , 0.12022753, 0.32707967]), array([1.        , 0.26800002, 0.93421293]), array([1.        , 0.08697559, 1.12119195]), array([ 1.        ,  0.03263929, -0.46009295]), array([ 1.        ,  0.59634987, -0.21149713]), array([ 1.        ,  0.50062591, -0.19634241]), array([ 1.        ,  0.01278344, -0.11661343]), array([1.        , 0.9201459 , 0.65362312]), array([1.        , 0.08614797, 0.67725558]), array([1.        , 0.61135639, 0.69386625]), array([1.        , 0.40117857, 1.0500862 ]), array([1.        , 0.49731047, 1.02272705]), array([1.        , 0.88596972, 2.26464099]), array([ 1.        ,  0.33415252, -0.24684369]), array([ 1.        , -1.38067454,  1.86383781]), array([1.        , 0.26194653, 0.50047253]), array([1.        , 0.5995053 , 0.64893891]), array([1.        , 0.43259124, 0.73471377]), array([1.        , 1.2375987 , 0.39655221]), array([1.        , 0.56120348, 0.61724594]), array([1.        , 0.64076708, 1.13022574]), array([0.80575861, 0.49911816, 0.87423345]), array([ 1.        ,  1.02671858, -0.27878861]), array([ 1.        ,  0.49581467, -0.21109829]), array([1.        , 0.92328443, 0.11059398]), array([1.        , 0.82275452, 0.19531227]), array([1.        , 0.24441813, 0.93095305]), array([1.03216695, 0.89004669, 0.9446258 ]), array([1.        , 0.75222592, 0.26327567]), array([ 1.        , -0.16196844, -4.34032194]), array([1.        , 0.82998381, 0.89746708]), array([ 1.        ,  4.34267421, -3.16126361]), array([1.        , 0.70553699, 0.17207083]), array([1.        , 0.62034911, 0.37816498]), array([1.        , 0.93039962, 0.18334241]), array([1.        , 0.86525672, 0.48666724]), array([1.        , 0.83335722, 0.33172702]), array([1.        , 0.45363415, 0.85436412]), array([ 1.        ,  1.20538937, -1.04855728]), array([1.        , 2.33049181, 0.9684172 ]), array([ 1.        ,  1.16381697, -0.35618121]), array([ 1.        ,  3.31516325, -1.01575505]), array([1.        , 0.83883147, 0.95756136]), array([1.        , 0.99677852, 0.30400187]), array([1.        , 0.95421358, 0.50691401]), array([1.        , 0.93811134, 1.15952815]), array([1.        , 0.99915833, 1.22641675]), array([ 1.        , -7.36830582,  1.42659845]), array([1.        , 0.89180829, 1.64850777]), array([1.        , 1.18293976, 0.15003855]), array([ 1.        , -2.60566995,  0.02793725]), array([1.        , 1.04067967, 0.14192398]), array([ 1.        , -0.15927697,  0.00472735]), array([1.        , 1.35252723, 0.28975448]), array([1.        , 1.21468183, 0.3106135 ]), array([ 1.        ,  0.11202471, -1.88325826]), array([  1.        ,  16.43375282, -16.19686905]), array([ 1.        , -0.9443754 ,  0.50849825]), array([ 1.        , -0.10954587,  0.11287822]), array([ 1.        , -0.12664118, -0.27004432]), array([1.        , 0.65641955, 1.5096081 ]), array([1.        , 0.82877105, 0.73060139]), array([ 1.        , -0.80292994,  0.46300847]), array([ 1.        , -0.14602384,  0.50713979]), array([ 1.        , -0.08199278,  0.63939616]), array([1.        , 0.06791651, 1.24448214]), array([ 1.        , -0.01776065,  1.09068738]), array([ 1.        , -0.33688196, -0.54484718]), array([ 1.        ,  0.08409397, -0.04291118]), array([ 1.        , -0.09113391,  0.37355036]), array([1.        , 0.18534693, 0.49056618]), array([1.        , 0.47056409, 0.20293788]), array([ 1.        , 13.62817294,  7.47928389]), array([1.        , 0.05874164, 0.80749685]), array([ 1.        , -0.05456613,  0.72512333]), array([1.        , 0.13266782, 1.13946703]), array([ 1.        ,  0.15020373, -0.05299898]), array([ 1.        ,  0.09991   , -0.02762147]), array([1.        , 0.31149543, 0.11669568]), array([ 1.       ,  0.0094413, -0.0071363]), array([1.        , 4.83222839, 6.5601977 ]), array([ 1.        ,  0.70433169, -1.39592365]), array([1.        , 0.28496176, 0.67604021]), array([1.        , 0.35964343, 0.93390543]), array([1.        , 0.80846834, 0.93923348]), array([ 1.        ,  0.37260929, -0.02046959]), array([1.        , 0.23332076, 0.05187303]), array([1.        , 0.27403971, 0.35774509]), array([ 1.        ,  0.30694296, -0.12473707]), array([  1.        ,  27.40506844, -17.4551299 ]), array([ 1.        , -0.55939522, -0.5923273 ]), array([1.        , 0.20593409, 0.81389289]), array([1.        , 0.12858767, 0.03771941]), array([1.        , 0.87629072, 1.91543732]), array([1.        , 0.38847808, 0.06818032]), array([1.        , 0.66835994, 0.18227478]), array([1.        , 0.51801425, 0.31855003]), array([1.        , 0.46222475, 0.49437187]), array([1.        , 0.25041695, 0.17976122]), array([1.        , 0.32982702, 0.87998464]), array([ 1.        , -0.06686374,  0.61724081]), array([1.        , 0.50461473, 0.83849029]), array([1.        , 0.65352429, 2.45746197]), array([1.        , 0.70750313, 0.05060491]), array([1.        , 0.67836154, 0.19900201]), array([1.        , 0.6925537 , 0.23890198]), array([1.        , 0.62292669, 0.21226513]), array([1.        , 0.97215078, 0.653731  ]), array([1.        , 0.49137341, 0.75789655]), array([1.        , 0.66774957, 0.91345203]), array([1.        , 0.66271486, 0.79321049]), array([ 1.        , -0.24759136, -0.41724736]), array([ 1.        ,  1.04539408, -0.05277561]), array([1.        , 0.72714496, 0.13061322]), array([1.        , 1.52449939, 0.62727133]), array([1.        , 0.63012851, 0.51110863]), array([1.        , 1.11397365, 0.8128457 ]), array([1.        , 0.43611077, 0.04900406]), array([1.        , 3.0013913 , 0.51568195]), array([1.        , 0.61233944, 0.8686142 ]), array([1.        , 1.00641195, 1.54409732]), array([1.        , 0.91025141, 0.21422958]), array([1.        , 0.67888341, 0.80239038]), array([ 1.        ,  0.24743353, -0.27716837]), array([1.        , 0.85883719, 0.07075055]), array([1.        , 0.95775891, 0.45784473]), array([1.        , 1.089     , 0.50940732]), array([1.        , 1.30829869, 0.65809112]), array([1.        , 1.33668853, 1.82431572]), array([1.        , 0.29061245, 4.49643612]), array([ 1.        , -0.09373278,  0.28601731]), array([1.        , 1.32283511, 0.34777087]), array([1.        , 1.12062685, 0.00968835]), array([ 1.        ,  1.46952912, -0.08082135]), array([1.        , 0.90701991, 0.58334666]), array([1.        , 1.16588183, 0.78518891]), array([1.        , 1.28620966, 1.11404196]), array([1.        , 0.97680565, 0.88297258]), array([1.        , 1.81684071, 1.39496756]), array([ 1.        ,  0.62020219, -0.11077743]), array([1.        , 0.09008872, 0.13474704]), array([ 1.        ,  0.10297272, -0.45291982]), array([ 1.        , -0.19587075,  0.57392508]), array([1.        , 0.08402751, 0.64909212]), array([ 1.        , -0.00604381,  0.55431528]), array([1.        , 0.00678169, 0.71049817]), array([ 1.        , -0.12189012,  0.83503862]), array([1.        , 0.02096501, 1.4349575 ]), array([ 1.        ,  0.19818837, -0.04377283]), array([1.        , 0.17718697, 0.20027244]), array([1.        , 1.33355763, 1.01239874]), array([ 1.        , -0.00969482,  0.48298044]), array([1.        , 0.20536103, 0.69435753]), array([1.        , 0.16104   , 0.58717249]), array([1.        , 0.27063805, 0.62859693]), array([1.        , 0.13910544, 0.80218818]), array([1.        , 0.21106361, 0.95325174]), array([ 1.        ,  0.2818838 , -0.01316302]), array([ 1.        ,  0.20713661, -0.22885909]), array([1.        , 0.14434386, 0.31694914]), array([1.        , 0.32992033, 0.39284669]), array([ 1.        ,  0.29728822, -0.77845311]), array([1.        , 0.19044058, 0.66238257]), array([1.        , 0.28507025, 0.76974497]), array([1.        , 0.11841125, 2.108702  ]), array([1.        , 0.27685947, 0.91327117]), array([1.        , 0.43166823, 0.10595015]), array([1.        , 0.31792065, 0.17558503]), array([1.        , 0.37955842, 0.64154664]), array([1.        , 0.24139089, 0.47866976]), array([1.        , 0.37265624, 0.47952719]), array([1.        , 0.39455529, 0.65707049]), array([1.        , 0.41984354, 0.74043343]), array([1.        , 0.40394661, 1.12195816]), array([1.        , 0.4143266 , 1.03478107]), array([1.        , 0.49396743, 0.04146449]), array([1.        , 0.56804165, 0.12476679]), array([1.        , 0.58402245, 0.32426141]), array([ 1.        , -0.92780631, -0.27081202]), array([1.        , 0.54207803, 0.45928736]), array([1.        , 0.54184103, 0.59588353]), array([1.        , 0.44809565, 0.76729203]), array([1.        , 0.42410148, 0.86341659]), array([1.        , 0.58122338, 0.9497873 ]), array([ 1.        ,  0.85257563, -0.25582361]), array([1.        , 0.70607437, 0.12177651]), array([1.        , 0.64473004, 0.27443283]), array([ 1.        , -0.17787872, -0.12706219]), array([1.        , 0.56552483, 0.76385078]), array([1.        , 0.91470325, 0.8402264 ]), array([1.        , 0.6656163 , 0.61453351]), array([1.        , 0.61855281, 0.82176135]), array([1.        , 0.55963919, 1.05674249]), array([1.        , 0.73502021, 0.26962501]), array([1.        , 0.83966039, 0.06247229]), array([1.        , 0.82611518, 0.16484015]), array([1.        , 0.66521736, 0.26529951]), array([1.        , 0.91559254, 0.40799108]), array([1.        , 0.70247953, 0.66505996]), array([1.        , 1.08666993, 1.28563306]), array([1.        , 0.4897334 , 0.88356841]), array([1.        , 0.67230454, 1.00672895]), array([ 1.        ,  0.88984262, -0.03611426]), array([1.        , 0.80422433, 0.12782708]), array([1.        , 0.79414705, 0.24544785]), array([1.        , 0.82192551, 0.40145418]), array([1.        , 0.77468726, 0.66133906]), array([1.        , 1.01895937, 0.66605158]), array([1.        , 0.84293729, 0.82278286]), array([ 1.        , 62.82449528, 81.58137106]), array([1.        , 0.7988689 , 1.01552192]), array([1.        , 1.11255807, 0.11219657]), array([1.        , 1.16575551, 0.06620174]), array([1.        , 0.98198565, 0.19429924]), array([1.        , 1.48621341, 0.14125921]), array([1.        , 0.93910379, 0.46380158]), array([1.        , 0.91825496, 0.55092666]), array([1.        , 0.96721374, 0.77747856]), array([1.        , 0.90030759, 0.94339069]), array([1.        , 1.60175698, 1.7170953 ]), array([ 0.        ,  0.12830536, -0.09953319]), array([ 0.        , -0.41096358, -0.33624695]), array([0.        , 0.06729948, 0.23606361]), array([ 0.        , -0.05415698,  0.43206533]), array([0.        , 0.02191379, 0.86676601]), array([ 0.        , -0.04115348,  0.81864734]), array([0.        , 0.06892245, 0.77361153]), array([ 0.        , -0.01337675,  0.75893442]), array([0.        , 0.03932738, 0.93687543]), array([0.        , 0.04964408, 0.05963497]), array([ 0.        ,  0.2780632 , -0.02104741]), array([0.        , 0.12313576, 0.2520967 ]), array([0.        , 0.15033736, 0.36748477]), array([0.        , 0.16492089, 0.58960964]), array([0.        , 0.19082998, 0.56854022]), array([0.        , 0.1362526 , 0.76137487]), array([0.        , 0.17245647, 0.81931074]), array([0.        , 0.13020036, 0.96794738]), array([0.        , 0.15262511, 0.92903123]), array([ 0.        ,  0.10713306, -0.64790852]), array([0.        , 1.64987882, 0.6516437 ]), array([0.        , 0.37358702, 0.4237786 ]), array([0.        , 0.26114087, 0.56660075]), array([0.        , 0.8135959 , 0.08495462]), array([0.        , 0.23583845, 0.8076958 ]), array([0.        , 0.27641162, 0.80722145]), array([0.        , 0.27987535, 1.05988276]), array([0.        , 0.35626865, 0.08890528]), array([0.        , 0.38432082, 0.22109281]), array([0.        , 0.37026549, 0.23048077]), array([0.        , 0.5090643 , 0.47208785]), array([0.        , 0.44322528, 0.47265261]), array([0.        , 0.35938187, 0.55349312]), array([0.        , 0.33778656, 0.81668405]), array([0.        , 0.44075387, 0.80548556]), array([0.       , 0.3446607, 1.1053925]), array([ 0.        ,  0.07737238, -0.05011667]), array([0.        , 0.54501418, 0.2078185 ]), array([ 0.        ,  0.06682054, -0.55234704]), array([0.        , 0.57137617, 0.42508261]), array([0.        , 0.49657497, 0.6826064 ]), array([ 0.        , -0.36220961,  0.28244025]), array([ 0.        , -0.19201614,  0.8303995 ]), array([0.        , 0.59439877, 0.74167824]), array([0.        , 0.3293781 , 1.21993527]), array([ 0.        ,  0.93387014, -0.18724169]), array([0.        , 0.12981197, 1.12527649]), array([0.        , 0.64138595, 0.31840592]), array([0.        , 0.58424979, 0.42643351]), array([0.        , 0.57541866, 0.50246129]), array([0.        , 0.62993176, 0.84096911]), array([0.        , 0.53927595, 0.73344033]), array([0.        , 0.69786927, 0.81842632]), array([0.        , 0.45798837, 0.97415175]), array([0.        , 0.47624918, 0.2980305 ]), array([0.        , 0.19543492, 0.7272091 ]), array([0.        , 0.71278534, 0.2239612 ]), array([0.        , 0.80729493, 0.42542309]), array([0.        , 0.81286146, 0.45844155]), array([0.        , 0.95922799, 0.8636447 ]), array([0.        , 0.66922378, 0.74958742]), array([0.        , 0.75095975, 0.80291247]), array([0.        , 0.67688634, 1.06456092]), array([0.        , 0.88018508, 0.02346799]), array([0.        , 0.90942599, 0.04741806]), array([0.        , 0.75853019, 0.32343396]), array([0.        , 0.90050974, 0.45545321]), array([0.        , 0.80912465, 0.48530558]), array([0.        , 0.80234434, 0.57851834]), array([0.        , 0.8946941 , 0.78658529]), array([0.        , 0.82110084, 0.80843533]), array([0.        , 0.85013996, 0.94836891]), array([    0.        , -9406.14666146,  6930.2686875 ]), array([0.        , 1.06448683, 0.14530637]), array([0.        , 1.00521955, 0.34367434]), array([0.        , 1.03728166, 0.31366539]), array([0.        , 0.87218988, 0.41142337]), array([0.        , 0.98155148, 0.66236403]), array([0.        , 1.01497272, 0.73509308]), array([0.        , 1.51871848, 1.92076667]), array([0.        , 1.18789446, 1.05825363]), array([0.        , 0.07780869, 0.14383965]), array([ 0.        , -0.20559986,  0.01070727]), array([ 0.        , -0.30479781,  0.06361664]), array([ 0.        , -0.05755944,  0.36356958]), array([0.        , 1.43819781, 0.45857413]), array([ 0.        , -0.26868681,  0.61707841]), array([ 0.        , -0.91220023,  0.86664648]), array([ 0.        , -0.08013568,  0.99927304]), array([0.        , 0.0996173 , 1.07029914]), array([0.        , 0.13594844, 1.82858274]), array([0.        , 0.18522344, 0.13660865]), array([0.        , 0.15238467, 0.4629713 ]), array([0.        , 0.82537616, 1.14939972]), array([0.        , 4.62021528, 0.02471231]), array([0.        , 0.0810756 , 0.74333389]), array([0.        , 0.16696908, 0.90356473]), array([0.        , 0.09867938, 0.87047253]), array([ 0.        ,  1.04956996, -0.41078706]), array([0.        , 0.07274329, 0.1035677 ]), array([0.        , 0.36879585, 0.19324336]), array([0.        , 0.19676357, 0.24906521]), array([0.        , 0.21820182, 0.4051148 ]), array([0.        , 0.70052353, 0.84968458]), array([0.       , 0.2795497, 0.5399668]), array([0.        , 0.02622494, 0.04275418]), array([0.        , 0.38486788, 0.80374606]), array([0.        , 0.33628721, 0.90532262]), array([ 0.        ,  3.98222483, -3.46757098]), array([ 0.        ,  0.48610601, -0.17097271]), array([0.        , 0.34607198, 0.13651378]), array([0.        , 0.45175978, 0.47640911]), array([0.        , 0.21910272, 0.74238901]), array([0.        , 0.54602922, 0.93807633]), array([0.        , 0.71329071, 0.29280055]), array([ 0.        , -0.82215328,  1.97862948]), array([-0.27921757,  0.78353249,  0.16841801]), array([0.        , 0.60469691, 0.12522374]), array([ 0.        ,  0.43761286, -0.03724752]), array([ 0.        , -0.31949219, -1.76903322]), array([0.        , 0.10186954, 0.14740422]), array([0.        , 0.5492349 , 0.80160035]), array([0.        , 0.62926745, 0.47868322]), array([0.        , 0.70594803, 0.7206821 ]), array([0.        , 6.49737254, 2.20011279]), array([0.        , 0.45529398, 1.51057043]), array([0.        , 0.03454971, 0.66928587]), array([0.        , 0.82599889, 0.30734386]), array([0.        , 0.89693947, 0.44702869]), array([0.        , 0.86265804, 0.42820086]), array([0.        , 0.73959915, 0.47082543]), array([0.        , 1.22939609, 0.16670812]), array([-0.10072835,  0.79606496,  0.5478966 ]), array([0.        , 0.21952206, 1.15107857]), array([ 0.        ,  0.59861272, -0.60839365]), array([0.        , 0.69722845, 0.06082038]), array([ 0.        ,  0.61601701, -0.08618762]), array([ 0.        ,  1.56457497, -0.41469107]), array([0.        , 0.5577464 , 1.43067255]), array([0.        , 1.00257232, 0.72265297]), array([0.        , 0.67757803, 1.00144745]), array([0.        , 0.85512047, 0.81169251]), array([0.        , 0.82608051, 0.77586723]), array([0.        , 0.91271338, 1.06796445]), array([0.        , 0.79594975, 0.06531152]), array([0.        , 0.94419215, 0.05935575]), array([0.        , 0.76996071, 0.18474545]), array([0.        , 0.59948374, 0.31764659]), array([0.        , 0.98723005, 0.55101204]), array([0.        , 1.67716495, 0.65422732]), array([0.        , 0.64750192, 0.50441988]), array([0.        , 0.8096446 , 0.94518441]), array([ 0.        , -2.18520701, -2.90302162]), array([0.        , 0.60462045, 1.18130377]), array([0.        , 0.94029754, 0.12526011]), array([0.        , 1.26017643, 0.05459041]), array([ 0.        ,  0.12567279, -1.20131924]), array([0.        , 5.38155642, 2.90139055]), array([0.        , 1.81916431, 1.00485958]), array([0.        , 0.62582697, 0.58570591]), array([0.        , 0.87443743, 1.06707788]), array([0.        , 1.11040934, 1.57527969]), array([ 0.        ,  8.27610673, 12.33542616]), array([ 0.        , -0.15709607,  0.32208037]), array([0.        , 5.54997903, 4.64605506]), array([0.00000000e+00, 9.51287876e-05, 1.78662038e-01]), array([0.        , 0.87305219, 0.48048413]), array([0.        , 0.06851768, 0.64956622]), array([0.        , 0.16467049, 0.94197746]), array([0.        , 0.02949644, 1.3952517 ]), array([ 0.        , -0.2280924 ,  0.94576766]), array([ 0.        , -1.27093951, -1.50471276]), array([ 0.        , -0.08286546,  0.85105403]), array([ 0.        ,  0.19270029, -0.5504626 ]), array([ 0.        , -0.67154084, -0.37145446]), array([0.0593353 , 0.10490475, 0.22920898]), array([0.        , 0.50909322, 0.07378165]), array([ 0.        , -0.16359328,  0.70056407]), array([0.        , 0.1329742 , 1.03126235]), array([ 0.        ,  2.55382668, -2.17966584]), array([0.08598842, 0.39346521, 0.08448326]), array([0.        , 2.40114267, 1.58679544]), array([0.        , 0.07071541, 0.41376715]), array([0.        , 0.41833139, 0.51742962]), array([0.19505937, 0.15172246, 0.30075209]), array([0.        , 0.07357841, 0.93167207]), array([0.        , 0.26137159, 0.13693405]), array([0.        , 0.9225753 , 0.81741897]), array([0.        , 0.42033779, 1.32909094]), array([ 0.00000000e+00,  8.43771448e-01, -6.19762372e-04]), array([ 0.        , -0.64435605, -2.11045197]), array([ 0.        , -1.8983441 ,  0.93592239]), array([ 0.        , -0.26757336,  1.49983107]), array([0.        , 0.37009109, 0.69920518]), array([0.        , 0.45044649, 1.05489024]), array([0.        , 0.82962317, 1.27954904]), array([0.        , 0.3904173 , 1.30119923]), array([ 0.        , -0.61989708,  2.65814953]), array([-0.02375862,  0.18532219,  0.55963591]), array([ 0.        ,  0.72550031, -0.17102896]), array([ 0.        ,  0.68005202, -7.05343808]), array([0.098053  , 0.29168353, 0.26891756]), array([0.        , 0.58495139, 0.39476481]), array([0.       , 0.3125152, 0.6099895]), array([0.        , 0.29139609, 0.64048475]), array([0.        , 1.44470893, 0.79256198]), array([ 0.        , -0.48778841, -0.83034906]), array([0.        , 0.27994561, 2.28462589]), array([ 0.        , 28.72374106, 65.52917238]), array([-0.18630522,  0.4097662 ,  0.24127547]), array([ 0.        ,  0.79908168, -0.4506787 ]), array([0.        , 0.37877668, 0.88376403]), array([0.        , 0.69387904, 0.6725541 ]), array([ 0.        , -2.88878481, -0.78916577]), array([0.        , 1.35365273, 1.72174687]), array([0.        , 0.78348919, 0.99494611]), array([0.        , 0.62571535, 0.15216845]), array([-0.03095481,  0.60131965, -0.12904822]), array([0.        , 0.92136014, 0.31730398]), array([-0.06522237,  0.3226756 ,  0.07985568]), array([ 0.        , -8.07126639, -4.9753556 ]), array([0.        , 1.38471922, 1.13338105]), array([0.        , 0.82617854, 0.83954301]), array([ 0.        ,  1.00854659, -0.18553711]), array([0.        , 0.78943805, 1.61798939]), array([ 0.        , -1.07031249, -0.14674595]), array([0.        , 0.67074118, 0.29053292]), array([0.        , 0.92444179, 0.45361088]), array([0.        , 0.47548657, 0.67906077]), array([ 0.        ,  0.02438032, -0.45810894]), array([ 0.        , -2.44292269,  1.66422915]), array([0.        , 0.83122347, 1.53127851]), array([0.        , 0.70809774, 0.82916868]), array([ 0.        , -1.63665621, -0.59586342]), array([  0.        , -31.29673512,  17.33550247]), array([ 0.        ,  1.62958443, -0.35211973]), array([ 0.        ,  3.38621189, -2.2699256 ]), array([0.        , 1.21409735, 0.13891688]), array([0.        , 2.21249501, 1.14937183]), array([0.        , 1.01682008, 0.61368372]), array([0.        , 0.94882854, 1.07850488]), array([ 0.        ,  0.52706427, -0.03112432]), array([ 0.        , -0.30704974,  0.33310266]), array([0.        , 1.47998196, 0.62795295]), array([ 0.        , -1.98978425, -1.19092383]), array([0.        , 0.22745006, 0.18479442]), array([  0.        , -40.14551608,   1.54365074]), array([ 0.        , -1.09572898,  1.09869479]), array([0.        , 1.25324867, 0.64910034]), array([ 0.        ,  1.76244446, -2.27668674]), array([0.        , 0.15304082, 0.06617582]), array([ 0.        , -0.00878751,  0.09414556]), array([0.3179774 , 0.28514534, 0.30832918]), array([ 0.        , -0.25771636, -0.32874259]), array([ 0.        , -5.28751053, -5.80382673]), array([ 0.        , -0.30738635,  0.23287062]), array([0.        , 0.29593116, 0.32741504]), array([ 0.        , -0.13083407,  0.85790663]), array([0.05817341, 0.36693977, 0.18704851]), array([0.        , 1.79142011, 0.23692881]), array([ 0.        , -1.42187272,  0.9218311 ]), array([ 0.        ,  0.54543411, -0.04446121]), array([ 0.        , -0.02202612,  0.40738711]), array([ 0.        ,  0.78926386, -0.92621498]), array([ 0.        ,  0.60790332, -0.07059367]), array([ 0.        ,  0.56980046, -2.03138483]), array([0.        , 0.94556605, 0.1856859 ]), array([ 0.        , -0.05660954,  0.48155557]), array([ 0.        , -0.02700181,  0.6705982 ]), array([ 0.        ,  3.75217128, -1.58449148]), array([ 0.        , -0.04214712,  1.26391875]), array([ 0.        , -0.03557999,  1.14354192]), array([ 0.        ,  0.44009362, -1.36646638]), array([0.11135277, 0.60263851, 0.19827916]), array([0.        , 0.77750423, 0.44390301]), array([0.        , 0.18035131, 1.02711334]), array([0.09800181, 0.15617734, 0.35071434]), array([0.        , 0.02811574, 0.53730008]), array([0.        , 1.25545455, 0.39009479]), array([ 0.        ,  0.52648656, -0.40410045]), array([ 0.        ,  0.53041342, -0.19715435]), array([ 0.        , -0.84114775,  0.29252875]), array([0.        , 0.80352352, 1.02628768]), array([0.21922127, 0.3146963 , 0.31727215]), array([ 0.04464015, -0.00554713,  0.41326119]), array([0.09062479, 0.4110774 , 0.1188202 ]), array([0.        , 0.67526073, 0.91222922]), array([0.        , 0.77793672, 1.55036962]), array([ 0.        , -0.25398566,  2.41116739]), array([ 0.        ,  1.22850377, -0.12481861]), array([ 0.        , -0.46093742,  3.81047843]), array([ 0.        , 32.23296598, 14.93926275]), array([0.31306986, 0.44713999, 0.29601525]), array([0.        , 0.83212174, 0.99759283]), array([-0.20378548,  0.39196577,  0.23410988]), array([ 0.        , -0.12955573,  0.52297067]), array([0.05491266, 0.20817114, 0.21545538]), array([0.        , 0.35294615, 1.00758281]), array([ 0.        ,  0.71685674, -0.44223481]), array([0.        , 0.64500343, 0.9730812 ]), array([ 0.        ,  0.95935097, -0.13066407]), array([0.18951184, 0.51243843, 0.1468736 ]), array([0.        , 0.68087238, 0.84302386]), array([0.        , 2.52486046, 2.35946713]), array([0.23543202, 0.20806592, 0.53952551]), array([0.        , 0.70442836, 1.60380024]), array([ 0.        ,  0.51510585, -1.29280381]), array([0.        , 0.76714747, 0.17168791]), array([ 0.        ,  5.23977931, -2.27798983]), array([0.        , 1.18110504, 0.2588069 ]), array([ 0.        , -0.93562399, -1.31399844]), array([0.        , 1.32605522, 0.59563901]), array([0.        , 0.71915789, 2.8547688 ]), array([0.        , 1.46514498, 0.69745065]), array([0.        , 0.53393841, 0.74168144]), array([ 0.        ,  2.09880924, -1.18445263]), array([ 0.        , -0.77045341, -0.7088713 ]), array([ 0.        ,  2.27180881, -0.75438344]), array([0.        , 1.73142853, 0.0979806 ]), array([0.        , 3.17549013, 0.43137763]), array([ 0.        ,  1.09438012, -0.08886009]), array([0.        , 0.8397169 , 0.35839506]), array([0.        , 0.2045002 , 0.87893103]), array([0.        , 5.34184328, 1.98434509]), array([0.16614395, 0.18393731, 0.11165868]), array([ 0.        , -2.08808808,  0.54861427]), array([0.25350672, 0.18564556, 0.24798562]), array([ 0.        , -1.71332793, -2.0528603 ]), array([0.        , 1.80873674, 1.3346454 ]), array([0.19057969, 0.09207171, 0.40253736]), array([0.03590636, 0.00072175, 0.16407759]), array([ 0.        , -0.63880138,  0.37686693]), array([-0.08353888,  0.5087962 ,  0.53608296]), array([ 0.        , -3.72854402, -1.10372924]), array([0.        , 0.24563837, 0.63494316]), array([0.51617917, 0.34380259, 0.37318062]), array([0.        , 0.02257028, 0.6093197 ]), array([ 0.        , -0.99390673,  0.26528651]), array([ 0.        ,  4.32958254, -3.92420005]), array([0.        , 0.34490309, 1.86138979]), array([0.        , 0.28297579, 0.65420201]), array([ 0.        , -1.53474526,  1.53939207]), array([ 0.        ,  1.06012048, -3.19134796]), array([ 0.        , -1.87914308, -2.50614863]), array([0.        , 0.20420311, 0.33446636]), array([ 0.        , -0.35086413,  0.19479291]), array([0.2725658 , 0.35314892, 0.33702613]), array([0.        , 0.72442906, 0.35692822]), array([0.        , 1.18730972, 0.56271357]), array([0.        , 0.75174226, 1.82265823]), array([ 0.        , -0.32483092,  2.26504712]), array([ 0.        ,  2.14205017, -3.6526283 ]), array([ 0.        ,  5.93313499, -3.23277613]), array([ 0.        ,  0.7756118 , -0.01421871]), array([0.        , 0.98632966, 2.90404827]), array([0.        , 0.62889846, 0.5158515 ]), array([0.        , 0.02798874, 1.13463378]), array([0.        , 0.86418842, 0.38521098]), array([0.        , 0.42959385, 1.28009682]), array([0.19243811, 0.30028798, 0.30893903]), array([0.        , 1.08444004, 0.10449902]), array([  0.        , -16.18445727,  -1.86010633]), array([0.41298968, 0.36132399, 0.46743573]), array([0.        , 0.3315283 , 0.71087778]), array([ 0.        ,  0.91534822, -1.22834074]), array([0.        , 0.54184147, 1.58550308]), array([0.22456689, 0.12293247, 0.48712761]), array([0.        , 0.99357356, 1.12446601]), array([0.20422303, 0.30662782, 0.3126057 ]), array([ 0.        ,  1.67724244, -1.05350362]), array([ 0.        ,  7.47002008, -8.78722926]), array([0.        , 3.07856213, 1.50922954]), array([0.        , 1.06365471, 1.69854169]), array([ 0.        ,  0.5189111 , -1.87788368]), array([0.19265722, 0.14384691, 0.28682264]), array([ 0.        ,  0.94711691, -1.11261708]), array([0.        , 2.88660515, 0.7312341 ]), array([0.        , 0.91967425, 1.7246853 ]), array([0.22768409, 0.14735723, 0.50216055]), array([0.49957973, 0.85817423, 0.37150804]), array([0.35220751, 0.53348872, 0.22795779]), array([ 0.        , -0.32032933, -0.16068052]), array([0.        , 0.86377117, 0.19875183]), array([0.31992587, 0.47104299, 0.34871521]), array([  0.        , -13.35394697,  23.94925385]), array([0.        , 1.41153063, 2.01147509]), array([0.        , 2.48008056, 2.22039658]), array([ 0.        , -1.22205927,  0.1604795 ]), array([0.        , 2.31274933, 0.22974968]), array([0.20821199, 0.31062811, 0.35745157]), array([0.        , 1.31970225, 0.49734874]), array([0.17649678, 0.87046224, 0.51007823]), array([0.        , 0.55018106, 1.12981583]), array([0.05122913, 0.04709699, 0.43257063]), array([0.        , 0.55453223, 1.04237738]), array([0.        , 2.39931402, 2.04894882]), array([0.        , 0.73087268, 1.79977308]), array([0.        , 0.24956753, 0.63297544]), array([ 0.        ,  2.6181401 , -0.28425768]), array([ 0.        ,  2.64019739, -0.15829563]), array([0.        , 3.3707208 , 2.75103673]), array([0.20022609, 0.87459583, 0.53047145]), array([0.        , 1.80506467, 1.33420795]), array([ 0.        ,  0.5017826 , -0.36100509]), array([-0.2215985 ,  0.26263451,  0.43225892]), array([ 0.        , -2.78578734, -2.10079167]), array([0.40911986, 0.4754702 , 0.56845647]), array([0.2938118 , 0.15432612, 0.05552276]), array([0.26263006, 0.46658025, 0.41502902]), array([0.        , 1.91559642, 0.47126305]), array([0.34204998, 0.49112006, 0.18348835]), array([ 0.        , -0.93120219,  1.1756058 ]), array([ 0.        , -1.03794888,  1.84081574]), array([0.        , 0.66274167, 0.15513486]), array([0.56983841, 0.40633329, 0.44614882]), array([0.        , 1.11646671, 3.93788509]), array([  0.        , -11.67875195,  -7.80127564]), array([ 0.        ,  1.03619237, -1.78997599]), array([ 0.        ,  0.36573503, -0.1022488 ]), array([0.        , 0.01025351, 0.80043573]), array([ 0.        ,  4.29649984, -0.34786995]), array([0.        , 1.14799581, 0.25957672]), array([ 0.        ,  0.3452989 , -0.30680439]), array([0.        , 0.0342448 , 1.34253343]), array([ 0.        , -3.16616876, -0.65307618]), array([0.51418265, 0.54672391, 0.20565629]), array([ 0.        , -0.14667757, -0.15123658]), array([ 0.        , -0.12754687,  0.85738506]), array([ 0.        , -8.63461683, -2.32051221]), array([ 0.        ,  1.96426933, -1.41821182]), array([0.12030219, 0.15336549, 0.55130806]), array([ 0.1563442 ,  0.69001342, -0.0504883 ]), array([0.49992355, 0.61537299, 0.45075679]), array([0.        , 0.99749844, 0.79704789]), array([ 0.        , -1.01234902,  0.61781791]), array([0.        , 0.80563228, 0.89125693]), array([0.07208074, 0.39669978, 0.28113138]), array([0.        , 0.45608097, 0.21545577]), array([0.37774595, 0.43891863, 0.55841596]), array([0.        , 0.16203813, 2.16653321]), array([0.40162917, 0.50297025, 0.60938616]), array([  0.        ,   4.06727615, -15.55509837]), array([0.        , 0.07160131, 0.33681129]), array([ 0.        ,  3.11265048, -1.89579344]), array([ 0.        , -2.39018083,  3.47919322]), array([  0.        , -13.3365158 ,  -4.65158994]), array([ 0.        , -2.74373197,  1.23599573]), array([0.5208251 , 0.51012399, 0.71204451]), array([ 0.        ,  4.23460398, -5.52276118]), array([0.        , 0.82549117, 1.39362231]), array([ 0.        ,  0.37960435, -1.20050775]), array([ 0.        ,  6.78503182, -3.82225911]), array([0.        , 0.47449619, 0.00696121]), array([ 0.        ,  1.28395459, -0.92745309]), array([0.        , 2.39253429, 0.97888172]), array([0.        , 0.80549727, 0.83939984]), array([ 0.        ,  1.1110785 , -0.12553364]), array([0.        , 0.59128912, 0.69537743]), array([ 0.        ,  7.88334855, -0.57744765]), array([0.        , 0.56593688, 1.13382314]), array([0.        , 0.8174455 , 2.08995712]), array([0.        , 0.15524347, 0.45488594]), array([0.        , 0.31763802, 0.14863396]), array([ 0.        ,  0.15454298, -0.34459269]), array([0.        , 0.57425025, 0.71261987]), array([0.        , 0.36226392, 0.30150532]), array([ 0.        , -0.14836875, -1.51914093]), array([0.        , 3.56370654, 3.30130816]), array([ 0.        , -3.10862952,  0.84508414]), array([ 0.        , -0.83858291,  1.85560904]), array([0.        , 0.86030372, 4.11654614]), array([ 0.        ,  1.2738802 , -0.66787768]), array([0.        , 2.68705888, 0.24942283]), array([0.28286757, 0.73284023, 0.30585556]), array([0.21133739, 0.56864672, 0.21654775]), array([ 0.        , -1.43342559, -3.3054237 ]), array([0.        , 0.58048418, 5.02055067]), array([ 0.        ,  0.21678578, -6.75022715]), array([ 0.        ,  3.51180126, -1.81856395]), array([0.        , 0.80604744, 0.41524411]), array([ 0.        ,  1.12999166, -0.35448048]), array([0.29861874, 0.75780138, 0.32728529]), array([0.        , 0.02540648, 0.44511313]), array([0.        , 2.25345203, 1.88086756]), array([0.        , 3.09003203, 0.40922359]), array([0.        , 4.19285436, 8.75573474]), array([ 0.        , -0.32719186,  4.30586195]), array([0.53016805, 0.67088462, 0.38405533]), array([0.55649338, 0.35010433, 0.41625447]), array([ 0.        , -0.07349654,  0.68549528]), array([0.03352744, 0.56744142, 0.05691195]), array([0.        , 0.30255174, 0.85173214]), array([  0.        , -19.15519259,  10.48607862]), array([  0.        ,   1.73063706, -12.69282561]), array([0.80574127, 0.4027731 , 0.76075328]), array([0.        , 2.45279012, 0.73411149]), array([0.        , 0.43425024, 2.6893467 ]), array([ 0.        ,  2.71887378, -1.84580744]), array([ 0.        , -1.54475364,  0.71990374]), array([0.36302166, 0.25852991, 0.68607469]), array([0.        , 0.32724076, 0.30938199]), array([ 0.        ,  3.13090294, -1.93529918]), array([0.        , 0.22615688, 0.17586553]), array([0.32927072, 0.79647481, 0.70749098]), array([0.36208289, 0.32143485, 0.50676656]), array([ 0.        ,  0.47152385, -0.02423639]), array([ 0.        , -0.08608688,  1.39216077]), array([ 0.        , -3.57592901,  3.0190272 ]), array([0.        , 1.83397481, 2.42130306]), array([0.        , 0.74740293, 1.29818789]), array([0.        , 0.39652422, 0.15261797]), array([0.        , 0.87099358, 0.05197794]), array([0.        , 1.03122168, 4.33753413]), array([0.54043938, 0.3639913 , 0.45294142]), array([0.68086877, 0.50676967, 0.46575077]), array([0.37185346, 0.80109426, 0.70904928]), array([ 0.        , -0.90176351, -0.20675436]), array([0.        , 1.16638451, 0.08934156]), array([ 0.        , -0.10539321,  0.38118901]), array([0.30626617, 0.4234194 , 0.24371747]), array([0.        , 0.20520323, 0.31762573]), array([ 0.        , -0.77654787, -2.56032476]), array([0.68861154, 0.62308866, 0.33518992]), array([ 0.        ,  5.11831554, -3.84604682]), array([ 0.        ,  1.06058985, -0.10957066]), array([ 0.        ,  0.27417234, -0.23366927]), array([ 0.        ,  0.66962848, -0.05782766]), array([ 0.        , -1.2387081 ,  1.13089488]), array([0.        , 0.38012165, 0.91711371]), array([0.        , 0.21113527, 0.01968646]), array([0.49617135, 0.50906505, 0.46889632]), array([ 0.        , -0.41105652,  1.0593762 ]), array([0.35482507, 0.53316502, 0.54408294]), array([ 0.        , -0.11597129,  0.77537678]), array([0.13784057, 0.28425011, 0.39443794]), array([ 0.        ,  1.46619366, -0.46069216]), array([ 0.        ,  0.55595469, -0.41784863]), array([0.38401499, 0.42665227, 0.63006207]), array([ 0.        ,  2.44540659, 13.21239818]), array([0.        , 0.0696753 , 0.92745276]), array([ 0.        , -7.90514199,  7.87498061]), array([0.        , 0.69589739, 0.26852679]), array([0.56826929, 0.67993177, 0.2538469 ]), array([0.        , 0.18848234, 0.985036  ]), array([0.51302438, 0.61790084, 0.51178033]), array([0.53940208, 0.61571458, 0.57684935]), array([0.64318588, 0.67414246, 0.61197508]), array([ 0.        , -0.47018563,  5.3719809 ]), array([ 0.       , -3.0831695,  0.8021056]), array([0.55703887, 0.64144023, 0.48795807]), array([ 0.        , -5.57694068,  3.35438717]), array([ 0.        ,  0.76810907, -1.28475695]), array([-0.14451084,  0.14945224,  0.42873705]), array([0.        , 0.46667717, 0.49899551]), array([ 0.        ,  0.5283606 , -0.74617083]), array([ 0.        ,  0.31349905, -0.55187175]), array([ 0.        , 20.74964494, -0.66962965]), array([ 0.        ,  0.61475065, -0.97244986]), array([ 0.        ,  0.24719069, -0.00297942]), array([0.        , 9.39536971, 0.67278768]), array([0.36198575, 0.57427527, 0.27076667]), array([0.        , 3.62858622, 1.4565151 ]), array([0.32436164, 0.38616074, 0.49530681]), array([0.40190462, 0.43064191, 0.60335295]), array([0.        , 2.98767423, 7.36090873]), array([  0.        , -38.58300057,  44.03198561]), array([0.        , 5.46404077, 1.90691163]), array([0.        , 1.00380395, 0.05923967]), array([0.        , 1.12000796, 2.29282899]), array([ 0.        , -2.10736006, -4.54767229]), array([ 0.        , -2.9430885 , -0.97504002]), array([0.        , 4.0201744 , 0.72303738]), array([0.55776934, 0.36813169, 0.71090962]), array([0.        , 1.0090027 , 0.98807923]), array([0.78502516, 0.52848122, 0.46866361]), array([0.        , 0.69766307, 0.15090324]), array([0.58854062, 0.63593113, 0.61932963]), array([0.84403949, 0.63599999, 0.51943821]), array([0.51442947, 0.55475533, 0.20664791]), array([0.81718117, 0.58706656, 0.49632992]), array([ 0.        , -0.92367724,  1.91824776]), array([  0.        , -57.20749886, -28.79836949]), array([0.        , 0.83414448, 0.33046952]), array([0.58219287, 0.38786575, 0.7780575 ]), array([0.38335565, 0.38604578, 0.3695222 ]), array([0.76400766, 0.52272921, 0.44065014]), array([0.27822194, 0.23677289, 0.55234053]), array([0.65102949, 0.52215118, 0.57328248]), array([0.57690308, 0.61861772, 0.80205424]), array([  0.        , -18.24687675, -24.26557189]), array([ 0.        , -1.3511548 ,  8.96011835]), array([0.45916536, 0.34199235, 0.74024012]), array([0.46186327, 0.22564711, 0.52560315]), array([ 0.        , -1.72762351,  0.7715981 ]), array([0.        , 0.59645752, 0.65978619]), array([0.        , 1.16906548, 0.91718131]), array([0.        , 1.01774671, 0.10966077]), array([0.40515619, 0.44100251, 0.67781946]), array([   0.        , -117.83242349,   79.1384295 ]), array([0.        , 4.41186847, 5.60275732]), array([0.32723988, 0.77139072, 0.36691464]), array([0.        , 1.56646891, 4.03399557]), array([ 0.        , -1.61746886, -3.49736972]), array([0.58420438, 0.76773877, 0.36714685]), array([ 0.        , -0.21545667,  0.3017603 ]), array([0.        , 0.43703472, 0.23271065]), array([0.03282878, 0.21603622, 0.00240551]), array([0.        , 1.52809205, 1.82150756]), array([ 0.        ,  1.20210072, -0.46050313]), array([0.        , 3.11652985, 1.11563879]), array([0.69201178, 0.50844554, 0.78819272]), array([ 0.        , -0.14139452, -5.69336811]), array([0.6698972 , 0.49211198, 0.35758748]), array([0.55308415, 0.62670963, 0.4610268 ]), array([0.        , 0.29410323, 0.57995442]), array([0.        , 0.55751235, 1.26560478]), array([ 0.        , -0.85830359, -0.08631917]), array([0.67498528, 0.74303843, 0.51297705]), array([0.50783137, 0.57099349, 0.45417297]), array([0.        , 0.38890319, 0.92995339]), array([0.        , 4.52745415, 6.889113  ]), array([ 0.        , -0.85146049,  0.79140739]), array([0.35432476, 0.54232009, 0.23722704]), array([ 0.        , -2.75459241, -1.08813998]), array([ 0.        ,  1.03993363, -0.04307021]), array([ 0.        , -0.97810643, -0.8335427 ]), array([0.        , 1.96532328, 3.14493164]), array([ 0.        , -8.49370365,  1.49904547]), array([0.22694611, 0.54119875, 0.42864988]), array([ 0.        , -0.53915196, -1.83847811]), array([0.4765076 , 0.6108548 , 0.34168544]), array([ 0.        ,  1.45137934, -2.42271103]), array([0.        , 3.3009491 , 2.86329276]), array([0.        , 0.57833616, 1.83008897]), array([0.52526418, 0.59568771, 0.52191278]), array([0.58579001, 0.50467023, 0.64195857]), array([0.63256477, 0.38997998, 0.76364798]), array([ 0.        , -1.44018475, -3.46366934]), array([  0.        ,   3.0387028 , -14.51473036]), array([ 0.        ,  5.35908391, -0.70648627]), array([ 0.        , -0.70451618, -0.37353582]), array([0.        , 0.08537511, 1.48972285]), array([ 0.        , -1.49843044,  2.55704435]), array([0.        , 0.95943003, 0.17025548]), array([0.28854942, 0.20579049, 0.2494014 ]), array([ 0.        , -0.68843595, -1.01689654]), array([0.        , 0.60734797, 0.52911509]), array([ 0.        , -3.04804627, -1.21042761]), array([ 0.        , -5.19740196,  2.45372865]), array([0.81735137, 0.61686149, 0.92991316]), array([0.        , 0.10100919, 7.96669828]), array([ 0.        ,  3.47297508, -1.29148403]), array([ 0.        ,  0.49219733, -0.97348301]), array([0.23426448, 0.1989177 , 0.53389497]), array([0.38399511, 0.3815269 , 0.55991399]), array([0.        , 2.46223999, 0.64264687]), array([ 0.        ,  0.95485111, -3.99757654]), array([0.        , 0.34594906, 1.43579309]), array([0.37099836, 0.39571104, 0.30393849]), array([  0.        , -12.62848053,  -7.62743931]), array([ 0.        ,  2.3320944 , -0.19547343]), array([ 0.        ,  0.30912678, -1.49710506]), array([0.        , 0.79204782, 0.73572682]), array([ 0.        , -0.40134493,  1.61508672]), array([0.20370885, 0.46222447, 0.43610339]), array([0.19388846, 0.14788359, 0.2939623 ]), array([0.87288107, 0.68854672, 0.54425286]), array([0.        , 1.3596999 , 4.72097766]), array([0.87606171, 0.50797779, 0.30791124]), array([ 0.05189626, -0.01155834,  0.51079289]), array([ 0.        ,  0.40968524, 14.24310613]), array([0.        , 1.57938459, 0.02297372]), array([ 0.        ,  0.5390832 , -0.07869047]), array([  0.        ,   2.52648384, -13.47067112]), array([0.57915586, 0.46210612, 0.87214236]), array([0.        , 0.09249234, 0.04424786]), array([0.        , 1.24490131, 0.46727379]), array([ 0.        ,  0.6270503 , -2.65842722]), array([ 0.        ,  1.96481928, -0.15976645]), array([0.86848796, 0.74382398, 0.5950492 ]), array([0.73640153, 0.57671968, 0.48785054]), array([0.        , 0.37712993, 0.48257814]), array([ 0.        ,  0.46774357, -2.28982772]), array([0.81829577, 0.48111996, 0.63705148]), array([0.78319728, 0.73465484, 0.67786616]), array([0.        , 0.03465681, 1.0272664 ]), array([0.63019862, 0.3226226 , 0.43909266]), array([ 0.        , 15.85832509,  7.68163777]), array([0.63377259, 0.58939122, 0.71096372]), array([0.55320397, 0.472823  , 0.71787787]), array([0.57478776, 0.77399999, 0.53404452]), array([0.63039257, 0.73325969, 0.63926986]), array([ 0.        , -0.34850645,  0.49807629]), array([ 0.        , -1.88218603,  3.92697701]), array([0.79208929, 0.48825216, 0.7628171 ]), array([0.858364  , 0.55342644, 0.52392981]), array([0.        , 8.60065337, 5.84948706]), array([ 0.        ,  1.05400694, -1.82298947]), array([ 0.        , -2.12520892, -1.4850402 ]), array([0.        , 0.0512864 , 1.72264955]), array([0.21826872, 0.43972552, 0.69537682]), array([ 0.        ,  1.18945742, -0.43510208]), array([ 0.        ,  0.38994652, -1.49772869]), array([0.87199144, 0.73783696, 0.55899825]), array([ 0.        , -0.510863  ,  1.99970824]), array([0.36650282, 0.76738372, 0.89102323]), array([ 0.        , -1.39802078,  1.64410946]), array([0.        , 0.78428173, 0.15440672]), array([ 0.        , -3.15946751, -5.09093393]), array([0.81851183, 0.92783484, 0.72730183]), array([0.6793359 , 0.77929201, 0.65618346]), array([0.70312964, 0.61152086, 0.59620289]), array([0.87579911, 0.83004769, 0.67906344]), array([0.50519875, 0.85768154, 0.47562974]), array([0.        , 0.63742823, 0.482653  ]), array([ 0.        ,  1.08111507, -0.98080466]), array([0.8647734 , 0.49730719, 0.65882712]), array([0.56668544, 0.66764468, 0.25006313]), array([   0.        , -610.91833072, -796.05327858]), array([ 0.        ,  0.94801854, -0.02721254]), array([ 0.        , -1.02698668,  0.05428892]), array([0.93262428, 0.57758576, 0.58947254]), array([0.        , 0.05518537, 1.62936572]), array([ 0.        , -5.37887142,  3.35781544]), array([0.77886451, 0.87818508, 0.59997976]), array([0.        , 0.89161763, 1.09069029]), array([0.79324657, 0.87794491, 0.67346311]), array([0.64064809, 0.91987683, 0.50743802]), array([ 0.        , -6.3705554 , -7.58614756])] )
			#self.ridge_vertices =[[729, 0], [730, 1], [731, 2], [732, 3], [733, 4], [734, 5], [735, 6], [736, 7], [737, 8], [738, 9], [739, 10], [740, 1118], [1118, 11], [741, 78], [78, 12], [742, 13], [927, 743], [743, 14], [744, 15], [745, 16], [746, 17], [747, 18], [748, 19], [749, 20], [750, 21], [751, 22], [752, 23], [753, 24], [754, 25], [755, 26], [756, 27], [757, 28], [758, 1165], [1165, 1137], [1137, 29], [759, 30], [760, 31], [1133, 761], [761, 32], [762, 33], [763, 34], [764, 35], [765, 1335], [1335, 36], [766, 37], [767, 909], [909, 1014], [1014, 38], [768, 39], [769, 40], [770, 1013], [1013, 1122], [1122, 41], [771, 42], [772, 43], [773, 44], [774, 45], [775, 46], [776, 47], [777, 48], [778, 49], [779, 50], [780, 51], [781, 52], [53, 782], [783, 54], [784, 55], [785, 56], [786, 57], [787, 58], [788, 59], [789, 60], [790, 61], [791, 62], [792, 63], [793, 64], [794, 65], [795, 66], [796, 67], [861, 797], [797, 68], [69, 798], [799, 70], [800, 71], [801, 72], [802, 73], [803, 74], [804, 75], [845, 805], [805, 76], [806, 77], [807, 78], [808, 79], [809, 80], [810, 81], [811, 1055], [1055, 981], [981, 1313], [1313, 24], [24, 82], [812, 83], [813, 1058], [1058, 84], [814, 85], [815, 86], [816, 87], [817, 88], [818, 89], [819, 90], [820, 91], [821, 92], [822, 93], [823, 94], [824, 95], [825, 96], [826, 97], [827, 98], [828, 99], [829, 100], [830, 1088], [1088, 1097], [1097, 1012], [1012, 101], [831, 102], [832, 103], [833, 104], [834, 904], [904, 1005], [1005, 43], [43, 105], [835, 106], [836, 107], [837, 108], [838, 109], [948, 839], [839, 987], [987, 110], [840, 111], [841, 112], [842, 113], [843, 114], [844, 115], [845, 116], [846, 117], [847, 47], [47, 118], [848, 119], [849, 1025], [1025, 930], [930, 1357], [1357, 120], [850, 121], [851, 122], [852, 123], [853, 124], [854, 125], [855, 126], [856, 127], [857, 128], [858, 129], [859, 130], [860, 131], [861, 132], [862, 133], [863, 1002], [1002, 134], [864, 135], [946, 865], [865, 1210], [1210, 1329], [1329, 136], [866, 137], [867, 138], [868, 139], [869, 140], [870, 141], [871, 142], [872, 143], [873, 144], [874, 145], [875, 146], [1023, 938], [938, 876], [876, 147], [877, 148], [878, 149], [69, 53], [53, 879], [879, 150], [880, 151], [881, 1053], [1053, 152], [882, 153], [883, 154], [884, 155], [885, 1160], [1160, 156], [886, 157], [887, 158], [1061, 888], [888, 159], [889, 160], [890, 161], [891, 162], [892, 163], [893, 164], [894, 165], [895, 166], [896, 167], [897, 168], [898, 169], [899, 170], [900, 1136], [1136, 3], [3, 171], [901, 172], [902, 1373], [1373, 173], [903, 1233], [1233, 1383], [1383, 204], [204, 1292], [1292, 174], [904, 175], [905, 176], [906, 14], [14, 104], [104, 1418], [1418, 177], [907, 178], [908, 179], [909, 180], [910, 181], [911, 182], [912, 183], [913, 184], [914, 185], [915, 1263], [1263, 16], [16, 186], [916, 187], [917, 188], [918, 189], [919, 1248], [1248, 1021], [1021, 1112], [1112, 190], [920, 191], [921, 192], [922, 193], [923, 194], [924, 195], [925, 196], [926, 197], [927, 198], [928, 199], [929, 200], [930, 201], [931, 202], [932, 203], [933, 32], [32, 1429], [1429, 204], [934, 205], [935, 1103], [1103, 1394], [1394, 913], [913, 206], [936, 207], [937, 208], [938, 209], [939, 210], [940, 211], [941, 212], [942, 122], [122, 1315], [1315, 213], [943, 214], [944, 215], [945, 216], [946, 217], [947, 218], [948, 161], [161, 219], [949, 220], [950, 221], [951, 222], [952, 223], [953, 224], [954, 225], [955, 226], [956, 227], [957, 228], [958, 1030], [1030, 1203], [1203, 76], [76, 229], [959, 230], [960, 231], [961, 232], [962, 1095], [1095, 1107], [1107, 1382], [1382, 1033], [1033, 233], [963, 234], [964, 235], [965, 236], [966, 237], [967, 238], [968, 239], [969, 240], [970, 1219], [1219, 241], [971, 1159], [1159, 242], [972, 243], [973, 1259], [1259, 237], [237, 244], [974, 245], [975, 246], [976, 247], [977, 248], [978, 249], [979, 250], [980, 1059], [1059, 137], [137, 251], [981, 252], [982, 163], [163, 1064], [1064, 100], [100, 1143], [1143, 253], [983, 1091], [1091, 75], [75, 254], [984, 1242], [1242, 253], [253, 15], [15, 309], [309, 255], [985, 256], [986, 257], [987, 258], [988, 259], [989, 260], [990, 261], [991, 262], [992, 263], [993, 264], [994, 265], [995, 266], [996, 1398], [1398, 88], [88, 1228], [1228, 1320], [1320, 267], [997, 1403], [1403, 268], [998, 269], [999, 270], [1000, 271], [1001, 272], [1002, 273], [1003, 274], [1004, 275], [1005, 276], [1006, 52], [52, 277], [1007, 278], [1008, 279], [1009, 280], [1010, 1387], [1387, 281], [1011, 282], [1012, 39], [39, 1075], [1075, 22], [22, 11], [11, 110], [110, 203], [203, 56], [56, 283], [1013, 1169], [1169, 155], [155, 284], [1014, 285], [1015, 286], [1016, 287], [1017, 288], [1018, 289], [1019, 290], [1020, 291], [1021, 292], [1022, 293], [1023, 294], [1024, 61], [61, 295], [1025, 35], [35, 143], [143, 142], [142, 1393], [1393, 296], [1026, 297], [1027, 298], [1028, 299], [1029, 300], [1030, 301], [1031, 302], [1032, 303], [1033, 304], [1034, 305], [1035, 307], [307, 1108], [1108, 240], [240, 146], [146, 306], [1036, 307], [1037, 308], [1038, 309], [1039, 1139], [1139, 28], [28, 1109], [1109, 1351], [1351, 1289], [1289, 120], [120, 310], [1040, 311], [1041, 312], [1042, 313], [1043, 1434], [1434, 314], [1044, 315], [1045, 216], [216, 316], [1046, 317], [1047, 318], [1048, 319], [1049, 320], [1050, 1120], [1120, 1130], [1130, 321], [1051, 322], [1052, 323], [1053, 238], [238, 324], [1054, 325], [1055, 326], [1056, 264], [264, 1312], [1312, 327], [1057, 328], [1058, 123], [123, 329], [1059, 263], [263, 330], [1060, 106], [106, 34], [34, 331], [1061, 332], [1062, 1446], [1446, 1271], [1271, 333], [1063, 334], [1064, 335], [1065, 336], [1066, 169], [169, 285], [285, 199], [199, 1316], [1316, 318], [318, 337], [1067, 338], [1068, 339], [1069, 340], [1070, 341], [1071, 342], [1072, 261], [261, 102], [102, 283], [283, 1305], [1305, 343], [1073, 344], [1074, 1217], [1217, 129], [129, 1243], [1243, 1409], [1409, 345], [1075, 346], [1076, 347], [1077, 348], [1078, 349], [1079, 350], [1080, 351], [1081, 352], [1082, 353], [1083, 354], [1084, 355], [1085, 356], [1086, 357], [1087, 358], [1088, 359], [1089, 360], [1090, 1154], [1154, 1307], [1307, 361], [1091, 362], [1092, 363], [1093, 364], [1094, 365], [1095, 1176], [1176, 366], [1096, 367], [1097, 26], [26, 86], [86, 1291], [1291, 277], [277, 1266], [1266, 1326], [1326, 368], [1098, 369], [1099, 370], [1100, 371], [1101, 372], [1102, 1359], [1359, 1161], [1161, 218], [218, 1273], [1273, 373], [1103, 1167], [1167, 1135], [1135, 1420], [1420, 67], [67, 374], [1104, 375], [1105, 376], [1106, 377], [1107, 378], [1108, 379], [1109, 380], [1110, 63], [63, 381], [1111, 382], [1112, 383], [1113, 384], [1114, 385], [1115, 386], [1116, 387], [1117, 388], [1118, 389], [1119, 390], [1120, 391], [1121, 392], [1122, 1232], [1232, 1244], [1244, 393], [1123, 394], [1124, 395], [1125, 396], [1126, 397], [1127, 398], [1128, 399], [1129, 400], [1130, 401], [1131, 402], [1132, 62], [62, 1202], [1202, 403], [1133, 1281], [1281, 404], [1134, 101], [101, 337], [337, 232], [232, 405], [1135, 406], [1136, 189], [189, 304], [304, 1223], [1223, 335], [335, 395], [395, 467], [467, 286], [286, 407], [1137, 408], [1138, 409], [1139, 410], [1140, 411], [1141, 412], [1142, 413], [1143, 414], [1144, 415], [1145, 1275], [1275, 416], [1146, 417], [1147, 418], [1148, 419], [1149, 420], [1150, 421], [1151, 422], [1152, 423], [1153, 424], [1154, 425], [1155, 1415], [1415, 183], [183, 1445], [1445, 426], [1156, 427], [1157, 1417], [1417, 300], [300, 428], [1158, 429], [1159, 430], [1160, 431], [1161, 432], [1162, 433], [1163, 434], [1164, 435], [1165, 436], [1166, 437], [1167, 438], [1168, 439], [1169, 440], [1170, 441], [1171, 442], [1172, 443], [1173, 444], [1174, 445], [1175, 446], [1176, 447], [1177, 448], [1178, 449], [1179, 450], [1180, 451], [1181, 27], [27, 247], [247, 1251], [1251, 452], [1182, 453], [1183, 454], [1184, 455], [1185, 456], [1186, 1317], [1317, 457], [1187, 458], [1188, 459], [1189, 460], [1190, 461], [1191, 66], [66, 1216], [1216, 462], [1192, 224], [224, 483], [483, 463], [1193, 464], [1194, 222], [222, 465], [1195, 121], [121, 1347], [1347, 1363], [1363, 1274], [1274, 466], [1196, 467], [1197, 468], [1198, 469], [1199, 470], [1200, 471], [1201, 472], [1202, 473], [1203, 474], [1204, 182], [182, 1342], [1342, 279], [279, 77], [77, 1278], [1278, 82], [82, 475], [1205, 476], [1206, 297], [297, 306], [306, 477], [1207, 478], [1208, 479], [1209, 480], [1210, 1442], [1442, 1456], [1456, 481], [1211, 1261], [1261, 151], [151, 482], [1212, 483], [1213, 484], [1214, 485], [1215, 486], [1216, 487], [1217, 488], [1218, 489], [1219, 490], [1220, 491], [1221, 492], [1222, 493], [1223, 494], [1224, 495], [1225, 496], [1226, 497], [1227, 498], [1228, 499], [1229, 500], [1230, 501], [1231, 502], [1232, 503], [1233, 504], [1234, 505], [1235, 506], [1236, 507], [1237, 508], [1238, 509], [1239, 510], [1240, 511], [1241, 512], [1242, 513], [1243, 514], [1244, 515], [1245, 1303], [1303, 1308], [1308, 91], [91, 1306], [1306, 1395], [1395, 516], [1246, 517], [1247, 1321], [1321, 399], [399, 1412], [1412, 311], [311, 518], [518, 280], [1248, 519], [1249, 97], [97, 1301], [1301, 1365], [1365, 215], [215, 520], [1250, 221], [221, 521], [1251, 522], [1252, 523], [1253, 524], [1254, 1341], [1341, 1314], [1314, 1424], [1424, 1450], [1450, 525], [1255, 526], [1256, 527], [1257, 528], [1258, 130], [130, 202], [202, 1339], [1339, 529], [1259, 1364], [1364, 530], [1260, 531], [1261, 532], [1262, 533], [1263, 534], [1264, 535], [1265, 1440], [1440, 71], [71, 536], [1266, 537], [1267, 538], [1268, 539], [1269, 540], [1270, 541], [1271, 542], [1272, 543], [1273, 544], [1274, 545], [1275, 546], [1276, 547], [1277, 548], [1278, 549], [1279, 550], [1280, 1441], [1441, 551], [1281, 552], [1282, 553], [1283, 1332], [1332, 48], [48, 473], [473, 393], [393, 554], [1284, 74], [74, 302], [302, 49], [49, 1439], [1439, 555], [1285, 556], [1286, 1419], [1419, 557], [1287, 558], [1288, 559], [1289, 560], [1290, 561], [1291, 562], [1292, 563], [1293, 564], [1294, 565], [1295, 566], [1296, 567], [1297, 568], [1298, 489], [489, 569], [1299, 57], [57, 570], [1300, 571], [1301, 572], [1302, 573], [1303, 574], [1304, 575], [1305, 576], [1306, 577], [1307, 578], [1308, 579], [1309, 580], [1310, 581], [1311, 582], [1312, 583], [1313, 584], [1314, 585], [1315, 586], [1316, 587], [1317, 588], [1318, 59], [59, 1378], [1378, 589], [1319, 590], [1320, 591], [1321, 592], [1322, 593], [593, 471], [471, 536], [1323, 594], [1324, 595], [1325, 596], [1326, 597], [1327, 598], [1328, 599], [1329, 600], [1330, 601], [1331, 602], [1332, 603], [1333, 604], [1334, 95], [95, 1397], [1397, 439], [439, 276], [276, 605], [1335, 606], [1336, 607], [1337, 608], [1338, 609], [1339, 610], [1340, 1423], [1423, 530], [530, 611], [1341, 612], [1342, 613], [1343, 614], [1344, 615], [1345, 1432], [1432, 616], [616, 255], [1346, 617], [1347, 618], [1348, 619], [1349, 620], [1350, 621], [1351, 622], [1352, 139], [139, 623], [1353, 624], [1354, 113], [113, 184], [184, 625], [1355, 626], [1356, 627], [1357, 628], [1358, 37], [37, 196], [196, 293], [293, 29], [29, 70], [70, 629], [1359, 630], [1360, 631], [1361, 632], [1362, 633], [1363, 634], [1364, 635], [1365, 636], [1366, 1413], [1413, 350], [350, 17], [17, 637], [1367, 638], [1368, 639], [1369, 640], [1370, 641], [1371, 642], [1372, 643], [1373, 1346], [1346, 1453], [1453, 374], [374, 644], [1374, 117], [117, 1455], [1455, 217], [217, 1438], [1438, 645], [1375, 646], [1376, 647], [1377, 648], [1378, 649], [1379, 650], [1380, 651], [1381, 652], [1382, 653], [1383, 654], [1384, 655], [1385, 656], [1386, 657], [1387, 658], [1388, 659], [1389, 660], [1390, 661], [1391, 662], [1392, 663], [1393, 664], [1394, 665], [1395, 666], [1396, 667], [1397, 668], [1398, 669], [1399, 670], [1400, 671], [1401, 672], [1402, 673], [1403, 674], [1404, 675], [1405, 676], [1406, 677], [1407, 678], [1408, 679], [1409, 680], [1410, 681], [1411, 682], [1412, 683], [1413, 684], [1414, 685], [1415, 686], [1416, 687], [1417, 688], [1418, 689], [1419, 690], [1420, 691], [1421, 453], [453, 692], [1422, 693], [1423, 694], [1424, 695], [1425, 696], [1426, 697], [1427, 698], [1428, 699], [1429, 700], [1430, 701], [1431, 702], [702, 305], [1432, 703], [1433, 704], [1434, 705], [1435, 706], [1436, 707], [1437, 385], [385, 411], [411, 450], [450, 708], [1438, 709], [1439, 710], [1440, 711], [1441, 712], [1442, 713], [1443, 714], [1444, 715], [212, 1445], [1445, 1408], [1408, 455], [455, 716], [1446, 717], [1447, 718], [1448, 719], [1449, 720], [1450, 721], [1451, 722], [1452, 723], [1453, 724], [1454, 725], [1455, 480], [480, 726], [1456, 727], [1457, 344], [344, 728]]
			self.save_network(1, path, name = 'creation')
			#print(list(self.vertices), self.ridge_vertices)
			print('0 ',len(self.ridge_vertices), len(self.vertices))
			self = self.sort_nodes()
			self = self.delete_boundary_to_boundary_fibre()
			print('1 ',len(self.ridge_vertices), len(self.vertices))
			self = self.delete_doubles()
			plot_geometry(self)
			print('2 ',len(self.ridge_vertices), len(self.vertices))
			self = self.cut_network()
			print('3 ',len(self.ridge_vertices), len(self.vertices))
			plot_geometry(self)
			self = self.delete_alone_points()
			print('4 ',len(self.ridge_vertices), len(self.vertices))
			self = self.delete_single_ridge_points()
			print('5 ',len(self.ridge_vertices), len(self.vertices))
			self.minimum_distance = 0.005
			self = self.merge_nodes()
			print('6 ',len(self.ridge_vertices), len(self.vertices))
			self = self.delete_single_ridge_points()
			print('7 ',len(self.ridge_vertices), len(self.vertices))
			self=self.delete_points_with_two_ridges()
			print('8 ',len(self.ridge_vertices), len(self.vertices))
			self = self.delete_single_ridge_points()
			print('9 ',len(self.ridge_vertices), len(self.vertices))
			self = self.sort_nodes()
			self = self.delete_boundary_to_boundary_fibre()
			print('10 ',len(self.ridge_vertices), len(self.vertices))
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
					#if self.dimension == 2: self = self.delete_crossings()
					#if self.creation == 'growth_network': self = self.add_boundary_nodes_adaptedtoGN()
					#if self.dimension == 2: self = self.delete_doubles_growth_network()
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
			if self.dimension == 2: 
				print('done')
				self = self.delete_doubles_growth_network()
			#from Plotting.network_plotting import plot_geometry
			#import matplotlib.pyplot as plt"""
		if self.dimension == 2: 
			print('done')
			self = self.delete_doubles_growth_network()
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
