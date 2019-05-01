import numpy as np
from scipy.spatial import Voronoi, voronoi_plot_2d
import network_types

# Class of the network, with different initial positions and all the corrections needed

class Network:
	def __init__(self, dimension, complexity_network, length_domain, min_distance, Ef, A, B, creation):
		self.complexity = complexity_network
		self.length=length_domain
		self.dimension = dimension
		self.Ef = Ef
		self.A = A
		self.B = B
		self.min_distance = min_distance
		
	def create_network(self, creation):
		self.vertices, self.ridge_vertices = network_types.select_network(self, creation)
		return self

# Deletion of the first point because it creates intersections in the network

	def delete_first_point(self):
		self.vertices = np.delete(self.vertices, -1, axis=0)
		ridges_to_delete = []
		for k in range(len(self.ridge_vertices)):
			if self.ridge_vertices[k][0]==-1 or self.ridge_vertices[k][1]==-1 or self.ridge_vertices[k][0]==len(self.vertices) or self.ridge_vertices[k][1]==len(self.vertices):
				ridges_to_delete= np.append(ridges_to_delete, k)
		ridges_to_delete=np.array(sorted(ridges_to_delete, reverse=True))
		if self.dimension == 2:
			for k in ridges_to_delete:
				self.ridge_vertices=np.delete(self.ridge_vertices, k, axis=0)
		if self.dimension == 3:
			for k in ridges_to_delete:
				self.ridge_vertices=np.delete(self.ridge_vertices, k, axis=0)
		return self

	def delete_point(self, i):
		self.vertices = np.delete(self.vertices, i, axis=0)
		for ridge in self.ridge_vertices:
			for i in range(2):
				if i < ridge[i]:
					ridge[i] -= 1
		return self

# Sort the nodes in three parts: two boundary lists and the interior nodes list

	def sort_nodes(self):
		nodes = self.vertices
		boundary_nodes_left =[]
		boundary_nodes_right=[]
		interior_nodes=[]
		for i in range(len(nodes)):
        		if nodes[i,0] <= 0.0:
        		        boundary_nodes_left.append(i)
        		elif nodes[i,0] >= 1.0:
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
			nodes=np.delete(nodes,point,0)
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
			self.ridge_vertices=np.delete(self.ridge_vertices,ridge,0)
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
        			if nodes[i,coord] < 0.0 or nodes[i,coord] > 1.0:
					nodes_to_delete = np.append(nodes_to_delete,i)
		nodes_to_delete = set(nodes_to_delete)
		# delete ridges that are completely out of the domain
		for i in range(len(self.ridge_vertices)):
			if self.ridge_vertices[i][0] in nodes_to_delete and self.ridge_vertices[i][1] in nodes_to_delete:
				ridges_to_delete = np.append(ridges_to_delete, i)
		ridges_to_delete=np.array(sorted(ridges_to_delete, reverse=True))
		for ridge in ridges_to_delete:
			self.ridge_vertices=np.delete(self.ridge_vertices,ridge, axis=0)
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
			nodes=np.delete(nodes,point,0)
		# Renumber points in the ridge_vertices after deleting some
		for ridge in self.ridge_vertices:
			for i in range(2):
				r=0
				for node in nodes_to_delete:
					if node < ridge[i]:
						r+=1
				ridge[i] = ridge[i] - r
		self.vertices=nodes
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
				nodes=np.delete(nodes,point,0)
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
			if nodes[i,0] > 1.0:
				nodes_to_delete_right = np.append(nodes_to_delete_right,i)
		nodes_to_delete = np.append(nodes_to_delete_left, nodes_to_delete_right)
		for i in range(len(self.ridge_vertices)):
			if self.ridge_vertices[i][0] in nodes_to_delete and self.ridge_vertices[i][1] in nodes_to_delete:
				ridges_to_delete = np.append(ridges_to_delete, i)
		ridges_to_delete=np.array(sorted(ridges_to_delete, reverse=True))
		for ridge in ridges_to_delete:
			self.ridge_vertices=np.delete(self.ridge_vertices,ridge, axis=0)
		for node in nodes_to_delete_left:
			node=int(node)
                	for k in range(len(self.ridge_vertices)):
                        	if self.ridge_vertices[k][0]==node:
					if self.dimension == 2:
			                	y = nodes[node,1]+(nodes[self.ridge_vertices[k][1],1]-nodes[node,1])/(nodes[self.ridge_vertices[k][1],0]-nodes[node,0])*(0.0-nodes[node,0])
        			                nodes=np.append(nodes,[[0.0,y]], axis=0)
					if self.dimension == 3:
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
		                		y = nodes[node,1]+(nodes[self.ridge_vertices[k][1],1]-nodes[node,1])/(nodes[self.ridge_vertices[k][1],0]-nodes[node,0])*(1.0-nodes[node,0])
        		                	nodes=np.append(nodes,[[1.0,y]], axis=0)
					if self.dimension == 3:
						y = nodes[node,1]+(nodes[self.ridge_vertices[k][1],1]-nodes[node,1])/(nodes[self.ridge_vertices[k][1],0]-nodes[node,0])*(1.0-nodes[node,0])
						z = nodes[node,2]+(nodes[self.ridge_vertices[k][1],2]-nodes[node,2])/(nodes[self.ridge_vertices[k][1],0]-nodes[node,0])*(1.0-nodes[node,0])
        			                nodes=np.append(nodes,[[1.0,y,z]], axis=0)
					self.ridge_vertices[k][0]=len(nodes)-1
		                elif self.ridge_vertices[k][1]==node:
					if self.dimension == 2:
						y = nodes[node,1]+(nodes[self.ridge_vertices[k][0],1]-nodes[node,1])/(nodes[self.ridge_vertices[k][0],0]-nodes[node,0])*(1.0-nodes[node,0])
						nodes=np.append(nodes,[[1.0,y]], axis=0)
					if self.dimension == 3:
						y = nodes[node,1]+(nodes[self.ridge_vertices[k][0],1]-nodes[node,1])/(nodes[self.ridge_vertices[k][0],0]-nodes[node,0])*(1.0-nodes[node,0])
						z = nodes[node,2]+(nodes[self.ridge_vertices[k][0],2]-nodes[node,2])/(nodes[self.ridge_vertices[k][0],0]-nodes[node,0])*(1.0-nodes[node,0])
        			                nodes=np.append(nodes,[[1.0,y,z]], axis=0)
					self.ridge_vertices[k][1]=len(nodes)-1
		nodes_to_delete=np.array(sorted(nodes_to_delete, reverse=True))
		for point in nodes_to_delete:
			nodes=np.delete(nodes,point,0)
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
		nodes_to_delete = []
		for i in range(len(self.vertices)):
			if len(self.list_nodes_ridges[i])==0:
				nodes_to_delete = np.append(nodes_to_delete,i)
		nodes_to_delete=np.array(sorted(nodes_to_delete, reverse=True))
		for point in nodes_to_delete:
			self.vertices=np.delete(self.vertices,point,0)
		# Renumber points after deleting some
		for ridge in self.ridge_vertices:
			for i in range(2):
				r=0
				for node in nodes_to_delete:
					if node < ridge[i]:
						r+=1
				ridge[i] = ridge[i] - r
		return self

# Global function that applies all the corrections to the network

	def set_fibers(self, creation):
		self = self.create_network(creation)
		if creation == 'Voronoi' or creation == 'random network 3':
			self = self.delete_first_point()
			self = self.cut_network()
		self = self.merge_nodes()
		self = self.sort_nodes()
		self = self.delete_doubles()
		self = self.create_ridge_node_list()
		print self.list_nodes_ridges
		self = self.delete_alone_points()
		self = self.create_ridge_node_list()
		self.vertices_ini = np.array(self.vertices.tolist())
		return self

# Plots the network

	def plot_network(self, **kw):
		import matplotlib.pyplot as plt
		fig = plt.figure()
		line_segments = []
		if self.dimension == 2:
			ax = fig.gca()
			from matplotlib.collections import LineCollection
			if kw.get('show_vertices', True):
				ax.scatter(self.vertices[:,0],self.vertices[:,1])
			for simplex in self.ridge_vertices:
			        simplex = np.asarray(simplex)
		            	line_segments.append([(x, y) for x, y in self.vertices[simplex]])
			lc = LineCollection(line_segments,linestyle='solid')
		if self.dimension ==3:
			from mpl_toolkits.mplot3d import Axes3D
			from mpl_toolkits.mplot3d.art3d import Line3DCollection
			ax = fig.add_subplot(111, projection='3d')
			if kw.get('show_vertices', True):
				ax.scatter(self.vertices[:,0],self.vertices[:,1],self.vertices[:,2])
			for simplex in self.ridge_vertices:
			        simplex = np.asarray(simplex)
		            	line_segments.append([(x, y, z) for x, y, z in self.vertices[simplex]])
			ax.set_xlim3d([0.0,1.0])
			ax.set_ylim3d([0.0,1.0])
			ax.set_zlim3d([0.0,1.0])		
  			lc = Line3DCollection(line_segments,linestyle='solid')
		ax.add_collection(lc)
		return ax.figure		

	def plot_network_extension(self, **kw):
		import matplotlib.pyplot as plt
		fig = plt.figure()
		ax = fig.gca()
		plt.xlim([-0.5,1.5])
		plt.ylim([-0.5,1.0])
		from matplotlib.collections import LineCollection
		ax.scatter(self.vertices_ini[1:,0],self.vertices_ini[1:,1], color='grey')
		if kw.get('show_vertices', True):
			ax.scatter(self.vertices[:,0],self.vertices[:,1], color='red')
		line_segments_ini = []
		for simplex_ini in self.ridge_vertices:
		        simplex_ini = np.asarray(simplex_ini)
	            	line_segments_ini.append([(x, y) for x, y in self.vertices_ini[simplex_ini]])
		line_segments = []
		for simplex in self.ridge_vertices:
		        simplex = np.asarray(simplex)
	            	line_segments.append([(x, y) for x, y in self.vertices[simplex]])
		lc_ini = LineCollection(line_segments_ini,linestyle='solid', color='grey', label='initial')
  		lc = LineCollection(line_segments,linestyle='solid', label='after tensile test', color='red')
		ax.add_collection(lc_ini)
		ax.add_collection(lc)
		ax.legend()
		return ax.figure
