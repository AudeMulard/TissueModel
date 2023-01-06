import numpy as np
import csv, os, fnmatch,sys

def length_square(x):
	if len(x) ==2: return x[0]**2+x[1]**2
	if len(x) ==3: return x[0]**2+x[1]**2+x[2]**2

class Network:
	def __init__(self, dimension, complexity_network, length_domain, min_distance, beam_Young, beam_poisson, beam_profile, connector_coeff, hyperstatic_param, creation, generation, path, **kw):
		self.complexity = complexity_network
		self.length=length_domain
		self.dimension = dimension
		self.beam_young = beam_Young
		self.beam_poisson = beam_poisson
		self.beam_profile = beam_profile
		self.connector_coeff = connector_coeff
		self.min_distance = min_distance
		self.creation = creation
		self.generation = generation
		self.lengths_ini = []
		self. hyperstatic_param =  hyperstatic_param
		if kw.get('name')!=None: self.kw = kw
		else: self.kw = kw
		self.interior_nodes = []
		self.boundary_nodes_left=[]
		self.boundary_nodes_right=[]
		
	def create_network(self):
		import Network_generation.network_types
		self.vertices, self.ridge_vertices = Network_generation.network_types.select_network(self)
		return self

	def delete_list_points(self,nodes_to_delete):
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
		self = self.create_ridge_node_list()
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
			self = self.delete_alone_points()
			self = self.sort_nodes()
			length_lists=[]
			for node in self.interior_nodes:
				length_lists.append(len(self.list_nodes_ridges[node]))
		self= self.delete_single_ridge_points()
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
		self = self.delete_list_points(nodes_to_delete)
		ridges_to_delete = []
		for i in range(len(self.ridge_vertices)):
			if self.ridge_vertices[i][0]==self.ridge_vertices[i][1]:
				ridges_to_delete = np.append(ridges_to_delete,i)
		ridges_to_delete=np.array(sorted(ridges_to_delete, reverse=True))
		for ridge in ridges_to_delete:
			self.ridge_vertices=np.delete(self.ridge_vertices,int(ridge),0)
		self = self.delete_doubles()
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
		self = self.delete_list_points(nodes_to_delete)
		# Delete ridges that are alone beacuse they were only attached to ridges that went outside
		self = self.delete_single_ridge_points()
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
		self.vertices=nodes
		self = self.delete_list_points(nodes_to_delete)
		return self

	def cut_network(self):
		if self.creation == 'Voronoi': self=self.cut_network_side()
		self=self.cut_network_updown()
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
		
	def delete_boundary_ridges(self):
		self = self.create_ridge_node_list()
		self = self.sort_nodes()
		ridges_to_delete = []
		for k in range(len(self.ridge_vertices)):
			if self.ridge_vertices[k][0] in self.boundary_nodes_left and self.ridge_vertices[k][1] in self.boundary_nodes_left:
				ridges_to_delete = np.append(ridges_to_delete,k)
			if self.ridge_vertices[k][0] in self.boundary_nodes_right and self.ridge_vertices[k][1] in self.boundary_nodes_right:
				ridges_to_delete = np.append(ridges_to_delete,k)
		ridges_to_delete=np.array(sorted(ridges_to_delete, reverse=True))
		for ridge in ridges_to_delete:
			self.ridge_vertices=np.delete(self.ridge_vertices,int(ridge), axis=0)
		self = self.create_ridge_node_list()
		self = self.delete_alone_points()
		return self

	def delete_alone_points(self):
		self = self.create_ridge_node_list()
		nodes_to_delete = []
		for i in range(len(self.vertices)):
			if len(self.list_nodes_ridges[i])==0:
				nodes_to_delete = np.append(nodes_to_delete,i)
		self = self.delete_list_points(nodes_to_delete)
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

	def save_network(self, step, path,**kw):
		filenames=sorted(fnmatch.filter(os.listdir(path), 'network_vertices_initial_*.csv'))
		number_network = len(filenames)
		filename = 'parameters_%02d_%05d.csv' % (number_network,len(self.vertices))
		with open(os.path.join(path,filename), 'w') as writeFile:
			writer = csv.writer(writeFile)
			writer.writerow(["dimension",self.dimension])
			writer.writerow(["complexity",self.complexity])
			writer.writerow(["length",self.length])
			writer.writerow(["merge_distance",self.min_distance])
			writer.writerow(["number of nodes",len(self.vertices)])
			writer.writerow(["number of springs",len(self.ridge_vertices)])
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
		last_network = 'network_vertices_initial_%02d_%05d.csv' % (number_network,len(self.vertices))
		with open(os.path.join(path,last_network), 'w') as writeFile:
			writer = csv.writer(writeFile)
			writer.writerows(self.vertices)
		writeFile.close()

# Global function that applies all the corrections to the network
	def generate_network(self, path):
		self = self.create_network()
		from Plotting.network_plotting import plot_geometry
		import matplotlib.pyplot as plt
		if self.creation == 'Voronoi': self = self.delete_first_point()
		length = []
		self = self.create_ridge_node_list()
		for i in range(len(self.list_nodes_ridges)):
			length.append(len(self.list_nodes_ridges[i]))
		print('before cut',np.mean(length))
		self = self.cut_network()
		length = []
		self = self.create_ridge_node_list()
		for i in range(len(self.list_nodes_ridges)):
			if i not in self.boundary_nodes_left and i not in self.boundary_nodes_right:
				length.append(len(self.list_nodes_ridges[i]))
		print('after cut',np.mean(length))
		distances = []
		for ridge in self.ridge_vertices:
			nodes = self.vertices
			i= ridge[0]
			j = ridge[1]
			if self.dimension == 2:
				distance = np.sqrt((nodes[i,0]-nodes[j,0])**2+(nodes[i,1]-nodes[j,1])**2)
			elif self.dimension == 3:
				distance = np.sqrt((nodes[i,0]-nodes[j,0])**2+(nodes[i,1]-nodes[j,1])**2+(nodes[i,2]-nodes[j,2])**2)
			distances.append(distance)
		if self.dimension == 2: self.min_distance = max(distances)/10.
		if self.dimension == 3: self.min_distance = max(distances)/50.
		self = self.merge_nodes()
		self = self.merge_nodes()
		self = self.sort_nodes()
		length = []
		self = self.create_ridge_node_list()
		for i in range(len(self.list_nodes_ridges)):
			length.append(len(self.list_nodes_ridges[i]))
		print(np.mean(length))
		while len(self.boundary_nodes_right) < min(self.complexity/10,2) or len(self.boundary_nodes_left) < min(self.complexity/10,2): #create a network with enough boundary points
			print(self.min_distance, len(self.boundary_nodes_right),len(self.boundary_nodes_left))
			self = self.create_network()
			if self.creation == 'Voronoi': self = self.delete_first_point()
			self = self.cut_network()
			self = self.merge_nodes()
			self = self.sort_nodes()
			print(self.min_distance, len(self.boundary_nodes_right),len(self.boundary_nodes_left))
		while len(self.ridge_vertices)-float(self.hyperstatic_param)*len(self.interior_nodes)<=0:  #create a network with enough constraints: Maxwell criterion
			self.min_distance +=0.001
			self = self.merge_nodes()
			self = self.sort_nodes()
		length = []
		self = self.create_ridge_node_list()
		for i in range(len(self.list_nodes_ridges)):
			length.append(len(self.list_nodes_ridges[i]))
		print(np.mean(length))
		self = self.delete_points_with_two_ridges()
		#nodes = self.vertices
		self = self.delete_doubles()
		self = self.delete_boundary_ridges()
		ridges_to_delete = []
		for i in range(len(self.ridge_vertices)):
			if self.ridge_vertices[i][0]==self.ridge_vertices[i][1]:
				ridges_to_delete.append(i)
		ridges_to_delete=np.array(sorted(ridges_to_delete, reverse=True))
		for k in ridges_to_delete:
			self.ridge_vertices=np.delete(self.ridge_vertices, int(k), axis=0)
		self.save_network('initial', path)
		distances = []
		for ridge in self.ridge_vertices:
			nodes = self.vertices
			i= ridge[0]
			j = ridge[1]
			if self.dimension == 2:
				distance = np.sqrt((nodes[i,0]-nodes[j,0])**2+(nodes[i,1]-nodes[j,1])**2)
			elif self.dimension == 3:
				distance = np.sqrt((nodes[i,0]-nodes[j,0])**2+(nodes[i,1]-nodes[j,1])**2+(nodes[i,2]-nodes[j,2])**2)
			distances.append(distance)
		self.min_distance = min(distances)
		print(max(distances))
		return self

	def set_fibers(self,path):
		self = self.generate_network(path)
		while len(self.boundary_nodes_left)==0 or len(self.boundary_nodes_right)==0: 
			print('No boundary nodes,starting again')
			self = self.generate_network(path)
		return self
