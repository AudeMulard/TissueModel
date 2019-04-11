import numpy as np
from scipy.spatial import Voronoi, voronoi_plot_2d

# Class of the network, with different initial positions and all the corrections needed

class Network:
	def __init__(self, dimension, complexity_network, length_domain, min_distance, Ef, A, B, creation):
		self.complexity = complexity_network
		self.length=length_domain
		self.dimension = dimension
		self.Ef = Ef
		self.A = A
		self.B = B
		
		if creation == "Voronoi":
			Seeds=np.random.rand(self.complexity,self.dimension)*self.length
			self.voronoi = Voronoi(Seeds)
			self.vertices = Voronoi(Seeds).vertices
			self.ridge_vertices = Voronoi(Seeds).ridge_vertices
			self.min_distance = min_distance
			#self.vertices_ini = self.set_fibers().vertices
		
		if creation == "big grid":
			self.vertices = np.array([[0.0,0.25],[0.0,0.5],[0.0,0.75],[0.25,0.0],[0.25,0.25],[0.25,0.5],[0.25,0.75],[0.25,1.0],[0.5,0.0],[0.5,0.25],[0.5,0.5],[0.5,0.75],[0.5,1.0],[0.75,0.0],[0.75,0.25],[0.75,0.5],[0.75,0.75],[0.75,1.0],[1.0,0.25],[1.0,0.5],[1.0,0.75]])
			self.ridge_vertices = np.array([[0,4],[4,9],[9,14],[14,18],[1,5],[5,10],[10,15],[15,19],[2,6],[6,11],[11,16],[16,20],[3,4],[4,5],[5,6],[6,7],[8,9],[9,10],[10,11],[11,12],[13,14],[14,15],[15,16],[16,17]])
			self.min_distance = min_distance

		if creation == "small grid":
			self.vertices = np.array([[0.0,0.33],[0.0,0.66],[0.33,0.0],[0.33,0.33],[0.33,0.66],[0.33,1.0],[0.66,0.0],[0.66,0.33],[0.66,0.66],[0.66,1.0],[1.0,0.33],[1.0,0.66]])
			self.ridge_vertices = np.array([[0,3],[3,7],[7,10],[1,4],[4,8],[8,11],[2,3],[3,4],[4,5],[6,7],[7,8],[8,9]])
			self.min_distance = min_distance

		if creation == "1 straight line":
			self.vertices = np.array([[0.0,0.0],[0.5,0.0],[1.0,0.0]])
			self.ridge_vertices = np.array([[0,1],[1,2]])
			self.min_distance = min_distance

		if creation == "1 long straight line":
			self.vertices = np.array([[0.0,0.0],[0.25,0.0],[0.5,0.0],[0.75,0.0],[1.0,0.0]])
			self.ridge_vertices = np.array([[0,1],[1,2],[2,3],[3,4]])
			self.min_distance = min_distance

		if creation == "1 long cross":
			self.vertices = np.array([[0.0,0.5],[0.5,0.5],[1.0,0.5],[0.5,0.0],[0.5,0.25],[0.5,0.75],[0.5,1.0]])
			self.ridge_vertices = np.array([[0,1],[1,2],[3,4],[4,1],[1,5],[5,6]])
			self.min_distance = min_distance

		if creation == "1 half long cross":
			self.vertices = np.array([[0.0,0.5],[0.25,0.5],[0.5,0.5],[0.75,0.5],[1.0,0.5],[0.5,0.25],[0.5,0.75]])
			self.ridge_vertices = np.array([[0,1],[1,2],[2,3],[3,4],[5,6],[5,2],[2,6]])
			self.min_distance = min_distance

		if creation == "manual":
			self.vertices = np.array([[0.0,0.0],[0.3,0.1],[0.6,0.3],[1.0,0.0]])
			self.ridge_vertices = np.array([[0,1],[1,2],[2,3]])
			self.min_distance = min_distance

		if creation == "1 straight line":
			self.vertices = np.array([[0.0,0.0],[0.5,0.0],[1.0,0.0]])
			self.ridge_vertices = np.array([[0,1],[1,2]])
			self.min_distance = min_distance

		if creation == "symmetric cross":
			self.vertices = np.array([[0.0,0.0],[0.5,0.0],[1.0,0.5],[1.0,-0.5]])
			self.ridge_vertices = np.array([[0,1],[1,2],[1,3]])
			self.min_distance = min_distance
#			self.vertices_ini = np.array(self.vertices.tolist())
	
		if creation == "asymmetric cross":
			self.vertices = np.array([[0.0,0.0],[0.7,-0.1],[1.0,0.2],[1.0,-0.3]])
			self.ridge_vertices = np.array([[0,1],[1,2],[1,3]])
			self.min_distance = min_distance

		if creation == "symmetric double point cross":
			self.vertices = np.array([[0.0,0.0],[0.0,1.0],[0.61,0.5],[0.69,0.5],[1.0,0.0],[1.0,1.0]])
			self.ridge_vertices = np.array([[0,2],[1,2],[2,3],[3,5],[3,4]])
			self.min_distance = min_distance

		if creation == "asymmetric double point cross":
			self.vertices = np.array([[-1.0,0.0],[-0.6,1.0],[0.19, 0.45],[0.22,0.51], [0.20,0.55],[0.69,0.62],[1.3,0.0],[2.0,1.0]])
			self.ridge_vertices = np.array([[0,2],[1,3],[2,3],[2,4],[3,4],[4,5],[4,6],[5,6],[5,7]])
			self.min_distance = min_distance

		if creation == "random network":
			self.vertices = np.array([[6.96810003e-01,5.85972462e-01], [  6.12222074e-01,   5.70042014e-01], [  1.53524456e-02,   8.24030940e-01], [  2.02945404e-02,   3.95001038e-01], [  5.91361447e-02,   6.07631016e-01], [  2.17734065e-01,   4.04153611e-01], [  1.00918445e-01,   6.11192022e-01], [  7.60619750e-01,  -4.27877180e-02], [  6.42935772e-01,   1.57907165e-01], [  8.73596877e-01,   3.18524365e-01], [  8.51906402e-01,   3.31465928e-01], [  4.47723467e-01,   3.51768688e-02], [  5.85307175e-01,   2.02853812e-01], [  4.64801092e-01,   3.75003275e-01], [  4.15171979e-01,   3.63331353e-01], [  3.16687341e-01,   1.87921226e-01], [  2.91696122e-01,   2.93138477e-01], [  3.09155556e-01,   3.21897785e-01], [  5.91396931e-01,   5.94012634e-01], [  8.14387805e-01,   1.72239083e+01], [  3.15164596e-01,   9.25079572e-01], [  1.13885589e-01,   6.27208867e-01], [  3.00439561e-01,   5.75803663e-01], [  4.17769233e-01,   7.61100709e-01], [  3.17801428e-01,   8.99124444e-01], [  2.97685323e-01,   6.19842149e-01], [  1.43974960e-01,   6.85770872e-01], [  1.24611979e-01,   6.43587311e-01]])
			self.ridge_vertices=np.array([[-1, 4], [-1, 3], [3, 4], [3, 5], [4, 6], [5, 6], [7, 9], [7, 8], [8, 10], [9, 10], [-1, 9], [-1, 7], [-1, 0], [0, 10], [11, 15], [11, 12], [12, 13], [13, 14], [14, 17], [15, 16], [16, 17], [-1, 15], [-1, 11], [8, 12], [0, 1], [1, 13], [-1, 16], [5, 17], [-1, 19], [1, 18], [18, 19], [-1, 2], [2, 20], [19, 20], [6, 21], [14, 22], [21, 22], [23, 24], [23, 25], [24, 26], [25, 27], [26, 27], [18, 23], [20, 24], [22, 25], [2, 26], [21, 27]])
			self.min_distance = min_distance

		if creation == "1 straight line 2":
			self.vertices = np.array([[-1.0,0.0],[0.5,0.0],[1.5,0.0]])
			self.ridge_vertices = np.array([[0,1],[1,2]])
			self.min_distance = min_distance

		if creation == "symmetric cross 2":
			self.vertices = np.array([[-1.0,0.0],[0.5,0.0],[1.5,0.5],[1.5,-0.5]])
			self.ridge_vertices = np.array([[0,1],[1,2],[1,3]])
			self.min_distance = min_distance

		if creation == "random network 2":
			self.vertices = np.array([[ 0.6886007,   1.22144112], [ 0.69731848,  1.1186779 ], [ 0.51616784,  0.63207458], [-0.25276059,  0.45293642], [ 0.49692315, -1.5152983 ], [ 0.34634344,  0.13602371], [ 0.33809343,  0.80025459], [ 0.19167735,  0.81100432], [-0.08463533,  0.51910635], [ 1.42874148,  0.32477703], [ 0.44490991,  0.64721725], [ 0.3765596,   0.67967657], [ 0.22545142,  0.52383797], [ 0.71710902,  1.04705891], [ 1.28244566,  0.44924527], [ 0.75903549,  0.8663037 ], [ 0.52775066,  0.61855448], [ 0.67122359, 0.02915664], [ 0.571157,    0.21613333], [ 0.36876754,  0.43831976], [ 0.35649777,  0.43378642],[ 0.54799455,  0.46092557], [ 0.56521616,  0.44491184], [ 0.51447422,  0.24044857], [ 0.34884099,  0.19190268], [ 0.3510838,   0.18885809], [ 0.61471224,  0.44134212], [ 0.84182501,  0.33648427], [ 0.75687701,  0.64930453], [ 0.74807286,  0.64506609], [ 1.05464284,  0.36505829], [ 0.91089201,  0.3296533 ]])
			self.ridge_vertices = np.array([[-1, 4], [-1, 5], [4, 5], [0, 6], [-1, 0], [-1, 7], [6, 7], [-1, 3], [3, 8], [7, 8], [0, 1], [1, 2], [2, 10], [6, 11], [10, 11], [8, 12], [11, 12], [13, 15], [13, 14], [14, 15], [-1, 9], [1, 13], [9, 14], [19, 20], [19, 21], [20, 24], [21, 22], [22, 23], [23, 25], [24, 25], [10, 19], [12, 20], [2, 16], [16, 21], [3, 24], [4, 17], [5, 25], [17, 18], [18, 23], [26, 29], [26, 27], [27, 31], [28, 30], [28, 29], [30, 31], [16, 29], [22, 26], [18, 27], [17, 31], [9, 30], [15, 28]])
			self.min_distance = min_distance


# Deletion of the first point because it creates intersections in the network

	def delete_first_point(self):
		self.vertices = np.delete(self.vertices, -1, axis=0)
		ridges_to_delete = []
		for k in range(len(self.ridge_vertices)):
			if self.ridge_vertices[k][0]==-1 or self.ridge_vertices[k][1]==-1 or self.ridge_vertices[k][0]==len(self.vertices) or self.ridge_vertices[k][1]==len(self.vertices):
				ridges_to_delete= np.append(ridges_to_delete, k)
		ridges_to_delete=np.array(sorted(ridges_to_delete, reverse=True))
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
		nodes_to_delete_down = []
		nodes_to_delete_up = []
		ridges_to_delete = []
		for i in range(len(nodes)):
        		if nodes[i,1] < 0.0:
				nodes_to_delete_down = np.append(nodes_to_delete_down,i)
			if nodes[i,1] > 1.0:
				nodes_to_delete_up = np.append(nodes_to_delete_up,i)
		nodes_to_delete = np.append(nodes_to_delete_down, nodes_to_delete_up)
		for i in range(len(self.ridge_vertices)):
			if self.ridge_vertices[i][0] in nodes_to_delete and self.ridge_vertices[i][1] in nodes_to_delete:
				ridges_to_delete = np.append(ridges_to_delete, i)
		ridges_to_delete=np.array(sorted(ridges_to_delete, reverse=True))
		for i in range(len(self.ridge_vertices)):
			if self.ridge_vertices[i][0] in nodes_to_delete and self.ridge_vertices[i][1] in nodes_to_delete:
				ridges_to_delete = np.append(ridges_to_delete, i)
		for ridge in ridges_to_delete:
			self.ridge_vertices=np.delete(self.ridge_vertices,ridge, axis=0)
		for node in nodes_to_delete_down:
			node= int(node)
                	for k in range(len(self.ridge_vertices)):
                        	if self.ridge_vertices[k][0]==node:
		                	x = nodes[node,0]+1.0/(nodes[self.ridge_vertices[k][1],1]-nodes[node,1])*(nodes[self.ridge_vertices[k][1],0]-nodes[node,0])*(0.0-nodes[node,1])
        		                nodes=np.append(nodes,[[x,0.0]], axis=0)
        		                self.ridge_vertices[k][0]=len(nodes)-1
		                elif self.ridge_vertices[k][1]==node:
                		        x = nodes[node,0]+1.0/(nodes[self.ridge_vertices[k][0],1]-nodes[node,1])*(nodes[self.ridge_vertices[k][0],0]-nodes[node,0])*(0.0-nodes[node,1])
                		        nodes=np.append(nodes,[[x,0.0]], axis=0)
					self.ridge_vertices[k][1]=len(nodes)-1
		for node in nodes_to_delete_up:
			node= int(node)
                	for k in range(len(self.ridge_vertices)):
                        	if self.ridge_vertices[k][0]==node:
		                	x = nodes[node,0]+1.0/(nodes[self.ridge_vertices[k][1],1]-nodes[node,1])*(nodes[self.ridge_vertices[k][1],0]-nodes[node,0])*(1.0-nodes[node,1])
        		                nodes=np.append(nodes,[[x,1.0]], axis=0)
					self.ridge_vertices[k][0]=len(nodes)-1
		                elif self.ridge_vertices[k][1]==node:
					x = nodes[node,0]+1.0/(nodes[self.ridge_vertices[k][0],1]-nodes[node,1])*(nodes[self.ridge_vertices[k][0],0]-nodes[node,0])*(1.0-nodes[node,1])
					nodes=np.append(nodes,[[x,1.0]], axis=0)
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
		                	y = nodes[node,1]+(nodes[self.ridge_vertices[k][1],1]-nodes[node,1])/(nodes[self.ridge_vertices[k][1],0]-nodes[node,0])*(0.0-nodes[node,0])
        		                nodes=np.append(nodes,[[0.0,y]], axis=0)
        		                self.ridge_vertices[k][0]=len(nodes)-1
		                elif self.ridge_vertices[k][1]==node:
                		        y = nodes[node,1]+(nodes[self.ridge_vertices[k][0],1]-nodes[node,1])/(nodes[self.ridge_vertices[k][0],0]-nodes[node,0])*(0.0-nodes[node,0])
                		        nodes=np.append(nodes,[[0.0,y]], axis=0)
					self.ridge_vertices[k][1]=len(nodes)-1
		for node in nodes_to_delete_right:
			node= int(node)
                	for k in range(len(self.ridge_vertices)):
                        	if self.ridge_vertices[k][0]==node:
		                	y = nodes[node,1]+(nodes[self.ridge_vertices[k][1],1]-nodes[node,1])/(nodes[self.ridge_vertices[k][1],0]-nodes[node,0])*(1.0-nodes[node,0])
        		                nodes=np.append(nodes,[[1.0,y]], axis=0)
					self.ridge_vertices[k][0]=len(nodes)-1
		                elif self.ridge_vertices[k][1]==node:
					y = nodes[node,1]+(nodes[self.ridge_vertices[k][0],1]-nodes[node,1])/(nodes[self.ridge_vertices[k][0],0]-nodes[node,0])*(1.0-nodes[node,0])
					nodes=np.append(nodes,[[1.0,y]], axis=0)
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

# Global function that applies all the corrections to the network

	def set_fibers(self, creation):
		if creation == 'Voronoi':
			self.delete_first_point()
		self = self.cut_network()
		self = self.merge_nodes()
		self = self.sort_nodes()
		self = self.create_ridge_node_list()
		self.vertices_ini = np.array(self.vertices.tolist())
		return self

# Plots the network

	def plot_network(self, **kw):
		import matplotlib.pyplot as plt
		fig = plt.figure()
		ax = fig.gca()
		plt.xlim([0.0,1.0])
		plt.ylim([0.0,1.0])
		from matplotlib.collections import LineCollection
		if kw.get('show_vertices', True):
			ax.scatter(self.vertices[:,0],self.vertices[:,1])
		line_segments = []
		for simplex in self.ridge_vertices:
		        simplex = np.asarray(simplex)
	            	line_segments.append([(x, y) for x, y in self.vertices[simplex]])
  		lc = LineCollection(line_segments,linestyle='solid')
		ax.add_collection(lc)
		return ax.figure

	def plot_network_extension(self, **kw):
		import matplotlib.pyplot as plt
		fig = plt.figure()
		ax = fig.gca()
#		plt.xlim([-0.5,1.5])
#		plt.ylim([0.0,1.0])
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
