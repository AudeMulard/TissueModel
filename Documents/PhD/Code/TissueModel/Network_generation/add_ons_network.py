from matplotlib.patches import Circle
import sympy.geometry as sg
import sympy
import numpy as np


def length_square(x):
	if len(x)==2:
		return x[0]**2+x[1]**2
	if len(x) ==3:
		return x[0]**2+x[1]**2+x[2]**2

class Cell:
	
	def __init__(self,centre,radius):
		self.centre_x = centre[0]
		self.centre_y = centre[1]
		if len(centre) == 3: self.centre_z = centre[2]
		self.radius = radius
		self.boundary_cell = []

	def add_cell(self,network):
		cell = Circle((self.centre_x,self.centre_y),self.radius)
		nodes = network.vertices
		ridges = network.ridge_vertices
		nodes_to_delete=[]
		ridges_to_delete=[]
		for i in range(len(nodes)):
			if cell.contains_point(nodes[i])==True:
				nodes_to_delete = np.append(nodes_to_delete,i)
		for i in range(len(ridges)):
			if ridges[i][0] in nodes_to_delete and ridges[i][1] in nodes_to_delete:
				ridges_to_delete = np.append(ridges_to_delete,i)
		ridges_to_delete = np.array(sorted(ridges_to_delete,reverse=True))
		for ridge in ridges_to_delete:
			network.ridge_vertices = np.delete(network.ridge_vertices,int(ridge),axis=0)
		center = sg.Point(self.centre_x,self.centre_y)
		circ = sg.Circle(center, self.radius)
		k=0
		for node in nodes_to_delete:
			node = int(node)
			for ridge in network.ridge_vertices:
				if ridge[0] == node:
					if len(nodes[node])==3:
						In_point = sg.Point(nodes[node][0],nodes[node][1],nodes[node][2])
						Out_point = sg.Point(nodes[ridge[1]][0],nodes[ridge[1]][1],nodes[ridge[1]][1])
					elif len(nodes[node])==2:
						In_point = sg.Point(nodes[node][0],nodes[node][1])
						Out_point = sg.Point(nodes[ridge[1]][0],nodes[ridge[1]][1])
					line = sg.Line(In_point,Out_point)
					intersec = sg.intersection(circ, line)
					if length_square(nodes[ridge[1]]-intersec[0]) >= length_square(nodes[ridge[1]]-intersec[1]):
						nodes=np.append(nodes,[[float(intersec[1].x),float(intersec[1].y)]],axis=0)
					else:
						nodes=np.append(nodes,[[float(intersec[0].x),float(intersec[0].y)]],axis=0)
					ridge[0]=len(nodes)-1
					k+=1
				if ridge[1] == node:
					if len(nodes[node])==3:
						In_point = sg.Point(nodes[node][0],nodes[node][1],nodes[node][2])
						Out_point = sg.Point(nodes[ridge[0]][0],nodes[ridge[0]][1],nodes[ridge[0]][1])
					elif len(nodes[node])==2:
						In_point = sg.Point(nodes[node][0],nodes[node][1])
						Out_point = sg.Point(nodes[ridge[0]][0],nodes[ridge[0]][1])
					line = sg.Line(In_point,Out_point)
					intersec = sg.intersection(circ, line)
					if length_square(nodes[ridge[0]]-intersec[0]) >= length_square(nodes[ridge[0]]-intersec[1]):
						nodes=np.append(nodes,[[float(intersec[1].x),float(intersec[1].y)]],axis=0)
					else:
						nodes=np.append(nodes,[[float(intersec[0].x),float(intersec[0].y)]],axis=0)
					ridge[1]=len(nodes)-1
					k+=1
		nodes_to_delete = np.array(sorted(nodes_to_delete, reverse=True))
		for point in nodes_to_delete:
			nodes=np.delete(nodes,int(point),0)
		# Renumber points after deleting some
		for ridge in network.ridge_vertices:
			for i in range(2):
				r=0
				for node in nodes_to_delete:
					if node < ridge[i]:
						r+=1
				ridge[i] = ridge[i] - r 
		network.vertices=nodes
		network.vertices_ini=np.array(network.vertices.tolist())
		network = network.create_ridge_node_list()
		network = network.sort_nodes()
		network.interior_nodes = network.interior_nodes[:-k]
		self.boundary_cell = list(range(len(nodes)-k,len(nodes)))
		return network

	def contraction_bc(self,network,distance):
		print(self.boundary_cell)
		for node in self.boundary_cell:
			network.vertices[node][0]=distance*(self.centre_x-network.vertices[node][0])+network.vertices[node][0]
			network.vertices[node][1]=distance*(self.centre_y-network.vertices[node][1])+network.vertices[node][1]
		return network
