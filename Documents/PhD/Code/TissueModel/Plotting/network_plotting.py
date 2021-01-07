import numpy
import os, sys, fnmatch, csv
sys.path.append('/home/aude/Documents/PhD/Code/TissueModel/')
import Network_generation.creation_network 
import matplotlib.pyplot as plt

import matplotlib.patches as patches

def length_square(network,x):
	if int(network.dimension) == 2:
		return x[0]**2+x[1]**2
	if int(network.dimension) == 3:
		return x[0]**2+x[1]**2+x[2]**2


####################################### PLOTTING THE GEOMETRY OF NETWORKS ###########################################


def plot_geometry(network, **kw):
	fig = plt.figure()
	line_segments = []
	if int(network.dimension) == 2:
		ax = fig.gca()
		import numpy as np
		from numpy import array
		from matplotlib.collections import LineCollection
		ax.scatter(network.vertices[:,0],network.vertices[:,1], s =2.)
		for i in range(len(network.vertices[:,0])):
			ax.annotate(i, (network.vertices[i,0],network.vertices[i,1]),fontsize=10)
		for simplex in network.ridge_vertices:
		        simplex = numpy.asarray(simplex)
		        line_segments.append([(x, y) for x, y in network.vertices[simplex]])
		lc = LineCollection(line_segments,linestyle='solid')
		ax.set_axis_off()
	if int(network.dimension) ==3:
		from mpl_toolkits.mplot3d import Axes3D
		from mpl_toolkits.mplot3d.art3d import Line3DCollection
		ax = fig.add_subplot(111, projection='3d')
		if kw.get('show_vertices', True):
			ax.scatter(network.vertices[:,0],network.vertices[:,1],network.vertices[:,2])
		for i in range(len(network.vertices[:, 0])):
			ax.text(network.vertices[i, 0], network.vertices[i, 1],network.vertices[i,2],i)
		for simplex in network.ridge_vertices:
		        simplex = numpy.asarray(simplex)
		        line_segments.append([(x, y, z) for x, y, z in network.vertices[simplex]])
		ax.set_xlim3d([0.0,network.length[0]])
		ax.set_ylim3d([0.0,network.length[0]])
		ax.set_zlim3d([0.0,network.length[0]])
		ax.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
		ax.set_axis_off()
		lc = Line3DCollection(line_segments,linestyle='solid')
	ax.add_collection(lc)
	return ax.figure		

def plot_network_extension(network, **kw):
	import matplotlib.pyplot as plt
	fig = plt.figure()
	ax = fig.gca()
	plt.xlim([-0.1,network.length*1.5])
	plt.ylim([-0.1,network.length*1.1])
	from matplotlib.collections import LineCollection
	if kw.get('show_vertices', True):
		ax.scatter(network.vertices[:,0],network.vertices[:,1], color='red')
	line_segments_ini = []
	for simplex_ini in network.ridge_vertices:
	        simplex_ini = numpy.asarray(simplex_ini)
	        line_segments_ini.append([(x, y) for x, y in network.vertices_ini[simplex_ini]])
	line_segments = []
	for simplex in network.ridge_vertices:
	        simplex = numpy.asarray(simplex)
	        line_segments.append([(x, y) for x, y in network.vertices[simplex]])
	lc_ini = LineCollection(line_segments_ini,linestyle='dashed', color='grey', label='initial')
	lc = LineCollection(line_segments,linestyle='solid', label='after tensile test', color='red')
	ax.add_collection(lc_ini)
	ax.add_collection(lc)
	ax.legend()
	return ax.figure




######################## MAIN ##################################

def sorted_ls(path):
    mtime = lambda f: os.stat(os.path.join(path, f)).st_mtime
    return list(sorted(os.listdir(path), key=mtime))
