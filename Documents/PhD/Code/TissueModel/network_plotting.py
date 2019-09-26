import numpy as np
from force_balance import length_square
import matplotlib.pyplot as plt

# Plot geometry

def plot_network(network,  **kw):
	import matplotlib.pyplot as plt
	fig = plt.figure()
	line_segments = []
	if network.dimension == 2:
		ax = fig.gca()
		from matplotlib.collections import LineCollection
		if kw.get('show_vertices', True):
			ax.scatter(network.vertices[:,0],network.vertices[:,1], s =2.)
		for simplex in network.ridge_vertices:
		        simplex = np.asarray(simplex)
		        line_segments.append([(x, y) for x, y in network.vertices[simplex]])
		lc = LineCollection(line_segments,linestyle='solid')
	if network.dimension ==3:
		from mpl_toolkits.mplot3d import Axes3D
		from mpl_toolkits.mplot3d.art3d import Line3DCollection
		ax = fig.add_subplot(111, projection='3d')
		if kw.get('show_vertices', True):
			ax.scatter(network.vertices[:,0],network.vertices[:,1],network.vertices[:,2])
		for simplex in network.ridge_vertices:
		        simplex = np.asarray(simplex)
		        line_segments.append([(x, y, z) for x, y, z in network.vertices[simplex]])
		ax.set_xlim3d([0.0,network.length])
		ax.set_ylim3d([0.0,network.length])
		ax.set_zlim3d([0.0,network.length])		
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
	#ax.scatter(network.vertices_ini[1:,0],network.vertices_ini[1:,1], color='grey')
	if kw.get('show_vertices', True):
		ax.scatter(network.vertices[:,0],network.vertices[:,1], color='red')
	line_segments_ini = []
	for simplex_ini in network.ridge_vertices:
	        simplex_ini = np.asarray(simplex_ini)
	        line_segments_ini.append([(x, y) for x, y in network.vertices_ini[simplex_ini]])
	line_segments = []
	for simplex in network.ridge_vertices:
	        simplex = np.asarray(simplex)
	        line_segments.append([(x, y) for x, y in network.vertices[simplex]])
	lc_ini = LineCollection(line_segments_ini,linestyle='dashed', color='grey', label='initial')
	lc = LineCollection(line_segments,linestyle='solid', label='after tensile test', color='red')
	ax.add_collection(lc_ini)
	ax.add_collection(lc)
	ax.legend()
	return ax.figure


# Constraints

def plot_constraints(network):
	from matplotlib import cm
	import matplotlib.colors
	ridge_constraints=[]
	for ridge in network.ridge_vertices:
		ridge_constraints.append(network.Ef*(np.sqrt(length_square(network,network.vertices[ridge[0]]-network.vertices[ridge[1]]))-np.sqrt(length_square(network,network.vertices_ini[ridge[0]]-network.vertices_ini[ridge[1]]))))
	max_constraint = max(ridge_constraints)
	cmap = plt.cm.rainbow
	norm = matplotlib.colors.Normalize(vmin=0., vmax=max_constraint)
	fig = plt.figure()
	line_segments = []
	ax = fig.gca()
	from matplotlib.collections import LineCollection
	ax.scatter(network.vertices[:,0],network.vertices[:,1])
	for simplex in network.ridge_vertices:
		simplex = np.asarray(simplex)
		line_segments.append([(x, y) for x, y in network.vertices[simplex]])
	lc = LineCollection(line_segments,linestyle='solid',cmap=cmap,norm=norm)
	sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
	lc.set_array(np.array(ridge_constraints))
	ax.add_collection(lc)
	sm.set_array([])
	fig.colorbar(sm, orientation='vertical')
	return ax.figure
