import numpy as np
import os, sys, fnmatch, csv
sys.path.append('/home/aude/Documents/PhD/Code/TissueModel/')
from Core_calculation.force_balance import length_square
import Network_generation.creation_network 
import matplotlib.pyplot as plt

import matplotlib.patches as patches

def length_square(network,x):
	if int(network.dimension) == 2:
		return x[0]**2+x[1]**2


####################################### PLOTTING THE GEOMETRY OF NETWORKS ###########################################


def plot_geometry(network,  **kw):
	fig = plt.figure()
	line_segments = []
	if int(network.dimension) == 2:
		ax = fig.gca()
		#ax.axis('equal')
		#ax.set(xlim=(0., 1.0), ylim=(0., 1.0))
		from matplotlib.collections import LineCollection
		ax.scatter(network.vertices[:,0],network.vertices[:,1], s =2.)
		for i in range(len(network.vertices[:,0])):
			ax.annotate(i, (network.vertices[i,0],network.vertices[i,1]),fontsize=10)
		#circ = patches.Circle((0.5, 0.5), 0.1, alpha=0.8, fc='yellow')
		#ax.add_patch(circ)
		for simplex in network.ridge_vertices:
		        simplex = np.asarray(simplex)
		        line_segments.append([(x, y) for x, y in network.vertices[simplex]])
		lc = LineCollection(line_segments,linestyle='solid')
	if int(network.dimension) ==3:
		from mpl_toolkits.mplot3d import Axes3D
		from mpl_toolkits.mplot3d.art3d import Line3DCollection
		ax = fig.add_subplot(111, projection='3d')
		if kw.get('show_vertices', True):
			ax.scatter(network.vertices[:,0],network.vertices[:,1],network.vertices[:,2])
		for simplex in network.ridge_vertices:
		        simplex = np.asarray(simplex)
		        line_segments.append([(x, y, z) for x, y, z in network.vertices[simplex]])
		#ax.set_xlim3d([0.0,network.length])
		#ax.set_ylim3d([0.0,network.length])
		#ax.set_zlim3d([0.0,network.length])		
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

def plot_constraints(network, **kw):
	from matplotlib import cm
	import matplotlib.colors
	ridge_constraints=[]
	### ATTENTION: TOOK NETWORK.EF OUT, SO NO COEFFICIENT
	for ridge in network.ridge_vertices:
		ridge_constraints.append((np.sqrt(length_square(network,network.vertices[ridge[0]]-network.vertices[ridge[1]]))-np.sqrt(length_square(network,network.vertices_ini[ridge[0]]-network.vertices_ini[ridge[1]]))))
	max_constraint = max(ridge_constraints)
	cmap = plt.cm.rainbow
	norm = matplotlib.colors.Normalize(vmin=0., vmax=max_constraint)
	fig = plt.figure()
	line_segments = []
	ax = fig.gca()
	ax.set(xlim=(0., 2.0), ylim=(0., 1.1))
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




######################## MAIN ##################################

def sorted_ls(path):
    mtime = lambda f: os.stat(os.path.join(path, f)).st_mtime
    return list(sorted(os.listdir(path), key=mtime))


if __name__ == '__main__':
	os.chdir('../Data/growth_network/')
	if len(sys.argv) != 1:
		os.chdir(sys.argv[1])
	else:
		os.chdir(sorted_ls('.')[-1])
	filenames=fnmatch.filter(os.listdir('.'), 'parameters_*.csv')
	for filename in filenames:
		network = load_network_info(int(filename[-7:-4]),step=0)
		plot_geometry(network)
		plt.savefig('network_%03d.pdf' % int(filename[-7:-4]))
		plt.show()
	"""
	type_plot = str(input("What graph do you want?\n enter 'geo' for geometry, 'const' for constraints"))

	if type_plot == 'geo':
		print 'Your geometry'
		plot_geometry(network)
		#plot_network_geometry(len(os.listdir('.'))-5)
		plt.show()
	elif type_plot == 'const':
		print 'Constraints'
		plot_constraints(network)
		plt.show()"""
