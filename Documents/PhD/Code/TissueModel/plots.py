import numpy as np
import matplotlib.pyplot as plt
from matplotlib_scalebar.scalebar import ScaleBar
import csv
import matplotlib.font_manager as fm
import sys
sys.path.insert(0,'../../../TissueModel/')
import os
from mpl_toolkits.axes_grid1 import make_axes_locatable
import seaborn as sns
#from creation_network import Network
from matplotlib.collections import LineCollection
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes, mark_inset
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
from matplotlib import cm
import matplotlib.colors 




def sorted_ls(path):
    mtime = lambda f: os.stat(os.path.join(path, f)).st_mtime
    return list(sorted(os.listdir(path), key=mtime))



def length_square(x):	
	return x[0]**2+x[1]**2

with open('network_vertices_initial.csv','r') as readFile:
	reader = csv.reader(readFile)
	list_vertices = np.array(list(reader))
	vertices_ini=list_vertices.astype(float)

def plot_network_geometry(step,**kw):
	with open('network_vertices_initial_59.csv','r') as readFile:
		reader = csv.reader(readFile)
		list_vertices = np.array(list(reader))
		vertices=list_vertices.astype(float)
	with open('network_ridge_vertices_59.csv','r') as readFile:
		reader = csv.reader(readFile)
		list_ridge_vertices=np.array(list(reader))
		ridge_vertices=list_ridge_vertices.astype(int)
	fig = plt.figure()
	ax = fig.add_subplot(111)
	try:
		with open('network_vertices_%03d.csv' % int(step) ,'r') as readFile:
			reader = csv.reader(readFile)
			list_vertices = np.array(list(reader))
			vertices=list_vertices.astype(float)
	except:
		with open('network_vertices_initial.csv' ,'r') as readFile:
			reader = csv.reader(readFile)
			list_vertices = np.array(list(reader))
			vertices=list_vertices.astype(float)
	line_segments = []
	#ax.set(xlim=(0., max(vertices[:,0])*1.1), ylim=(0., max(vertices[:,1])*1.1))
	if kw.get('show_vertices', True):
		ax.scatter(vertices[:,0],vertices[:,1], s =2.)
	for i in range(len(vertices[:,0])):
		ax.annotate(i, (vertices[i,0],vertices[i,1]),fontsize=5)
	for simplex in ridge_vertices:
	        simplex = np.asarray(simplex)
            	line_segments.append([(x, y) for x, y in vertices[simplex]])
	lc = LineCollection(line_segments,linestyle='solid')
	ax.add_collection(lc)	
	return ax.figure

def plot_network_constraints(step,**kw):
	with open('network_vertices_initial_59.csv','r') as readFile:
		reader = csv.reader(readFile)
		list_vertices = np.array(list(reader))
		vertices_ini=list_vertices.astype(float)
	with open('network_ridge_vertices_59.csv','r') as readFile:
		reader = csv.reader(readFile)
		list_ridge_vertices=np.array(list(reader))
		ridge_vertices=list_ridge_vertices.astype(int)
	fig = plt.figure()
	ax = fig.add_subplot(111)
	try:
		with open('network_vertices_%03d.csv' % int(step) ,'r') as readFile:
			reader = csv.reader(readFile)
			list_vertices = np.array(list(reader))
			vertices=list_vertices.astype(float)
	except:
		with open('network_vertices_initial.csv' ,'r') as readFile:
			reader = csv.reader(readFile)
			list_vertices = np.array(list(reader))
			vertices=list_vertices.astype(float)
	# Creation color map
	ridge_constraints = []
	for ridge in ridge_vertices:
		ridge_constraints.append((np.sqrt(length_square(vertices[ridge[0]]-vertices[ridge[1]]))-np.sqrt(length_square(vertices_ini[ridge[0]]-vertices_ini[ridge[1]]))))
	cmap = plt.cm.rainbow
	norm = matplotlib.colors.Normalize(vmin=0., vmax=max(ridge_constraints))
	line_segments = []
	for i in range(len(vertices[:,0])):
		ax.annotate(i, (vertices[i,0],vertices[i,1]),fontsize=5)
	ax.set(xlim=(0., max(vertices[:,0])*1.1), ylim=(0., max(vertices[:,1])*1.1))
	for simplex in ridge_vertices:
	        simplex = np.asarray(simplex)
            	line_segments.append([(x, y) for x, y in vertices[simplex]])
	lc = LineCollection(line_segments,linestyle='solid',cmap=cmap,norm=norm)
	sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
	lc.set_array(np.array(ridge_constraints))
	ax.add_collection(lc)
	sm.set_array([])
	divider = make_axes_locatable(ax)
	cax = divider.append_axes('right', size='5%', pad=0.05)
	fig.colorbar(sm, cax=cax, orientation='vertical')


if __name__ == '__main__':
	os.chdir('../Data/default/')
	if len(sys.argv) != 1:
		os.chdir(sys.argv[1])
	else:
		os.chdir(sorted_ls('.')[-1])
	folder = str(input("What graph do you want?\n enter 'geo' for geometry, 'const' for constraints"))
	if folder == 'geo':
		print 'Your geometry'
		print len(os.listdir('.'))-5
		plot_network_geometry(len(os.listdir('.'))-5)
		plt.show()
	elif folder == 'const':
		print 'Constraints'
		plot_network_constraints(len(os.listdir('.'))-5)
		plt.show()




