import numpy as np
import matplotlib.pyplot as plt
from matplotlib_scalebar.scalebar import ScaleBar
import csv
import matplotlib.font_manager as fm
import sys
sys.path.insert(0,'../../../TissueModel/')
import os

import seaborn as sns
#from creation_network import Network
from matplotlib.collections import LineCollection
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes, mark_inset
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar

os.chdir('../Data/default/')
os.chdir(os.listdir('.')[0])

def plot_network_geometry(step,**kw):
	fig = plt.figure()
	ax = fig.add_subplot(111)
	with open('network_vertices_%s.csv' % step ,'r') as readFile:
		reader = csv.reader(readFile)
		list_vertices = np.array(list(reader))
		vertices=list_vertices.astype(float)
	line_segments = []
	ax.set(xlim=(0., max(vertices[:,0])*1.1), ylim=(0., max(vertices[:,1])*1.1))
	if kw.get('show_vertices', True):
		ax.scatter(vertices[:,0],vertices[:,1], s =2.)
	for simplex in ridge_vertices:
	        simplex = np.asarray(simplex)
            	line_segments.append([(x, y) for x, y in vertices[simplex]])
	lc = LineCollection(line_segments,linestyle='solid')
	ax.add_collection(lc)	
	return ax.figure

with open('network_vertices_initial.csv','r') as readFile:
	reader = csv.reader(readFile)
	list_vertices = np.array(list(reader))
	vertices=list_vertices.astype(float)
with open('network_ridge_vertices.csv','r') as readFile:
	reader = csv.reader(readFile)
	list_ridge_vertices=np.array(list(reader))
	ridge_vertices=list_ridge_vertices.astype(int)

"""
def distribution_length_fiber():
	lengths = []
	for ridge in ridge_vertices:
		lengths.append(np.linalg.norm(vertices[ridge[0]]-vertices[ridge[1]]))
	mean_length = np.mean(lengths)
	return mean_length

mean_length = distribution_length_fiber()
"""
#print os.listdir('.')
#plot_network_geometry(len(os.listdir('.'))-5)
#plt.show()

"""
fig.tight_layout()

mean_length = distribution_length_fiber(ax2)
plot_network_SumSch('initial',mean_length, ax3)


#plt.savefig('initial.pdf',bbox_inches='tight')

plot_network_SumSch(5,mean_length, ax4)

#plt.savefig('step_5.pdf',bbox_inches='tight')

plot_network_SumSch(10,mean_length, ax5)

#plt.savefig('step_10.pdf',bbox_inches='tight')

plot_network_SumSch(15,mean_length, ax6)
#plt.savefig('step_15.pdf',bbox_inches='tight')

plot_stress_strain(ax1)
plt.subplots_adjust(hspace = 0.6)


#plt.show()
"""

