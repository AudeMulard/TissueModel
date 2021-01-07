import matplotlib.pyplot as plt
import sys
sys.path.append('/home/am2548/Documents/TissueModel/')
sys.path.append('/users/am2548/TissueModel/Documents/PhD/Code/TissueModel')
#from Plotting.network_plotting import *
from Plotting.information_network import *


def sorted_ls(path):
    mtime = lambda f: os.stat(os.path.join(path, f)).st_mtime
    return list(sorted(os.listdir(path), key=mtime))

os.chdir('../Data/Study_networks/Oct-13-2020_0021')
"""if len(sys.argv) == 1:
	os.chdir(sorted_ls('.')[-1])
else:
	os.chdir(sys.argv[1])
print(os.getcwd())
"""
list_modes = [['Voronoi','random'],['Voronoi','grid'],['growth_network','grid']]#,['Voronoi','regular']]

fig = plt.figure()
ax1 = fig.add_subplot(221)
ax2 = fig.add_subplot(222, sharex=ax1, sharey=ax1)
ax3 = fig.add_subplot(223, sharex=ax1, sharey=ax1)
ax4 = fig.add_subplot(224, sharex=ax1, sharey=ax1)

line_segments=[]
#test_number = fnmatch.filter(sorted_ls('.'), 'network_vertices_01_00_*.rpt')[i][-13:-4]
#creation , generation = list_modes[0]
filename = fnmatch.filter(os.listdir('.'), 'network_ridge_vertices_*.csv')

network = load_network_info(int(filename[0][-7:-4]))
from matplotlib.collections import LineCollection
ax1.scatter(network.vertices[:,0],network.vertices[:,1], s =2.)
for simplex in network.ridge_vertices:
        simplex = np.asarray(simplex)
        line_segments.append([(x, y) for x, y in network.vertices[simplex]])
lc = LineCollection(line_segments,linestyle='solid')
ax1.add_collection(lc)
ax1.set_title('%s, %s' % (network.creation, network.generation))

line_segments=[]

network = load_network_info(int(filename[1][-7:-4]))
from matplotlib.collections import LineCollection
ax2.scatter(network.vertices[:,0],network.vertices[:,1], s =2.)
for simplex in network.ridge_vertices:
        simplex = np.asarray(simplex)
        line_segments.append([(x, y) for x, y in network.vertices[simplex]])
lc = LineCollection(line_segments,linestyle='solid')
ax2.add_collection(lc)
ax2.set_title('%s, %s' % (network.creation, network.generation))

line_segments=[]

network = load_network_info(int(filename[2][-7:-4]))
from matplotlib.collections import LineCollection
ax3.scatter(network.vertices[:,0],network.vertices[:,1], s =2.)
for simplex in network.ridge_vertices:
        simplex = np.asarray(simplex)
        line_segments.append([(x, y) for x, y in network.vertices[simplex]])
lc = LineCollection(line_segments,linestyle='solid')
ax3.add_collection(lc)
ax3.set_title('%s, %s' % (network.creation, network.generation))


line_segments=[]

network = load_network_info(int(filename[3][-7:-4]))
from matplotlib.collections import LineCollection
ax4.scatter(network.vertices[:,0],network.vertices[:,1], s =2.)
for simplex in network.ridge_vertices:
        simplex = np.asarray(simplex)
        line_segments.append([(x, y) for x, y in network.vertices[simplex]])
lc = LineCollection(line_segments,linestyle='solid')
ax4.add_collection(lc)
ax4.set_title('%s, %s' % (network.creation, network.generation))


plt.subplots_adjust(hspace=0.5)
plt.savefig('Geometries.png')
