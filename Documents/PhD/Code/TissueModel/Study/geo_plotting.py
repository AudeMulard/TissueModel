import matplotlib.pyplot as plt
import sys
sys.path.append('/home/aude/Documents/PhD/Code/TissueModel/')
from Plotting.network_plotting import *
from Plotting.information_network import *


def sorted_ls(path):
    mtime = lambda f: os.stat(os.path.join(path, f)).st_mtime
    return list(sorted(os.listdir(path), key=mtime))

os.chdir('../Data/Study_networks/')
if len(sys.argv) == 1:
	os.chdir(sorted_ls('.')[-1])
else:
	os.chdir(sys.argv[1])
print(os.getcwd())

list_modes = [['Voronoi','random'],['Voronoi','regular'],['Voronoi','grid'],['growth_network','grid']]

fig = plt.figure()
ax1 = fig.add_subplot(221)
ax2 = fig.add_subplot(222, sharex=ax1, sharey=ax1)
ax3 = fig.add_subplot(223, sharex=ax1, sharey=ax1)
ax4 = fig.add_subplot(224, sharex=ax1, sharey=ax1)

line_segments=[]
creation , generation = list_modes[0]
filename = fnmatch.filter(os.listdir('.'), 'network_vertices_*_%s_%s.csv' % (creation,generation))

network = load_network_info(int(filename[0][25:28]))
from matplotlib.collections import LineCollection
ax1.scatter(network.vertices[:,0],network.vertices[:,1], s =2.)
for simplex in network.ridge_vertices:
        simplex = np.asarray(simplex)
        line_segments.append([(x, y) for x, y in network.vertices[simplex]])
lc = LineCollection(line_segments,linestyle='solid')
ax1.add_collection(lc)
ax1.set_title('Random Voronoi')

line_segments=[]
creation , generation = list_modes[1]
filename = fnmatch.filter(os.listdir('.'), 'network_vertices_*_%s_%s.csv' % (creation,generation))

network = load_network_info(int(filename[0][25:28]))
from matplotlib.collections import LineCollection
ax2.scatter(network.vertices[:,0],network.vertices[:,1], s =2.)
for simplex in network.ridge_vertices:
        simplex = np.asarray(simplex)
        line_segments.append([(x, y) for x, y in network.vertices[simplex]])
lc = LineCollection(line_segments,linestyle='solid')
ax2.add_collection(lc)
ax2.set_title('Hexagonal network')
line_segments=[]
creation , generation = list_modes[2]
filename = fnmatch.filter(os.listdir('.'), 'network_vertices_*_%s_%s.csv' % (creation,generation))

network = load_network_info(int(filename[0][25:28]))
from matplotlib.collections import LineCollection
ax3.scatter(network.vertices[:,0],network.vertices[:,1], s =2.)
for simplex in network.ridge_vertices:
        simplex = np.asarray(simplex)
        line_segments.append([(x, y) for x, y in network.vertices[simplex]])
lc = LineCollection(line_segments,linestyle='solid')
ax3.add_collection(lc)
ax3.set_title('Disturbed grid')
line_segments=[]
creation , generation = list_modes[3]
filename = fnmatch.filter(os.listdir('.'), 'network_vertices_*_%s_%s.csv' % (creation,generation))

network = load_network_info(int(filename[0][25:28]))
from matplotlib.collections import LineCollection
ax4.scatter(network.vertices[:,0],network.vertices[:,1], s =2.)
for simplex in network.ridge_vertices:
        simplex = np.asarray(simplex)
        line_segments.append([(x, y) for x, y in network.vertices[simplex]])
lc = LineCollection(line_segments,linestyle='solid')
ax4.add_collection(lc)
ax4.set_title('Growth network')


plt.subplots_adjust(hspace=0.5)
plt.savefig('Geometries.pdf')
