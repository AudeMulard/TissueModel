from creation_network import Network
from force_balance import *
from tensile_test import *
from scipy.spatial import Voronoi, voronoi_plot_2d
import matplotlib.pyplot as plt
from Several_networks_test import *
from network_plotting import *
import os
from datetime import date
import datetime
from optimize_geometry import energy_minimum
import sys
from stress_strain import plot_second_der
import random


## PARAMETERS
dimension=2 #dimension of the problem
complexity_network=150 #number of random seed points
length_domain=1.0
min_distance = 0.02*length_domain
space_discretization = 0.0025*length_domain
Ef=1.0
A=1.4E-8
B=3.8
traction_distance = 0.2*length_domain
#iteration = 15


"""
data_path = '../Data/influence_points_curves/'

today = date.today()
new_dir = data_path+today.strftime("%b-%d-%Y")+'_'+'%04d' % len(os.listdir(data_path))
os.mkdir(new_dir)
path = new_dir

"""
## EXPERIMENT
creation="Voronoi"
constitutive = 'linear2'
scheme='nonlinear'
side = 'right'
plot = True
video = True
phase = 'traction'
stress_rep = True

"""

for k in [10,30,50,75,90,100, 110,120,130,140,150]:
	print 'Complexity = ' , k
	try:
		network = Network(dimension, k, length_domain, min_distance, Ef, A, B, creation, path)
		network = network.set_fibers(creation, path)
		test_1 = Tensile_test(constitutive, scheme, side, space_discretization, traction_distance, plot, video, phase, path)
		network = test_1.full_test(network, path)
	except (numpy.linalg.linalg.LinAlgError, ValueError):
			print('Singular Matrix')
			pass

"""
def sorted_ls(path):
    mtime = lambda f: os.stat(os.path.join(path, f)).st_mtime
    return list(sorted(os.listdir(path), key=mtime))

os.chdir('../Data/influence_points_curves/')

def random_color():
    return numpy.random.rand(3,)

if len(sys.argv) != 1:
	os.chdir(sys.argv[1])
else:
	os.chdir(sorted_ls('.')[-1])

fig, ax1 = plt.subplots()
ax2 =ax1.twinx()
labels=[]
print sorted_ls('.')
for filename in sorted_ls('.'):
	if 'stress_strain' in filename:
		with open(filename, 'r') as readFile:
			reader = csv.reader(readFile)
			curve = np.array(list(reader))
			strain = [float(i) for i in curve[0]]
			stress = [float(i) for i in curve[1]]
		color = random_color()
		ax1.plot(strain, stress, marker='o',linestyle='dashed', markersize=5.,color=color)
		plot_second_der(ax2,strain, stress,color)
		labels.append(filename[14:18])
ax1.set_ylabel(r'$\frac{\sigma}{k/\bar{l_0}}$', fontsize=30)

ax2.set_ylabel('second derivative')
plt.xlabel(r'$\frac{L-L_0}{L_0}$',fontsize=30)
plt.legend(labels)
plt.tight_layout()
plt.savefig('Curves_number_points.pdf')


fig = plt.figure()
labels=[]
for filename in sorted_ls('.'):
	if 'stress_strain' in filename:
		with open(filename, 'r') as readFile:
			reader = csv.reader(readFile)
			curve = np.array(list(reader))
			strain = [float(i) for i in curve[0]]
			stress = []
			for i in range(len(curve[1])):
				stress.append(float(curve[1][i])/(1.+strain[i]))
		plt.plot(strain, stress, marker='o',linestyle='dashed', markersize=5.)
		labels.append(filename[14:18])


ax1.set_ylabel(r'$\frac{\sigma}{k/\bar{l_0}}$', fontsize=30)

ax2.set_ylabel('second derivative')
plt.xlabel(r'$\frac{L-L_0}{L_0}$',fontsize=30)
plt.legend(labels)
plt.tight_layout()
plt.savefig('Curves_number_points_lagrangian.pdf')

