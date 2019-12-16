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
space_discretization = 0.002*length_domain
Ef=1.0
A=1.4E-8
B=3.8
traction_distance = 0.1*length_domain
#iteration = 15

current_dir = os.getcwd()

data_path = '../Data/study_96/'

today = date.today()
new_dir = data_path+today.strftime("%b-%d-%Y")+'_'+'%04d' % len(os.listdir(data_path))
os.mkdir(new_dir)
path = new_dir


## EXPERIMENT
creation="Voronoi"
constitutive = 'linear2'
scheme='nonlinear'
side = 'right'
plot = True
video = True
phase = 'traction'
stress_rep = True


def sorted_ls(path):
    mtime = lambda f: os.stat(os.path.join(path, f)).st_mtime
    return list(sorted(os.listdir(path), key=mtime))

os.chdir('../Data/default/Dec-03-2019_0096/')

network = Network(dimension, complexity_network, length_domain, min_distance, Ef, A, B, creation, path)

with open('network_vertices_initial.csv','r') as readFile:
	reader = csv.reader(readFile)
	list_vertices = np.array(list(reader))
	network.vertices=list_vertices.astype(float)
with open('network_ridge_vertices.csv','r') as readFile:
	reader = csv.reader(readFile)
	list_ridge_vertices=np.array(list(reader))
	network.ridge_vertices=list_ridge_vertices.astype(int)

os.chdir(current_dir)

network = network.sort_nodes()
network.vertices_ini = np.array(network.vertices.tolist())
network = network.distribution_length_fiber(path)
network = network.create_ridge_node_list()
network.save_network('initial', path)

test_1 = Tensile_test(constitutive, scheme, side, space_discretization, traction_distance, plot, video, phase, path)

network = test_1.full_test(network, path)

test_2 = Tensile_test(constitutive, scheme, side, space_discretization, -traction_distance, plot, video, phase, path)

network = test_2.full_test(network, path)
