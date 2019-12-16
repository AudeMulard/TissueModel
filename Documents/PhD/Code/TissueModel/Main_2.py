from creation_network import Network
from force_balance import *
from tensile_test import *
from scipy.spatial import Voronoi, voronoi_plot_2d
import matplotlib.pyplot as plt
from Several_networks_test import *
from network_plotting import *
import os
from datetime import date

## PARAMETERS
dimension=2 #dimension of the problem
complexity_network=30 #number of random seed points
length_domain=1.0
min_distance = 0.02*length_domain
space_discretization = 0.0025*length_domain
Ef=1.0
A=1.4E-8
B=3.8
traction_distance = 0.2*length_domain
#iteration = 15



data_path = '../Data/default/'

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

'''
while True:
	folder = str(input("Did you change the folder? enter 'y' for yes: "))
	if folder != "y":
		print('Exit and change the folder')
	else:
		break
'''

## One tensile_test



network = Network(dimension, complexity_network, length_domain, min_distance, Ef, A, B, creation, path)


network = network.set_fibers(creation, path)

print len(network.vertices)
test_1 = Tensile_test(constitutive, scheme, side, space_discretization, traction_distance, plot, video, phase, path)

network = test_1.full_test(network, path)



