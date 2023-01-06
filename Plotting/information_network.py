import os
import csv
import numpy as np
import sys
import fnmatch
import math
sys.path.append('C:/Users/am2548/TissueModel/Documents/PhD/Code/TissueModel/')
import Network_generation.creation_network
import Core_calculation.tensile_test


def sorted_ls(path):
    mtime = lambda f: os.stat(os.path.join(path, f)).st_ctime
    return list(sorted(os.listdir(path), key=mtime))

def load_network_info(path,**kw):
	filename = fnmatch.filter(os.listdir('.'), 'parameters_??_%05d.csv' % path)
	with open(filename[-1],'r') as readFile:
		reader = csv.reader(readFile)
		params = list(reader)
		params = dict([param for param in params if param !=[]])
	network = Network_generation.creation_network.Network(int(params['dimension']), 
															params['complexity'], 
															eval(params['length']), 
															eval(params['merge_distance']), 
															eval(params["beam Young"]),
															eval(params["beam Poisson"]),
															eval(params["beam Profile"]),
															eval(params["connector coeff"]),
															params['hyperstatic_param'], 
															params['creation'],
															params['generation'],
															path)
	filename = fnmatch.filter(os.listdir('.'),'network_vertices_initial_*_%05d.csv' % path)
	with open(filename[-1], 'r') as readFile:
		reader = csv.reader(readFile)
		list_vertices = list(reader)
		list_vertices = np.array([vertice for vertice in list_vertices if vertice !=[]],dtype=object)
		network.vertices=list_vertices.astype(float)
	filename = fnmatch.filter(os.listdir('.'), 'network_ridge_vertices_*_%05d.csv' % path)
	with open(filename[-1], 'r') as readFile:
		reader = csv.reader(readFile)
		list_ridge_vertices=list(reader)
		list_ridge_vertices=np.array([ridge for ridge in list_ridge_vertices if ridge !=[]],dtype=object)
		network.ridge_vertices=list_ridge_vertices.astype(int)
	network.ridge_vertices = [list(l) for l in network.ridge_vertices]
	network=network.sort_nodes()
	network=network.create_ridge_node_list()
	return network

def write_node_label(test_number,network):
	import pandas as pd
	list_node_label = []
	node_label = pd.read_csv('node_label_%09d.csv' % int(test_number),header=None).to_numpy()[0]
	number_elements = pd.read_csv('number_elements_%09d.csv' % int(test_number),header=None).to_numpy()[0]
	number_elements=[number+1 for number in number_elements]
	network = network.create_ridge_node_list()
	for k in range(len(node_label)):
		ridge=network.list_nodes_ridges[k]
		if node_label[k]==0:
			list_node_label.append(sum(number_elements[:(ridge[0]+1)-1]))
		elif node_label[k]==1:
			list_node_label.append(sum(number_elements[:(ridge[0]+1)])-1)
	return list_node_label

def calculate_strain(network,test_number):
	import pandas as pd
	list_node_label = write_node_label(test_number,network)
	strain=[]
	for k in range(51):
		filename = 'network_vertices_01_%02d_%09d.csv' % (k,test_number)
		if network.dimension ==2:
			df = pd.read_csv(filename,usecols=[4,12,13])
			network.vertices_end=df.loc[:,['   COORD-COOR1','   COORD-COOR2']].to_numpy()
		elif network.dimension==3:
			df = pd.read_csv(filename,usecols=[4,12,13,14])
			network.vertices_end=df.loc[:,['   COORD-COOR1','   COORD-COOR2','   COORD-COOR3']].to_numpy()
		strain.append(np.nanmax(network.vertices_end[:,0]-1.0))
	return strain

def load_info_test(path):
	filename = fnmatch.filter(sorted_ls('.'), 'parameters_test_*_%05d.csv' % path)
	print(filename)
	with open(filename[-1],'r') as readFile:
		reader = csv.reader(readFile)
		params = list(reader)
		params = dict([param for param in params if param !=[]])
	test = Core_calculation.tensile_test.Tensile_test(params["side"],eval(params["space discretization"]),
		eval(params["traction_distance"]),eval(params["element size"]),'.')
	return test

def stress_strain_curve(test_number,network):
	import pandas as pd
	filenames = fnmatch.filter(os.listdir('.'), 'stress_data_*_%09d.rpt' % int(test_number))
	stress=[]
	lengths = []
	for ridge in network.ridge_vertices:
		length = np.sqrt(length_square(network.vertices[ridge[0]]-network.vertices[ridge[1]]))
		lengths.append(length)
	df_1 = pd.read_fwf(filenames[0],skiprows=2,skipfooter=4,colspecs=([0,30],[31,50]), sep='\t').to_numpy()
	stress_1 = df_1[:,1]
	for i in range(len(stress_1)):
		stress.append(stress_1[i]/(np.mean(lengths)*len(network.ridge_vertices))/(network.beam_young*network.beam_profile**2*np.pi))
	return stress


def plot_second_der(ax,strain, stress,color):
	der_sec = []
	for i in range(len(strain)-1):
		der_sec.append(2*stress[i-1]/[(strain[i]-strain[i-1])*(strain[i+1]-strain[i-1])]-2*stress[i]/[(strain[i+1]-strain[i])*(strain[i]-strain[i-1])]+2*stress[i+1]/[(strain[i+1]-strain[i])*(strain[i+1]-strain[i-1])])
	ax.plot(strain[:-1],der_sec, color = color,linestyle='dashed')

def plot_first_der(ax,strain, stress,color):
	der_sec = []
	for i in range(1,len(stress)-1):
		der_sec.append((stress[i+1]-stress[i-1])/(strain[i+1]-strain[i-1]))
	ax.plot(strain[1:-1],der_sec, color = color,linestyle='dashed')

def length_square(x):
	if len(x) ==2: return x[0]**2+x[1]**2
	if len(x) ==3: return x[0]**2+x[1]**2+x[2]**2
