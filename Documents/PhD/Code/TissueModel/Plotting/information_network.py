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
    mtime = lambda f: os.stat(os.path.join(path, f)).st_mtime
    return list(sorted(os.listdir(path), key=mtime))

def load_network_info(path,**kw):
	filename = fnmatch.filter(os.listdir('.'), 'parameters_??_%05d.csv' % path)
	#print(path,filename)
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
															params['disturbance'],
															params['hyperstatic_param'], 
															params['creation'],
															params['generation'],
															path)
	filename = fnmatch.filter(os.listdir('.'),'network_vertices_initial_*_%05d.csv' % path)
	with open(filename[0], 'r') as readFile:
		reader = csv.reader(readFile)
		list_vertices = list(reader)
		list_vertices = np.array([vertice for vertice in list_vertices if vertice !=[]],dtype=object)
		network.vertices_ini=list_vertices.astype(float)
	filename = fnmatch.filter(os.listdir('.'), 'network_ridge_vertices_*_%05d.csv' % path)
	with open(filename[0], 'r') as readFile:
		reader = csv.reader(readFile)
		list_ridge_vertices=list(reader)
		list_ridge_vertices=np.array([ridge for ridge in list_ridge_vertices if ridge !=[]],dtype=object)
		network.ridge_vertices=list_ridge_vertices.astype(int)
	network.ridge_vertices = [list(l) for l in network.ridge_vertices]
	if kw.get('step'): 
		try:
			with open('network_vertices_%03d_%03d.csv' % (int(kw['step']),int(params['number of nodes'])) ,'r') as readFile:
				reader = csv.reader(readFile)
				list_vertices = np.array(list(reader))
				network.vertices=list_vertices.astype(float)
		except:
			print('No info on step, Initial step')
			network.vertices = network.vertices_ini
			pass
	else:
		network.vertices = network.vertices_ini
	network=network.sort_nodes()
	network=network.create_ridge_node_list()
	return network

def write_node_label(test_number,network):
	import pandas as pd
	filename = 'nodes_label_%09d.csv' % int(test_number)
	list_node_label = []
	node_label = pd.read_csv('node_label_%09d.csv' % int(test_number),header=None).to_numpy()[0]
	number_elements = pd.read_csv('number_elements_%s.csv' % test_number,header=None).to_numpy()[0]
	number_elements=[number+1 for number in number_elements]
	for k in range(len(node_label)):
		ridge=network.list_nodes_ridges[k]
		if node_label[k]==0:
			list_node_label.append(sum(number_elements[:(ridge[0]+1)-1]))
		elif node_label[k]==1:
			list_node_label.append(sum(number_elements[:(ridge[0]+1)])-1)
	return list_node_label


def load_info_test(path):
	filename = fnmatch.filter(os.listdir('.'), 'parameters_test_%05d.csv' % path)
	#print(filename)
	with open(filename[0],'r') as readFile:
		reader = csv.reader(readFile)
		params = list(reader)
		params = dict([param for param in params if param !=[]])
	test = Core_calculation.tensile_test.Tensile_test(params["constitutive"],params["side"],eval(params["space discretization"]),
		eval(params["traction_distance"]),eval(params["element size"]),params["plot"], params["video"],'.',params["details"])
	return test

def stress_strain_curve(test_number,network):
	import pandas as pd
	filenames = fnmatch.filter(os.listdir('.'), 'stress_data_*_%09d.rpt' % int(test_number))
	stress=[]
	lengths = []
	for ridge in network.ridge_vertices:
		# Minimum length
		length = np.sqrt(length_square(network.vertices[ridge[0]]-network.vertices[ridge[1]]))
		lengths.append(length)
	df_1 = pd.read_fwf(filenames[0],skiprows=2,skipfooter=4,colspecs=([0,30],[31,50]), sep='\t').to_numpy()
	stress_1 = df_1[:,1]
	strain=[]
	for i in range(len(stress_1)):
		#stress.append(stress_1[i]/(len(network.boundary_nodes_right)+len(network.boundary_nodes_left))/(network.beam_young*network.beam_profile**2*np.pi))
		stress.append(stress_1[i]/(np.mean(lengths)*len(network.ridge_vertices))/(network.beam_young*network.beam_profile**2*np.pi))#(len(network.boundary_nodes_right)+len(network.boundary_nodes_left)))#np.mean(lengths)))#*len(network.ridge_vertices)))#/(len(network.ridge_vertices)*np.mean(lengths)))
		strain.append(df_1[i,0]*0.5)
	#print(test_number,strain,stress)
	return strain, stress


def plot_second_der(ax,strain, stress,color):
	der_sec = []
	#print(strain,stress)
	for i in range(len(strain)-1):
		der_sec.append(2*stress[i-1]/[(strain[i]-strain[i-1])*(strain[i+1]-strain[i-1])]-2*stress[i]/[(strain[i+1]-strain[i])*(strain[i]-strain[i-1])]+2*stress[i+1]/[(strain[i+1]-strain[i])*(strain[i+1]-strain[i-1])])
	#print(der_sec)
	ax.plot(strain[:-1],der_sec, color = color,linestyle='dashed')

def plot_first_der(ax,strain, stress,color):
	der_sec = []
	for i in range(1,len(stress)-1):
		der_sec.append((stress[i+1]-stress[i-1])/(strain[i+1]-strain[i-1]))
	ax.plot(strain[1:-1],der_sec, color = color,linestyle='dashed')

def plot_stress_strain(k,ax1,ax2,color):
	filename = fnmatch.filter(os.listdir('.'), 'stress_strain_%03d.csv' % k)
	with open(filename[0], 'r') as readFile:
		reader = csv.reader(readFile)
		curve = np.array(list(reader))
		stress = [float(i) for i in curve[1]]
		strain=[float(i) for i in curve[0]]

	ax1.plot(strain, stress, marker='o',linestyle='dashed', markersize=5., color = color,label=k)
	return strain

def length_square(x):
	if len(x) ==2: return x[0]**2+x[1]**2
	if len(x) ==3: return x[0]**2+x[1]**2+x[2]**2


def test_isotropy(network):
	import matplotlib.pyplot as plt
	fig = plt.figure()
	ax1 = fig.add_subplot()
	angles,lengths=[],[]
	for ridge in network.ridge_vertices:
		if 0.1*network.length[0]<network.vertices[ridge[0]][0]<0.9*network.length[0] and 0.1*network.length[0]<network.vertices[ridge[1]][0]<0.9*network.length[0]:
			if 0.1*network.length[1]<network.vertices[ridge[0]][1]<0.9*network.length[1] and 0.1*network.length[1]<network.vertices[ridge[1]][1]<0.9*network.length[1]:
				lengths.append(np.sqrt(length_square(network.vertices[ridge[0]]-network.vertices[ridge[1]])))
				angles.append(np.arccos(np.dot(network.vertices[ridge[0]]-network.vertices[ridge[1]], [1.,0.])/lengths[-1])*180/np.pi)
	import seaborn as sns
	plt.hist(angles,weights=lengths)
	fig1 = plt.figure()
	plt.hist(angles)
	#ax1.scatter(angles,lengths)
	#print(angles)
				

def calculate_network_data(network, length_ini, **kw):
	import pandas as pd
	if kw.get('step')!=None:
		if kw.get('test_number')!=None:
			filename = 'network_vertices_01_%02d_%s.csv' % (int(kw['step']),kw['test_number'])
		else:
			filename = 'network_vertices_%03d_%03d.csv' % (int(kw['step']),len(network.vertices))
		df = pd.read_csv(filename,usecols=[4,12,13])
		if network.dimension==3:
			network.vertices=df.reindex(list(kw['node_label'].transpose()[0])).loc[:,['   COORD-COOR1','   COORD-COOR2','   COORD-COOR3']].to_numpy()
		elif network.dimension ==2:
			network.vertices=df.reindex(list(kw['node_label'])).loc[:,['   COORD-COOR1','   COORD-COOR2']].to_numpy()
		from Plotting.network_plotting import plot_geometry
		import matplotlib.pyplot as plt
	strain_omega_step=(np.max(network.vertices)-1.0)/1.0
	stretch_ratios,cos_theta_square = [],[]
	list_angle = []
	lengths=[]
	#print(len(network.vertices))
	for ridge in network.ridge_vertices:
		# Minimum length
		length = np.sqrt(length_square(network.vertices[ridge[0]]-network.vertices[ridge[1]]))
		lengths.append(length)		
		stretch_ratios.append(length / length_ini[list(network.ridge_vertices).index(ridge)])
		# Minimum & max angle
		try:
			if network.dimension==2: cos_theta_square.append(np.dot(network.vertices[ridge[0]]-network.vertices[ridge[1]], [1.,0.])**2/length**2)
			if network.dimension==3: cos_theta_square.append(np.dot(network.vertices[ridge[0]]-network.vertices[ridge[1]], [1.,0.,0])**2/length**2)
		except RuntimeWarning:
			continue
		list_angle.append(cos_theta_square[-1]*length)
	omega_xx = sum(list_angle)/sum(lengths)
	return lengths, stretch_ratios, cos_theta_square, omega_xx,strain_omega_step


