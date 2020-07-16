import os
import csv
import numpy as np
import matplotlib.pyplot as plt
import sys
import fnmatch
import math
from statistics import mean
sys.path.append('/home/aude/Documents/PhD/Code/TissueModel/')
import Network_generation.creation_network

def sorted_ls(path):
    mtime = lambda f: os.stat(os.path.join(path, f)).st_mtime
    return list(sorted(os.listdir(path), key=mtime))

def load_network_info(path,**kw):
	with open('parameters_%03d.csv' % path, 'r') as readFile:
		reader = csv.reader(readFile)
		params = dict(list(reader))
	network = Network_generation.creation_network.Network(int(params['dimension']), params['complexity'], eval(params['length']), params['merge_distance'], eval(params['k_tension']), eval(params['k_compression']), params['A'], params['disturbance'], params['creation'],params['generation'],path)
	filename = fnmatch.filter(os.listdir('.'),'network_vertices_initial_%03d*.csv' % path)
	with open(filename[0], 'r') as readFile:
		reader = csv.reader(readFile)
		list_vertices = np.array(list(reader))
		network.vertices_ini=list_vertices.astype(float)
	filename = fnmatch.filter(os.listdir('.'), 'network_ridge_vertices_*.csv')
	with open('network_ridge_vertices_%03d.csv' % path, 'r') as readFile:
		reader = csv.reader(readFile)
		list_ridge_vertices=np.array(list(reader))
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


def plot_second_der(ax,strain, stress,color):
	der_sec = [0]
	for i in range(1,len(stress)-1):
		der_sec.append((stress[i+1]+stress[i-1]-2*stress[i])/(strain[i+1]-strain[i])**2)
		#print stress[i+1],stress[i],stress[i-1]
		#der_sec.append((stress[i]))
	ax.plot(strain[0:len(der_sec)],der_sec, color = color)

def plot_first_der(ax,strain, stress,color):
	der_sec = []
	for i in range(0,len(stress)-1):
		der_sec.append((stress[i+1]-stress[i])/(strain[i+1]-strain[i]))
	ax.plot(strain[0:len(der_sec)],der_sec, color = color)

def plot_stress_strain(k,ax1,ax2,color):
	filename = fnmatch.filter(os.listdir('.'), 'stress_strain_%03d.csv' % k)
	with open(filename[0], 'r') as readFile:
		reader = csv.reader(readFile)
		curve = np.array(list(reader))
		stress = [float(i) for i in curve[1]]
		strain=[float(i) for i in curve[0]]
	#for i in range(int(len(curve[1]))):
	#	strain.append(i*0.0025)
	
	
	ax1.plot(strain, stress, marker='o',linestyle='dashed', markersize=5., color = color,label=k)
	
	
	#plot_second_der(ax2,strain, stress, color='tab:blue')
	#plot_first_der(ax2,strain, stress, 'green')
	


	#plt.plot(stress)
	#plt.show()
	return strain

def length_square(x):	
	return x[0]**2+x[1]**2

def calculate_network_data(network, length_ini, **kw):
	if kw.get('step')!=None:
		if kw.get('name')!=None:
			filename = 'network_vertices_%03d_%03d_%s.csv' % (int(kw['step']),len(network.vertices),kw['name'])
		else:
			filename = 'network_vertices_%03d_%03d.csv' % (int(kw['step']),len(network.vertices))
		with open(filename  ,'r') as readFile:
			reader = csv.reader(readFile)
			list_vertices = np.array(list(reader))
			network.vertices=list_vertices.astype(float)
	stretch_ratios,cos_theta_square = [],[]
	list_angle = []
	lengths=[]
	for ridge in network.ridge_vertices:
		# Minimum length
		length = np.sqrt(length_square(network.vertices[ridge[0]]-network.vertices[ridge[1]]))
		lengths.append(length)		
		stretch_ratios.append(length / length_ini[list(network.ridge_vertices).index(ridge)])
		# Minimum & max angle
		cos_theta_square.append(np.dot(network.vertices[ridge[0]]-network.vertices[ridge[1]], [1.,0.])**2/length**2)
		list_angle.append(cos_theta_square[-1]*length)
	omega_xx = sum(list_angle)/sum(lengths)
	return lengths, stretch_ratios, cos_theta_square, omega_xx

def plot_avg_angle(strain,network,ax2,color):
	avg_angles = []
	min_lengths=[]
	min_angles=[]
	max_angles=[]
	for k in range(len(strain)-1):
		values = smallest_values(k,network.ridge_vertices)
		min_lengths.append(values[0])
		avg_angles.append(values[1])
		min_angles.append(values[2])
		max_angles.append(values[3])
	"""fig = plt.figure()
	ax = fig.add_subplot(111)
	ax2.plot(strain[:-1],avg_angles,color=color)
	#plt.plot(strain,min_angles)
	#plt.plot(strain,max_angles)
	#plt.plot(strain,min_lengths)
	#plt.show()"""
	return avg_angles

if __name__ == "__main__":
	os.chdir('../Data/3D_GN/')
	if len(sys.argv) == 1:
		os.chdir(sorted_ls('.')[-1])
	else:
		os.chdir(sys.argv[1])
	filenames=fnmatch.filter(os.listdir('.'), 'stress_strain_*.csv')
	print(filenames)
	fig, ax1 = plt.subplots()
	color = 'tab:red'
	ax1.set_xlabel('strain')
	ax1.set_ylabel('stress', color=color)
	ax1.tick_params(axis='y',labelcolor=color)
	color = 'tab:blue'
	ax2 =ax1.twinx()
	ax2.set_ylabel('second derivative', color=color)
	ax2.tick_params(axis='y',labelcolor=color)
	for filename in filenames: 
		color = np.random.rand(3,)
		network = load_network_info(int(filename[-7:-4]))
		#print smallest_values()
		strain = plot_stress_strain(len(network.vertices),ax1,ax2,color)
		#avg_angles = plot_avg_angle(strain, network,ax2,color)
	#handles, labels = ax1.get_legend_handles_labels()
	# sort both labels and handles by labels
	#labels, handles = zip(*sorted(zip(labels, handles), key=lambda t: t[0]))
	#ax1.legend(handles, labels)
	#ax1.legend(loc='upper left')
	plt.savefig('stress_strain.pdf')
	plt.show()

	os.chdir('../../../TissueModel/')
