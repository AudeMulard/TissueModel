import os
import csv
import numpy as np
import matplotlib.pyplot as plt
import sys
import fnmatch
import math
from statistics import mean 
from network_plotting import load_network_info

def sorted_ls(path):
    mtime = lambda f: os.stat(os.path.join(path, f)).st_mtime
    return list(sorted(os.listdir(path), key=mtime))


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

def smallest_values(step,ridge_vertices):
	with open('network_vertices_%03d_%03d.csv' % (int(step),len(network.vertices))  ,'r') as readFile:
		reader = csv.reader(readFile)
		list_vertices = np.array(list(reader))
		vertices=list_vertices.astype(float)
	min_length = 1.	
	list_angle = []
	lengths=[]
	for ridge in ridge_vertices:
		# Minimum length
		length = np.sqrt(length_square(vertices[ridge[0]]-vertices[ridge[1]]))
		lengths.append(length)
		if length < min_length:
			min_length=length
		# Minimum & max angle
		list_angle.append(abs(np.dot(vertices[ridge[0]]-vertices[ridge[1]], [1.,0.]))**2/length)
	avg_angle = sum(list_angle)/sum(lengths)
	max_angle = max(list_angle)

	min_angle = min(list_angle)
	return min_length, avg_angle, min_angle, max_angle

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
	ax = fig.add_subplot(111)"""
	ax2.plot(strain[:-1],avg_angles,color=color)
	#plt.plot(strain,min_angles)
	#plt.plot(strain,max_angles)
	#plt.plot(strain,min_lengths)
	#plt.show()"""
	return avg_angles

if __name__ == "__main__":
	os.chdir('../Data/compression_tension')
	if len(sys.argv) == 1:
		os.chdir(sorted_ls('.')[-1])
	else:
		os.chdir(sys.argv[1])
	filenames=fnmatch.filter(os.listdir('.'), 'stress_strain_*.csv')
	print filenames
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
		print filename[-7:-4]
		network = load_network_info(int(filename[-7:-4]))
		#print smallest_values()
		strain = plot_stress_strain(len(network.vertices),ax1,ax2,color)
		avg_angles = plot_avg_angle(strain, network,ax2,color)
		print avg_angles[-1]-avg_angles[0]
	handles, labels = ax1.get_legend_handles_labels()
	# sort both labels and handles by labels
	labels, handles = zip(*sorted(zip(labels, handles), key=lambda t: t[0]))
	ax1.legend(handles, labels)
	#ax1.legend(loc='upper left')
	plt.savefig('stress_strain.pdf')
	plt.show()

	os.chdir('../../../TissueModel/')
