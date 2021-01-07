

from Core_calculation.tensile_test import *
import fnmatch
import sys
sys.path.append('/home/aude/Documents/PhD/Code/TissueModel/')
from Plotting.information_network import stress_strain_curve

def length_square(x):
	if len(x) ==3:
		return x[0]**2+x[1]**2+x[2]**2
	elif len(x) ==2:
		return x[0]**2+x[1]**2

def idparams_check(path):
	current_path = os.getcwd()
	os.chdir(path)
	# Compare geometries
	filenames = fnmatch.filter(os.listdir('.'), 'network_vertices_01_*.rpt')
	import pandas as pd
	test_number_1 = int(filenames[0][-13:-4])
	test_number_2 = int(filenames[1][-13:-4])
	filenames = fnmatch.filter(os.listdir('.'), 'network_vertices_*_%09d.rpt' % test_number_2)
	step=0
	distance=[]
	node_label_1 = pd.read_csv('node_label_%09d.csv' % test_number_1,header=None)
	node_label_2 = pd.read_csv('node_label_%09d.csv' % test_number_2,header=None)
	for filename in filenames:
		dist=0
		df_1 = pd.read_fwf(filename,skiprows=19,sep='\t',header=None)
		df_1.columns=['Node_label','Magnitude','X','Y','Z']
		
		vertices_table_1=df_1.reindex(list(node_label_1.transpose()[0])).loc[:,['X','Y','Z']].to_numpy()
		file_2 = 'network_vertices_01_%d_%s.rpt' % (step, test_number_2)
		df_2 = pd.read_fwf(filename,skiprows=19,sep='\t',header=None)
		df_2.columns=['Node_label','Magnitude','X','Y','Z']
		
		vertices_table_2=df_2.reindex(list(node_label_2.transpose()[0])).loc[:,['X','Y','Z']].to_numpy()
		step=step+1
		for i in range(min(len(vertices_table_1),len(vertices_table_2))):
			dist += length_square(vertices_table_1[i]-vertices_table_2[i])
		distance.append(dist/len(vertices_table_1))
	with open('testing.txt', 'a') as writeFile:
		writeFile.write('Difference in network positions :\n')
		writeFile.write(str(distance)+'\n')
	
	stress_1 = stress_strain_curve(test_number_1)
	stress_2 = stress_strain_curve(test_number_2)
	diff = []
	for i in range(len(stress_1)):
		diff.append(stress_1[i]-stress_2[i])
	with open('testing.txt', 'a') as writeFile:
		writeFile.write('Difference in global stress :\n')
		writeFile.write(str(diff)+'\n')

	os.chdir(current_path)


def network_def_idparams(network,path,test_1):
	network_1 = network.set_fibers(path)
	test_1.save_parameters(network_1,path)
	os.system("abaqus cae noGUI=solver.py")
	print('First calculation done')
	
	network_2 = network.set_fibers(path)
	test_1.save_parameters(network_2,path)
	os.system("abaqus cae noGUI=solver.py")

