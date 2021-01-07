import sys
sys.path.append('C:/Users/am2548/TissueModel/Documents/PhD/Code/TissueModel/')
from Plotting.information_network import *
from Core_calculation.tensile_test import *
import fnmatch
import pandas as pd
import numpy as np
def length_square(x):
	if len(x) ==3:
		return x[0]**2+x[1]**2+x[2]**2
	elif len(x) ==2:
		return x[0]**2+x[1]**2

def identical_check(path):
	current_path = os.getcwd()
	os.chdir(path)
	filenames = fnmatch.filter(os.listdir('.'), 'network_vertices_01_*.csv')
	test_number_1 = int(filenames[0][-13:-4])
	test_number_2 = int(filenames[1][-13:-4])
	filenames = fnmatch.filter(os.listdir('.'), 'network_vertices_*_%09d.csv' % test_number_1)
	filename = fnmatch.filter(os.listdir('.'), 'parameters_*.csv')
	network_1 = load_network_info(int(filename[0][-9:-4]))
	network_2 = load_network_info(int(filename[0][-9:-4]))
	node_label_1 = write_node_label(test_number_1,network_1)
	#node_label_1 = pd.read_csv('node_label_%09d.csv' % test_number_1,header=None)
	#node_label_2 = pd.read_csv('node_label_%09d.csv' % test_number_2,header=None)
	node_label_2 = write_node_label(test_number_2,network_2)
	step=0
	criteria = True
	for filename in filenames:
		df_1 = pd.read_csv(filename,usecols=[4,12,13])
		#df_1 = pd.read_fwf(filename,skiprows=19,sep='\t',header=None)
		#df_1.columns=['Node_label','Magnitude','X','Y','Z']
		#vertices_table_1=df_1.reindex(list(node_label_1.transpose()[0])).loc[:,['X','Y','Z']].to_numpy()
		vertices_table_1=df_1.reindex(node_label_1).loc[:,['   COORD-COOR1','   COORD-COOR2']].to_numpy()
		file_2 = 'network_vertices_01_%02d_%09d.csv' % (step, test_number_2)
		df_2 = pd.read_csv(file_2,usecols=[4,12,13])
		#df_2.columns=['Node_label','Magnitude','X','Y','Z']
		#vertices_table_2=df_2.reindex(list(node_label_2.transpose()[0])).loc[:,['X','Y','Z']].to_numpy()
		vertices_table_2=df_2.reindex(node_label_2).loc[:,['   COORD-COOR1','   COORD-COOR2']].to_numpy()
		step=step+1
		if np.allclose(vertices_table_1,vertices_table_2,rtol = 1e-5,atol = 10e-8)==False:
			criteria = False
		# read file with pandas, convert to list and then use isclose## use df.to_numy()
	with open('testing.txt', 'a') as writeFile:
		if criteria==True:
			writeFile.write('Identical network test : all good\n')
		else:
			writeFile.write('Identical network test : not good\n')
	os.chdir(current_path)


def network_def_tests(network,path,test_1):
	network = network.set_fibers(path)
	test_1.save_parameters(network,path)
	os.system("abaqus cae noGUI=new_solver_1.py")
	#network = test_1.full_test(network, path,test_1.details,name='first')
	print('First calculation done')
	
	#network.vertices = np.array(network.vertices_ini)
	os.system("abaqus cae noGUI=new_solver_1.py")
	#network = test_1.full_test(network, path,test_1.details,name='second')

