

from Core_calculation.tensile_test import *
import sys
sys.path.append('C:/Users/am2548/TissueModel/Documents/PhD/Code/TissueModel/')
from Plotting.information_network import *
def length_square(x):
	if len(x) ==3:
		return x[0]**2+x[1]**2+x[2]**2
	elif len(x) ==2:
		return x[0]**2+x[1]**2

def discretization_check(path):
	import pandas as pd
	current_path = os.getcwd()
	os.chdir(path)
	filenames = fnmatch.filter(os.listdir('.'), 'network_vertices_01_*.csv')
	test_number_1 = int(filenames[0][-13:-4])
	test_number_2 = int(filenames[1][-13:-4])
	test_number_3 = int(filenames[2][-13:-4])
	test_number_4 = int(filenames[3][-13:-4])
	filenames = fnmatch.filter(os.listdir('.'), 'network_vertices_*_%09d.csv' % test_number_1)
	filename = fnmatch.filter(os.listdir('.'), 'parameters_*.csv')
	network = load_network_info(int(filename[0][-9:-4]))
	node_label_1 = write_node_label(test_number_1,network)
	node_label_2 = write_node_label(test_number_2,network)
	node_label_3 = write_node_label(test_number_3,network)
	node_label_4 = write_node_label(test_number_4,network)
	step=0
	criteria_1 = True
	criteria_2 = True
	for filename in filenames:
		df_1 = pd.read_csv(filename,usecols=[4,12,13])
		vertices_table_1=df_1.reindex(node_label_1).loc[:,['   COORD-COOR1','   COORD-COOR2']].to_numpy()
		file_2 = 'network_vertices_01_%02d_%09d.csv' % (step, test_number_2)
		df_2 = pd.read_csv(file_2,usecols=[4,12,13])
		vertices_table_2=df_2.reindex(node_label_2).loc[:,['   COORD-COOR1','   COORD-COOR2']].to_numpy()
		file_3 = 'network_vertices_01_%02d_%09d.csv' % (step, test_number_3)
		df_3 = pd.read_csv(file_3,usecols=[4,12,13])
		vertices_table_3=df_3.reindex(node_label_3).loc[:,['   COORD-COOR1','   COORD-COOR2']].to_numpy()
		file_4 = 'network_vertices_01_%02d_%09d.csv' % (step, test_number_4)
		df_4 = pd.read_csv(file_4,usecols=[4,12,13])
		vertices_table_4=df_4.reindex(node_label_4).loc[:,['   COORD-COOR1','   COORD-COOR2']].to_numpy()
		table = np.array([[1,len(df_1.index),vertices_table_1],[2,len(df_2.index),vertices_table_2],[3,len(df_3.index),vertices_table_3],[4,len(df_4.index),vertices_table_4]])
		sorted_table=table[np.argsort(table[:, 1])]
		step=step+1
		print(np.amax(sorted_table[0,2]-sorted_table[1,2]),np.amax(sorted_table[1,2]-sorted_table[2,2]),np.amax(sorted_table[2,2]-sorted_table[3,2]))
		print(np.mean(sorted_table[0,2]-sorted_table[1,2]),np.mean(sorted_table[1,2]-sorted_table[2,2]),np.mean(sorted_table[2,2]-sorted_table[3,2]))
		#print(np.isclose(vertices_table_1,vertices_table_2,rtol = 1e-3,atol=1e-3))
		#print(np.isclose(vertices_table_2,vertices_table_3,rtol=1e-3, atol=1e-3))
		if np.allclose(vertices_table_1,vertices_table_2,rtol = 1e-4,atol = 1e-5)==False:
			criteria_1= False
		if np.allclose(vertices_table_2,vertices_table_3,rtol = 1e-4,atol = 1e-5)==False:
			criteria_2 = False
		# read file with pandas, convert to list and then use isclose## use df.to_numy()
	with open('testing.txt', 'a') as writeFile:
		writeFile.write('Different step size test : %s, %s all good\n' % (criteria_1,criteria_2))
	os.chdir(current_path)
	
	"""
	filenames = fnmatch.filter(os.listdir('.'), 'network_vertices_*_first_step.csv')
	

	for filename in filenames: # change, will not work because of teh different step size, need to adapt
		with open(filename,'r') as readFile:
			reader = csv.reader(readFile)
			list_vertices = np.array(list(reader))
			vertices_first=list_vertices.astype(float)
		with open('network_vertices_%03d_%s_second_step.csv' % (3*(int(filename[17:20])+1)-1,filename[21:24])) as readFile:
			reader = csv.reader(readFile)
			list_vertices = np.array(list(reader))
			vertices_second=list_vertices.astype(float)
		with open('network_vertices_%03d_%s_third_step.csv' % (10*(int(filename[17:20])+1)-1,filename[21:24])) as readFile:
			reader = csv.reader(readFile)
			list_vertices = np.array(list(reader))
			vertices_third=list_vertices.astype(float)
		for i in range(len(vertices_first)):
			norm_1 = length_square(vertices_first[i]-vertices_second[i])
			norm_2 = length_square(vertices_first[i]-vertices_third[i])
			if norm_1 >10e-15:
				with open('testing.txt', 'a') as writeFile:
					writeFile.write(filename[17:20]+' ' + str(i)+' ' +str(norm_1)+'\n')
			if norm_2 >10e-15:
				with open('testing.txt', 'a') as writeFile:
					writeFile.write(filename[17:20]+' ' +str(i)+' ' +str(norm_2)+'\n')
	with open('testing.txt', 'a') as writeFile:
		writeFile.write('Different step size test : done\n')
	
	"""

def network_def_tests(network,path,test_1):
	network = network.set_fibers(path)
	print('strating')
	test_1.save_parameters(network,path)
	os.system("abaqus cae noGUI=new_solver_1.py")
	print('First calculation done')

	test_1.element_size  = test_1.element_size/2.
	test_1.save_parameters(network,path)
	os.system("abaqus cae noGUI=new_solver_1.py")
	print('Second calculation done')
	
	test_1.element_size  = test_1.element_size/2.
	test_1.save_parameters(network,path)
	os.system("abaqus cae noGUI=new_solver_1.py")
	print('Third calculation done')
	
	test_1.element_size  = test_1.element_size/2.
	test_1.save_parameters(network,path)
	os.system("abaqus cae noGUI=new_solver_1.py")
	print('Fourth calculation done')

