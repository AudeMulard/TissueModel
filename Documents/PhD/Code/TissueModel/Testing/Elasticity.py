
import sys
sys.path.append('C:/Users/am2548/TissueModel/Documents/PhD/Code/TissueModel/')
from Plotting.information_network import *
from Core_calculation.tensile_test import *
import sys
sys.path.append('/home/aude/Documents/PhD/Code/TissueModel/')
from Plotting.information_network import stress_strain_curve
import matplotlib.pyplot as plt

def length_square(x):
	if len(x) ==3:
		return x[0]**2+x[1]**2+x[2]**2
	elif len(x) ==2:
		return x[0]**2+x[1]**2

def elasticity_check(path):
	import pandas as pd
	current_path = os.getcwd()
	os.chdir(path)
	filenames = fnmatch.filter(os.listdir('.'), 'network_vertices_03_00_*.csv')
	test_number = int(filenames[0][-13:-4])
	filenames = fnmatch.filter(os.listdir('.'), 'network_vertices_01_*_%09d.csv' % test_number)
	step=0
	filename = fnmatch.filter(os.listdir('.'), 'parameters_*.csv')
	network = load_network_info(int(filename[0][-9:-4]))
	node_label = write_node_label(test_number,network)
	distance_tract_1, distance_compr_1 = [], []
	distance_tract_2, distance_compr_2 = [], []
	distance_compr_ten=[]
	for filename in filenames:
		dist= 0
		df_1 = pd.read_csv(filename,usecols=[4,12,13])
		vertices_table_1=df_1.reindex(node_label).loc[:,['   COORD-COOR1','   COORD-COOR2']].to_numpy()
		file_3 = 'network_vertices_03_%02d_%09d.csv' % (step, test_number)
		df_3 = pd.read_csv(file_3,usecols=[4,12,13])
		vertices_table_3=df_3.reindex(node_label).loc[:,['   COORD-COOR1','   COORD-COOR2']].to_numpy()
		file_5 = 'network_vertices_05_%02d_%09d.csv' % (step, test_number)
		df_5 = pd.read_csv(file_5,usecols=[4,12,13])
		vertices_table_5=df_5.reindex(node_label).loc[:,['   COORD-COOR1','   COORD-COOR2']].to_numpy()
		for i in range(len(vertices_table_1)):
			dist += length_square(vertices_table_1[i]-vertices_table_3[i])
		distance_tract_1.append(dist/len(vertices_table_1))
		dist=0
		for i in range(len(vertices_table_3)):
			dist += length_square(vertices_table_5[i]-vertices_table_3[i])
		distance_tract_2.append(dist/len(vertices_table_3))
		dist = 0
		file_4 = 'network_vertices_04_%02d_%09d.csv' % (step, test_number)
		df_4 = pd.read_csv(file_4,usecols=[4,12,13])
		vertices_table_4=df_4.reindex(node_label).loc[:,['   COORD-COOR1','   COORD-COOR2']].to_numpy()
		file_2 = 'network_vertices_02_%02d_%09d.csv' % (step, test_number)
		df_2 = pd.read_csv(file_2,usecols=[4,12,13])
		vertices_table_2=df_2.reindex(node_label).loc[:,['   COORD-COOR1','   COORD-COOR2']].to_numpy()
		file_6 = 'network_vertices_06_%02d_%09d.csv' % (step, test_number)
		df_6 = pd.read_csv(file_6,usecols=[4,12,13])
		vertices_table_6=df_6.reindex(node_label).loc[:,['   COORD-COOR1','   COORD-COOR2']].to_numpy()
		step=step+1
		for i in range(len(vertices_table_2)):
			dist += length_square(vertices_table_2[i]-vertices_table_4[i])
		distance_compr_1.append(dist/len(vertices_table_2))
		dist = 0
		for i in range(len(vertices_table_4)):
			dist += length_square(vertices_table_4[i]-vertices_table_6[i])
		distance_compr_2.append(dist/len(vertices_table_4))
		dist=0
		for i in range(len(vertices_table_4)):
			dist += length_square(vertices_table_5[i]-vertices_table_6[i])
		distance_compr_ten.append(dist/len(vertices_table_5))
	with open('testing.txt', 'a') as writeFile:
		writeFile.write('Difference in network positions in traction 1:\n')
		writeFile.write(str(distance_tract_1)+'\n')
	with open('testing.txt', 'a') as writeFile:
		writeFile.write('Difference in network positions in traction 2:\n')
		writeFile.write(str(distance_tract_2)+'\n')
	with open('testing.txt', 'a') as writeFile:
		writeFile.write('Difference in network positions in compression :\n')
		writeFile.write(str(distance_compr_2)+'\n')
	with open('testing.txt', 'a') as writeFile:
		writeFile.write('Difference in network positions in compression :\n')
		writeFile.write(str(distance_compr_2)+'\n')
	with open('testing.txt', 'a') as writeFile:
		writeFile.write('Difference in network positions in compression and tension :\n')
		writeFile.write(str(distance_compr_2)+'\n')
	
	fig_stress = plt.figure()
	strain, stress_curve = stress_strain_curve(test_number,network)
	print(strain,stress_curve)
	stress_1 = stress_curve[0:10]
	stress_2 = stress_curve[10:21]
	stress_3 = stress_curve[21:32]
	stress_4 = stress_curve[32:43]
	stress_5 = stress_curve[43:54]
	stress_6 = stress_curve[54:65]
	difftract,diffcomp = [],[]
	for i in range(len(stress_1)):
		difftract.append(stress_1[i]-stress_3[i])
		difftract.append(stress_3[i]-stress_5[i])
		diffcomp.append(stress_2[i]-stress_4[i])
		diffcomp.append(stress_4[i]-stress_6[i])
	with open('testing.txt', 'a') as writeFile:
		writeFile.write('Difference in global stress in traction:\n')
		writeFile.write(str(difftract)+'\n')
	with open('testing.txt', 'a') as writeFile:
		writeFile.write('Difference in global stress in compression:\n')
		writeFile.write(str(diffcomp)+'\n')
	fig_stress = plt.figure()
	plt.plot(strain[0:10],stress_1,label='Traction 1')
	plt.plot(strain[10:21],stress_2,label='Compression 1')
	plt.plot(strain[21:32],stress_3,label='Traction 2')
	plt.plot(strain[32:43],stress_4,label='Compression 2')
	plt.plot(strain[43:54],stress_5,label='Traction 3')
	plt.plot(strain[54:65],stress_6,label='Compression 3')
	plt.title('Comparison of stress for elasticity test')
	plt.savefig('elas_stress_strain.pdf')
	os.chdir(current_path)


def network_def_elasticity(network,path,test_1):
	network= network.set_fibers(path)
	test_1.save_parameters(network,path)
	os.system("abaqus cae script=elasticity.py")
	

