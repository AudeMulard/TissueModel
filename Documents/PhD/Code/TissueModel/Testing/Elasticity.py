

from Core_calculation.tensile_test import *

def length_square(x):
	if len(x) ==3:
		return x[0]**2+x[1]**2+x[2]**2
	elif len(x) ==2:
		return x[0]**2+x[1]**2

def elasticity_check(path):
	current_path = os.getcwd()
	os.chdir(path)
	# Compare geometries in traction & compression
	filenames = fnmatch.filter(os.listdir('.'), 'network_vertices_*_First_tract.csv')
	distance_tract, distance_compr = [], []
	for filename in filenames:
		dist = 0
		with open(filename,'r') as readFile:
			reader = csv.reader(readFile)
			list_vertices = np.array(list(reader))
			vertices_first=list_vertices.astype(float)
		with open('network_vertices_%s_%s_Second_tract.csv' % (filename[17:20],filename[21:24])) as readFile:
			reader = csv.reader(readFile)
			list_vertices = np.array(list(reader))
			vertices_second=list_vertices.astype(float)
		for i in range(len(vertices_first)):
			dist += length_square(vertices_first[i]-vertices_second[i])
		distance_tract.append(dist/len(vertices_first))
		dist = 0
		with open('network_vertices_%s_%s_First_compr.csv' % (filename[17:20],filename[21:24])) as readFile:
			reader = csv.reader(readFile)
			list_vertices = np.array(list(reader))
			vertices_first=list_vertices.astype(float)
		with open('network_vertices_%s_%s_Second_compr.csv' % (filename[17:20],filename[21:24])) as readFile:
			reader = csv.reader(readFile)
			list_vertices = np.array(list(reader))
			vertices_second=list_vertices.astype(float)
		for i in range(len(vertices_first)):
			dist += length_square(vertices_first[i]-vertices_second[i])
		distance_compr.append(dist/len(vertices_first))
	with open('testing.txt', 'a') as writeFile:
		writeFile.write('Difference in network positions in traction :')
		writeFile.write(str(distance_tract))
	with open(path,'testing.txt', 'a') as writeFile:
		writeFile.write('Difference in network positions in compression :')
		writeFile.write(str(distance_compr))
	
	# Compare stress_strain_curve
	filenames = fnmatch.filter(os.listdir('.'), 'stress_strain_*.csv')
	print filenames
	stress = []
	difftract,diffcomp = [],[]
	for filename in filenames:
		with open(filename,'r') as readFile:
			reader = csv.reader(readFile)
			curve = np.array(list(reader))
			stress.append([float(i) for i in curve[1]])
	print stress
	for i in range(len(curve[1])-1):
		print i, len(curve[1])
		difftract.append(stress[0][i]-stress[2][i-1])
		diffcomp.append(stress[1][i]-stress[3][i])
		
	with open('testing.txt', 'a') as writeFile:
		writeFile.write('Difference in global stress in traction:')
		writeFile.write(str(difftract))
	with open('testing.txt', 'a') as writeFile:
		writeFile.write('Difference in global stress in traction:')
		writeFile.write(str(diffcomp))
	
	os.chdir(current_path)


def network_def_elasticity(network,path,test_1):
	network= network.set_fibers(path)
	network_1 = test_1.full_test(network, path,test_1.details,name='First_tract')
	test_1.traction_distance = -test_1.traction_distance
	network_1 = test_1.full_test(network, path,test_1.details,name='First_compr')
	test_1.traction_distance = -test_1.traction_distance
	network_1 = test_1.full_test(network, path,test_1.details,name='Second_tract')
	test_1.traction_distance = -test_1.traction_distance
	network_1 = test_1.full_test(network, path,test_1.details,name='Second_compr')
	

