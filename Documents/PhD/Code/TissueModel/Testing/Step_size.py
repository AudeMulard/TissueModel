

from Core_calculation.tensile_test import *

def length_square(x):
	if len(x) ==3:
		return x[0]**2+x[1]**2+x[2]**2
	elif len(x) ==2:
		return x[0]**2+x[1]**2

def identical_check(path):
	current_path = os.getcwd()
	os.chdir(path)
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
			if norm_1 <10e-15:
				with open('testing.txt', 'a') as writeFile:
					writeFile.write(filename[17:20]+str(i)+str(norm_1))
			if norm_2 <10e-15:
				with open('testing.txt', 'a') as writeFile:
					writeFile.write(filename[17:20]+str(i)+str(norm_2))
	with open('testing.txt', 'a') as writeFile:
		writeFile.write('Different step size test : done')
	os.chdir(current_path)


def network_def_tests(network,path,test_1):
	network = network.set_fibers(path)
	
	network = test_1.full_test(network, path,test_1.details,name='first_step')

	network.vertices = np.array(network.vertices_ini)
	test_1.space_discretization  = test_1.space_discretization/3.
	network = test_1.full_test(network, path,test_1.details,name='second_step')
	
	test_1.space_discretization  = test_1.space_discretization*3./10.
	network.vertices = np.array(network.vertices_ini)
	network = test_1.full_test(network, path,test_1.details,name='third_step')

