

from Core_calculation.tensile_test import *

def length_square(x):
	if len(x) ==3:
		return x[0]**2+x[1]**2+x[2]**2
	elif len(x) ==2:
		return x[0]**2+x[1]**2

def identical_check(path):
	current_path = os.getcwd()
	os.chdir(path)
	filenames = fnmatch.filter(os.listdir('.'), 'network_vertices_*_first.csv')

	for filename in filenames:
		with open(filename,'r') as readFile:
			reader = csv.reader(readFile)
			list_vertices = np.array(list(reader))
			vertices_first=list_vertices.astype(float)
		with open('network_vertices_%s_%s_second.csv' % (filename[17:20],filename[21:24])) as readFile:
			reader = csv.reader(readFile)
			list_vertices = np.array(list(reader))
			vertices_second=list_vertices.astype(float)
		for i in range(len(vertices_first)):
			norm = length_square(vertices_first[i]-vertices_second[i])
			if norm >= 10e-5:
				with open('testing.txt', 'a') as writeFile:
					writeFile.write(filename[17:20]+' '+str(i),+'\n')
	with open('testing.txt', 'a') as writeFile:
		writeFile.write('Identical network test : all good\n')
	os.chdir(current_path)


def network_def_tests(network,path,test_1):
	network = network.set_fibers(path)
	
	network = test_1.full_test(network, path,test_1.details,name='first')
	print('First calculation done')
	
	network.vertices = np.array(network.vertices_ini)
	network = test_1.full_test(network, path,test_1.details,name='second')

