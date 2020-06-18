

from Core_calculation.tensile_test import *

def length_square(x):
	if len(x) ==3:
		return x[0]**2+x[1]**2+x[2]**2
	elif len(x) ==2:
		return x[0]**2+x[1]**2

def idparams_check(path):
	current_path = os.getcwd()
	os.chdir(path)
	# Compare geometries
	filenames = fnmatch.filter(os.listdir('.'), 'network_vertices_*_first.csv')
	len_second = fnmatch.filter(os.listdir('.'),'network_vertices_000_*_second.csv')
	distance = []
	for filename in filenames:
		dist = 0
		with open(filename,'r') as readFile:
			reader = csv.reader(readFile)
			list_vertices = np.array(list(reader))
			vertices_first=list_vertices.astype(float)
		with open('network_vertices_%s_%s_second.csv' % (filename[17:20],len_second[0][21:24])) as readFile:
			reader = csv.reader(readFile)
			list_vertices = np.array(list(reader))
			vertices_second=list_vertices.astype(float)
		for i in range(min(len(vertices_first),len(vertices_second))):
			dist += length_square(vertices_first[i]-vertices_second[i])
		distance.append(dist/len(vertices_first))
	with open('testing.txt', 'a') as writeFile:
		writeFile.write('Difference in network positions :\n')
		writeFile.write(str(distance)+'\n')
	
	# Compare stress_strain_curve
	filenames = fnmatch.filter(os.listdir('.'), 'stress_strain_*.csv')
	stress = []
	diff = []
	for filename in filenames:
		with open(filename,'r') as readFile:
			reader = csv.reader(readFile)
			curve = np.array(list(reader))
			stress.append([float(i) for i in curve[1]])
	for i in range(len(curve[1])):
		diff.append(stress[0][i]-stress[1][i]) # si pas meme implementer linearisation
		
	with open('testing.txt', 'a') as writeFile:
		writeFile.write('Difference in global stress :\n')
		writeFile.write(str(diff)+'\n')
	
	os.chdir(current_path)


def network_def_idparams(network,path,test_1):
	network_1 = network.set_fibers(path)
	network_1 = test_1.full_test(network, path,test_1.details,name='first')

	network_2 = network.set_fibers(path)
	network_2 = test_1.full_test(network, path,test_1.details,name='second')

