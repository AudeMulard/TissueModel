from Network_generation.creation_network import Network
from Network_generation.add_ons_network import Cell
from Core_calculation.tensile_test import *
from scipy.spatial import Voronoi, voronoi_plot_2d
import matplotlib.pyplot as plt
#from Several_networks_test import *
from Plotting.network_plotting import *
import os
from datetime import date
import matplotlib.patches as patches
import fnmatch,time

## PARAMETERS
dimension=2 #dimension of the problem
complexity_network=500 #number of random seed points
length_domain=(1.0,1.0,1.0)
min_distance = 0.0001*length_domain[0]
space_discretization = 0.1*length_domain[0]
beam_Young = 33000
beam_poisson = 0.3
beam_profile = 0.01
connector_coeff = 0.001
disturbance=0.02
traction_distance = 0.5*length_domain[0]
#iteration = 15



## EXPERIMENT
creation="growth_network"
generation = 'grid'
constitutive = 'spring'
side = 'right'
plot = True
video = True
phase = 'only_one'
stress_rep = True
details = True
element_size=0.001
hyperstatic_param=2
data_path = '../Data_1/default/' 

today = date.today()

new_dir = data_path+today.strftime("%b-%d-%Y")+'_'+'%04d' % len(os.listdir(data_path))
os.mkdir(new_dir)
path = new_dir


#list_modes = [['Voronoi','random'],['Voronoi','grid'],['growth_network','grid'],['Voronoi','regular']]#,['growth_network','grid']]
list_modes = [['growth_network','grid']]
test_1 = Tensile_test(constitutive, side, space_discretization, traction_distance, element_size, plot, video, path,details)


	
	#network = Network(dimension, complexity_network, length_domain, min_distance, k_tension, k_compression, A, disturbance,hyperstatic_param, creation,generation, path,name='%s_%s' % (creation, generation))
	#network = network.set_fibers(path)
	#### Changing disturbance


	
	#### Changing geometry
	# Changing complexity
"""
for k in range(6):
	for complexity in [200,400,600,800,1000,1200,1500,2000,2500]:
		complexity_network=complexity
		with open(os.path.join(path,'time_complexity.txt'), 'a') as writeFile:
			writeFile.write('%s \n' % complexity)
			writeFile.write('%d \n' % time.time())
			start = time.time()
			network = Network(dimension, complexity_network, length_domain, min_distance, beam_Young, beam_poisson, beam_profile, connector_coeff, disturbance, hyperstatic_param, creation, generation, path)
			network = network.set_fibers(path)
			test_1 = Tensile_test(constitutive, side, space_discretization, traction_distance,element_size, plot, video, path,details)
			test_1.save_parameters(network,path)
			os.system("abaqus cae noGUI=new_solver_1.py")
			writeFile.write('Test time: %d \n' % (time.time() - start))
			writeFile.close()
			
"""
#new_dir = data_path+today.strftime("%b-%d-%Y")+'_'+'%04d' % len(os.listdir(data_path))
#os.mkdir(new_dir)
#path = new_dir
creation = 'Voronoi'
generation = 'random'
#complexity = [250,250,300,350,400,550,700,1200,3000
#complexity = [275,275,350,400,500,600,700,1200,2500]#,5000]
#hyperstatic_param =[0,1.6,1.75,1.9,2,2.2,2.4,2.6,2.8]#,3.0]
complexity = [700,700,1200,3000]
hyperstatic_param =[2.4,2.4,2.6,2.8]
#hyperstatic_param =[2.6,2.8,3.0]

for i in range(len(hyperstatic_param)): 
	network = Network(dimension, complexity[i], length_domain, min_distance, beam_Young, beam_poisson, beam_profile, connector_coeff, disturbance, hyperstatic_param[i], creation, generation, path)
	network = network.set_fibers(path)
	test_1.save_parameters(network,path)	
	os.system("abaqus cae noGUI=new_solver_4.py")

"""	

#new_dir = data_path+today.strftime("%b-%d-%Y")+'_'+'%04d' % len(os.listdir(data_path))
#os.mkdir(new_dir)
#path = new_dir
creation = 'Voronoi'
generation='random'
network = Network(dimension, complexity_network, length_domain, min_distance, beam_Young, beam_poisson, beam_profile, connector_coeff, disturbance, hyperstatic_param, creation, generation, path)
network = network.set_fibers(path)
test_1.save_parameters(network,path)
#network.beam_profile = 0.01	
for connector_coeff in [0.001,0.005,0.01,0.05,0.1,0.5,1.0,1.1]:
	network.connector_coeff = connector_coeff
	network.save_network('initial',path)
	os.system("abaqus cae noGUI=new_solver_1.py")


	if mode[1] == 'grid':
		for disturbance in [0.,0.0001,0.001,0.01,0.1,1.]:
			network = Network(dimension, complexity_network, length_domain, min_distance, beam_Young, beam_poisson, beam_profile, connector_coeff, disturbance, hyperstatic_param, creation, generation, path)
			network = network.set_fibers(path)
			network = test_1.full_test(network, path,test_1.details,name='%s_%s_disturbance_%04d' % (creation,generation,disturbance))

for k in range(5):
	for mode in list_modes:
		print(mode)
		
		##################### Parameter study (comparison on same network) #########################
		creation = mode[0]
		generation = mode[1]	
		#### Time calculation
		with open(os.path.join(path,'time.txt'), 'a') as writeFile:
			if creation=='Voronoi':
				complexity_network=500
			elif creation=='growth_network':
				complexity_network=300
			
			writeFile.write('%s \n' % mode)
			writeFile.write('%d \n' % time.time())
			start = time.time()
			network = Network(dimension, complexity_network, length_domain, min_distance, beam_Young, beam_poisson, beam_profile, connector_coeff, disturbance, hyperstatic_param, creation, generation, path)
			network = network.set_fibers(path)
			#plot_geometry(network)
			#plt.show()
			writeFile.write('Generation time: %d \n' % (time.time() - start))
			test_1 = Tensile_test(constitutive, side, space_discretization, traction_distance,element_size, plot, video, path,details)
			test_1.save_parameters(network,path)
			os.system("abaqus cae script=new_solver_3.py")
			#network = test_1.full_test(network, path,test_1.details,name='%s_%s_time' % (creation,generation))
			writeFile.write('Test time: %d \n' % (time.time() - start))
			writeFile.close()


#### Without changing network
# Change of k_tension, k_compression
network = Network(dimension, complexity_network, length_domain, min_distance, beam_Young, beam_poisson, beam_profile, connector_coeff, disturbance, hyperstatic_param, creation, generation, path)
network = network.set_fibers(path)
test_1 = Tensile_test(constitutive, side, space_discretization, traction_distance, element_size, plot, video, path,details)
test_1.save_parameters(network,path)

	for truss_Young in [1.,100.,1000.,10000.,33000.,1e5,1e6]:
		network.truss_young = truss_Young
		network.save_network('initial',path)
		os.system("abaqus cae noGUI=solver.py")

for beam_profile in [0.01]:
	network.beam_profile = beam_profile
	network.save_network('initial',path)
	os.system("abaqus cae script=new_solver_4.py")


creation = 'Voronoi'
generation='grid'
for disturbance in [0.01,0.02,0.03,0.04,0.05]:#, 0.06,0.07,0.08,0.09,0.1]:
	network = Network(dimension, complexity_network, length_domain, min_distance, beam_Young, beam_poisson, beam_profile, connector_coeff, disturbance, hyperstatic_param, creation, generation, path)
	network = network.set_fibers(path)
	test_1 = Tensile_test(constitutive, side, space_discretization, traction_distance,element_size, plot, video, path,details)
	test_1.save_parameters(network,path)
	os.system("abaqus cae noGUI=new_solver_1.py")
"""