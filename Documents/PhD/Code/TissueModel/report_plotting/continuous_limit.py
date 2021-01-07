import sys, fnmatch, os
import matplotlib.pyplot as plt
sys.path.append('C:/Users/am2548/TissueModel/Documents/PhD/Code/TissueModel/')
from Plotting.information_network import *
#from Study.Study_plotting import *

current_dir=os.getcwd()
os.chdir('../Data_1/Study_networks/continuous_limit/')

def stress_strain_curve(test_number,network,type):
	import pandas as pd
	filenames = fnmatch.filter(os.listdir('.'), 'stress_data_*_%09d.rpt' % int(test_number))
	stress=[]
	lengths = []
	for ridge in network.ridge_vertices:
		# Minimum length
		length = np.sqrt(length_square(network.vertices[ridge[0]]-network.vertices[ridge[1]]))
		lengths.append(length)
	df_1 = pd.read_fwf(filenames[0],skiprows=2,skipfooter=4,colspecs=([0,30],[31,50]), sep='\t').to_numpy()
	stress_1 = df_1[:,1]
	strain=[]
	for i in range(len(stress_1)):
		if type == 'global_length':
			stress.append(stress_1[i]*df_1[i,0]/(1+df_1[i,0]))
		elif type == 'number_bc':
			stress.append(stress_1[i]*df_1[i,0]/(len(network.boundary_nodes_right)+len(network.boundary_nodes_left)))
		elif type == 'mean_length':
			stress.append(stress_1[i]*df_1[i,0]/(np.mean(lengths)))
		elif type == 'number_ridge':
			stress.append(stress_1[i]*df_1[i,0]/(np.mean(lengths)*len(network.ridge_vertices)))#(len(network.boundary_nodes_right)+len(network.boundary_nodes_left)))#np.mean(lengths)))#*len(network.ridge_vertices)))#/(len(network.ridge_vertices)*np.mean(lengths)))
		strain.append(df_1[i,0]*1.0)
	return strain,stress




test_numbers_file = fnmatch.filter(sorted_ls('.'), 'network_vertices_01_00_*.csv')
number_vertices_file = fnmatch.filter(sorted_ls('.'), 'parameters_??_*.csv')

#print(len(number_vertices_file),len(test_numbers_file))
test_numbers = []
number_vertices = []
for i in range(len(test_numbers_file)):
	test_numbers.append(test_numbers_file[i][-13:-4])
	number_vertices.append(int(number_vertices_file[i][-9:-4]))

test_numbers = [x for _,x in sorted(zip(number_vertices,test_numbers))]
number_vertices = sorted(number_vertices)
print(number_vertices)
test_1 = load_info_test(number_vertices[0])

complexity = [100,200,400,600,800,1000,1200,1500,2000]
for type in ['global_length','number_bc','mean_length','number_ridge']:
	end_stress = []
	end_std=[]
	mean_n_vertices=[]
	fig_stress_strain,axss = plt.subplots()
	#axssd =axss.twinx()
	#axssd.set_ylabel('second derivative')
	#axssd.tick_params(axis='y')
	for i in range(9):
		k = i*5
		strain=[]
		network_0 = load_network_info(number_vertices[k])
		strain_0,stress_0 = stress_strain_curve(test_numbers[k],network_0,type)
		#print(len(stress_0),len(strain_0))
		strain.extend(strain_0)
		network_1 = load_network_info(number_vertices[k+1])
		strain_1,stress_1 = stress_strain_curve(test_numbers[k+1],network_1,type)
		strain.extend(strain_1)
		network_2 = load_network_info(number_vertices[k+2])
		strain_2,stress_2 = stress_strain_curve(test_numbers[k+2],network_2,type)
		strain.extend(strain_2)
		network_3 = load_network_info(number_vertices[k+3])
		strain_3,stress_3 = stress_strain_curve(test_numbers[k+3],network_3,type)
		strain.extend(strain_3)
		network_4 = load_network_info(number_vertices[k+4])
		strain_4,stress_4 = stress_strain_curve(test_numbers[k+4],network_4,type)
		strain.extend(strain_4)
		strain = sorted(list(set(strain)))
		stress_0 = np.interp(strain,strain_0,stress_0)
		stress_1 = np.interp(strain,strain_1,stress_1)
		stress_2 = np.interp(strain,strain_2,stress_2)
		stress_3 = np.interp(strain,strain_3,stress_3)
		stress_4 = np.interp(strain,strain_4,stress_4)
		stress_list = [stress_0,stress_1,stress_2,stress_3,stress_4]
		mean_stress = np.mean(stress_list,axis=0)
		std_stress = np.std(stress_list,axis=0)
		#print(mean_stress)
		#color = np.random.rand(3,)
		#axss.errorbar(strain,mean_stress,std_stress,label = complexity[i])
		end_stress.append(mean_stress[-1])
		end_std.append(std_stress[-1])
		mean_vertices= np.mean(number_vertices[k:k+4])
		mean_n_vertices.append(mean_vertices)
		#plot_second_der(axssd,strain, mean_stress,color = color)
		"""
		color = np.random.rand(3,)
		axss.plot(strain, stress, label='%s, %s' % (creation, generation),color = color)
		plot_second_der(axssd,strain, stress,color = color)
		axom.plot(strain_omega,omega_xx, label='%s, %s' % (creation, generation))"""
	axss.errorbar(mean_n_vertices,end_stress,end_std)
	if type == 'global_length':
		axss.set_ylabel(r'$\frac{1}{L}*\sum_{bn} f_k x_f$ at 100% strain',fontsize=20)
	elif type == 'number_bc':
		axss.set_ylabel(r'$\frac{1}{n_{bn}}*\sum_{bn} f_k x_f$ at 100% strain',fontsize=20)
	elif type == 'mean_length':
		axss.set_ylabel(r'$\frac{1}{l_0}*\sum_{bn} f_k x_f$ at 100% strain',fontsize=20)
	elif type == 'number_ridge':
		axss.set_ylabel(r'$\frac{1}{l_0 * n_{ridges}}*\sum_{bn} f_k x_f$ at 100% strain',fontsize=20)
	axss.set_xlabel('Number of points in network',fontsize=20)
	#axss.set_xlabel(r'$\frac{L-L_0}{L_0}$',fontsize=20)
	fig_stress_strain.tight_layout()
	#axss.legend(loc = 'upper left')
	fig_stress_strain.savefig('comp_networks_figure_stress_strain_%s.pdf' % type)


plt.show()
os.chdir(current_dir)