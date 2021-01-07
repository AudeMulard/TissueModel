import matplotlib.pyplot as plt
from scipy import stats
import fnmatch
import os, sys, csv
import numpy as np
import scipy
import statistics
import seaborn as sns
sys.path.append('/users/am2548/TissueModel/Documents/PhD/Code/TissueModel')
from Plotting.information_network import *

list_modes = [['growth_network','grid'],['Voronoi','random'],['growth_network','grid'],['Voronoi','random'],['growth_network','grid'],['Voronoi','random']]

############################################### NUMERICAL OUTPUTS #################################################

def plot_initial_info(filename):
	print(filename)
	network = load_network_info(int(filename[0][-9:-4]))
	from Plotting.network_plotting import plot_geometry
	length_ini = []
	for ridge in network.ridge_vertices:
		length_ini.append(np.sqrt(length_square(network.vertices_ini[ridge[1]]-network.vertices_ini[ridge[0]])))
	lengths_ini, stretch_ratio_ini, cos_theta_square_ini, omega_xx_ini,strain_omega_step = calculate_network_data(network,length_ini)
	return network, length_ini, stretch_ratio_ini, cos_theta_square_ini, omega_xx_ini 

def angles_info_global(network,lengths_ini,filenames, cos_theta_square_ini,name,test_number):
	import pandas as pd
	omega_xx = []
	cos_theta_square = []
	corr_lambda_l0, corr_theta_0, corr_lamd_theta, theta= [],[],[],[]
	theta_0 = [np.arccos(np.sqrt(i)) for i in cos_theta_square_ini]
	node_label = write_node_label(test_number,network)
	strain_omega=[]
	if len(filenames)==0: lengths, stretch_ratio, cos_theta_square, omega_xx, strain_omega, corr_lambda_l0, corr_theta_0, corr_lamd_theta=[],[],[],[],[],[],[],[]
	for step in range(len(filenames)):
		lengths, stretch_ratio, cos_theta_square_step, omega_xx_step,strain_omega_step = calculate_network_data(network,lengths_ini,node_label=node_label,step = step,name = '%s_%s_%s' % (creation,generation,name),test_number = test_number)
		omega_xx.append(omega_xx_step)
		strain_omega.append(strain_omega_step)
		cos_theta_square.append(cos_theta_square_step)
		corr_lambda_l0.append(statistics.mean(scipy.stats.pearsonr(stretch_ratio, lengths_ini)))
		for i in cos_theta_square:
			theta =  list(np.arccos(np.sqrt(i)))
		corr_theta_0.append(statistics.mean(scipy.stats.pearsonr(theta, theta_0)))
		corr_lamd_theta.append(statistics.mean(scipy.stats.pearsonr(stretch_ratio, theta_0)))
	return lengths, stretch_ratio, cos_theta_square, omega_xx, strain_omega, corr_lambda_l0, corr_theta_0, corr_lamd_theta


def sorted_ls(path):
    mtime = lambda f: os.stat(os.path.join(path, f)).st_mtime
    return list(sorted(os.listdir(path), key=mtime))
current_dir=os.getcwd()
os.chdir('../Data/Study_networks/')

if len(sys.argv) == 1:
	os.chdir(sorted_ls('.')[-1])
else:
	os.chdir(sys.argv[1])
print(os.getcwd())


########################## GRAPHICALS OUTPUTS #################################

# Initial info: segment length and angles, alone then cumulative curves together
def graph_comp_ini_end_info(name,number, generation,test_number):
	filename = fnmatch.filter(os.listdir('.'), 'parameters_%02d_*.csv' % number) 
	
	network, lengths_ini, stretch_ratio_ini, cos_theta_square_ini, omega_xx_ini = plot_initial_info(filename)
	fig_ini_end, axs = plt.subplots(2,2,sharex='row',sharey='row')
	sns.histplot(lengths_ini,ax = axs[0,0],kde='True')
	axs[0,0].tick_params(labelsize=8)
	axs[0,0].set_xlabel('0% strain lengths')
	sns.histplot(stretch_ratio_ini, ax = axs[1,0],kde='True')
	axs[1,0].tick_params(labelsize=8)
	axs[1,0].set_xlabel('1% strain stretch ratios')
	
	filenames = fnmatch.filter(sorted_ls('.'), 'network_vertices_*_%s.csv' % test_number)
	lengths, stretch_ratio, cos_theta_square, omega_xx, strain_omega, corr_lambda_l0, corr_theta_0, corr_lamd_theta= angles_info_global(network,lengths_ini,filenames,cos_theta_square_ini,name,test_number)
	sns.histplot(lengths, ax=axs[0,1],kde='True')
	axs[0,1].tick_params(labelsize=8)
	axs[0,1].set_xlabel('50% strain lengths')
	sns.histplot(stretch_ratio,ax=axs[1,1],kde='True')
	axs[1,1].tick_params(labelsize=8)
	axs[1,1].set_xlabel('50% strain stretch ratios')
	plt.subplots_adjust(hspace=1.0)
	fig_ini_end.suptitle('%s %s' % (network.creation, network.generation))
	#if name=='complexity': label=len(network.vertices)
	fig_ini_end.savefig('ini_end_%s_%s.png' % (network.creation, number))
	plt.close()
	return network,lengths_ini,cos_theta_square_ini, lengths, stretch_ratio, cos_theta_square, omega_xx, strain_omega,corr_lambda_l0, corr_theta_0, corr_lamd_theta

def correlation_figures(name, creation, generation, corr_lambda_l0, corr_theta_0, corr_lamd_theta):
	fig_corr, axs = plt.subplots(2,2)
	axs[0,0].tick_params(labelsize=8)
	axs[0,0].scatter(strain_omega,corr_lambda_l0)
	axs[0,0].set_xlabel(r'corr($\lambda$,$\lambda_0$)')
	axs[1,0].scatter(strain_omega,corr_theta_0)
	axs[1,0].tick_params(labelsize=8)
	axs[1,0].set_xlabel(r'corr($\theta_0$,$\frac{L-L_0}{L_0}$)')
	axs[0,1].tick_params(labelsize=8)
	axs[0,1].scatter(strain_omega,corr_lamd_theta)
	axs[0,1].set_xlabel(r'corr($\lambda$,$\theta$)')
	fig_corr.suptitle('%s %s' % (creation, generation))
	fig_corr.savefig('corr_%s_%s_%s.png' % (creation, generation,name))
	plt.subplots_adjust(hspace=1.0)
	plt.close()


def micro_info(axct,axcl,axlf,axll,creation, generation,cos_theta_square_ini,cos_theta_square,stretch_ratio,lengths_ini,name):
	if name == 'time': title = '%s %s' % (creation, generation)
	else: title = name
	a = int(i/2.)
	b = i % 2
	
	# cos sq theta on cos sq theta0
	axct[a,b].tick_params(labelsize=8)
	axct[a,b].scatter(cos_theta_square_ini, cos_theta_square[-1])
	axct[a,b].set_xlabel(r'$cos^2\theta_i$',labelpad=-2,fontsize=8)
	axct[a,b].set_ylabel(r'$cos^2\theta$',labelpad=-2,fontsize=8)
	axct[a,b].set_title(title,fontsize=10)
	# stretch ration on cos sq theta0
	axcl[a,b].tick_params(labelsize=8)
	axcl[a,b].scatter(cos_theta_square_ini,stretch_ratio)
	axcl[a,b].set_xlabel(r'$cos^2\theta_i$',labelpad=-2,fontsize=8)
	axcl[a,b].set_ylabel(r'$\lambda$',labelpad=-2,fontsize=8)
	axcl[a,b].set_title(title,fontsize=10)
	# stretch ratio frequency
	sns.histplot(stretch_ratio,ax=axlf[a,b],kde='True')
	axlf[a,b].tick_params(labelsize=8)
	axlf[a,b].set_xlabel(r'$\lambda$',labelpad=-2,fontsize=8)
	axlf[a,b].set_title(title,fontsize=10)
	# stretch ratio on segment length initial
	axll[a,b].tick_params(labelsize=8)
	axll[a,b].scatter(lengths_ini, stretch_ratio)
	axll[a,b].set_xlabel('Initial length',fontsize=8)
	axll[a,b].set_ylabel(r'$\lambda$',labelpad=-2,fontsize=8)
	axll[a,b].set_title(title,fontsize=10)
	"""
		# cos sq theta on cos sq theta0
	axct[i].tick_params(labelsize=8)
	axct[i].scatter(cos_theta_square_ini, cos_theta_square[-1])
	axct[i].set_xlabel(r'$cos^2\theta_i$',labelpad=-2,fontsize=8)
	axct[i].set_ylabel(r'$cos^2\theta$',labelpad=-2,fontsize=8)
	axct[i].set_title(title,fontsize=10)
	# stretch ration on cos sq theta0
	axcl[i].tick_params(labelsize=8)
	axcl[i].scatter(cos_theta_square_ini,stretch_ratio)
	axcl[i].set_xlabel(r'$cos^2\theta_i$',labelpad=-2,fontsize=8)
	axcl[i].set_ylabel(r'$\lambda$',labelpad=-2,fontsize=8)
	axcl[i].set_title(title,fontsize=10)
	# stretch ratio frequency
	sns.histplot(stretch_ratio,ax=axlf[i],kde='True')
	axlf[i].tick_params(labelsize=8)
	axlf[i].set_xlabel(r'$\lambda$',labelpad=-2,fontsize=8)
	axlf[i].set_title(title,fontsize=10)
	# stretch ratio on segment length initial
	axll[i].tick_params(labelsize=8)
	axll[i].scatter(lengths_ini, stretch_ratio)
	axll[i].set_xlabel('Initial length',labelpad=-2,fontsize=8)
	axll[i].set_ylabel(r'$\lambda$',labelpad=0,fontsize=8)
	axll[i].set_title(title,fontsize=10)
"""

#################### COMPARISON OF DIFFERENT COMPLEXITY #########################

creation, generation = 'Voronoi', 'random'

#complexities = [10,50,100,200]
#parameter_values=[0.001,0.005,0.01,0.05,0.1,0.5,1.0,1.1]
parameter_values=range(5)#, 0.06,0.07,0.08,0.09]#,0.1]
#parameter_name = 'truss_area'
#parameter_values = [1,2,3,4,5,6]#,7]
#parameter_values=[0,1,2,3]
#parameter_values = [0.01,0.1,1.0]
parameter_name = 'complexity'


figure_omega_xx,axom = plt.subplots()
fig_stress_strain,axss = plt.subplots()
axssd =axss.twinx()
axssd.set_ylabel('second derivative')
axssd.tick_params(axis='y')
if len(parameter_values)==3:
	fig_cos_theta,axct = plt.subplots(1,3,figsize=(6,2.5))
	fig_costh_lambda,axcl = plt.subplots(1,3,figsize=(6,2.5))
	fig_lambda_freq,axlf = plt.subplots(1,3,figsize=(6,2.5))
	fig_lambda_length,axll = plt.subplots(1,3,figsize=(6,2.5))
else:
	fig_cos_theta,axct = plt.subplots(int(len(parameter_values)/2.)+len(parameter_values)%2,2,figsize=(8, 8))
	fig_costh_lambda,axcl = plt.subplots(int(len(parameter_values)/2.)+len(parameter_values)%2,2)
	fig_lambda_freq,axlf = plt.subplots(int(len(parameter_values)/2.)+len(parameter_values)%2,2)
	fig_lambda_length,axll = plt.subplots(int(len(parameter_values)/2.)+len(parameter_values)%2,2)

for i in range(5):
	parameter = parameter_values[i]
	test_number = fnmatch.filter(sorted_ls('.'), 'network_vertices_01_00_*.csv')[i][-13:-4]
	print(test_number)
	# Initial info: segment length and angles, alone then cumulative curves together
	network,lengths_ini,cos_theta_square_ini,lengths, stretch_ratio, cos_theta_square, omega_xx, strain_omega,corr_lambda_l0, corr_theta_0, corr_lamd_theta=graph_comp_ini_end_info('%s' % (parameter_name) ,i, generation,test_number)
	
	# Stress-strain responses and omega xx responses
	strain,stress = stress_strain_curve(test_number,network)
	test_1 = load_info_test(len(network.vertices))
	
	#strain = [i/test_1.iterations * test_1.traction_distance for i in range(test_1.iterations+1)]
	if parameter_name == 'hyperstatic_param': label_name = int(len(network.ridge_vertices)/len(network.vertices)*1e3)/1e3
	elif parameter_name == 'complexity': label_name = len(network.vertices)
	else: label_name = parameter
	color = np.random.rand(3,)
	axss.plot(strain, stress, label='%s' % label_name,color=color)
	print(parameter)
	#if parameter!=2 and parameter!=0.3: 
	#	print('I am printing',parameter)
	plot_second_der(axssd,strain, stress,color = color)
	axom.plot(strain_omega,omega_xx, label='%s' % label_name,color=color)
	# Information on the microscopic geometry
	micro_info(axct,axcl,axlf,axll,creation, generation,cos_theta_square_ini,cos_theta_square,stretch_ratio,lengths_ini,'%s' % str(label_name))

	# Correlation figures: # comparison: correlation theta/theta0; correlation lambdaf/L0; correlation lambdaf/theta0 all on strain
	correlation_figures('%s_%04d' % (parameter_name,label_name), creation, generation, corr_lambda_l0, corr_theta_0, corr_lamd_theta)

plt.subplots_adjust(wspace=0.5)
fig_lambda_length.savefig('%s_lambda_length_%s_%s.png' % (parameter_name,creation, generation))
plt.tight_layout()
plt.close()

plt.subplots_adjust(wspace=0.5)
fig_lambda_freq.savefig('%s_lambda_freq_%s_%s.png' % (parameter_name,creation, generation))
plt.tight_layout()
plt.close()

plt.subplots_adjust(wspace=1.0)
fig_costh_lambda.savefig('%s_costh_lambda_%s_%s.png' % (parameter_name,creation, generation))
plt.tight_layout()
plt.close()

plt.subplots_adjust(wspace=1.0)
fig_cos_theta.savefig('%s_cos_theta_%s_%s.png' % (parameter_name,creation, generation))
plt.tight_layout()
plt.close()

axom.set_xlabel(r'$\frac{L-L_0}{L_0}$',fontsize=20)
axom.set_ylabel(r'$\Omega_{xx}$',fontsize=20)
figure_omega_xx.tight_layout()
axom.legend()
figure_omega_xx.savefig('%s_figure_omega_xx_%s_%s.png' % (parameter_name,creation, generation))
plt.close()

axss.set_xlabel(r'$\frac{L-L_0}{L_0}$',fontsize=20)
axss.set_ylabel(r'$\frac{\sigma_{xx}}{E}$',fontsize=20)
fig_stress_strain.tight_layout()
axss.legend(loc = 'upper left')
fig_stress_strain.savefig('%s_figure_stress_strain_%s_%s.png' % (parameter_name,creation, generation))

"""
#################### COMPARISON OF DIFFERENT NETWORKS #########################
figure_omega_xx,axom = plt.subplots()
fig_stress_strain,axss = plt.subplots()
axssd =axss.twinx()
axssd.set_ylabel('second derivative')
axssd.tick_params(axis='y')
print(int(len(list_modes)/2.))
fig_cos_theta,axct = plt.subplots(int(len(list_modes)/2.+1),2)
fig_costh_lambda,axcl = plt.subplots(int(len(list_modes)/2.+1),2)
fig_lambda_freq,axlf = plt.subplots(int(len(list_modes)/2.+1),2)
fig_lambda_length,axll = plt.subplots(int(len(list_modes)/2.+1),2)


stress_V=np.zeros((3,20))
stress_G=np.zeros((3,20))
omega_V=np.zeros((3,11))
omega_G=np.zeros((3,11))
for i in range(len(list_modes)):
	creation, generation = list_modes[i]
	test_number = fnmatch.filter(sorted_ls('.'), 'node_label_*.csv')[i][-13:-4]
	print(test_number)
	# Initial info: segment length and angles, alone then cumulative curves together
	network,lengths_ini,cos_theta_square_ini,lengths, stretch_ratio, cos_theta_square, omega_xx,strain_omega, corr_lambda_l0, corr_theta_0, corr_lamd_theta=graph_comp_ini_end_info('time',i, generation,test_number)
	if creation=='Voronoi': omega_V[int(i/2)-1,:len(omega_xx)]=omega_xx
	if creation=='growth_network': omega_G[int(i/2),:len(omega_xx)]=omega_xx
	# Stress-strain responses and omega xx responses
	strain,stress = stress_strain_curve(test_number,network)
	test_1 = load_info_test(len(network.vertices))
	if creation=='Voronoi':  stress_V[int(i/2)-1,:len(stress)]=stress
	if creation=='growth_network': stress_G[int(i/2),:len(stress)]=stress
	color = np.random.rand(3,)
	axss.plot(strain, stress, label='%s, %s' % (creation, generation),color = color)
	plot_second_der(axssd,strain, stress,color = color)
	axom.plot(strain_omega,omega_xx, label='%s, %s' % (creation, generation))

	# Information on the microscopic geometry
	micro_info(axct,axcl,axlf,axll,creation, generation,cos_theta_square_ini,cos_theta_square,stretch_ratio,lengths_ini,'time')

	# Correlation figures: # comparison: correlation theta/theta0; correlation lambdaf/L0; correlation lambdaf/theta0 all on strain
	correlation_figures('time', creation, generation, corr_lambda_l0, corr_theta_0, corr_lamd_theta)


stress_G_mean=[]
stress_G_std=[]
for i in range(len(stress_G[0])):
	stress_G_mean.append(np.mean(stress_G[:,i]))
	stress_G_std.append(np.std(stress_G[:,i]))

axss.errorbar(strain,stress_G_mean,stress_G_std,label='Growth')


stress_V=np.array(stress_V)
stress_V_mean=[]
stress_V_std=[]
for i in range(len(stress_V[0])):
	stress_V_mean.append(np.mean(stress_V[:,i]))
	stress_V_std.append(np.std(stress_V[:,i]))
axss.errorbar(strain,stress_V_mean,stress_V_std,label='Voronoi')

omega_G=np.array(omega_G)
omega_G_mean=[]
omega_G_std=[]
for i in range(len(omega_G[0])):
	omega_G_mean.append(np.mean(omega_G[:,i]))
	omega_G_std.append(np.std(omega_G[:,i]))
print(len(strain_omega),len(omega_G_mean),len(omega_G_std))
axom.errorbar(strain_omega,omega_G_mean,omega_G_std,label='Growth')


omega_V=np.array(omega_V)
omega_V_mean=[]
omega_V_std=[]
for i in range(len(omega_G[0])):
	omega_V_mean.append(np.mean(omega_V[:,i]))
	omega_V_std.append(np.std(omega_V[:,i]))
axom.errorbar(strain_omega,omega_V_mean,omega_V_std,label='Voronoi')

plt.subplots_adjust(hspace=0.5,wspace=0.5)
fig_lambda_length.savefig('comp_networks_lambda_length.png')
plt.close()

plt.subplots_adjust(hspace=0.5,wspace=0.5)
fig_lambda_freq.savefig('comp_networks_lambda_freq.png')
plt.close()

plt.subplots_adjust(hspace=0.5,wspace=0.5)
fig_costh_lambda.savefig('comp_networks_costh_lambda.png')
plt.close()
	
plt.subplots_adjust(hspace=0.5,wspace=1.0)
fig_cos_theta.savefig('comp_networks_cos_theta.png')
plt.close()

axom.set_xlabel(r'$\frac{L-L_0}{L_0}$',fontsize=20)
axom.set_ylabel(r'$\Omega_{xx}$',fontsize=20)
axom.legend(loc = 'upper left')
figure_omega_xx.tight_layout()
figure_omega_xx.savefig('comp_networks_figure_omega_xx.png')
plt.close()

axss.set_xlabel(r'$\frac{L-L_0}{L_0}$',fontsize=20)
axss.set_ylabel(r'$\frac{\sigma}{k/\bar{l_0}}$',fontsize=20)
fig_stress_strain.tight_layout()
axss.legend(loc = 'upper left')
fig_stress_strain.savefig('comp_networks_figure_stress_strain.png')

os.chdir(current_dir)
"""