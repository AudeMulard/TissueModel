import matplotlib.pyplot as plt
from scipy import stats
import fnmatch
import os, sys, csv
import numpy as np
import scipy
import statistics
import seaborn as sns
sys.path.append('/users/am2548/TissueModel/Documents/PhD/Code/TissueModel')
from Plotting.information_network import load_network_info, length_square, calculate_network_data

list_modes = [['Voronoi','random'],['Voronoi','grid'],['growth_network','grid']]#,['Voronoi','regular'],['growth_network','grid']]#,['growth_network','grid']]

############################################### NUMERICAL OUTPUTS #################################################

def plot_initial_info(filename):
	print(filename)
	network = load_network_info(int(filename[0][21:24]))
	length_ini = []
	for ridge in network.ridge_vertices:
		length_ini.append(length_square(network.vertices_ini[ridge[1]]-network.vertices_ini[ridge[0]]))
	lengths_ini, stretch_ratio_ini, cos_theta_square_ini, omega_xx_ini = calculate_network_data(network,length_ini)
	return network, length_ini, stretch_ratio_ini, cos_theta_square_ini, omega_xx_ini 
	
def stress_strain_info(creation, generation,name):
	filename = fnmatch.filter(os.listdir('.'), 'stress_strain_*_%s_%s_%s.csv' % (creation,generation,name))
	with open(filename[0], 'r') as readFile:
		reader = csv.reader(readFile)
		curve = np.array(list(reader))
		stress = [float(i) for i in curve[1]]
		strain=[float(i) for i in curve[0]]
	print(strain)
	return stress, strain

def angles_info_global(network,lengths_ini,filenames, cos_theta_square_ini,name):
	omega_xx = []
	cos_theta_square = []
	corr_lambda_l0, corr_theta_0, corr_lamd_theta, theta= [],[],[],[]
	theta_0 = [np.arccos(scipy.sqrt(i)) for i in cos_theta_square_ini]
	for step in range(len(filenames)):
		lengths, stretch_ratio, cos_theta_square_step, omega_xx_step = calculate_network_data(network,lengths_ini,step = step,name = '%s_%s_%s' % (creation,generation,name))
		omega_xx.append(omega_xx_step)
		cos_theta_square.append(cos_theta_square_step)
		corr_lambda_l0.append(statistics.mean(scipy.stats.pearsonr(stretch_ratio, lengths_ini)))
		for i in cos_theta_square:
			theta =  list(np.arccos(scipy.sqrt(i)))
		corr_theta_0.append(statistics.mean(scipy.stats.pearsonr(theta, theta_0)))
		corr_lamd_theta.append(statistics.mean(scipy.stats.pearsonr(stretch_ratio, theta_0)))
	return lengths, stretch_ratio, cos_theta_square, omega_xx, corr_lambda_l0, corr_theta_0, corr_lamd_theta


def sorted_ls(path):
    mtime = lambda f: os.stat(os.path.join(path, f)).st_mtime
    return list(sorted(os.listdir(path), key=mtime))

os.chdir('../Data/Study_networks/')
if len(sys.argv) == 1:
	os.chdir(sorted_ls('.')[-1])
else:
	os.chdir(sys.argv[1])
print(os.getcwd())


########################## GRAPHICALS OUTPUTS #################################

# Initial info: segment length and angles, alone then cumulative curves together
def graph_comp_ini_end_info(name,creation, generation):
	filename = fnmatch.filter(os.listdir('.'), 'parameters_01_*.rpt' % (creation,generation,name))
	network, lengths_ini, stretch_ratio_ini, cos_theta_square_ini, omega_xx_ini = plot_initial_info(filename)
	fig_ini_end, axs = plt.subplots(2,2,sharex='row',sharey='row')
	sns.distplot(lengths_ini,ax = axs[0,0])
	axs[0,0].set_xlabel('0% strain lengths')
	sns.distplot(stretch_ratio_ini, ax = axs[1,0])
	axs[1,0].set_xlabel('1% strain stretch ratios')
	
	filenames = fnmatch.filter(os.listdir('.'), 'network_vertices_*_%s_%s_%s.csv' % (creation,generation,name))
	lengths, stretch_ratio, cos_theta_square, omega_xx, corr_lambda_l0, corr_theta_0, corr_lamd_theta= angles_info_global(network,lengths_ini,filenames,cos_theta_square_ini,name)
	sns.distplot(lengths, ax=axs[0,1])
	axs[0,1].set_xlabel('50% strain lengths')
	sns.distplot(stretch_ratio,ax=axs[1,1])
	axs[1,1].set_xlabel('50% strain stretch ratios')
	plt.subplots_adjust(hspace=0.5)
	fig_ini_end.suptitle('%s %s' % (creation, generation))
	fig_ini_end.savefig('ini_end_%s_%s.pdf' % (creation, generation))
	plt.close()
	return lengths_ini,cos_theta_square_ini, lengths, stretch_ratio, cos_theta_square, omega_xx, corr_lambda_l0, corr_theta_0, corr_lamd_theta

def correlation_figures(name, creation, generation, corr_lambda_l0, corr_theta_0, corr_lamd_theta):
	fig_corr, axs = plt.subplots(2,2)
	axs[0,0].scatter(strain,corr_lambda_l0)
	axs[0,0].set_xlabel(r'corr($\lambda$,$\lambda_0$)')
	axs[1,0].scatter(strain,corr_theta_0)
	axs[1,0].set_xlabel(r'corr($\theta_0$,$\frac{L-L_0}{L_0}$)')
	axs[0,1].scatter(strain,corr_lamd_theta)
	axs[0,1].set_xlabel(r'corr($\lambda$,$\theta$)')
	plt.subplots_adjust(hspace=0.5)
	fig_corr.suptitle('%s %s' % (creation, generation))
	fig_corr.savefig('corr_%s_%s_%s.pdf' % (creation, generation,name))
	plt.close()


def micro_info(axct,axcl,axlf,axll,creation, generation,cos_theta_square_ini,cos_theta_square,stretch_ratio,lengths_ini,name):
	if name == 'time': title = '%s %s' % (creation, generation)
	else: title = name
	a = int(i/2.)
	b = i % 2
	# cos sq theta on cos sq theta0
	print(len(cos_theta_square_ini), len(cos_theta_square[-1]))
	axct[a,b].scatter(cos_theta_square_ini, cos_theta_square[-1])
	axct[a,b].set_xlabel(r'$cos^2\theta_i$')
	axct[a,b].set_ylabel(r'$cos^2\theta$')
	axct[a,b].set_title(title)
	# stretch ration on cos sq theta0
	axcl[a,b].scatter(cos_theta_square_ini,stretch_ratio)
	axcl[a,b].set_xlabel(r'$cos^2\theta_i$')
	axcl[a,b].set_ylabel(r'$\lambda$')
	axcl[a,b].set_title(title)
	# stretch ratio frequency
	sns.distplot(stretch_ratio,ax=axlf[a,b])
	axlf[a,b].set_xlabel(r'$\lambda$')
	axlf[a,b].set_title(title)
	# stretch ratio on segment length initial
	axll[a,b].scatter(lengths_ini, stretch_ratio)
	axll[a,b].set_xlabel('Initial length')
	axll[a,b].set_ylabel(r'$\lambda$')
	axll[a,b].set_title(title)

"""
#################### COMPARISON OF DIFFERENT COMPLEXITY #########################

creation, generation = 'Voronoi', 'random'

complexities = [10,50,100,200]

figure_omega_xx,axom = plt.subplots()
fig_stress_strain,axss = plt.subplots()
fig_cos_theta,axct = plt.subplots(1,len(complexities))
fig_costh_lambda,axcl = plt.subplots(1,len(complexities))
fig_lambda_freq,axlf = plt.subplots(1,len(complexities))
fig_lambda_length,axll = plt.subplots(1,len(complexities))

for i in range(len(complexities)):
	complexity = complexities[i]
	
	# Initial info: segment length and angles, alone then cumulative curves together
	lengths_ini,cos_theta_square_ini,lengths, stretch_ratio, cos_theta_square, omega_xx, corr_lambda_l0, corr_theta_0, corr_lamd_theta=graph_comp_ini_end_info('complexity_%04d' % complexity ,creation, generation)
	
	# Stress-strain responses and omega xx responses
	stress, strain = stress_strain_info(creation, generation,'complexity_%04d' % complexity)
	axss.plot(strain, stress, label='%s' % complexity)
	axom.plot(strain,omega_xx, label='%s' % complexity)

	# Information on the microscopic geometry
	micro_info(axct,axcl,axlf,axll,creation, generation,cos_theta_square_ini,cos_theta_square,stretch_ratio,lengths_ini,'%04d' % complexity)

	# Correlation figures: # comparison: correlation theta/theta0; correlation lambdaf/L0; correlation lambdaf/theta0 all on strain
	correlation_figures('complexity_%04d' % complexity, creation, generation, corr_lambda_l0, corr_theta_0, corr_lamd_theta)

plt.subplots_adjust(hspace=0.2)
fig_lambda_length.savefig('complexity_lambda_length_%s_%s.pdf' % (creation, generation))
plt.close()

plt.subplots_adjust(hspace=0.5)
fig_lambda_freq.savefig('complexity_lambda_freq_%s_%s.pdf' % (creation, generation))
plt.close()

plt.subplots_adjust(hspace=0.5)
fig_costh_lambda.savefig('complexity_costh_lambda_%s_%s.pdf' % (creation, generation))
plt.close()
	
plt.subplots_adjust(hspace=0.5)
fig_cos_theta.savefig('complexity_cos_theta_%s_%s.pdf' % (creation, generation))
plt.close()

axom.set_xlabel(r'$\frac{L-L_0}{L_0}$',fontsize=20)
axom.set_ylabel(r'$\Omega_{xx}$',fontsize=20)
figure_omega_xx.tight_layout()
plt.legend()
figure_omega_xx.savefig('complexity_figure_omega_xx_%s_%s.pdf' % (creation, generation))
plt.close()

axss.set_xlabel(r'$\frac{L-L_0}{L_0}$',fontsize=20)
axss.set_ylabel(r'$\frac{\sigma}{k/\bar{l_0}}$',fontsize=20)
fig_stress_strain.tight_layout()
plt.legend(loc = 'upper left')
fig_stress_strain.savefig('complexity_figure_stress_strain_%s_%s.pdf' % (creation, generation))

"""
#################### COMPARISON OF DIFFERENT NETWORKS #########################
figure_omega_xx,axom = plt.subplots()
fig_stress_strain,axss = plt.subplots()
print(int(len(list_modes)/2.))
fig_cos_theta,axct = plt.subplots(int(len(list_modes)/2.)+1,2)
fig_costh_lambda,axcl = plt.subplots(int(len(list_modes)/2.)+1,2)
fig_lambda_freq,axlf = plt.subplots(int(len(list_modes)/2.)+1,2)
fig_lambda_length,axll = plt.subplots(int(len(list_modes)/2.)+1,2)

for i in range(len(list_modes)):
	creation, generation = list_modes[i]
	
	# Initial info: segment length and angles, alone then cumulative curves together
	lengths_ini,cos_theta_square_ini,lengths, stretch_ratio, cos_theta_square, omega_xx, corr_lambda_l0, corr_theta_0, corr_lamd_theta=graph_comp_ini_end_info('time',creation, generation)
	
	# Stress-strain responses and omega xx responses
	stress, strain = stress_strain_info(creation, generation,'time')
	axss.plot(strain, stress, label='%s, %s' % (creation, generation))
	axom.plot(strain,omega_xx, label='%s, %s' % (creation, generation))

	# Information on the microscopic geometry
	micro_info(axct,axcl,axlf,axll,creation, generation,cos_theta_square_ini,cos_theta_square,stretch_ratio,lengths_ini,'time')

	# Correlation figures: # comparison: correlation theta/theta0; correlation lambdaf/L0; correlation lambdaf/theta0 all on strain
	correlation_figures('time', creation, generation, corr_lambda_l0, corr_theta_0, corr_lamd_theta)


plt.subplots_adjust(hspace=0.5,wspace=0.5)
fig_lambda_length.savefig('comp_networks_lambda_length.pdf')
plt.close()

plt.subplots_adjust(hspace=0.5,wspace=0.5)
fig_lambda_freq.savefig('comp_networks_lambda_freq.pdf')
plt.close()

plt.subplots_adjust(hspace=0.5,wspace=0.5)
fig_costh_lambda.savefig('comp_networks_costh_lambda.pdf')
plt.close()
	
plt.subplots_adjust(hspace=0.5,wspace=1.0)
fig_cos_theta.savefig('comp_networks_cos_theta.pdf')
plt.close()

axom.set_xlabel(r'$\frac{L-L_0}{L_0}$',fontsize=20)
axom.set_ylabel(r'$\Omega_{xx}$',fontsize=20)
axom.legend(loc = 'upper left')
figure_omega_xx.tight_layout()
figure_omega_xx.savefig('comp_networks_figure_omega_xx.pdf')
plt.close()

axss.set_xlabel(r'$\frac{L-L_0}{L_0}$',fontsize=20)
axss.set_ylabel(r'$\frac{\sigma}{k/\bar{l_0}}$',fontsize=20)
fig_stress_strain.tight_layout()
axss.legend(loc = 'upper left')
fig_stress_strain.savefig('comp_networks_figure_stress_strain.pdf')


