import matplotlib.pyplot as plt
from scipy import stats
import fnmatch
import os, sys, csv
import numpy as np
import scipy
import statistics
sys.path.append('/home/aude/Documents/PhD/Code/TissueModel/')
from Plotting.information_network import load_network_info, length_square, calculate_network_data

list_modes = [['Voronoi','random'],['Voronoi','grid']]#,['Voronoi','regular'],['growth_network','grid']]#,['growth_network','grid']]

############################################### NUMERICAL OUTPUTS #################################################

def plot_initial_info(filename):
	print filename[0]
	network = load_network_info(int(filename[0][21:24]))
	length_ini = []
	for ridge in network.ridge_vertices:
		length_ini.append(length_square(network.vertices_ini[ridge[1]]-network.vertices_ini[ridge[0]]))
	lengths_ini, stretch_ratio_ini, cos_theta_square_ini, omega_xx_ini = calculate_network_data(network,length_ini)
	return network, length_ini, stretch_ratio_ini, cos_theta_square_ini, omega_xx_ini 
	
def stress_strain_info(creation, generation):
	filename = fnmatch.filter(os.listdir('.'), 'stress_strain_*_%s_%s_time.csv' % (creation,generation))
	with open(filename[0], 'r') as readFile:
		reader = csv.reader(readFile)
		curve = np.array(list(reader))
		stress = [float(i) for i in curve[1]]
		strain=[float(i) for i in curve[0]]
	return stress, strain

def angles_info_global(network,lengths_ini,filenames, cos_theta_square_ini):
	omega_xx = []
	cos_theta_square = []
	corr_lambda_l0, corr_theta_0, corr_lamd_theta, theta= [],[],[],[]
	theta_0 = [np.arccos(scipy.sqrt(i)) for i in cos_theta_square_ini]
	for step in range(len(filenames)):
		lengths, stretch_ratio, cos_theta_square_step, omega_xx_step = calculate_network_data(network,lengths_ini,step = step,name = '%s_%s_time' % (creation,generation))
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
figure_omega_xx,axom = plt.subplots()
fig_stress_strain,axss = plt.subplots()




for mode in list_modes:
	creation, generation = mode
	
	# Initial info: segment length and angles, alone then cumulative curves together
	filename = fnmatch.filter(os.listdir('.'), 'network_vertices_000_*_%s_%s_time.csv' % (creation,generation))
	network, lengths_ini, stretch_ratio_ini, cos_theta_square_ini, omega_xx_ini = plot_initial_info(filename)
	fig_ini_end, axs = plt.subplots(2,2,sharex='row',sharey='row')
	import seaborn as sns
	sns.distplot(lengths_ini,ax = axs[0,0])
	axs[0,0].set_xlabel('0% strain lengths')
	sns.distplot(stretch_ratio_ini, ax = axs[1,0])
	axs[1,0].set_xlabel('1% strain stretch ratios')
	
	# Stress-strain responses and omega xx responses
	filenames = fnmatch.filter(os.listdir('.'), 'network_vertices_*_%s_%s_time.csv' % (creation,generation))
	lengths, stretch_ratio, cos_theta_square, omega_xx, corr_lambda_l0, corr_theta_0, corr_lamd_theta= angles_info_global(network,lengths_ini,filenames,cos_theta_square_ini)
	sns.distplot(lengths, ax=axs[0,1])
	axs[0,1].set_xlabel('50% strain lengths')
	sns.distplot(stretch_ratio,ax=axs[1,1])
	axs[1,1].set_xlabel('50% strain stretch ratios')
	plt.subplots_adjust(hspace=0.5)
	fig_ini_end.suptitle('%s %s' % (creation, generation))
	fig_ini_end.savefig('ini_end_%s_%s.pdf' % (creation, generation))
	plt.close()
	
	stress, strain = stress_strain_info(creation, generation)
	axss.plot(strain, stress, label='%s, %s' % (creation, generation))
	axss.set_xlabel(r'$\frac{L-L_0}{L_0}$',fontsize=20)
	axss.set_ylabel(r'$\frac{\sigma}{k/\bar{l_0}}$',fontsize=20)
	axom.plot(strain,omega_xx, label='%s, %s' % (creation, generation))
	axom.set_xlabel(r'$\frac{L-L_0}{L_0}$',fontsize=20)
	axom.set_ylabel(r'$\Omega_{xx}$',fontsize=20)

	# cos sq theta on cos sq theta0
	fig_cos_theta,axct = plt.subplots()
	axct.scatter(cos_theta_square_ini, cos_theta_square[-1])
	axct.set_xlabel(r'$cos^2\theta_i$')
	axct.set_ylabel(r'$cos^2\theta$')
	plt.subplots_adjust(hspace=0.5)
	fig_cos_theta.suptitle('%s %s' % (creation, generation))
	fig_cos_theta.savefig('cos_theta_%s_%s.pdf' % (creation, generation))
	plt.close()
	# stretch ration on cos sq theta0
	fig_costh_lambda,axcl = plt.subplots()
	axcl.scatter(cos_theta_square_ini,stretch_ratio)
	axcl.set_xlabel(r'$cos^2\theta_i$')
	axcl.set_ylabel(r'$\lambda$')
	plt.subplots_adjust(hspace=0.5)
	fig_costh_lambda.suptitle('%s %s' % (creation, generation))
	fig_costh_lambda.savefig('costh_lambda_%s_%s.pdf' % (creation, generation))
	plt.close()
	# stretch ratio frequency
	fig_lambda_freq,axlf = plt.subplots()
	axlf = sns.distplot(stretch_ratio)
	axlf.set_xlabel(r'$\lambda$')
	plt.subplots_adjust(hspace=0.5)
	fig_lambda_freq.suptitle('%s %s' % (creation, generation))
	fig_lambda_freq.savefig('lambda_freq_%s_%s.pdf' % (creation, generation))
	plt.close()
	# stretch ration on segment length initial
	fig_lambda_length,axll = plt.subplots()
	axll.scatter(lengths_ini, stretch_ratio)
	axll.set_xlabel('Initial length')
	axll.set_ylabel(r'$\lambda$')
	plt.subplots_adjust(hspace=0.2)
	fig_lambda_length.suptitle('%s %s' % (creation, generation))
	fig_lambda_length.savefig('lambda_length_%s_%s.pdf' % (creation, generation))
	plt.close()
	# Correlation figures: # comparison: correlation theta/theta0; correlation lambdaf/L0; correlation lambdaf/theta0 all on strain
	fig_corr, axs = plt.subplots(2,2)
	axs[0,1].scatter(corr_lambda_l0,strain)
	axs[0,1].set_xlabel(r'corr($\lambda$,$\lambda_0$)')
	axs[1,0].scatter(corr_theta_0,strain)
	axs[1,0].set_xlabel(r'corr($\theta_0$,$\frac{L-L_0}{L_0}$)')
	axs[1,1].scatter(corr_lamd_theta,strain)
	axs[1,1].set_xlabel(r'corr($\lambda$,$\theta$)')
	plt.subplots_adjust(hspace=0.5)
	fig_corr.suptitle('%s %s' % (creation, generation))
	fig_corr.savefig('corr_%s_%s.pdf' % (creation, generation))
	plt.close()
	

figure_omega_xx.tight_layout()
plt.legend()
figure_omega_xx.savefig('figure_omega_xx.pdf')
plt.close()
fig_stress_strain.tight_layout()
plt.legend(loc = 'upper left')
fig_stress_strain.savefig('figure_stress_strain.pdf')
