
from Plotting.information_network import *
from Core_calculation.tensile_test import *
from scipy.interpolate import interp1d

def length_square(x):
	if len(x) ==3:
		return x[0]**2+x[1]**2+x[2]**2
	elif len(x) ==2:
		return x[0]**2+x[1]**2

def check_bio(path):
	current_path = os.getcwd()
	os.chdir(path)
	# Evolution of geometry: average + distribution of angles decreases
	filenames_ini = fnmatch.filter(os.listdir('.'), 'network_vertices_initial_*.csv')
	filenames = fnmatch.filter(os.listdir('.'), 'network_vertices_*_00.csv')
	angle_global, distr_global = [],[]
	for filename in filenames_ini:
		network = load_network_info(int(filename[-6:-4]))
		filenames = fnmatch.filter(os.listdir('.'), 'network_vertices_*_%s_*.csv' % (filename[26:28]))
		length_ini = []
		for ridge in network.ridge_vertices:
			length_ini.append(length_square(network.vertices_ini[ridge[1]]-network.vertices_ini[ridge[0]]))
		angle ,dist_angle = [],[]
		for step in range(len(filenames)):
			data_network = calculate_network_data(step,network.ridge_vertices,length_ini)
			omega_xx = data_network[3]
			distr = std_dev(data_network[2])
			angle.append(omega_xx)
			dist_angle.append(distr)
		angle_global.append(angle)
		distr_global.append(dist_angle)
	with open('testing.txt', 'a') as writeFile:
		writeFile.write('Evolution of the geometry :'+ ' Orientation :' + str(angle_global)+ ' Distribution (stdev):'+str( distr_global))
	
	# Compare stress_strain_curve to experimental
	stress_exp, strain_exp = get_exp_data()
	exp_curve = interp1d(strain_exp,stress_exp,fill_value='extrapolate')
	filenames = fnmatch.filter(os.listdir('.'), 'stress_strain_*_0*.csv')
	print filenames
	stress,strain,diff_global = [],[],[]
	diff = []
	for filename in filenames:
		with open(filename,'r') as readFile:
			reader = csv.reader(readFile)
			curve = np.array(list(reader))
			stress = [float(i) for i in curve[1]]
			strain = [float(i) for i in curve[0]]
		calc_curve = interp1d(strain, stress,fill_value='extrapolate')
		for i in range(1,10):
			diff.append(exp_curve(i/10.)-calc_curve(i/10.)) 
		diff_global.append(diff)		
	with open('testing.txt', 'a') as writeFile:
		writeFile.write('Difference in global stress compared to experiment:'+ str(diff_global))
	
	os.chdir(current_path)


def get_exp_data():
	strain_exp=[1.388,2.775,4.097,5.286,6.476,7.401,8.722,9.912,11.233,12.291,13.348,14.273,15.066,16.123,17.181,18.238,19.295,20.22,21.542,22.731,23.656,24.78]
	stress_exp=[2.798,2.752,2.752,2.798,2.798,2.89,3.028,3.257,3.486,3.716,3.991,4.266,4.495,4.862,5.183,5.505,5.826,6.101,6.514,6.927,7.248,7.523]

	offset_x=strain_exp[0]
	offset_y=stress_exp[0]
	
	for i in range(len(strain_exp)):
		strain_exp[i] = strain_exp[i]-offset_x
		stress_exp[i] = stress_exp[i]-offset_y
		stress_exp[i] = 2*stress_exp[i]/(4.862-offset_y)
		strain_exp[i] = 0.1*strain_exp[i]/(7.401-offset_x)
	return stress_exp, strain_exp

	

