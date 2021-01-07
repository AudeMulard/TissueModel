
import sys
sys.path.append('C:/Users/am2548/TissueModel/Documents/PhD/Code/TissueModel/')

import os,fnmatch
from Plotting.information_network import load_network_info, sorted_ls, load_info_test

current_dir = os.getcwd()
os.chdir('C:/Users/am2548/TissueModel/Documents/PhD/Code/Data/validity/number_points/')


filenames = fnmatch.filter(sorted_ls('.'),'stress_data_*.rpt')
print(len(filenames))
filename_param = fnmatch.filter(sorted_ls('.'), 'parameters_??_*.csv')
print(len(filename_param))
print(filenames)
import pandas as pd
end_strain=[]
number_points=[]
for i in range(len(filename_param)):
	df_1 = pd.read_fwf(filenames[i],skiprows=2,skipfooter=4,sep='\t').to_numpy()
	stress = df_1[1]
	strain = df_1[0]
	
	network = load_network_info(int(filename_param[i][-9:-4]))
	number_points.append(len(network.vertices))
	end_strain.append(float(stress[-1])/len(network.ridge_vertices))


import matplotlib.pyplot as plt
fig = plt.figure()
plt.scatter(number_points,end_strain)
plt.xlabel('number of nodes in the network')
plt.ylabel('30% stress')
plt.show()
os.chdir(current_dir)