import os
import csv
import numpy as np
import matplotlib.pyplot as plt
import sys

def sorted_ls(path):
    mtime = lambda f: os.stat(os.path.join(path, f)).st_mtime
    return list(sorted(os.listdir(path), key=mtime))





def plot_second_der(ax,strain, stress,color):
	der_sec = [0]
	for i in range(1,len(stress)-1):
		der_sec.append((stress[i+1]+stress[i-1]-2*stress[i])/(strain[i+1]-strain[i])**2)
		#print stress[i+1],stress[i],stress[i-1]
		#der_sec.append((stress[i]))
	ax.plot(strain[0:len(der_sec)],der_sec, color = color)

def plot_first_der(ax,strain, stress,color):
	der_sec = []
	for i in range(0,len(stress)-1):
		der_sec.append((stress[i+1]-stress[i])/(strain[i+1]-strain[i]))
	ax.plot(strain[0:len(der_sec)],der_sec, color = color)





if __name__ == "__main__":
	os.chdir('../Data/study_96/')
	if len(sys.argv) == 1:
		os.chdir(sorted_ls('.')[-1])
	else:
		os.chdir(sys.argv[1])

	fig, ax1 = plt.subplots()

	with open("stress_strain_260.csv", 'r') as readFile:
		reader = csv.reader(readFile)
		curve = np.array(list(reader))
		stress = [float(i) for i in curve[1]]
	strain=[]
	for i in range(int(len(curve[1]))):
		strain.append(i*0.0025)
	print stress

	color = 'tab:red'
	ax1.set_xlabel('strain')
	ax1.set_ylabel('stress', color=color)
	ax1.plot(strain, stress, marker='o',linestyle='dashed', markersize=5., color = color)
	ax1.tick_params(axis='y',labelcolor=color)

	ax2 =ax1.twinx()
	color = 'tab:blue'
	ax2.set_ylabel('second derivative', color=color)
	plot_second_der(ax2,strain, stress, color)
	#plot_first_der(ax2,strain, stress, 'green')
	ax2.tick_params(axis='y',labelcolor=color)


	#plt.plot(stress)
	#plt.show()
	plt.savefig('stress_strain.pdf')
	plt.show()

	os.chdir('../../../TissueModel/')
