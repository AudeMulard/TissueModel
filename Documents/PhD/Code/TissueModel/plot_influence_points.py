import csv
import matplotlib.pyplot as plt
import numpy as np
import os
import sys

def sorted_ls(path):
    mtime = lambda f: os.stat(os.path.join(path, f)).st_mtime
    return list(sorted(os.listdir(path), key=mtime))

os.chdir('../Data/influence_points/')
if len(sys.argv) != 1:
	os.chdir(sys.argv[1])
else:
	os.chdir(sorted_ls('.')[-1])



fig = plt.figure()
ax = fig.gca()



with open('strain_points.csv', 'r') as readFile:
	reader = csv.reader(readFile)
	curve = np.array(list(reader))
	end_stress = [float(i) for i in curve[1]]
	number_points = [float(i) for i in curve[0]]
	ax.plot(number_points, end_stress, marker='o',linestyle='none', markersize=5.)
	ax.set_xlabel('Number of nodes')
	ax.set_ylabel('Strain at 0.2\%')


#plt.show()

category = []
for k in range(len(end_stress)):
	category.append(int(number_points[k]/20.))

print category

mean=[]
st_dev=[]
set_10 =[]
k=0
for i in list(set(category)):
	values = []
	print 'i ',i
	for k in range(len(category)):
		if category[k] == i:
			values.append(end_stress[k])
	print values
	mean.append(np.mean(values))
	st_dev.append(np.std(values))
	print 'done'
	set_10.append(i*20+10)

plt.plot(set_10,mean,marker='o',linestyle='none', markersize=5.)
plt.errorbar(set_10,mean,st_dev,marker='o',linestyle='none', markersize=5.)
plt.ylabel(r'$\frac{\sigma}{k/\bar{l_0}}$ at 20\% strain', fontsize=10)
plt.xlabel('number of points in the network',fontsize=10)
plt.tight_layout()
plt.savefig('mean_number_points.pdf')
