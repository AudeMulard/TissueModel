import os
import csv
import numpy as np
import matplotlib.pyplot as plt

os.chdir('../Data/default/')
os.chdir(sorted(os.listdir('.'))[-1])

with open('stress_strain.csv', 'r') as readFile:
	reader = csv.reader(readFile)
	curve = np.array(list(reader))
	strain = curve[0]
	stress = curve[1]
	plt.plot(strain, stress, marker='o',linestyle='dashed', markersize=5.)

plt.show()

os.chdir('../../../TissueModel/')
