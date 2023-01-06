import numpy as np
from Network_generation.creation_network import Network
import csv, os,fnmatch


class Tensile_test:
	def __init__(self, side, space_discretization, traction_distance, element_size, path):
		self.side = side
		self.space_discretization = space_discretization
		self.traction_distance = traction_distance
		self.element_size = element_size

	# Adds the parameters of the test to the parameters file

	def save_parameters(self,network,path):		
		filenames=sorted(fnmatch.filter(os.listdir(path), 'parameters_test_*_%05d.csv' % len(network.vertices)))
		number_file = len(filenames)
		filename = 'parameters_test_%02d_%05d.csv' % (number_file,len(network.vertices))
		print(number_file)
		with open(os.path.join(path,filename), 'w') as writeFile:
			writer = csv.writer(writeFile)
			writer.writerow(["side",self.side])
			writer.writerow(["space discretization",self.space_discretization])
			writer.writerow(["traction_distance",self.traction_distance])
			writer.writerow(["element size",self.element_size])
		writeFile.close()


