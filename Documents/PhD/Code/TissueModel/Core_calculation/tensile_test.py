import numpy as np
from Network_generation.creation_network import Network
import csv, os,fnmatch


class Tensile_test:
	def __init__(self, constitutive, side, space_discretization, traction_distance, element_size, _plot, video, path,details):
		self.constitutive = constitutive
		self.side = side
		self.space_discretization = space_discretization
		self.traction_distance = traction_distance
		self.iterations = abs(int(traction_distance / space_discretization))
		self.element_size = element_size
		self._plot = _plot
		self.video = video
		self.details = details

	# Adds the parameters of the test to the parameters file

	def save_parameters(self,network,path):
		filename = 'parameters_test_%05d.csv' % len(network.vertices)
		with open(os.path.join(path,filename), 'w') as writeFile:
			writer = csv.writer(writeFile)
			writer.writerow(["constitutive",self.constitutive])
			writer.writerow(["side",self.side])
			writer.writerow(["space discretization",self.space_discretization])
			writer.writerow(["traction_distance",self.traction_distance])
			writer.writerow(["element size",self.element_size])
			writer.writerow(["plot",self._plot])
			writer.writerow(["video",self.video])
			writer.writerow(["details",self.details])
		writeFile.close()


