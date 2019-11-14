from plots import *
import argparse

parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('strings', type=str)
args = parser.parse_args()
print (args.strings)


for k in range(int(args.strings[-5])):
	plot_network_geometry(k)
	plt.savefig('step_%d.png' % k,bbox_inches='tight')


