from network_plotting import *
import argparse
import sys
from Plotting.information_network import *


'''
parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('strings', type=str)

args = parser.parse_args()
print 'here'
print (args.strings)
'''
def sorted_ls(path):
    mtime = lambda f: os.stat(os.path.join(path, f)).st_mtime
    return list(sorted(os.listdir(path), key=mtime))


os.chdir('../Data/3D_GN/')

if len(sys.argv) != 2:
	os.chdir(sys.argv[2])
else:
	os.chdir(sorted_ls('.')[-1])

print sys.argv[1][-11:-8]
filenames=fnmatch.filter(os.listdir('.'), 'parameters_*.csv')

for k in range(int(sys.argv[1][-11:-8])):
	network = load_network_info(int(filenames[0][-7:-4]),step=k)
	plot_geometry(network)
	#plt.axis('equal')
	#plt.xlim(0.0,2.0)
	#plt.ylim([-0.1,1.1])
	#plt.xlim([0.0,2.0])
	plt.savefig('step_%03d.png' % k,bbox_inches='tight',dpi=500)
	plt.close()


for k in range(int(sys.argv[1][-11:-8])):
	network = load_network_info(int(filenames[0][-7:-4]),step=k)
	plot_constraints(network)
	#plt.xlim(0.0,2.0)
	plt.savefig('constraints_step_%03d.png' % k,bbox_inches='tight',dpi=200)
	plt.close()
