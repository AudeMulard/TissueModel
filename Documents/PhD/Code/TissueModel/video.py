from plots import *
import argparse
import sys


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


os.chdir('../Data/default/')

if len(sys.argv) != 2:
	os.chdir(sys.argv[2])
else:
	os.chdir(sorted_ls('.')[-1])

print sys.argv[1][-7:-4]

for k in range(int(sys.argv[1][-7:-4])):
	plot_network_geometry(k)
	plt.savefig('step_%03d.png' % k,bbox_inches='tight',dpi=200)
	plt.close()


for k in range(int(sys.argv[1][-7:-4])):
	plot_network_constraints(k)
	plt.savefig('constraints_step_%03d.png' % k,bbox_inches='tight',dpi=200)
	plt.close()