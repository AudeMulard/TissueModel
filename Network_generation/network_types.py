import numpy as np
from scipy.spatial import Voronoi, voronoi_plot_2d
import matplotlib.pyplot as plt
from sympy import Symbol
from sympy.solvers import solve
import time
import cProfile, pstats
import scipy
#profile = cProfile.Profile()

def length_square(x):
	return x[0]**2+x[1]**2


from scipy.special import gammainc
###Taken from here: https://github.com/diregoblin/poisson_disc_sampling/blob/main/poisson_disc.py
# Uniform sampling in a hyperspere
# Based on Matlab implementation by Roger Stafford
# Can be optimized for Bridson algorithm by excluding all points within the r/2 sphere
def hypersphere_volume_sample(center,radius,k=1):
    ndim = center.size
    x = np.random.normal(size=(k, ndim))
    ssq = np.sum(x**2,axis=1)
    fr = radius*gammainc(ndim/2,ssq/2)**(1/ndim)/np.sqrt(ssq)
    frtiled = np.tile(fr.reshape(k,1),(1,ndim))
    p = center + np.multiply(x,frtiled)
    return p


# Uniform sampling on the sphere's surface
def hypersphere_surface_sample(center,radius,k=1):
    ndim = center.size
    vec = np.random.standard_normal(size=(k, ndim))
    #vec /= np.linalg.norm(vec, axis=1)[:,None]
    vec /= np.linalg.norm(vec, axis=1)[:,None]
    p = center + np.multiply(vec, radius)
    return p


def squared_distance(p0, p1):
    return np.sum(np.square(p0-p1))

def Bridson_sampling(network, dims=np.array([1.0,1.0]), radius=0.05, k=30, hypersphere_sample=hypersphere_volume_sample):
    # References: Fast Poisson Disk Sampling in Arbitrary Dimensions
    #             Robert Bridson, SIGGRAPH, 2007

    ndim=dims.size

    # size of the sphere from which the samples are drawn relative to the size of a disc (radius)
    sample_factor = 2
    if hypersphere_sample == hypersphere_volume_sample:
        sample_factor = 2
        
    # for the surface sampler, all new points are almost exactly 1 radius away from at least one existing sample
    # eps to avoid rejection
    if hypersphere_sample == hypersphere_surface_sample:
        eps = 0.001
        sample_factor = 1 + eps
    
    def in_limits(p):
        #return np.all(np.zeros(ndim) <= p) and np.all(p < dims)
        if len(dims) == 2:
            if network.creation =='Voronoi': return np.all(np.array((-0.1,-0.1)) <= p) and np.all(p < dims)
            if network.creation =='growth_network': return np.all(np.array((0.,0.)) <= p) and np.all(p < dims)
        if len(dims) == 3:
            if network.creation =='Voronoi': return np.all(np.array((-0.1,-0.1,-0.1)) <= p) and np.all(p < dims)
            if network.creation =='growth_network': return np.all(np.array((0.,0.,0.)) <= p) and np.all(p < dims)
        

    # Check if there are samples closer than "squared_radius" to the candidate "p"
    def in_neighborhood(p, n=2):
        indices = (p / cellsize).astype(int)
        indmin = np.maximum(indices - n, np.zeros(ndim, dtype=int))
        indmax = np.minimum(indices + n + 1, gridsize)
        
        # Check if the center cell is empty
        if not np.isnan(P[tuple(indices)][0]):
            return True
        a = []
        for i in range(ndim):
            a.append(slice(indmin[i], indmax[i]))
        with np.errstate(invalid='ignore'):
            if np.any(np.sum(np.square(p - P[tuple(a)]), axis=ndim) < squared_radius):
                return True

    def add_point(p):
        points.append(p)
        indices = (p/cellsize).astype(int)
        P[tuple(indices)] = p

    cellsize = radius/np.sqrt(ndim)
    gridsize = (np.ceil(dims/cellsize)).astype(int)

    # Squared radius because we'll compare squared distance
    squared_radius = radius*radius

    # Positions of cells
    P = np.empty(np.append(gridsize, ndim), dtype=np.float32) #n-dim value for each grid cell
    # Initialise empty cells with NaNs
    P.fill(np.nan)

    points = []
    initial = np.random.uniform(np.zeros(ndim), dims)
    add_point(initial)
    while len(points):
        i = np.random.randint(len(points))
        p = points[i]
        del points[i]
        Q = hypersphere_sample(np.array(p), radius * sample_factor, k)
        for q in Q:
            if in_limits(q) and not in_neighborhood(q):
                add_point(q)
    return P[~np.isnan(P).any(axis=ndim)]


def Seeds_generation(network):
	Seeds = []
	if network.generation == 'random':
		for node in range(network.complexity):
			if network.dimension==2:
				Seeds.append([np.random.rand()*network.length[0],np.random.rand()*network.length[1]])
			elif network.dimension ==3:
				Seeds.append([np.random.rand()*network.length[0],np.random.rand()*network.length[1],np.random.rand()*network.length[2]])
	elif network.generation == 'Bridson_sampling':
		if network.dimension==2:
			if network.creation == 'Voronoi': dims2d = np.array([network.length[0]+0.2,network.length[1]+0.2])
			elif network.creation == 'growth_network': dims2d = np.array([network.length[0],network.length[1]])
			Seeds = Bridson_sampling(network,dims=dims2d, radius=1./np.sqrt(network.complexity), k=30, hypersphere_sample=hypersphere_surface_sample)
		elif network.dimension==3:
			dims3d = np.array([network.length[0],network.length[1],network.length[2]])
			Seeds = Bridson_sampling(network,dims3d, radius=1./network.complexity**(1./3.), k=30)
	elif network.generation == 'grid':
		number_points = int(np.sqrt(network.complexity))
		for x in range(-2,number_points+2):
			for y in range(-2,number_points+2):
				if network.dimension == 2:
					#Seeds.append([x/float(number_points)+np.random.rand()*0.02,y/float(number_points)+np.random.rand()*0.02])
					Seeds.append([x/float(number_points),y/float(number_points)])
				if network.dimension == 3:
					for z in range(1,number_points):
						Seeds.append([x/float(number_points)+np.random.rand()*network.disturbance,y/float(number_points)+np.random.rand()*network.disturbance,z/float(number_points)+np.random.rand()*network.disturbance])
	elif network.generation == 'regular':
		for x in range(int(np.sqrt(network.complexity/2.))):
			for y in range(-2,int(np.sqrt(network.complexity))+2):
				Seeds.append([x/np.sqrt(network.complexity)*2,y/np.sqrt(network.complexity)+1./np.sqrt(network.complexity)])
		for x in range(int(np.sqrt(network.complexity/2.))+2):
			for y in range(-2,int(np.sqrt(network.complexity))+1):
				Seeds.append([x/np.sqrt(network.complexity)*2-1./np.sqrt(network.complexity),y/np.sqrt(network.complexity)-1./(2*np.sqrt(network.complexity))])
	return Seeds

def select_network(network):
		creation = network.creation
		
		if creation == "Voronoi":
			Seeds = Seeds_generation(network)
			voronoi = Voronoi(Seeds)
			vertices = Voronoi(Seeds).vertices
			ridge_vertices = Voronoi(Seeds).ridge_vertices
			if network.dimension == 3:
				new_ridge_vertices = []
				for ridge in ridge_vertices:
					for i in range(len(ridge)-1):
						new_ridge_vertices.append([ridge[i],ridge[i+1]])
				ridge_vertices = new_ridge_vertices
			return vertices, ridge_vertices

		if creation == "growth_network":
			div_space = network.complexity
			#div_space = network.complexity**(3/2)
			Seeds = Seeds_generation(network)
			#for l in range(len(Seeds)):
			#	plt.annotate(l, (Seeds[l][0], Seeds[l][1]))
			Seeds = np.array(Seeds)
			#print(Seeds)
			len_seeds = len(Seeds)
			from numpy import array
			start = time.time()
			vertices = np.zeros((len_seeds*2,network.dimension))
			lines = []
			if network.dimension == 2:
				import random
				for i in range(len_seeds):
					a_coeff = random.uniform(-1,1)*np.pi/2. # angle
					b_coeff = float(Seeds[i][1]-np.tan(a_coeff)*Seeds[i][0])
					lines.append([np.cos(a_coeff), np.sin(a_coeff),b_coeff,Seeds[i][0],Seeds[i][1],Seeds[i][0],Seeds[i][1],'empty','empty'])
				#print(lines)
				#lines = [[0.39869299504387395, 0.9170844539642712, -0.09151351750767646, 0.04250977, 0.0062686056, 0.04250977, 0.0062686056, 'empty', 'empty'], [0.5127428530023312, -0.8585422334952602, 0.13634455215889688, 0.0050142473, 0.12794864, 0.0050142473, 0.12794864, 'empty', 'empty'], [0.8295804960487898, 0.5583871422010394, 0.17303193881642806, 0.030896451, 0.19382821, 0.030896451, 0.19382821, 'empty', 'empty'], [0.3841398404195805, 0.9232749227626726, 0.3133182593193253, 0.0020711431, 0.31829622, 0.0020711431, 0.31829622, 'empty', 'empty'], [0.88928013031534, -0.4573629300963643, 0.4221227828997279, 0.031482268, 0.40593123, 0.031482268, 0.40593123, 'empty', 'empty'], [0.999728663860634, 0.02329374713161411, 0.4834512912869571, 0.01352423, 0.4837664, 0.01352423, 0.4837664, 'empty', 'empty'], [0.9818928897257755, -0.18943693701590025, 0.5688045043637219, 0.031602185, 0.5627075, 0.031602185, 0.5627075, 'empty', 'empty'], [0.5925441364322345, -0.8055379856839635, 0.6835261923844239, 0.0030708476, 0.6793515, 0.0030708476, 0.6793515, 'empty', 'empty'], [0.7129442880635516, 0.7012206800412802, 0.7644368458013648, 0.04164039, 0.8053925, 0.04164039, 0.8053925, 'empty', 'empty'], [0.994358812886389, -0.10606861569366871, 0.8830201954030386, 0.008798207, 0.8820817, 0.008798207, 0.8820817, 'empty', 'empty'], [0.30582295332027354, 0.9520883998991195, 0.844510873099726, 0.035346977, 0.9545531, 0.035346977, 0.9545531, 'empty', 'empty'], [0.04171038617067018, 0.9991297431692711, -1.6044826985057585, 0.06996708, 0.0715074, 0.06996708, 0.0715074, 'empty', 'empty'], [0.9857897509314467, -0.16798382945544577, 0.2338290293816117, 0.09771564, 0.21717776, 0.09771564, 0.21717776, 'empty', 'empty'], [0.8848203502233192, -0.4659323425463003, 0.30782585748988667, 0.059276253, 0.27661192, 0.059276253, 0.27661192, 'empty', 'empty'], [0.9912443326941571, 0.13204042146901476, 0.46411792402321206, 0.08379503, 0.47528, 0.08379503, 0.47528, 'empty', 'empty'], [0.9990621131562645, 0.043300046832991934, 0.6276061322602505, 0.053746954, 0.62993556, 0.053746954, 0.62993556, 'empty', 'empty'], [0.7077940072228208, 0.7064188865959499, 0.6782317723019244, 0.058540177, 0.7366582, 0.058540177, 0.7366582, 'empty', 'empty'], [0.8161836205203942, -0.5777926077704881, 0.9529745495732671, 0.0777878, 0.897907, 0.0777878, 0.897907, 'empty', 'empty'], [0.9946785728959251, 0.1030268732987949, -0.009060485457668481, 0.11319934, 0.002664482, 0.11319934, 0.002664482, 'empty', 'empty'], [0.7784644582823433, 0.6276886865247594, 0.04377796163373475, 0.108236715, 0.131051, 0.108236715, 0.131051, 'empty', 'empty'], [0.007869632450477513, 0.9999690339630994, -15.012039015685138, 0.120699905, 0.32491198, 0.120699905, 0.32491198, 'empty', 'empty'], [0.8853953276645852, 0.4648388040490186, 0.35326326384499523, 0.10225634, 0.40694857, 0.10225634, 0.40694857, 'empty', 'empty'], [0.8400543849492624, -0.5425021938458097, 0.6130281938947248, 0.100821316, 0.5479184, 0.100821316, 0.5479184, 'empty', 'empty'], [0.8029206876138728, -0.5960858741839684, 0.7159208649974287, 0.13272116, 0.6173891, 0.13272116, 0.6173891, 'empty', 'empty'], [0.10768426474112848, -0.994185143284269, 1.633837411922759, 0.10313113, 0.68168867, 0.10313113, 0.68168867, 'empty', 'empty'], [0.9448529437715762, 0.3274949078171852, 0.687557284804819, 0.1458637, 0.738115, 0.1458637, 0.738115, 'empty', 'empty'], [0.6814568230347414, -0.7318583184875315, 0.9687879199266654, 0.14579044, 0.8122146, 0.14579044, 0.8122146, 'empty', 'empty'], [0.6246267885200426, 0.7809234117780937, 0.7039832271454778, 0.14801021, 0.88902915, 0.14801021, 0.88902915, 'empty', 'empty'], [0.9325062032387652, -0.3611539573661386, 1.0013739923426674, 0.12107148, 0.95448375, 0.12107148, 0.95448375, 'empty', 'empty'], [0.6801721431987183, 0.7330524235117583, -0.1279695141386168, 0.16739658, 0.05244138, 0.16739658, 0.05244138, 'empty', 'empty'], [0.735059506872314, -0.6780025968656983, 0.2986654065482337, 0.19748582, 0.116508864, 0.19748582, 0.116508864, 'empty', 'empty'], [0.4332742349048367, -0.9012621357671853, 0.510576002210073, 0.15849093, 0.18089595, 0.15849093, 0.18089595, 'empty', 'empty'], [0.3769398697576807, 0.9262377311398314, -0.20653313837343684, 0.18571052, 0.24980514, 0.18571052, 0.24980514, 'empty', 'empty'], [0.5731794175908976, 0.8194298964830118, 0.04682606891902974, 0.19133495, 0.32036272, 0.19133495, 0.32036272, 'empty', 'empty'], [0.818432496930947, 0.5746026870519972, 0.26832274042791066, 0.17045248, 0.3879935, 0.17045248, 0.3879935, 'empty', 'empty'], [0.9887801202078574, 0.14937829119967555, 0.4763126694625303, 0.15228026, 0.49931815, 0.15228026, 0.49931815, 'empty', 'empty'], [0.7511626007023852, 0.6601172223976808, 0.41070921166668506, 0.1930376, 0.5803495, 0.1930376, 0.5803495, 'empty', 'empty'], [0.8618198064169538, 0.5072145712294199, 0.5347357646595894, 0.1975274, 0.65098834, 0.1975274, 0.65098834, 'empty', 'empty'], [0.9368900284946262, 0.34962419039211023, 0.8801358984284241, 0.19283587, 0.9520975, 0.19283587, 0.9520975, 'empty', 'empty'], [0.8794078989737983, 0.47606905719915216, 0.06965820381903638, 0.22819601, 0.19319253, 0.22819601, 0.19319253, 'empty', 'empty'], [0.6180749181383581, -0.7861191993382822, 0.6867649135918656, 0.24083133, 0.38045555, 0.24083133, 0.38045555, 'empty', 'empty'], [0.988007067741044, 0.15440865938717319, 0.4187305196780725, 0.20357631, 0.45054603, 0.20357631, 0.45054603, 'empty', 'empty'], [0.9880667562087924, 0.15402624865598338, 0.48246841002891816, 0.22596489, 0.5176933, 0.22596489, 0.5176933, 'empty', 'empty'], [0.4938199136570473, -0.8695641970986078, 1.1628119766656666, 0.23914851, 0.74169695, 0.23914851, 0.74169695, 'empty', 'empty'], [0.5514317853741946, -0.8342199866216513, 1.1362598671125215, 0.21648712, 0.80875266, 0.21648712, 0.80875266, 'empty', 'empty'], [0.7800937739471314, -0.6256626118355817, 1.0544778375687978, 0.2181492, 0.8795145, 0.2181492, 0.8795145, 'empty', 'empty'], [0.9684403959901694, -0.2492452595625518, 0.13115055523498725, 0.26596582, 0.06269955, 0.26596582, 0.06269955, 'empty', 'empty'], [0.15925863165695878, -0.9872368957057638, 1.783656860121931, 0.2662024, 0.13348055, 0.2662024, 0.13348055, 'empty', 'empty'], [0.9998874656562229, -0.01500186733628896, 0.2005855780442063, 0.29920626, 0.19609642, 0.29920626, 0.19609642, 'empty', 'empty'], [0.7461739516098228, -0.6657510299946834, 0.527728690445187, 0.29269782, 0.26657796, 0.29269782, 0.26657796, 'empty', 'empty'], [0.8312752008467915, -0.555861080178426, 0.5364785892748857, 0.29808444, 0.33715406, 0.29808444, 0.33715406, 'empty', 'empty'], [0.10872689944466775, 0.9940716580494331, -2.0511518472355186, 0.27417693, 0.45560142, 0.27417693, 0.45560142, 'empty', 'empty'], [0.642274482684399, 0.7664747151031713, 0.2697922718423478, 0.2636982, 0.58448327, 0.2636982, 0.58448327, 'empty', 'empty'], [0.04262058084911855, -0.9990913302036425, 7.62313302263164, 0.2953908, 0.6987224, 0.2953908, 0.6987224, 'empty', 'empty'], [0.9534215677794248, -0.3016410351610398, 0.8871953388536127, 0.29131606, 0.7950295, 0.29131606, 0.7950295, 'empty', 'empty'], [0.46589456388028927, -0.8848402428386692, 1.43468276228067, 0.26200485, 0.9370757, 0.26200485, 0.9370757, 'empty', 'empty'], [0.28487525060077923, 0.9585645995941762, -0.00380124854084829, 0.29777104, 0.99815583, 0.29777104, 0.99815583, 'empty', 'empty'], [0.5236274921909274, 0.8519473278447678, -0.4991127792897466, 0.34077555, 0.055332553, 0.34077555, 0.055332553, 'empty', 'empty'], [0.9476585409498481, -0.319285592792417, 0.23939269375508024, 0.3365863, 0.12598987, 0.3365863, 0.12598987, 'empty', 'empty'], [0.6830545933692942, 0.7303673202417453, 0.19018893952506882, 0.30691114, 0.5183587, 0.30691114, 0.5183587, 'empty', 'empty'], [0.2645944734274348, -0.9643597692934202, 1.7859601766851798, 0.31681815, 0.6312623, 0.31681815, 0.6312623, 'empty', 'empty'], [0.7122477172623539, -0.7019281938023048, 1.0892108072939921, 0.34514424, 0.74906725, 0.34514424, 0.74906725, 'empty', 'empty'], [0.4482621213156914, -0.8939021594075922, 1.5135129874321465, 0.3134453, 0.88845587, 0.3134453, 0.88845587, 'empty', 'empty'], [0.9584730359421302, -0.2851831681075796, 0.31557660130663384, 0.36933517, 0.20568496, 0.36933517, 0.20568496, 'empty', 'empty'], [0.7905222068803828, -0.6124333763185751, 0.6239093147429304, 0.37471294, 0.3336117, 0.37471294, 0.3336117, 'empty', 'empty'], [0.5837695788503486, 0.8119193795007523, -0.10201417001854585, 0.36345917, 0.40349272, 0.36345917, 0.40349272, 'empty', 'empty'], [0.8298872565925318, -0.5579311259871789, 0.7837947664396301, 0.3948579, 0.5183328, 0.3948579, 0.5183328, 'empty', 'empty'], [0.7549409718752149, 0.6557927484991779, 0.522613619017402, 0.35632128, 0.83213836, 0.35632128, 0.83213836, 'empty', 'empty'], [0.9457382318979187, 0.3249295257845898, 0.8285858198078059, 0.35024357, 0.94891983, 0.35024357, 0.94891983, 'empty', 'empty'], [0.19981460210475438, -0.9798337230294325, 1.997072930415991, 0.40289313, 0.021400144, 0.40289313, 0.021400144, 'empty', 'empty'], [0.39609418988645595, -0.9182098849055113, 1.1172659755964305, 0.4441531, 0.087647825, 0.4441531, 0.087647825, 'empty', 'empty'], [0.9979280614130999, -0.06433960090249541, 0.17954745181424156, 0.4162769, 0.15270875, 0.4162769, 0.15270875, 'empty', 'empty'], [0.9928222157374172, 0.11959953151348597, 0.2199310770151594, 0.40218785, 0.2683803, 0.40218785, 0.2683803, 'empty', 'empty'], [0.7944278612337583, -0.6073585212175395, 0.7707382660469091, 0.41775292, 0.45135647, 0.41775292, 0.45135647, 'empty', 'empty'], [0.6629618295920655, 0.7486531990874954, 0.1746208906765862, 0.41017506, 0.63781327, 0.41017506, 0.63781327, 'empty', 'empty'], [0.9983149007328482, -0.05802894945424543, 0.7309656800513096, 0.4290739, 0.70602494, 0.4290739, 0.70602494, 'empty', 'empty'], [0.8766637041320894, -0.4811036789065372, 1.0064985938114248, 0.40726343, 0.7829968, 0.40726343, 0.7829968, 'empty', 'empty'], [0.8370765953970495, 0.5470857094080271, 0.6194156154855146, 0.4233213, 0.8960845, 0.4233213, 0.8960845, 'empty', 'empty'], [0.9091720095189071, 0.4164207693035403, 0.79286665078642, 0.41284177, 0.98195726, 0.41284177, 0.98195726, 'empty', 'empty'], [0.9831342329414052, 0.18288542866700688, -0.06495714944660357, 0.4736528, 0.023153095, 0.4736528, 0.023153095, 'empty', 'empty'], [0.9974972870908707, 0.0707047540576523, 0.12940469494630585, 0.48617345, 0.16386572, 0.48617345, 0.16386572, 'empty', 'empty'], [0.8943702968174194, 0.4473273657744639, -0.0005154218552404777, 0.4620243, 0.2305702, 0.4620243, 0.2305702, 'empty', 'empty'], [0.3502562801091451, 0.9366539052628265, -0.8846226826104095, 0.450502, 0.32010794, 0.450502, 0.32010794, 'empty', 'empty'], [0.5885243580368864, -0.8084794864418458, 1.0149717545687051, 0.45437145, 0.3907835, 0.45437145, 0.3907835, 'empty', 'empty'], [0.8180125741521529, -0.5752003377337053, 0.8772980895112582, 0.49529076, 0.52902544, 0.49529076, 0.52902544, 'empty', 'empty'], [0.44740397088593326, -0.8943319779788147, 1.5223489826172605, 0.46493763, 0.59296834, 0.46493763, 0.59296834, 'empty', 'empty'], [0.8051511346525885, -0.5930696842426271, 1.1700576939637333, 0.4566902, 0.83366233, 0.4566902, 0.83366233, 'empty', 'empty'], [0.7059755997954859, 0.7082361558783934, -0.42104469713331527, 0.5224724, 0.10310066, 0.5224724, 0.10310066, 'empty', 'empty'], [0.6777750007172855, 0.7352693713209355, -0.39853295952542334, 0.5488096, 0.1968311, 0.5488096, 0.1968311, 'empty', 'empty'], [0.9983459910929746, 0.05749158258898731, 0.36389017986131883, 0.5250738, 0.39412752, 0.5250738, 0.39412752, 'empty', 'empty'], [0.0724438456138518, -0.9973724927190826, 7.6966413501489015, 0.5252747, 0.4649086, 0.5252747, 0.4649086, 'empty', 'empty'], [0.9999941755840509, 0.0034130335442782806, 0.6845202332539814, 0.5095775, 0.68625945, 0.5095775, 0.68625945, 'empty', 'empty'], [0.848192851495189, 0.5296875368294597, 0.4416752658506675, 0.5047314, 0.75687474, 0.5047314, 0.75687474, 'empty', 'empty'], [0.6683967002794835, 0.743804981870583, 0.2793287355371863, 0.5222066, 0.8604505, 0.5222066, 0.8604505, 'empty', 'empty'], [0.8198783258459488, -0.5725377985845513, 1.2832924954693508, 0.5066362, 0.9294981, 0.5066362, 0.9294981, 'empty', 'empty'], [0.9695077596850918, -0.24506061272753352, 1.1269973073777666, 0.54432946, 0.9894082, 0.54432946, 0.9894082, 'empty', 'empty'], [0.8750363986349904, -0.4840571258270103, 0.3336101202752574, 0.5613049, 0.02310458, 0.5613049, 0.02310458, 'empty', 'empty'], [0.16809764879681718, -0.9857703487471013, 3.555482596295195, 0.5914339, 0.08715339, 0.5914339, 0.08715339, 'empty', 'empty'], [0.14886118070065282, 0.9888581035115238, -3.464633756423785, 0.5616693, 0.26643452, 0.5616693, 0.26643452, 'empty', 'empty'], [0.7704213303838763, -0.6375350764385739, 0.8749701490995446, 0.5949208, 0.38266435, 0.5949208, 0.38266435, 'empty', 'empty'], [0.4512763031802965, 0.8923842771967271, -0.7307600261864322, 0.59879595, 0.45333958, 0.59879595, 0.45333958, 'empty', 'empty'], [0.9891100091533679, -0.14717808869741617, 0.6469115318414383, 0.5567925, 0.56406164, 0.5567925, 0.56406164, 'empty', 'empty'], [0.7694099209579369, -0.6387553315092573, 1.0982466754989355, 0.5582078, 0.6348289, 0.5582078, 0.6348289, 'empty', 'empty'], [0.8348363427622001, -0.5504982114443555, 1.099505110321353, 0.5694581, 0.72399956, 0.5694581, 0.72399956, 'empty', 'empty'], [0.9988532255795384, -0.04787727800692036, 0.9026940926944177, 0.59161246, 0.8743368, 0.59161246, 0.8743368, 'empty', 'empty'], [0.8604141257648064, -0.5095954593443548, 0.4048536293311056, 0.6330651, 0.029909642, 0.6330651, 0.029909642, 'empty', 'empty'], [0.4148951675162303, 0.9098692213563877, -1.1773142300397514, 0.61316526, 0.16736327, 0.61316526, 0.16736327, 'empty', 'empty'], [0.9431715645410511, 0.33230618387442906, 0.06908744293670446, 0.6282557, 0.2904398, 0.6282557, 0.2904398, 'empty', 'empty'], [0.40927440703657975, -0.9124113434985648, 1.9041421713588453, 0.6205662, 0.52068985, 0.6205662, 0.52068985, 'empty', 'empty'], [0.9724821948706125, -0.232977210601458, 0.8218169697962232, 0.62136954, 0.6729557, 0.62136954, 0.6729557, 'empty', 'empty'], [0.32225358778462104, -0.9466533817400852, 2.581770307439581, 0.60773927, 0.7964731, 0.60773927, 0.7964731, 'empty', 'empty'], [0.9921972042903633, -0.12467841753241427, 0.23013986450762275, 0.6802087, 0.14466558, 0.6802087, 0.14466558, 'empty', 'empty'], [0.9901339857030135, 0.14012383935599543, 0.12326915876533256, 0.66346455, 0.21716271, 0.66346455, 0.21716271, 'empty', 'empty'], [0.4765848211263355, -0.8791284936071512, 1.5680651435405908, 0.69826525, 0.28001556, 0.69826525, 0.28001556, 'empty', 'empty'], [0.21603357579257523, -0.9763859350330041, 3.381804852664206, 0.6652843, 0.37498423, 0.6652843, 0.37498423, 'empty', 'empty'], [0.17073233045282413, -0.9853174469875927, 4.323608154378896, 0.6686832, 0.46455467, 0.6686832, 0.46455467, 'empty', 'empty'], [0.7392285502746551, 0.6734546387536667, -0.04073630153011576, 0.67089224, 0.57046235, 0.67089224, 0.57046235, 'empty', 'empty'], [0.14442653549650952, -0.9895155258228521, 5.326464663285492, 0.68401915, 0.64001584, 0.68401915, 0.64001584, 'empty', 'empty'], [0.5908994906681242, -0.8067451840129889, 1.6444343203379108, 0.68384147, 0.710797, 0.68384147, 0.710797, 'empty', 'empty'], [0.9825616017437349, 0.18593735175802134, 0.6531460225763568, 0.6768605, 0.7812333, 0.6768605, 0.7812333, 'empty', 'empty'], [0.8582918438092283, 0.5131618758740317, 0.44419988337393124, 0.68180394, 0.85184187, 0.68180394, 0.85184187, 'empty', 'empty'], [0.6550229612017018, 0.7556089731458685, 0.16112743810089714, 0.65585357, 0.91769457, 0.65585357, 0.91769457, 'empty', 'empty'], [0.022641868794425195, -0.9997436400285306, 30.54042779263141, 0.66931224, 0.98718464, 0.66931224, 0.98718464, 'empty', 'empty'], [0.8167088011001459, 0.5770500274721095, -0.4574532767910997, 0.7032072, 0.039401576, 0.7032072, 0.039401576, 'empty', 'empty'], [0.3865216133240924, 0.9222803491527621, -1.5330401294364302, 0.73431015, 0.21909946, 0.73431015, 0.21909946, 'empty', 'empty'], [0.6343550298437561, 0.7730418462877202, -0.5059080256450517, 0.74078524, 0.3968325, 0.74078524, 0.3968325, 'empty', 'empty'], [0.23375185667550646, -0.9722962868903458, 3.637902404552033, 0.74051553, 0.5577108, 0.74051553, 0.5577108, 'empty', 'empty'], [0.7691838818131423, -0.6390275079829238, 1.4142544078368462, 0.74646187, 0.79410404, 0.74646187, 0.79410404, 'empty', 'empty'], [0.98203351964351, -0.18870656135010164, 1.0289012648549924, 0.74381375, 0.8859708, 0.74381375, 0.8859708, 'empty', 'empty'], [0.6995137141269903, 0.7146191739299074, 0.20174354466640898, 0.738881, 0.9565801, 0.738881, 0.9565801, 'empty', 'empty'], [0.013801310140240764, -0.9999047573836285, 55.685733036107884, 0.7676617, 0.06865257, 0.7676617, 0.06865257, 'empty', 'empty'], [0.8312641845746424, 0.5558775543619789, -0.1869871395339493, 0.75713253, 0.31931758, 0.75713253, 0.31931758, 'empty', 'empty'], [0.9437806022745293, 0.330572495483709, 0.2226903952487605, 0.75952536, 0.4887249, 0.75952536, 0.4887249, 'empty', 'empty'], [0.994144576598458, 0.10805813629649809, 0.5527929846122931, 0.75460917, 0.6348149, 0.75460917, 0.6348149, 'empty', 'empty'], [0.037665565203539746, -0.9992904008334603, 20.730611563659778, 0.7546194, 0.71009845, 0.7546194, 0.71009845, 'empty', 'empty'], [0.2667142588841023, -0.9637756502982964, 3.877398944868939, 0.79720676, 0.99668133, 0.79720676, 0.99668133, 'empty', 'empty'], [0.21248242095608152, 0.9771648892508585, -3.7902624785610497, 0.8334159, 0.04245355, 0.8334159, 0.04245355, 'empty', 'empty'], [0.9413178343901911, 0.33752145807186934, -0.17202129810301003, 0.81583863, 0.120508, 0.81583863, 0.120508, 'empty', 'empty'], [0.7912534545036894, 0.6114883242842646, -0.4169712454265776, 0.8034497, 0.20394245, 0.8034497, 0.20394245, 'empty', 'empty'], [0.011259269593238684, -0.9999366124151204, 73.6383250829539, 0.82576567, 0.30201146, 0.82576567, 0.30201146, 'empty', 'empty'], [0.9989994713905291, -0.04472198744961249, 0.40627783253383787, 0.8063569, 0.37017983, 0.8063569, 0.37017983, 'empty', 'empty'], [0.969811215888111, 0.24385693661986244, 0.2366967954072132, 0.8115843, 0.4407679, 0.8115843, 0.4407679, 'empty', 'empty'], [0.9697948796561572, -0.24392189609114567, 0.7175094446910881, 0.8356952, 0.5073162, 0.8356952, 0.5073162, 'empty', 'empty'], [0.6233008792024828, -0.7819821059240498, 1.6218581851247627, 0.83203435, 0.5780028, 0.83203435, 0.5780028, 'empty', 'empty'], [0.8649026623237818, 0.5019396225695221, 0.19676017701434545, 0.81588155, 0.6702508, 0.81588155, 0.6702508, 'empty', 'empty'], [0.6254365141200338, 0.7802750584283602, -0.3267922316008838, 0.84923184, 0.7326829, 0.84923184, 0.7326829, 'empty', 'empty'], [0.5531543612295944, 0.8330787793796451, -0.4256134294500049, 0.8192801, 0.8082643, 0.8192801, 0.8082643, 'empty', 'empty'], [0.9885750874780084, 0.1507292155418054, 0.802274552468264, 0.8030488, 0.92471635, 0.8030488, 0.92471635, 'empty', 'empty'], [0.17444016711336635, 0.9846677754945883, -4.9028845649141815, 0.8847343, 0.09120219, 0.8847343, 0.09120219, 'empty', 'empty'], [0.6252198848695095, -0.7804486501774203, 1.27055561307086, 0.8788313, 0.17352921, 0.8788313, 0.17352921, 'empty', 'empty'], [0.30123855641608205, -0.9535488095154621, 2.979921174972065, 0.8646675, 0.24287897, 0.8646675, 0.24287897, 'empty', 'empty'], [0.9123412556444567, -0.40943062080051623, 0.7596267972157519, 0.8770176, 0.36604834, 0.8770176, 0.36604834, 'empty', 'empty'], [0.5488270304036318, -0.835935936958288, 1.9699956888944008, 0.879325, 0.6306678, 0.879325, 0.6306678, 'empty', 'empty'], [0.20983289682505696, 0.9777372629750822, -3.3464857622931232, 0.888125, 0.79182106, 0.888125, 0.79182106, 'empty', 'empty'], [0.017597639109900098, 0.9998451495595494, -47.98982417539897, 0.8601554, 0.8816297, 0.8601554, 0.8816297, 'empty', 'empty'], [0.6775432772585215, -0.7354829076476097, 1.9051257923539557, 0.88038075, 0.9494599, 0.88038075, 0.9494599, 'empty', 'empty'], [0.923149285857491, -0.3844416679039934, 0.3959374407062577, 0.9008042, 0.020801282, 0.9008042, 0.020801282, 'empty', 'empty'], [0.8190752352454078, 0.5736861154043036, -0.428954042773664, 0.93303454, 0.22454998, 0.93303454, 0.22454998, 'empty', 'empty'], [0.9954591069822868, 0.09519015876669285, 0.20975109501892486, 0.91050065, 0.29681715, 0.91050065, 0.29681715, 'empty', 'empty'], [0.9473174147209321, 0.32029629371325785, 0.04983340367644651, 0.94767416, 0.3702503, 0.94767416, 0.3702503, 'empty', 'empty'], [0.942757206999805, 0.3334799074156143, 0.1092481658783358, 0.90923375, 0.43086988, 0.90923375, 0.43086988, 'empty', 'empty'], [0.9786881866931699, -0.20535197400374547, 0.6917396688876147, 0.90624446, 0.5015881, 0.90624446, 0.5015881, 'empty', 'empty'], [0.709218835765964, 0.7049883991916254, -0.35712435040501866, 0.93069, 0.56801414, 0.93069, 0.56801414, 'empty', 'empty'], [0.5738349581611237, -0.8189709645599296, 2.0235546419706383, 0.9261505, 0.7017628, 0.9261505, 0.7017628, 'empty', 'empty'], [0.8910167149812729, -0.453970498627368, 1.3235575985980974, 0.92464405, 0.85245407, 0.92464405, 0.85245407, 'empty', 'empty'], [0.5333411799130342, -0.8459002221355498, 2.475619291341845, 0.93114704, 0.99878323, 0.93114704, 0.99878323, 'empty', 'empty'], [0.7654170555140755, -0.643534560943049, 0.8397285327052333, 0.97927785, 0.016387668, 0.97927785, 0.016387668, 'empty', 'empty'], [0.9907478953031172, -0.1357151721527239, 0.21370618755120446, 0.95502526, 0.08288439, 0.95502526, 0.08288439, 'empty', 'empty'], [0.3652059436828888, -0.9309267525958692, 2.664605663882084, 0.9785185, 0.17031702, 0.9785185, 0.17031702, 'empty', 'empty'], [0.5803029759412476, -0.8144006729575631, 1.655600633298722, 0.98741037, 0.269863, 0.98741037, 0.269863, 'empty', 'empty'], [0.9330456838728287, -0.3597579072185702, 0.8580437701164458, 0.9743627, 0.48235512, 0.9743627, 0.48235512, 'empty', 'empty'], [0.9879077315958443, 0.1550429419649707, 0.4255952957581143, 0.9999668, 0.5825308, 0.9999668, 0.5825308, 'empty', 'empty'], [0.9954493992537509, 0.09529162358437604, 0.5445881694054757, 0.9514868, 0.6356714, 0.9514868, 0.6356714, 'empty', 'empty'], [0.5409877712929447, 0.841030458016529, -0.7871350746768122, 0.9839788, 0.74257815, 0.9839788, 0.74257815, 'empty', 'empty'], [0.9722053008110456, 0.23412999183125657, 0.5748584049885257, 0.9894701, 0.8131462, 0.9894701, 0.8131462, 'empty', 'empty'], [0.675362818809183, -0.7374856357720566, 1.9606512713679696, 0.9573744, 0.9152134, 0.9573744, 0.9152134, 'empty', 'empty'], [0.7887505288738617, 0.6147134317722394, 0.19976888696146955, 0.9988849, 0.97825074, 0.9988849, 0.97825074, 'empty', 'empty']]
				lines_done_up = []
				lines_done_down = []
				len_seeds = len(Seeds)
				intersection_points = np.zeros((len_seeds,len_seeds),dtype=object)
				for i in range(len_seeds):
					for j in range(len_seeds):
						if i!=j:
							x_inter = (lines[j][2]-lines[i][2])/(lines[i][1]/lines[i][0]-lines[j][1]/lines[j][0])
							y_inter = lines[j][1]/lines[j][0]*x_inter+lines[j][2]
							intersection_points[i,j] = [x_inter,y_inter]
				#plt.scatter(intersection_points[5,20][0],intersection_points[5,20][1], color='red')
				for k in range(int(div_space*1.5)):
					#print(k)
					for i in range(len_seeds):
						if i not in lines_done_up: # look at one side of each initial node
							xk = Seeds[i][0]+k/div_space * lines[i][0]
							yk = Seeds[i][1]+k/div_space * lines[i][1]
							#if i==5: plt.scatter(xk,yk, color='green')
							#if i==5: plt.annotate(k,(xk,yk), color='green')
							#if i==20: plt.scatter(xk,yk, color='blue')
							#if i==20: plt.annotate(k,(xk,yk), color='blue')
							if xk >= 1.0: # if the line goes out of the domain, stop the line
								vertices[i][0] = 1.0; lines[i][5] = 1.0
								vertices[i][1] = lines[i][1]/lines[i][0]*1.0+lines[i][2]; lines[i][6]=lines[i][1]/lines[i][0]*1.0+lines[i][2]
								lines_done_up.append(i)
							else:
								lines[i][5]=xk
								lines[i][6]=yk
								#if i==13 and k ==7 : print(xk,yk)
								for j in range(len_seeds):
									if i !=j:
										#if i==13 and k ==7 and j==19: 
										#	print('j', j, intersection_points[13,19], xk-1./div_space*lines[i][0],yk-1./div_space*lines[i][1])
										#	print(lines[j][3]-1./div_space*lines[j][0],lines[j][5]+1./div_space*lines[j][0])
										#	print(lines[j][4]-1./div_space*lines[j][1],yk, lines[j][6]+1./div_space*lines[j][1])
										#	print(lines[j][4]-1./div_space*lines[j][1]>= yk >=lines[j][6]+1./div_space*lines[j][1])
										x_inter = intersection_points[i,j][0]
										y_inter = intersection_points[i,j][1]
										if lines[j][3]-1./div_space*lines[j][0]<= x_inter <=lines[j][5]+1./div_space*lines[j][0] or lines[j][3]-1./div_space*lines[j][0]>= x_inter >=lines[j][5]+1./div_space*lines[j][0]: 
											if lines[j][4]-1./div_space*lines[j][1]<= y_inter <=lines[j][6]+1./div_space*lines[j][1] or lines[j][4]-1./div_space*lines[j][1]>= y_inter >=lines[j][6]+1./div_space*lines[j][1]:
												#print('hi')
												x_inter = intersection_points[i,j][0]
												y_inter = intersection_points[i,j][1]
												if xk-1./div_space*lines[i][0]<=x_inter<=xk or xk<=x_inter<=xk-1./div_space*lines[i][0]:
													if yk-1./div_space*lines[i][1]<=y_inter<=yk or yk<=y_inter<=yk-1./div_space*lines[i][1]:
														#print(i,j)
														if i in lines_done_up: 
															if length_square(intersection_points[i,j]-Seeds[i])<length_square(intersection_points[i,lines[i][7]]-Seeds[i]):
																vertices[i][0] = x_inter ; lines[i][5] = x_inter
																vertices[i][1] = y_inter ; lines[i][6] = y_inter
																lines[i][7] = j
																#lines_done_up.append(i)
														else:
																vertices[i][0] = x_inter ; lines[i][5] = x_inter
																vertices[i][1] = y_inter ; lines[i][6] = y_inter
																lines[i][7] = j
																lines_done_up.append(i)
									if k == int(div_space*1.5)-1 and i not in lines_done_up:
										vertices[i][0] = xk
										vertices[i][1] = yk
						if i not in lines_done_down: # look at the other side of each initial node
							xk = Seeds[i][0]-k/div_space * lines[i][0]
							yk = Seeds[i][1]-k/div_space * lines[i][1]
							#if i==5: plt.scatter(xk,yk, color='green')
							#if i==5: plt.annotate(k,(xk,yk), color='green')
							#if i==20: plt.scatter(xk,yk, color='blue')
							#if i==20: plt.annotate(k,(xk,yk), color='blue')
							if xk <=0.0:
								vertices[len_seeds+i][0] = 0.0; lines[i][3] = 0.0
								vertices[len_seeds+i][1] = lines[i][2]; lines[i][4] = lines[i][2]
								lines_done_down.append(i)
							else:
								lines[i][3]=xk
								lines[i][4]=yk
								for j in range(len_seeds):
									if i !=j:
										x_inter = intersection_points[i,j][0]
										y_inter = intersection_points[i,j][1]
										if lines[j][3]-1./div_space*lines[j][0]<= x_inter <=lines[j][5]+1./div_space*lines[j][0] or lines[j][3]-1./div_space*lines[j][0]>= x_inter >=lines[j][5]+1./div_space*lines[j][0]:
											if lines[j][4]-1./div_space*lines[j][1]<= y_inter <=lines[j][6]+1./div_space*lines[j][1] or lines[j][4]-1./div_space*lines[j][1]>= y_inter >=lines[j][6]+1./div_space*lines[j][1]:
												x_inter = intersection_points[i,j][0]
												y_inter = intersection_points[i,j][1]
												if xk+1./div_space*lines[i][0]<=x_inter<=xk or xk<=x_inter<=xk+1./div_space*lines[i][0]:
													if yk+1./div_space*lines[i][1]<=y_inter<=yk or yk<=y_inter<=yk+1./div_space*lines[i][1]:
														#print(i,j)
														if i in lines_done_down: 
															#print('in', lines[i][8], length_square(intersection_points[i,j]-Seeds[i]),length_square(intersection_points[i,lines[i][8]]-Seeds[i]))
															if length_square(intersection_points[i,j]-Seeds[i])<length_square(intersection_points[i,lines[i][8]]-Seeds[i]):#
																vertices[len_seeds+i][0] = x_inter; lines[i][3] = x_inter
																vertices[len_seeds+i][1] = y_inter; lines[i][4] = y_inter
																lines[i][8]=j
																#lines_done_down.append(i)
														else:
															vertices[len_seeds+i][0] = x_inter; lines[i][3] = x_inter
															vertices[len_seeds+i][1] = y_inter; lines[i][4] = y_inter
															lines[i][8]=j
															lines_done_down.append(i)
									if k == int(div_space*1.5)-1 and i not in lines_done_down:
										vertices[len_seeds+i][0] = xk
										vertices[len_seeds+i][1] = yk
					if len(lines_done_down)==len_seeds and len(lines_done_up)==len_seeds: break
				ridge_vertices = []
				list_points = []
				import operator
				for i in range(len_seeds):
					list_point =[[i, vertices[i][0]], [len_seeds + i, vertices[len_seeds + i][0]]]
					for j in range(len_seeds):
						if lines[j][7] == i:
							list_point.append([j, vertices[j][0]])
						if lines[j][8] == i:
							list_point.append([len_seeds+j, vertices[len_seeds+j][0]])
					list_point = sorted(list_point,key=operator.itemgetter(1))
					list_points.append(list_point)
				for list_point in list_points:
					for i in range(len(list_point)-1):
						ridge_vertices.append([list_point[i][0],list_point[i+1][0]])
			elif network.dimension == 3:
				len_seeds = len(Seeds)
				vertices = np.zeros((len_seeds*2,network.dimension))
				import random
				from mpl_toolkits.mplot3d import Axes3D
				diameter_fiber = 20.
				for i in range(len_seeds):
					phi = random.uniform(-1,1)*np.pi/2.
					theta = random.uniform(-1,1)*np.pi/2.
					lines.append([np.cos(phi)*np.sin(theta),np.sin(phi)*np.sin(theta),np.cos(theta), Seeds[i][0],Seeds[i][1],Seeds[i][2],Seeds[i][0],Seeds[i][1],Seeds[i][2],'empty','empty'])
				lines_done_up = []
				lines_done_down = []
				points = []
				for k in range(int(div_space*np.sqrt(3)*1.5)):
					print(k)
					for i in range(len_seeds):
						if i not in lines_done_up:
							xk = Seeds[i][0]+k/div_space*lines[i][0]
							if xk >= network.length[0]:
								vertices[i][0] = 1.0; lines[i][6]=1.0
								vertices[i][1] = Seeds[i][1]+(1.0-Seeds[i][0])*lines[i][1]/lines[i][0]; lines[i][7] = Seeds[i][1]+(1.0-Seeds[i][0])*lines[i][1]/lines[i][0]
								vertices[i][2] = Seeds[i][2]+(1.0-Seeds[i][0])*lines[i][2]/lines[i][0]; lines[i][8] = Seeds[i][2]+(1.0-Seeds[i][0])*lines[i][2]/lines[i][0]
								lines_done_up.append(i)
							elif xk <=0.0:
								vertices[i][0] = 0.0; lines[i][6] = 0.0
								vertices[i][1] = Seeds[i][1]+(-Seeds[i][0])*lines[i][1]/lines[i][0]; lines[i][7] = Seeds[i][1]+(-Seeds[i][0])*lines[i][1]/lines[i][0]
								vertices[i][2] = Seeds[i][2]+(-Seeds[i][0])*lines[i][2]/lines[i][0]; lines[i][8] = Seeds[i][2]+(-Seeds[i][0])*lines[i][2]/lines[i][0]
								lines_done_up.append(i)
							else:
								yk = Seeds[i][1]+k/div_space*lines[i][1]
								zk = Seeds[i][2]+k/div_space*lines[i][2]
								lines[i][6]=xk
								lines[i][7]=yk
								lines[i][8]=zk
								for j in range(len_seeds):
									if i !=j:
										distance_x = diameter_fiber/2.*lines[i][0]
										distance_y = diameter_fiber/2.*lines[i][1]
										distance_z = diameter_fiber/2.*lines[i][2]
										if lines[j][3]-1./div_space*lines[j][0]<= xk<= lines[j][6]+1./div_space*lines[j][0] or  lines[j][3]-1./div_space*lines[j][0]>= xk>= lines[j][6]+1./div_space*lines[j][0]:#see to readd the distance buffer so that the diameter is included
											if lines[j][4]-1./div_space*lines[j][1]<= yk<= lines[j][7]+1./div_space*lines[j][1] or lines[j][4]-1./div_space*lines[j][1]>= yk>= lines[j][7]+1./div_space*lines[j][1]:
												if lines[j][5]-1./div_space*lines[j][2]<= zk<=lines[j][8] +1./div_space*lines[j][2] or lines[j][5]-1./div_space*lines[j][2]>= zk>=lines[j][8]+1./div_space*lines[j][2]:
													vector1 = np.array([lines[j][6]-lines[j][3], lines[j][7]-lines[j][4],lines[j][8]-lines[j][5]]).astype(float)
													vector2 = np.array([lines[j][6]-xk, lines[j][7]-yk,lines[j][8]-zk])
													d = np.linalg.norm(np.cross(vector1,vector2).astype(float))/np.linalg.norm(vector1)
													if d<diameter_fiber/2:
														if i in lines_done_up:
															if np.linalg.norm(vertices[i]-Seeds[i])>d:
																new_vertices = np.dot(-vector1,vector2)/np.linalg.norm(vector1)**2*vector1 + np.array([lines[j][6],lines[j][7],lines[j][8]])
																vertices[i][0] = new_vertices[0]; lines[i][6] = new_vertices[0]
																vertices[i][1] = new_vertices[1]; lines[i][7] = new_vertices[1]
																vertices[i][2] = new_vertices[2]; lines[i][8] = new_vertices[2]
																lines_done_up.append(i)
																lines[i][9] = j
														else:
															new_vertices = np.dot(-vector1,vector2)/np.linalg.norm(vector1)**2*vector1 + np.array([lines[j][6],lines[j][7],lines[j][8]])
															vertices[i][0] = new_vertices[0]; lines[i][6] = new_vertices[0]
															vertices[i][1] = new_vertices[1]; lines[i][7] = new_vertices[1]
															vertices[i][2] = new_vertices[2]; lines[i][8] = new_vertices[2]
															lines_done_up.append(i)
															lines[i][9] = j
									if k == int(div_space*1.5*np.sqrt(3))-1 and i not in lines_done_up:
										vertices[i][0] = xk
										vertices[i][1] = yk
										vertices[i][2] = zk
						if i not in lines_done_down:
							xk = Seeds[i][0]-k/div_space*lines[i][0]
							if xk <=0.0:
								vertices[len_seeds+i][0] = 0.0; lines[i][3] = 0.0
								vertices[len_seeds+i][1] = Seeds[i][1]+(-Seeds[i][0])*lines[i][1]/lines[i][0]; lines[i][4] = Seeds[i][1]+(-Seeds[i][0])*lines[i][1]/lines[i][0]
								vertices[len_seeds+i][2] = Seeds[i][2]+(-Seeds[i][0])*lines[i][2]/lines[i][0]; lines[i][5] = Seeds[i][2]+(-Seeds[i][0])*lines[i][2]/lines[i][0]
								lines_done_down.append(i)
							elif xk >= network.length[0]:
								vertices[len_seeds+i][0] = 1.0; lines[i][3]=1.0
								vertices[len_seeds+i][1] = Seeds[i][1]+(1.0-Seeds[i][0])*lines[i][1]/lines[i][0] ; lines[i][4]=Seeds[i][1]+(1.0-Seeds[i][0])*lines[i][1]/lines[i][0]
								vertices[len_seeds+i][2] = Seeds[i][2]+(1.0-Seeds[i][0])*lines[i][2]/lines[i][0]; lines[i][5] = Seeds[i][2]+(1.0-Seeds[i][0])*lines[i][2]/lines[i][0]
								lines_done_down.append(i)
							else:
								yk = Seeds[i][1]-k/div_space*lines[i][1]
								zk = Seeds[i][2]-k/div_space*lines[i][2]
								lines[i][3]=xk
								lines[i][4]=yk
								lines[i][5]=zk
								for j in range(len_seeds):
									if i !=j:
										distance_x = diameter_fiber/2.*lines[i][0]
										distance_y = diameter_fiber/2.*lines[i][1]
										distance_z = diameter_fiber/2.*lines[i][2]
										if lines[j][3]-1./div_space*lines[j][0]<= xk<= lines[j][6]+1./div_space*lines[j][0] or  lines[j][3]-1./div_space*lines[j][0]>= xk>= lines[j][6]+1./div_space*lines[j][0]:#see to readd the distance buffer so that the diameter is included
											if lines[j][4]-1./div_space*lines[j][1]<= yk<= lines[j][7]+1./div_space*lines[j][1] or lines[j][4]-1./div_space*lines[j][1]>= yk>= lines[j][7]+1./div_space*lines[j][1]:
												if lines[j][5]-1./div_space*lines[j][2]<= zk<=lines[j][8] +1./div_space*lines[j][2] or lines[j][5]-1./div_space*lines[j][2]>= zk>=lines[j][8]+1./div_space*lines[j][2]:
													vector1 = np.array([lines[j][6]-lines[j][3], lines[j][7]-lines[j][4],lines[j][8]-lines[j][5]]).astype(float)
													vector2 = np.array([lines[j][6]-xk, lines[j][7]-yk,lines[j][8]-zk])
													d = np.linalg.norm(np.cross(vector1,vector2).astype(float))/np.linalg.norm(vector1)
													if d<diameter_fiber/2:
														if i in lines_done_down:
															if np.linalg.norm(vertices[len_seeds+i]-Seeds[i])>d:
																new_vertices = np.dot(-vector1,vector2)/np.linalg.norm(vector1)**2*vector1 + np.array([lines[j][6],lines[j][7],lines[j][8]])
																vertices[len_seeds+i][0] = new_vertices[0]; lines[i][3] = new_vertices[0]
																vertices[len_seeds+i][1] = new_vertices[1]; lines[i][4] = new_vertices[1]
																vertices[len_seeds+i][2] = new_vertices[2]; lines[i][5] = new_vertices[2]
																lines[i][10]=j
																lines_done_down.append(i)
														else:
															new_vertices = np.dot(-vector1,vector2)/np.linalg.norm(vector1)**2*vector1 + np.array([lines[j][6],lines[j][7],lines[j][8]])
															vertices[len_seeds+i][0] = new_vertices[0]; lines[i][3] = new_vertices[0]
															vertices[len_seeds+i][1] = new_vertices[1]; lines[i][4] = new_vertices[1]
															vertices[len_seeds+i][2] = new_vertices[2]; lines[i][5] = new_vertices[2]
															lines[i][10]=j
															lines_done_down.append(i)
									if k == int(div_space*np.sqrt(3)*1.5)-1 and i not in lines_done_down:
										vertices[i][0] = xk
										vertices[i][1] = yk
										vertices[i][2] = zk
				ridge_vertices = []
				list_points = []
				import operator
				for i in range(len_seeds):
					list_point =[[i, vertices[i][0]], [len_seeds + i, vertices[len_seeds + i][0]]]
					for j in range(len_seeds):
						if lines[j][9] == i:
							list_point.append([j, vertices[j][0]])
						if lines[j][10] == i:
							list_point.append([len_seeds+j, vertices[len_seeds+j][0]])
					list_point = sorted(list_point,key=operator.itemgetter(1))
					list_points.append(list_point)
				for list_point in list_points:
					for i in range(len(list_point)-1):
						ridge_vertices.append([list_point[i][0],list_point[i+1][0]])
			return vertices,ridge_vertices



