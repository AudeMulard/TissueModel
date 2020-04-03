import numpy as np
from scipy.spatial import Voronoi, voronoi_plot_2d
import matplotlib.pyplot as plt

def length_square(x):
	return x[0]**2+x[1]**2

def select_network(network,creation):
		
		if creation == "Voronoi":
			Seeds = []
			for node in range(network.complexity):
				if network.dimension==2:
					Seeds.append([np.random.rand()*network.length[0],np.random.rand()*network.length[1]])
				elif network.dimension ==3:
					Seeds.append([np.random.rand()*network.length[0],np.random.rand()*network.length[1],np.random.rand()*network.length[2]])
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
			Seeds = []
			for node in range(network.complexity):
				if network.dimension==2:
					Seeds.append([np.random.rand()*network.length[0],np.random.rand()*network.length[1]])
				elif network.dimension ==3:
					Seeds.append([np.random.rand()*network.length[0],np.random.rand()*network.length[1],np.random.rand()*network.length[2]])
			Seeds = np.array(Seeds)
			#Seeds = np.array([[ 0.92374815,  0.90291537], [ 0.87529856,  0.91491711], [ 0.59413968,  0.65649066], [ 0.70452016,  0.04474947], [ 0.79686153,  0.75654823]])
			#plt.scatter(Seeds[:,0],Seeds[:,1],label='Seeds',color='blue')
			orientations = np.random.rand(network.complexity,1)*180
			#orientations = [[  44.12501178], [  99.64676771], [  38.26756236], [ 142.78837107], [  42.27702724]]
			lines = []
			for i in range(len(Seeds)):
				a = float(np.tan(orientations[i]))
				b = float(Seeds[i][1]-np.tan(orientations[i])*Seeds[i][0])
				lines.append([a,b,0.0,b,i,1.0, a+b, i])
				#plt.plot([[0.0,Seeds[i][1]-np.tan(orientations[i])*Seeds[i][0]],[1.0, np.tan(orientations[i])+Seeds[i][1]-np.tan(orientations[i])*Seeds[i][0]]])
			inter_points = []
			for line1 in lines:
				for line2 in lines:
					if line1 != line2:
						#print line1,line2
						x_inter = (line2[1]-line1[1])/(line1[0]-line2[0])
						y_inter = line1[0]*x_inter+line1[1]
						#inter_points.append([x_inter,y_inter])
						if min(1.0,max(line2[5],line2[2])) > x_inter > max(Seeds[lines.index(line1)][0],min(line2[5],line2[2])) and max(0.0,min(line2[6],line2[3])) < y_inter < min(1.0,max(line2[6],line2[3])) and length_square(Seeds[lines.index(line1)]-[x_inter,y_inter])<=length_square(Seeds[lines.index(line1)]-[line1[5],line1[6]]):
							print y_inter
							line1[5] = x_inter
							line1[6] = y_inter
							line1[7] = lines.index(line2)
						elif min(Seeds[lines.index(line1)][0],max(line2[5],line2[2])) > x_inter > max(0.0,min(line2[5],line2[2])) and max(0.0,min(line2[6],line2[3])) < y_inter < min(1.0,max(line2[6],line2[3])) and length_square(Seeds[lines.index(line1)]-[x_inter,y_inter])<=length_square(Seeds[lines.index(line1)]-[line1[2],line1[3]]):
							print y_inter
							line1[2] = x_inter
							line1[3] = y_inter
							line1[4] = lines.index(line2)
						plt.scatter(x_inter,y_inter)
				plt.scatter(Seeds[lines.index(line1)][0],Seeds[lines.index(line1)][1],color='blue')
				plt.scatter(line1[2],line1[3], color='red')
				plt.scatter(line1[5],line1[6], color='red')
				plt.show()
			#inter_points = np.array([list(l) for l in inter_points])
			#plt.scatter(inter_points[:,0],inter_points[:,1],label = 'inter_points')
			lines = np.array(lines)
			vertices = []
			ridge_vertices=[]
			for k in range(len(lines)):
				vertices.append([lines[k,2], lines[k,3]])
				vertices.append([lines[k,5], lines[k,6]])
			for line in lines:
				ridge_vertices.append([vertices.index([line[2],line[3]]), vertices.index([lines[int(line[4])][2],lines[int(line[4])][3]])])
				ridge_vertices.append([vertices.index([line[2],line[3]]), vertices.index([lines[int(line[4])][5],lines[int(line[4])][6]])])
				ridge_vertices.append([vertices.index([line[5],line[6]]), vertices.index([lines[int(line[7])][2],lines[int(line[7])][3]])])
				ridge_vertices.append([vertices.index([line[5],line[6]]), vertices.index([lines[int(line[7])][5],lines[int(line[7])][6]])])
			vertices = np.array(vertices)
			print vertices
			for i in range(len(vertices[:,0])):
				plt.annotate(i, (vertices[i,0],vertices[i,1]),fontsize=10)
			plt.scatter(vertices[:,0],vertices[:,1],label='vertices')
			plt.legend()
			#plt.show()
			print ridge_vertices
			return vertices,ridge_vertices
				

		if creation == "reg_Voronoi":
			Seeds = []
			for x in range(50):
				for y in range(50):
					Seeds.append([x/10.,y/20.-0.05])
			for x in range(50):
				for y in range(50):
					Seeds.append([x/10.-0.05,y/20.-0.025])
			#print Seeds
			#Seeds=np.random.rand(network.complexity,network.dimension)*network.length
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
		
		if creation == "disturbed_grid":
			Seeds = []
			for x in range(network.complexity):
				for y in range(network.complexity):
					Seeds.append([x/20.+np.random.rand()*network.disturbance,y/20.+np.random.rand()*network.disturbance])
			voronoi = Voronoi(Seeds)
			vertices = Voronoi(Seeds).vertices
			ridge_vertices = Voronoi(Seeds).ridge_vertices
			return vertices, ridge_vertices

		if creation == "big grid":
			vertices = np.array([[0.0,0.25],[0.0,0.5],[0.0,0.75],[0.25,0.0],[0.25,0.25],[0.25,0.5],[0.25,0.75],[0.25,1.0],[0.5,0.0],[0.5,0.25],[0.5,0.5],[0.5,0.75],[0.5,1.0],[0.75,0.0],[0.75,0.25],[0.75,0.5],[0.75,0.75],[0.75,1.0],[1.0,0.25],[1.0,0.5],[1.0,0.75]])
			ridge_vertices = np.array([[0,4],[4,9],[9,14],[14,18],[1,5],[5,10],[10,15],[15,19],[2,6],[6,11],[11,16],[16,20],[3,4],[4,5],[5,6],[6,7],[8,9],[9,10],[10,11],[11,12],[13,14],[14,15],[15,16],[16,17]])
			return vertices, ridge_vertices

		if creation == "cut big grid":
			vertices = np.array([[0.0,0.25],[0.0,0.5],[0.0,0.75],[0.25,0.0],[0.25,0.25],[0.25,0.5],[0.25,0.75],[0.25,1.0],[0.5,0.0],[0.5,0.25],[0.5,0.5],[0.5,0.75],[0.5,1.0],[0.75,0.0],[0.75,0.25],[0.75,0.5],[0.75,0.75],[0.75,1.0],[1.0,0.25],[1.0,0.5],[1.0,0.75]])
			ridge_vertices = np.array([[0,4],[4,9],[9,14],[14,18],[1,5],[5,10],[10,15],[15,19],[2,6],[6,11],[11,16],[16,20],[3,4],[4,5],[5,6],[6,7],[8,9],[9,10],[10,11],[11,12],[13,14],[14,15],[15,16],[16,17]])
			return vertices, ridge_vertices

		if creation == "vertical line":
			vertices = np.array([[0.0,0.0],[0.0,1.0],[0.5,0.0],[0.5,1.0],[1.0,0.0],[1.0,1.0]])
			ridge_vertices = np.array([[0,2],[2,4],[1,3],[3,5],[2,3]])
			return vertices, ridge_vertices

		if creation == "small grid":
			vertices = np.array([[0.0,0.33],[0.0,0.66],[0.33,0.33],[0.33,0.66],[0.66,0.33],[0.66,0.66],[1.0,0.33],[1.0,0.66]])
			ridge_vertices = np.array([[0,2],[2,4],[4,6],[1,3],[3,5],[5,7],[2,3],[4,5]])
			return vertices, ridge_vertices

		if creation == "1 straight line":
			vertices = np.array([[0.0,0.0],[0.5,0.0],[1.0,0.0]])
			ridge_vertices = np.array([[0,1],[1,2]])
			return vertices, ridge_vertices

		if creation == "1 long straight line":
			vertices = np.array([[0.0,0.0],[0.25,0.0],[0.5,0.0],[0.75,0.0],[1.0,0.0]])
			ridge_vertices = np.array([[0,1],[1,2],[2,3],[3,4]])
			return vertices, ridge_vertices

		if creation == "1 long cross":
			vertices = np.array([[0.0,0.5],[0.5,0.5],[1.0,0.5],[0.5,0.0],[0.5,0.25],[0.5,0.75],[0.5,1.0]])
			ridge_vertices = np.array([[0,1],[1,2],[3,4],[4,1],[1,5],[5,6]])
			return vertices, ridge_vertices

		if creation == "1 half long cross":
			vertices = np.array([[0.0,0.5],[0.25,0.5],[0.5,0.5],[0.75,0.5],[1.0,0.5],[0.5,0.25],[0.5,0.75]])
			ridge_vertices = np.array([[0,1],[1,2],[2,3],[3,4],[5,6],[5,2],[2,6]])
			return vertices, ridge_vertices

		if creation == "2 manual":
			vertices = np.array([[0.0,0.0],[0.3,0.5],[1.0,0.0]])
			ridge_vertices = np.array([[0,1],[1,2]])
			return vertices, ridge_vertices

		if creation == "3 manual":
			vertices = np.array([[0.0,0.0],[0.3,0.1],[0.6,0.3],[1.0,0.0]])
			ridge_vertices = np.array([[0,1],[1,2],[2,3]])
			return vertices, ridge_vertices

		if creation == "1 straight line":
			vertices = np.array([[0.0,0.0],[0.5,0.0],[1.0,0.0]])
			ridge_vertices = np.array([[0,1],[1,2]])
			return vertices, ridge_vertices

		if creation == "symmetric cross":
			vertices = np.array([[0.0,0.0],[0.5,0.0],[1.0,0.5],[1.0,-0.5]])
			ridge_vertices = np.array([[0,1],[1,2],[1,3]])
#			self.vertices_ini = np.array(self.vertices.tolist())
			return vertices, ridge_vertices
	
		if creation == "asymmetric cross":
			vertices = np.array([[0.0,0.0],[0.7,-0.1],[1.0,0.2],[1.0,-0.3]])
			ridge_vertices = np.array([[0,1],[1,2],[1,3]])
			return vertices, ridge_vertices

		if creation == "symmetric double point cross":
			vertices = np.array([[0.0,0.0],[0.0,1.0],[0.61,0.5],[0.69,0.5],[1.0,0.0],[1.0,1.0]])
			ridge_vertices = np.array([[0,2],[1,2],[2,3],[3,5],[3,4]])
			return vertices, ridge_vertices

		if creation == "asymmetric double point cross":
			vertices = np.array([[-1.0,0.0],[-0.6,1.0],[0.19, 0.45],[0.22,0.51], [0.20,0.55],[0.69,0.62],[1.3,0.0],[2.0,1.0]])
			ridge_vertices = np.array([[0,2],[1,3],[2,3],[2,4],[3,4],[4,5],[4,6],[5,6],[5,7]])
			return vertices, ridge_vertices

		if creation == "random network":
			vertices = np.array([[6.96810003e-01,5.85972462e-01], [  6.12222074e-01,   5.70042014e-01], [  1.53524456e-02,   8.24030940e-01], [  2.02945404e-02,   3.95001038e-01], [  5.91361447e-02,   6.07631016e-01], [  2.17734065e-01,   4.04153611e-01], [  1.00918445e-01,   6.11192022e-01], [  7.60619750e-01,  -4.27877180e-02], [  6.42935772e-01,   1.57907165e-01], [  8.73596877e-01,   3.18524365e-01], [  8.51906402e-01,   3.31465928e-01], [  4.47723467e-01,   3.51768688e-02], [  5.85307175e-01,   2.02853812e-01], [  4.64801092e-01,   3.75003275e-01], [  4.15171979e-01,   3.63331353e-01], [  3.16687341e-01,   1.87921226e-01], [  2.91696122e-01,   2.93138477e-01], [  3.09155556e-01,   3.21897785e-01], [  5.91396931e-01,   5.94012634e-01], [  8.14387805e-01,   8.72239083e-01], [  3.15164596e-01,   9.25079572e-01], [  1.13885589e-01,   6.27208867e-01], [  3.00439561e-01,   5.75803663e-01], [  4.17769233e-01,   7.61100709e-01], [  3.17801428e-01,   8.99124444e-01], [  2.97685323e-01,   6.19842149e-01], [  1.43974960e-01,   6.85770872e-01], [  1.24611979e-01,   6.43587311e-01]])
			ridge_vertices=np.array([[-1, 4], [-1, 3], [3, 4], [3, 5], [4, 6], [5, 6], [7, 9], [7, 8], [8, 10], [9, 10], [-1, 9], [-1, 7], [-1, 0], [0, 10], [11, 15], [11, 12], [12, 13], [13, 14], [14, 17], [15, 16], [16, 17], [-1, 15], [-1, 11], [8, 12], [0, 1], [1, 13], [-1, 16], [5, 17], [-1, 19], [1, 18], [18, 19], [-1, 2], [2, 20], [19, 20], [6, 21], [14, 22], [21, 22], [23, 24], [23, 25], [24, 26], [25, 27], [26, 27], [18, 23], [20, 24], [22, 25], [2, 26], [21, 27]])
			return vertices, ridge_vertices

		if creation == "1 straight line 2":
			vertices = np.array([[-1.0,0.0],[0.5,0.0],[1.5,0.0]])
			ridge_vertices = np.array([[0,1],[1,2]])
			return vertices, ridge_vertices

		if creation == "symmetric cross 2":
			vertices = np.array([[-1.0,0.0],[0.5,0.0],[1.5,0.5],[1.5,-0.5]])
			ridge_vertices = np.array([[0,1],[1,2],[1,3]])
			return vertices, ridge_vertices

		if creation == "random network 2":
			vertices = np.array([[ 0.6886007,   1.22144112], [ 0.69731848,  1.1186779 ], [ 0.51616784,  0.63207458], [-0.25276059,  0.45293642], [ 0.49692315, -1.5152983 ], [ 0.34634344,  0.13602371], [ 0.33809343,  0.80025459], [ 0.19167735,  0.81100432], [-0.08463533,  0.51910635], [ 1.42874148,  0.32477703], [ 0.44490991,  0.64721725], [ 0.3765596,   0.67967657], [ 0.22545142,  0.52383797], [ 0.71710902,  1.04705891], [ 1.28244566,  0.44924527], [ 0.75903549,  0.8663037 ], [ 0.52775066,  0.61855448], [ 0.67122359, 0.02915664], [ 0.571157,    0.21613333], [ 0.36876754,  0.43831976], [ 0.35649777,  0.43378642],[ 0.54799455,  0.46092557], [ 0.56521616,  0.44491184], [ 0.51447422,  0.24044857], [ 0.34884099,  0.19190268], [ 0.3510838,   0.18885809], [ 0.61471224,  0.44134212], [ 0.84182501,  0.33648427], [ 0.75687701,  0.64930453], [ 0.74807286,  0.64506609], [ 1.05464284,  0.36505829], [ 0.91089201,  0.3296533 ]])
			ridge_vertices = np.array([[-1, 4], [-1, 5], [4, 5], [0, 6], [-1, 0], [-1, 7], [6, 7], [-1, 3], [3, 8], [7, 8], [0, 1], [1, 2], [2, 10], [6, 11], [10, 11], [8, 12], [11, 12], [13, 15], [13, 14], [14, 15], [-1, 9], [1, 13], [9, 14], [19, 20], [19, 21], [20, 24], [21, 22], [22, 23], [23, 25], [24, 25], [10, 19], [12, 20], [2, 16], [16, 21], [3, 24], [4, 17], [5, 25], [17, 18], [18, 23], [26, 29], [26, 27], [27, 31], [28, 30], [28, 29], [30, 31], [16, 29], [22, 26], [18, 27], [17, 31], [9, 30], [15, 28]])
			return vertices, ridge_vertices

		if creation == "random network 3":
			vertices = np.array([[ 0.96747413,  1.30350344], [ 0.71732007, -1.68553583], [ 0.66949702,  0.53769709], [ 0.26746628,  1.05725235], [ 0.22433102,  0.92974385], [ 0.54132183,  0.63446214], [ 0.37429159,  0.5678483 ], [ 1.25975916,  0.35115262], [ 0.77160248,  0.67884074], [ 0.68232243,  0.54247464], [ 0.66374711,  0.52136863], [ 0.58848258, -0.16118395], [ 0.16895845,  0.69422255], [-0.17027436,  0.82316349], [ 0.51283584,  0.69352915], [ 0.57329291,  0.90709476], [ 0.40196155,  0.94183993], [ 0.4250409,   0.32231426], [ 0.40522241,  0.19204137], [ 0.25371474,  0.4674381 ], [ 0.14688212,  0.47416062], [ 0.07251124,  0.23379788], [-0.043192,    0.33432779], [ 1.14379205,  0.31863161], [ 0.67418869,  0.0693474 ], [ 0.71738841,  0.23097817], [ 0.41781785,  0.77486727], [ 0.50973872,  0.69527462], [ 0.23891247,  0.82375975], [ 0.1903611,   0.71261965], [ 0.25607053,  0.68906206]])
			ridge_vertices = np.array([[-1, 7], [-1, 8], [7, 9], [8, 9], [-1, 0], [0, 8], [-1, 3], [3, 4], [4, 13], [-1, 13], [0, 15], [2, 9], [2, 5], [5, 14], [14, 15], [3, 16], [15, 16], [17, 18], [17, 19], [18, 21], [19, 20], [20, 22], [21, 22], [2, 10], [5, 6], [6, 19], [10, 17], [1, 11], [1, 21], [11, 18], [12, 13], [12, 20], [-1, 22], [-1, 1], [23, 25], [23, 24], [24, 25], [7, 23], [10, 25], [11, 24], [26, 27], [26, 28], [27, 30], [28, 29], [29, 30], [14, 27], [16, 26], [4, 28], [6, 30], [12, 29]])
			return vertices, ridge_vertices

		if creation == "3d random":
			vertices = np.array([[0.6278694211103966, 0.24678507302282998, 0.7325883165906247], [0.40377824097623816, 0.47412042047443537, 0.7766300136696769], [0.06338111085118947, 0.15596552277517112, 0.6614135118426991], [0.6622292581241507, 0.4784196253580352, 3.6821043305610774], [0.3518731102795726, 0.7331321179162598, 3.7509605040903953], [0.4068731477447932, 0.4793301715293799, 0.7809808208750539], [0.37729780988585193, 0.49030084334726604, 0.8260182582776893], [0.3430887236649174, 0.5201820202200497, 0.8573530707695022], [-0.12846213610137913, 1.4083234432598106, 0.9384908038250677], [0.24793075366807435, 0.6961328709250224, 0.8324216857940325], [0.5541003142096788, 0.4320136642365377, 0.8126522394347525], [1.1399969470543374, 0.2657605817261536, 1.71181057968245], [0.09886819229282612, 1.2428927323509478, 0.6684154915242873], [0.5805601626583469, 0.6025622217185238, 0.7245037139746031], [0.7313701162354891, 0.6917975346803149, 0.86047511714591], [0.3340405830887804, 0.7645844455547355, 0.6280802355348203], [0.3510330591546621, 0.6441693042285355, 0.6919825514323923], [0.3028720067687286, 0.5137715691941943, 0.2410244322397597], [0.3027978726453227, 0.5899642209532896, 0.6230789727107579], [0.2990295049650231, 0.5602321888958395, 0.6187300681606293], [-1.7030528083655483, -1.3834174288744845, -0.010761337861453446], [0.9253016003257195, 0.31184430336736335, 0.40117702217489726], [0.7269990020216389, 0.574391800865046, 0.4804238345776466], [0.6669699675064145, 0.20169133763894195, 0.6768388613113479], [0.538647853544071, -0.5217159003140682, -0.34124082410202916], [0.3744545726460393, 0.4407985342239835, 0.223819314675487], [1.0565794923918652, -0.9735171858490231, 0.4327540396155169], [1.128062731977535, -0.6677630242667772, 0.5913293361063472], [0.7480562856385427, 0.008708921374270373, 0.6770859871276028], [0.740504712888167, 0.5935173466140313, 0.14195387533239887], [0.7537255478684216, 0.6062565717524786, 0.16512333192313194], [1.47637849787721, 0.40482460061793685, 0.4668427569647079], [1.3279940326565134, 0.9266966297583683, 0.08449744659796876], [1.2984803644477583, 1.0621217974286354, 1.051582965992494], [0.7135446726931893, 0.654972827262529, 0.20148340820727062], [-1.2199199959266949, 2.6968619375767084, -1.3484684521644723], [0.3145578710762814, 0.6066274764778224, 0.2577738608457551], [0.7007226698167393, 0.6503840534723225, 0.1884972070617002], [0.5082709776854598, 0.5999302134797301, 0.24682868389781157], [0.3298847950446264, 0.7726023927195362, 0.6185396138202078], [0.32659746191879296, 0.7702879441813375, 0.61345176109884], [0.6662751660645184, 0.6797293491851745, 0.558149757560152], [0.6990430702377806, 0.6573511547536317, 0.5075643356197435], [0.566784974613209, 0.7331221246872901, 0.5371957705923112], [0.7369550076669578, 0.7278986776654902, 0.5634044517682504], [0.7277998677711747, 0.6845997157740601, 0.5187748530383501], [0.685185014508221, 0.6547293915945818, 0.383996608328494], [0.6785649878549804, 0.6576841117742704, 0.3754296651024585], [0.6866870621647759, 0.6586311594608923, 0.35671464852153606], [0.9488786584389965, 0.9274838476002214, 0.2294148248276664], [0.7247515767760557, 0.7081623983572901, 0.22017820296006765]])
			ridge_vertices = np.array([[3, 6], [6, 7], [7, 4], [3, 11], [11, 14], [14, 12], [12, 8], [8, 4], [3, 11], [11, 10], [10, 5], [5, 6], [4, 8], [8, 9], [9, 7], [5, 6], [6, 7], [7, 9], [9, 16], [5, 10], [10, 13], [13, 15], [15, 16], [8, 12], [12, 15], [15, 16], [16, 9], [10, 11], [11, 14], [14, 13], [12, 15], [15, 13], [13, 14], [-1, 8], [8, 4], [-1, 9], [9, 8], [-1, 7], [7, 4], [-1, 9], [9, 7], [-1, 3], [3, 4], [-1, 2], [2, 6], [6, 3], [2, 6], [6, 7], [7, -1], [-1, 20], [2, -1], [-1, 20], [1, 2], [2, 20], [20, 19], [1, 5], [5, 16], [16, 18], [18, 19], [1, 5], [5, 6], [6, 2], [-1, 9], [9, 16], [16, 18], [-1, 18], [18, 19], [19, 20], [-1, 26], [26, 28], [28, 0], [0, 1], [1, 2], [-1, 26], [26, 24], [0, 28], [28, 23], [0, 1], [1, 19], [19, 17], [17, 25], [25, 23], [-1, 24], [24, 25], [25, 17], [-1, 20], [20, 19], [19, 17], [23, 25], [25, 24], [24, 26], [26, 28], [-1, 27], [27, 26], [0, 10], [10, 11], [11, -1], [-1, 27], [27, 28], [0, 1], [1, 5], [5, 10], [-1, 3], [3, 11], [26, 27], [27, 28], [-1, 27], [27, 21], [21, 31], [-1, 24], [24, 29], [-1, 32], [32, 30], [30, 29], [21, 27], [27, 26], [26, 24], [24, 29], [29, 30], [21, 31], [31, 32], [32, 30], [-1, 32], [32, 31], [-1, 35], [35, 36], [36, 17], [-1, 18], [18, 40], [17, 19], [19, 18], [18, 40], [40, 36], [35, -1], [-1, 40], [40, 36], [-1, 29], [29, 37], [37, 35], [17, 25], [25, 38], [38, 36], [24, 29], [29, 37], [37, 38], [38, 25], [35, 37], [37, 38], [38, 36], [-1, 8], [8, 12], [12, -1], [-1, 39], [39, 15], [15, 39], [39, 40], [40, 18], [18, 16], [-1, 39], [39, 40], [-1, 11], [11, 14], [14, 33], [12, -1], [-1, 33], [33, 14], [41, 44], [44, 45], [45, 42], [41, 43], [43, 47], [47, 46], [46, 42], [41, 43], [43, -1], [-1, 44], [42, 46], [46, 48], [48, 50], [50, 49], [49, 45], [43, -1], [-1, 50], [50, 48], [48, 47], [44, -1], [-1, 49], [49, 45], [46, 47], [47, 48], [-1, 50], [50, 49], [0, 10], [10, 13], [13, 41], [41, 42], [42, 22], [22, 23], [13, 41], [41, 44], [44, 33], [33, 14], [21, 27], [27, 28], [28, 23], [23, 22], [21, 31], [31, 45], [45, 42], [42, 22], [-1, 31], [31, 45], [45, 44], [44, 33], [13, 41], [41, 43], [43, 39], [39, 15], [22, 23], [23, 25], [25, 38], [38, 47], [47, 46], [22, 42], [42, 46], [36, 38], [38, 47], [47, 43], [43, 39], [39, 40], [33, -1], [-1, 44], [39, -1], [-1, 43], [21, 30], [30, 34], [34, 48], [48, 46], [46, 22], [30, 32], [32, 49], [49, 50], [50, 34], [31, 32], [32, 49], [49, 45], [34, 50], [50, 48], [34, 37], [37, 38], [38, 47], [47, 48], [34, 50], [50, -1], [-1, 35], [35, 37], [-1, 32], [32, 49], [29, 30], [30, 34], [34, 37]])
			"""
			new_ridge_vertices = []
			for ridge in ridge_vertices:
				for i in range(len(ridge)-1):
				#	if [ridge[i],ridge[i+1]] not in new_ridge_vertices or [ridge[i+1],ridge[i]] not in new_ridge_vertices:
						new_ridge_vertices.append([ridge[i],ridge[i+1]])
			ridge_vertices = new_ridge_vertices
			"""
			return vertices, ridge_vertices
