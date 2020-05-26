import numpy as np
import matplotlib.pyplot as plt


Seeds = []
fig = plt.figure()
line_segments = []
line_segments_easy = []
ax = fig.gca()
for node in range(50):
	Seeds.append([np.random.rand(),np.random.rand()])
#print Seeds
#Seeds = [[0.7977318360336223, 0.41156693136196887], [0.5378664497295724, 0.9246717646577932], [0.12905048117007756, 0.08959250996676771], [0.025007378232321242, 0.6685836716280847], [0.7173221306608901, 0.00990468462747962], [0.8666020320145748, 0.10992230429330319], [0.5165236968890671, 0.10155483866910864], [0.7137668801911321, 0.991071591130373], [0.07289510266847343, 0.6671075395403172], [0.7201681962925923, 0.5364449265037047]]

Seeds = np.array(Seeds)
orientations = np.random.rand(50,1)*180
#orientations = [[  48.37602184], [  77.37429736], [ 110.63274126], [ 172.23345393], [  74.15807816], [ 100.34703817], [ 151.02175529], [  70.52035998], [  26.19793658], [ 165.73859278]]

print orientations
#ax.scatter(Seeds[:,0],Seeds[:,1],label='Seeds',color='green')
vertices = np.zeros((len(Seeds)*2,2))
lines = []
for i in range(len(Seeds)):
	a = float(np.tan(orientations[i]))
	b = float(Seeds[i][1]-np.tan(orientations[i])*Seeds[i][0])
	lines.append([a,b, Seeds[i][0],Seeds[i][1],Seeds[i][0],Seeds[i][1],'empty','empty','empty','empty' ])
	plt.annotate(i, (Seeds[i][0],Seeds[i][1]),fontsize=10)
	#ax.plot([[0.0,Seeds[i][1]-np.tan(orientations[i])*Seeds[i][0]],[1.0, np.tan(orientations[i])+Seeds[i][1]-np.tan(orientations[i])*Seeds[i][0]]], linestyle = 'dashed')


plt.xlim([0.0,1.0])
plt.ylim([0.0,1.0])
# espace divided par carres, si dans le nouveau carre, il y a deja un point alors couper la fibre sinon continuer
# pour savoi si qqch dans le square, faire toutes les droites et deux cas, soit leur x is still not pas encore passe par la. si x gone par la, calculer y et voir si y appartient au square defined by the 100th

lines_done_up = []
lines_done_down = []

for k in range(500):
	for i in range(len(Seeds)):
		if i not in lines_done_up:
			xk = Seeds[i][0]+k/500.
			#if i==2: print xk
			yk = lines[i][0]*xk+lines[i][1]
			if xk >= 1.0:
				vertices[i][0] = 1.0
				vertices[i][1] = lines[i][0]*1.0+lines[i][1]
				lines_done_up.append(i)
			else:
				lines[i][4]=xk+0.002
				lines[i][5]=lines[i][0]*(xk+0.002)+lines[i][1]
				for j in range(len(Seeds)):
					if i !=j and i not in lines_done_up:
						#if i==7 and j==4: print k, i,j, xk,  (lines[j][1]-lines[i][1])/(lines[i][0]-lines[j][0]), xk+0.002, lines[j][2], lines[j][4]
						if lines[j][2]<= xk <xk+0.002<=lines[j][4]:
							x_inter = (lines[j][1]-lines[i][1])/(lines[i][0]-lines[j][0])
							#if i==8: print i,j, xk, x_inter, xk+0.001
							#plt.scatter(x_inter, y_inter, color = 'blue')
							#print x_inter, y_inter
							if xk<=x_inter<=xk+0.002:
								vertices[i][0] = x_inter
								vertices[i][1] = lines[j][0]*x_inter+lines[j][1]
								lines[i][6] = j
								if x_inter <= Seeds[j][0]:
									lines[i][7]='down'
								else:
									lines[i][7]='up'
								lines_done_up.append(i)
				#lines[i][4]=xk
				#lines[i][5]=yk
		#if i ==2: print xk, lines[i][4]
		if i not in lines_done_down:
			xk = Seeds[i][0]-k/500.
			#if i==7: print 2, xk
			yk = lines[i][0]*xk+lines[i][1]
			if xk <=0.0:
				vertices[len(Seeds)+i][0] = 0.0
				vertices[len(Seeds)+i][1] = lines[i][1] 
				lines_done_down.append(i)
			#plt.scatter(xk, yk, color = 'red')
			else:
				lines[i][2]=xk-0.002
				lines[i][3]=lines[i][0]*(xk-0.002)+lines[i][1]
				for j in range(len(Seeds)):
					if i !=j and i not in lines_done_down:
						#if i==4 and j==7: print k, i,j, xk,  (lines[j][1]-lines[i][1])/(lines[i][0]-lines[j][0]), xk-0.002, lines[j][2], lines[j][4]
						#if i==7 and j==4: print k, i,j, xk,  (lines[j][1]-lines[i][1])/(lines[i][0]-lines[j][0]), xk-0.002, lines[j][2], lines[j][4]
						if lines[j][2]<= xk-0.002<= xk<= lines[j][4]:
							x_inter = (lines[j][1]-lines[i][1])/(lines[i][0]-lines[j][0])
							if xk-0.002<=x_inter<=xk:
								if i==5: print j
								vertices[len(Seeds)+i][0] = x_inter
								vertices[len(Seeds)+i][1] = lines[j][1] + x_inter*lines[j][0]
								lines[i][8]=j
								if x_inter <= Seeds[j][0]:
									lines[i][9]='down'
								else:
									lines[i][9]='up'
								lines_done_down.append(i)
				#lines[i][2]=xk
				#lines[i][3]=yk




print lines_done_down, lines_done_up
print vertices
ridge_vertices = []
ridge_vertices_easy = []

for i in range(int(len(vertices)/2.)):
	ridge_vertices_easy.append([i, len(Seeds)+i])

for i in range(len(vertices)):
	plt.annotate(i, (vertices[i][0],vertices[i][1]),fontsize=10, color = 'green')

"""
for i in range(len(lines)):
	if lines[i][6]!='empty':
		if i==9: print lines[i][6]
		ridge_vertices.append([i,lines[i][6]])
		ridge_vertices.append([i,len(Seeds)+lines[i][6]])
	if lines[i][8]!='empty':
		ridge_vertices.append([len(Seeds)+i,len(Seeds)+lines[i][8]])
		ridge_vertices.append([len(Seeds)+i,lines[i][8]])
"""
list_points = []
import operator
for i in range(len(lines)):
	list_point =[[i, vertices[i][0]], [len(Seeds) + i, vertices[len(Seeds) + i][0]]]
	for j in range(len(lines)):
		if lines[j][6] == i:
			list_point.append([j, vertices[j][0]])
		if lines[j][8] == i:
			list_point.append([len(Seeds)+j, vertices[j][0]])
	list_point = sorted(list_point,key=operator.itemgetter(1))
	list_points.append(list_point)
	



for list_point in list_points:
	for i in range(len(list_point)-1):
		ridge_vertices.append([list_point[i][0],list_point[i+1][0]])


from matplotlib.collections import LineCollection
ax.scatter(vertices[:,0],vertices[:,1], color = 'blue')

for simplex in ridge_vertices_easy:
        simplex = np.asarray(simplex)
        line_segments_easy.append([(x, y) for x, y in vertices[simplex]])

for simplex in ridge_vertices:
        simplex = np.asarray(simplex)
        line_segments.append([(x, y) for x, y in vertices[simplex]])

lc_easy = LineCollection(line_segments_easy,linestyle='dashed', color = 'red')
lc = LineCollection(line_segments,linestyle='solid', color = 'green')
ax.add_collection(lc)
ax.add_collection(lc_easy)

plt.show()
