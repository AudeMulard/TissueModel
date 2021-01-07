import matplotlib.pyplot as plt
import numpy as np


from shapely.geometry import LineString


def deriv(f,x):

    h = 0.000000001                 #step-size 
    return (f(x+h) - f(x))/h        #definition of derivative


def tangent_line(f,x_0,a,b):
    x = np.linspace(-1,6,200)
    y_0 = f(x_0)
    y_tan = deriv(f,x_0) * (x - x_0) + y_0 
    return y_tan



fig = plt.figure()

x = np.linspace(-1,6,200)
y = np.linspace(0,0,200)
# the function, which is y = x^2 here
def f(x):
	return (x+1)**2-3

ax = fig.add_subplot(1, 1, 1)
ax.spines['left'].set_position('zero')
ax.spines['bottom'].set_position('zero')
ax.spines['right'].set_color('none')
ax.spines['top'].set_color('none')
ax.xaxis.set_ticks_position('bottom')
#ax.yaxis.set_ticks_position('left')
ax.plot(x,tangent_line(f,2.5,-2,2),linestyle = 'dashed')
ax.plot(x,tangent_line(f,1.5,-2,2),linestyle = 'dashed')
ax.plot(x,tangent_line(f,1.0,-2,2),linestyle = 'dashed')

#### Intersections ####
first_line = LineString(np.column_stack((x,y)))
second_line = LineString(np.column_stack((x,tangent_line(f,2.5,-2,2))))
intersection = first_line.intersection(second_line)
#print intersection.xy[0],intersection.xy[1]
plt.scatter(intersection.xy[0], 'o')
plt.annotate(r'$t_1$',(intersection.xy[0][0],intersection.xy[1][0]+0.2),fontsize=15)
second_line = LineString(np.column_stack((x,tangent_line(f,1.5,-2,2))))
intersection = first_line.intersection(second_line)
#print intersection.xy[0],intersection.xy[1]
plt.scatter(intersection.xy[0], 'o')
plt.annotate(r'$t_2$',(intersection.xy[0][0]+0.05,intersection.xy[1][0]+0.2),fontsize=15)
second_line = LineString(np.column_stack((x,tangent_line(f,1.0,-2,2))))
intersection = first_line.intersection(second_line)
#print intersection.xy[0][0],intersection.xy[1]
plt.scatter(intersection.xy[0], 'o')
plt.annotate(r'$t_3$',(intersection.xy[0][0],intersection.xy[1][0]+0.2),fontsize=15)

second_line = LineString(np.column_stack((x,f(x))))
intersection = first_line.intersection(second_line)
#print intersection.xy[0][0],intersection.xy[1]
plt.scatter(intersection.xy[0], 'o')
plt.annotate(r'$t$',(intersection.xy[0][0]-0.1,intersection.xy[1][0]+0.2),fontsize=15)

plt.annotate(r'$f$',(-0.5,f(-0.5)+0.1),fontsize=15)

ax.plot(x,f(x))
plt.xlim([-1,3])
plt.ylim([-5,10])
plt.xlabel('x')     
plt.ylabel('y')     
plt.savefig('newton_graph.pdf')

plt.show()  
