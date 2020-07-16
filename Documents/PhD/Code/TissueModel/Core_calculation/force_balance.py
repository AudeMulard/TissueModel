import numpy as np
import matplotlib.pyplot as plt
#np.set_printoptions(precision=0)
import sympy as sp
import scipy
from scipy.sparse import csc_matrix
from scipy.sparse.linalg import spsolve, isolve,gmres,cgs,qmr
import warnings
import time
from sympy import Matrix
from sympy import *
from scipy.sparse.linalg.dsolve import linsolve
from scipy.optimize import newton_krylov,anderson,fsolve,root, newton

warnings.filterwarnings("error", category=RuntimeWarning)

#with warnings.catch_warnings():
#    warnings.simplefilter("error")

### Set boundary conditions on exterior nodes, chosse what side you are pulling from, then add the displacements to the corresponding points.

def new_bc(network, defo, side):
	new_positions = network.vertices
	if side == 'left':
		for k in network.boundary_nodes_left:
		        new_positions[k,0] = new_positions[k,0] - defo
	if side == 'right':
		for k in network.boundary_nodes_right:
		        new_positions[k,0] = new_positions[k,0] + defo
	if side == 'both':
		for k in network.boundary_nodes_left:
		        new_positions[k,0] = new_positions[k,0] - defo
		for k in network.boundary_nodes_right:
		        new_positions[k,0] = new_positions[k,0] + defo
#	for k in network.interior_nodes:
#		new_positions[k] = new_positions[k] * 1.05 # essential movement to not have a singular matrix at first step
	network.vertices=new_positions
	return network

### Calculate all the functions and jacobians
import numba

# Def of the square of the norm
#@numba.autojit(nopython=True)
def length_square(dimension,x):	
	if dimension == 2:
		return x[0]**2+x[1]**2
	elif dimension == 3:
		return x[0]**2+x[1]**2+x[2]**2

######################### LINEAR SCHEME 3D #########################################
# Write the contribution of the node j to the equilibrium of the point i

def write_system(network):
	matrix = np.zeros((len(network.interior_nodes)*network.dimension,len(network.interior_nodes)*network.dimension))
	rest = np.zeros(len(network.interior_nodes)*network.dimension)
	Ef = network.Ef
	if network.dimension == 2:
		coords = [0,1]
	elif network.dimension == 3:
		coords = [0,1,2]
	for i in range(len(network.interior_nodes)):
		for k in network.list_nodes_ridges[network.interior_nodes[i]]:
			for coord in coords:
				matrix[network.dimension*i+coord,network.dimension*i+coord] += Ef
			if network.interior_nodes[i]==network.ridge_vertices[k][0]:
				if network.ridge_vertices[k][1] in network.interior_nodes:
					r=network.ridge_vertices[k][1]
					j=network.interior_nodes.index(r)
					for coord in coords:
						matrix[network.dimension*i+coord,network.dimension*j+coord]=-Ef
						rest[network.dimension*i+coord] += Ef*(network.vertices_ini[network.interior_nodes[i]][coord]-network.vertices_ini[network.ridge_vertices[k][1]][coord])
				else:
					for coord in coords:
						rest[network.dimension*i+coord] = rest[network.dimension*i+coord] + Ef*network.vertices[network.ridge_vertices[k][1]][coord] + Ef*(network.vertices_ini[network.interior_nodes[i]][coord]-network.vertices_ini[network.ridge_vertices[k][1]][coord])
			elif network.interior_nodes[i]==network.ridge_vertices[k][1]:
				if network.ridge_vertices[k][0] in network.interior_nodes:
					r=network.ridge_vertices[k][0]
					j=network.interior_nodes.index(r)
					for coord in coords:
						matrix[network.dimension*i+coord,network.dimension*j+coord]=-Ef
						rest[network.dimension*i+coord] += Ef*(network.vertices_ini[network.interior_nodes[i]][coord]-network.vertices_ini[network.ridge_vertices[k][0]][coord])
				else:
					for coord in coords:
						rest[network.dimension*i+coord] = rest[network.dimension*i+coord] + Ef*network.vertices[network.ridge_vertices[k][0]][coord] + Ef*(network.vertices_ini[network.interior_nodes[i]][coord]-network.vertices_ini[network.ridge_vertices[k][0]][coord])
	return matrix, rest

def linear_scheme(network):
	matrix, rest = write_system(network)
	new_pos = np.linalg.solve(matrix,rest)
	if network.dimension == 2:
		coords = [0,1]
	elif network.dimension == 3:
		coords = [0,1,2]
	for i in range(len(network.interior_nodes)):
		for coord in coords:
			j=network.interior_nodes[i]
			network.vertices[j,coord] = new_pos[network.dimension*i+coord]
	return network


######################### NONLINEAR SCHEME #########################################
# calculate the whole jacobian at one point, differentiated by itself
# calculate the force at each point

# Calculate jacobian for the force contributed by j on i, differentiated by x
#@numba.autojit()
def find_jacobian_x(network, i, j, constitutive):
	if constitutive == 'exponential':
		lambda_f = length_square(network.dimension,network.vertices[i]-network.vertices[j])/length_square(network.dimension,network.vertices_ini[i]-network.vertices_ini[j])
		jac = network.Ef*network.A/length_square(network.dimension,network.vertices_ini[i]-network.vertices_ini[j])*(network.vertices[i,0]-network.vertices[j,0])*np.exp(network.B*0.5*(lambda_f-1.))
		force = network.Ef*network.A/network.B*(np.exp(network.B*0.5*(lambda_f-1.))-1.)
		return jac*(network.vertices[i]-network.vertices[j]) + force*np.array([1,0])
	elif constitutive == 'constant':
		return network.Ef*np.array([1,0])
	elif constitutive == 'spring':
		"""try:
			ridge = list(network.ridge_vertices).index([i,j])
		except:
			ridge = list(network.ridge_vertices).index([j,i])
		if network.state_ridge[ridge] == 'tension':
			Ef = network.k_tension
		else:
			Ef = network.k_compression"""
		Ef = network.k_compression
		length = sp.sqrt(length_square(network.dimension,network.vertices[i]-network.vertices[j]))
		diff_dist = Ef*(length-sp.sqrt(length_square(network.dimension,network.vertices_ini[i]-network.vertices_ini[j])))
		degree_3_x = -Ef*diff_dist*(network.vertices[i][0]-network.vertices[j][0])**2/length**3
		degree_2_x = Ef*(network.vertices[i][0]-network.vertices[j][0])**2/length**2
		degree_1 = Ef*diff_dist/length
		degree_3_y = -Ef*diff_dist*(network.vertices[i][0]-network.vertices[j][0])*(network.vertices[i][1]-network.vertices[j][1])/length**3
		degree_2_y = Ef*(network.vertices[i][0]-network.vertices[j][0])*(network.vertices[i][1]-network.vertices[j][1])/length**2
		return -degree_3_x-degree_2_x-degree_1,-degree_3_y-degree_2_y

# Calculate jacobian for the force contributed by j on i, differentiated by y
#@numba.autojit()
def find_jacobian_y(network, i, j, constitutive):
	if constitutive == 'exponential':
		lambda_f = length_square(network.dimension,network.vertices[i]-network.vertices[j])/length_square(network.dimension,network.vertices_ini[i]-network.vertices_ini[j])
		jac = network.Ef*network.A/length_square(network.dimension,network.vertices_ini[i]-network.vertices_ini[j])*(network.vertices[i,1]-network.vertices[j,1])*np.exp(network.B*0.5*(lambda_f-1.))
		force = network.Ef*network.A/network.B*(np.exp(network.B*0.5*(lambda_f-1.))-1.)
		return jac*(network.vertices[i]-network.vertices[j]) + force*np.array([0,1])
	elif constitutive == 'constant':
		return network.Ef*np.array([0,1])
	elif constitutive == 'spring':
		"""try:
			ridge = list(network.ridge_vertices).index([i,j])
		except:
			ridge = list(network.ridge_vertices).index([j,i])
		if network.state_ridge[ridge] == 'tension':
			Ef = network.k_tension
		else:
			Ef = network.k_compression"""
		Ef = network.k_compression
		length = sp.sqrt(length_square(network.dimension,network.vertices[i]-network.vertices[j]))
		diff_dist = Ef*(length-sp.sqrt(length_square(network.dimension,network.vertices_ini[i]-network.vertices_ini[j])))
		degree_3_x = -Ef*diff_dist*(network.vertices[i][0]-network.vertices[j][0])*(network.vertices[i][1]-network.vertices[j][1])/length**3
		degree_2_x = Ef*(network.vertices[i][0]-network.vertices[j][0])*(network.vertices[i][1]-network.vertices[j][1])/length**2
		degree_1 = Ef*diff_dist/length
		degree_3_y = -Ef*diff_dist*(network.vertices[i][1]-network.vertices[j][1])**2/length**3
		degree_2_y = Ef*(network.vertices[i][1]-network.vertices[j][1])**2/length**2
		return -degree_3_x-degree_2_x,-degree_3_y-degree_2_y-degree_1

## Calculate the jacobian at the position i, where it is needed to differentiate the whole sum of the force

def write_jacobian_point(network, i, constitutive):
	jacobian_x=0
	jacobian_y=0
	for j in network.list_nodes_ridges[network.interior_nodes[i]]:
		if network.interior_nodes[i]==network.ridge_vertices[j][0]:
			jacobian_x -= np.array(find_jacobian_x(network, network.interior_nodes[i], network.ridge_vertices[j][1], constitutive))
			jacobian_y -= np.array(find_jacobian_y(network, network.interior_nodes[i], network.ridge_vertices[j][1], constitutive))
		if network.interior_nodes[i]==network.ridge_vertices[j][1]:
			jacobian_x -= np.array(find_jacobian_x(network, network.interior_nodes[i], network.ridge_vertices[j][0], constitutive))
			jacobian_y -= np.array(find_jacobian_y(network, network.interior_nodes[i], network.ridge_vertices[j][0], constitutive))
	return np.array([jacobian_x, jacobian_y])

### Iterative Newton scheme to solve the nonlinear equations

# Calculate vector F

def write_vector_F(network, constitutive):
	F = np.zeros(len(network.interior_nodes)*2)
	for i in range(len(network.interior_nodes)):
		F[2*i] = np.array(write_force_eq_point(network,i, constitutive))[0]
		F[2*i+1] = np.array(write_force_eq_point(network,i, constitutive))[1]
	return F

# Writes the matrix J with the jacobian values calculated by the previous functions and placed in the right spots

def write_matrix_J(network, constitutive):
    J = np.zeros((len(network.interior_nodes)*2,len(network.interior_nodes)*2))
    for i in range(len(network.interior_nodes)):
        J[2*i,2*i]=write_jacobian_point(network, i, constitutive)[0,0]
        J[2*i,2*i+1]=write_jacobian_point(network, i, constitutive)[1,0]
        J[2*i+1, 2*i] = write_jacobian_point(network,i, constitutive)[0,1]
        J[2*i+1,2*i+1]=write_jacobian_point(network, i, constitutive)[1,1]
        for k in network.list_nodes_ridges[network.interior_nodes[i]]:
            if network.interior_nodes[i]==network.ridge_vertices[k][0] and network.ridge_vertices[k][1] in network.interior_nodes:
                r=network.ridge_vertices[k][1]
                j=network.interior_nodes.index(r)
                J[2*i,2*j],J[2*i+1,2*j]=find_jacobian_x(network, network.interior_nodes[i], r, constitutive)
                J[2*i,2*j+1],J[2*i+1,2*j+1]=find_jacobian_y(network, network.interior_nodes[i], r, constitutive)
            elif network.interior_nodes[i]==network.ridge_vertices[k][1] and network.ridge_vertices[k][0] in network.interior_nodes:
                r=network.ridge_vertices[k][0]
                j=network.interior_nodes.index(r)
                J[2*i,2*j],J[2*i+1,2*j]=find_jacobian_x(network, network.interior_nodes[i], r, constitutive)
                J[2*i,2*j+1],J[2*i+1,2*j+1]=find_jacobian_y(network, network.interior_nodes[i], r, constitutive)	
    return J

## Iterative scheme of the Newton method

def iterative_newton(network, constitutive,details):
	from scipy.optimize.nonlin import InverseJacobian
	from scipy.optimize.nonlin import BroydenFirst, KrylovJacobian
	max_iter = 500
	epsilon = 1e-8
	vertices=[]
	global vertices_ini
	global list_nodes_ridges
	global interior_nodes
	global ridge_vertices
	ridge_vertices = network.ridge_vertices
	interior_nodes = network.interior_nodes
	list_nodes_ridges=network.list_nodes_ridges
	vertices_ini=network.vertices_ini
	global vertices_all
	vertices_all = network.vertices
	#print vertices_all, vertices_ini
	for i in range(len(network.interior_nodes)):
		j = network.interior_nodes[i]
		vertices.append(vertices_ini[j][0])
		vertices.append(vertices_ini[j][1])
	start = time.time()
	jac=BroydenFirst()
	sol = newton(create_F, vertices,tol=epsilon,maxiter=500)#,method='df-sane')#,verbose=True,f_tol=epsilon, inner_M=KrylovJacobian(inner_M=InverseJacobian(jac)),method='gmres')
	print sol
	vertices_sol = np.array(network.vertices)
	for i in range(len(network.interior_nodes)):
		j = network.interior_nodes[i]
		vertices_sol[j]=[sol[2*i],sol[2*i+1]]
	print 'krylov_time', time.time()-start
	middle = time.time()
	for k in range(max_iter):
		F = write_vector_F(network, constitutive)
		#print F
		J = write_matrix_J(network, constitutive)
		#if np.linalg.cond(J) > 10e10:
		#	print('STOP SIMULATION: ERROR')
		#print np.linalg.cond(J)
		#conversion in sparse
		J_sparse = csc_matrix(J)
		"""#diff = spsolve(J_sparse,-F)#, use_umfpack=True)
		diff = isolve.bicg(J_sparse,-F)
		print 'isolve', time.time()-start
		middle = time.time()
		#diff_h = linsolve.spsolve(J_sparse,-F, use_umfpack=True)
		sol = newton_krylov(create_F,vertices)
		print 'krylov', time.time()-middle
		middle = time.time()
		diff_h = gmres(J_sparse,-F)
		print 'gmres', time.time() - middle
		#for i in range(len(diff)):
		#	if diff[i]-diff_h[i]>10e-10: print(diff[i], diff_h[i], i)
		middle = time.time()
		"""
		diff= spsolve(J_sparse, -F)
		"""
		print 'spsolve',time.time() - middle
		for i in range(len(diff)):
			if diff[i]-diff_h[i]>10e-10: print(diff[i], diff_h[i], i)
		middle = time.time()
		diff_h = qmr(J_sparse, -F)
		print time.time() - middle
		for i in range(len(diff)):
			if diff[i]-diff_h[i]>10e-10: print(diff[i], diff_h[i], i)
		middle = time.time()
		diff_h = cgs(J_sparse, -F)
		print time.time() - middle
		for i in range(len(diff)):
			if diff[i]-diff_h[i]>10e-10: print(diff[i], diff_h[i], i)
		"""
		#diff = scipy.linalg.solve(J,-F)
		#print(len(diff_h),len(diff))
		#print(np.allclose(diff,diff_h))
		#print diff_hey, 'not', np.linalg.norm(diff_hey)
		#print F, J, diff
		for i in range(len(network.interior_nodes)):
			j=network.interior_nodes[i]
			network.vertices[j,0] = network.vertices[j,0] + diff[2*i]
			network.vertices[j,1] = network.vertices[j,1] + diff[2*i+1]
		if details == True: # This is just information on what is happening during the simulation
			if k%25==0:
				print('Step: ', k, ', convergence: ',np.linalg.norm(diff))
		if np.linalg.norm(diff) < epsilon:
			if details == True:
				print(['convergence!, nre iter:', k])
			break
		elif k== max_iter-1:
			raise ValueError
	print 'mymethod', time.time()-middle
	for i in range(len(network.vertices)):
		if np.linalg.norm(network.vertices[i]-vertices_sol[i])>10e-10: print i, network.vertices[i],vertices_sol[i], network.vertices_ini[i]
	print vertices_sol
	return network



######################### NONLINEAR SCHEME 3D #########################################


## calculate the force at applied on node i, by the ridge linked to j.

def create_F(array):
	result = []
	for i in range(len(interior_nodes)):
		force_equation_x = 0.0
		force_equation_y = 0.0
		#print i, interior_nodes[i],list_nodes_ridges[interior_nodes[i]]
		for j in list_nodes_ridges[interior_nodes[i]]:
			if interior_nodes[i] == ridge_vertices[j][0]:
				k = interior_nodes[i]
				j = ridge_vertices[j][1]
			elif interior_nodes[i] == ridge_vertices[j][1]:
				k = interior_nodes[i]
				j = ridge_vertices[j][0]
				#print k,j
			force_equation_x += (float(sp.sqrt((array[2*i]-vertices_all[j][0])**2+(array[2*i+1]-vertices_all[j][1])**2))
								 -float(np.linalg.norm(vertices_ini[k]-vertices_ini[j])))*(array[2*i]-vertices_all[j][0])\
								 /float(sp.sqrt((array[2*i]-vertices_all[j][0])**2+(array[2*i+1]-vertices_all[j][1])**2))
			force_equation_y+= (float(sp.sqrt((array[2 * i] - vertices_all[j][0]) ** 2 + (array[2 * i + 1] - vertices_all[j][1]) ** 2))
								- float(np.linalg.norm(vertices_ini[k] - vertices_ini[j]))) * (array[2 * i + 1] -vertices_all[j][1])\
							    / float(sp.sqrt((array[2 * i] - vertices_all[j][0]) ** 2 + (array[2 * i + 1] - vertices_all[j][1]) ** 2))
			#print i,'force_eq', force_equation_x
		result.append(force_equation_x)
		result.append(force_equation_y)
	return result

def write_force(network, i, j, constitutive):
	if constitutive == 'spring':
		if network.dimension ==2:
			try:			
				ridge = list(network.ridge_vertices).index([i,j])
			except:
				ridge = list(network.ridge_vertices).index([j,i])
			if network.state_ridge[ridge] == 'tension':
				Ef = network.k_tension
			else:
				Ef = network.k_compression
		elif network.dimension == 3:
			Ef = network.k_tension
		return Ef*(float(sp.sqrt(length_square(network.dimension,network.vertices[i]-network.vertices[j])))-float(sp.sqrt(length_square(network.dimension,network.vertices_ini[i]-network.vertices_ini[j]))))*(network.vertices[i]-network.vertices[j])/float(sp.sqrt(length_square(network.dimension,network.vertices[i]-network.vertices[j])))

## Sums up all the forces applied on the node i.

def write_force_eq_point(network, i, constitutive):
	force_equation =0.
	for j in network.list_nodes_ridges[network.interior_nodes[i]]:
		if network.interior_nodes[i]==network.ridge_vertices[j][0]:
			force_equation += write_force(network, network.interior_nodes[i], network.ridge_vertices[j][1], constitutive)
		if network.interior_nodes[i]==network.ridge_vertices[j][1]:
			force_equation += write_force(network, network.interior_nodes[i], network.ridge_vertices[j][0], constitutive)
	return force_equation

# Calculate jacobian for the force contributed by j on i, differentiated by x
def find_jacobian_x_3d(network, i, j, constitutive):
	if constitutive == 'spring':
		length = sp.sqrt(length_square(network.dimension,network.vertices[i]-network.vertices[j]))
		jac = network.k_tension*(network.vertices[i][0]-network.vertices[j][0])/length
		force = network.k_tension*(1-sp.sqrt(length_square(network.dimension,network.vertices_ini[i]-network.vertices_ini[j]))/sp.sqrt(length_square(network.dimension,network.vertices[i]-network.vertices[j])))
		jac_x = jac*(network.vertices[i][0]-network.vertices[j][0])/length+force
		jac_y = jac*(network.vertices[i][1]-network.vertices[j][1])/length
		jac_z = jac*(network.vertices[i][2]-network.vertices[j][2])/length
		return -jac_x, -jac_y, -jac_z

# Calculate jacobian for the force contributed by j on i, differentiated by y
def find_jacobian_y_3d(network, i, j, constitutive):
	if constitutive == 'spring':
		length = sp.sqrt(length_square(network.dimension,network.vertices[i]-network.vertices[j]))
		jac = network.k_tension*(network.vertices[i][1]-network.vertices[j][1])/length
		force = network.k_tension*(1-sp.sqrt(length_square(network.dimension,network.vertices_ini[i]-network.vertices_ini[j]))/sp.sqrt(length_square(network.dimension,network.vertices[i]-network.vertices[j])))
		jac_x = jac*(network.vertices[i][0]-network.vertices[j][0])/length
		jac_y = jac*(network.vertices[i][1]-network.vertices[j][1])/length+force
		jac_z = jac*(network.vertices[i][2]-network.vertices[j][2])/length
		return -jac_x, -jac_y, -jac_z

# Calculate jacobian for the force contributed by j on i, differentiated by z
def find_jacobian_z_3d(network, i, j, constitutive):
	if constitutive == 'spring':
		length = sp.sqrt(length_square(network.dimension,network.vertices[i]-network.vertices[j]))
		jac = network.k_tension*(network.vertices[i][2]-network.vertices[j][2])/length
		force = network.k_tension*(1-sp.sqrt(length_square(network.dimension,network.vertices_ini[i]-network.vertices_ini[j]))/sp.sqrt(length_square(network.dimension,network.vertices[i]-network.vertices[j])))
		jac_x = jac*(network.vertices[i][0]-network.vertices[j][0])/length
		jac_y = jac*(network.vertices[i][1]-network.vertices[j][1])/length
		jac_z = jac*(network.vertices[i][2]-network.vertices[j][2])/length+force
		return -jac_x, -jac_y, -jac_z

## Calculate the jacobian at the node i, where it is needed to differentiate the whole sum of the force
def write_jacobian_point_3d(network, i, constitutive):
	jacobian_x=0
	jacobian_y=0
	jacobian_z=0
	for j in network.list_nodes_ridges[network.interior_nodes[i]]:
		if network.interior_nodes[i]==network.ridge_vertices[j][0]:
			jacobian_x -= np.array(find_jacobian_x_3d(network, network.interior_nodes[i], network.ridge_vertices[j][1], constitutive))
			jacobian_y -= np.array(find_jacobian_y_3d(network, network.interior_nodes[i], network.ridge_vertices[j][1], constitutive))
			jacobian_z -= np.array(find_jacobian_z_3d(network, network.interior_nodes[i], network.ridge_vertices[j][1], constitutive))
		if network.interior_nodes[i]==network.ridge_vertices[j][1]:
			jacobian_x -= np.array(find_jacobian_x_3d(network, network.interior_nodes[i], network.ridge_vertices[j][0], constitutive))
			jacobian_y -= np.array(find_jacobian_y_3d(network, network.interior_nodes[i], network.ridge_vertices[j][0], constitutive))
			jacobian_z -= np.array(find_jacobian_z_3d(network, network.interior_nodes[i], network.ridge_vertices[j][0], constitutive))
	return np.array([jacobian_x, jacobian_y, jacobian_z])

### Iterative Newton scheme to solve the nonlinear equations

# Calculate vector F

def write_vector_F_3d(network, constitutive):
	F = np.zeros(len(network.interior_nodes)*network.dimension)
	if network.dimension == 2:
		coords = [0,1]
	elif network.dimension == 3:
		coords = [0,1,2]
	for i in range(len(network.interior_nodes)):
		force = write_force_eq_point(network,i, constitutive)
		F[network.dimension*i],F[network.dimension*i+1],F[network.dimension*i+2] = force[0],force[1],force[2]
	return F

# Calculate matrix J

def write_matrix_J_3d(network, constitutive):
	J = np.zeros((len(network.interior_nodes)*network.dimension,len(network.interior_nodes)*network.dimension))
	if network.dimension == 2:
		coords = [0,1]
	elif network.dimension == 3:
		coords = [0,1,2]
	for i in range(len(network.interior_nodes)):
		for coord1 in coords:
			for coord2 in coords:
				J[network.dimension*i+coord1,network.dimension*i+coord2]=write_jacobian_point_3d(network, i, constitutive)[coord2,coord1]
		for k in network.list_nodes_ridges[network.interior_nodes[i]]:
			if network.interior_nodes[i]==network.ridge_vertices[k][0] and network.ridge_vertices[k][1] in network.interior_nodes:
				r=network.ridge_vertices[k][1]
				j=network.interior_nodes.index(r)
				J[network.dimension*i,network.dimension*j],J[network.dimension*i+1,network.dimension*j],J[network.dimension*i+2,network.dimension*j]=find_jacobian_x_3d(network, network.interior_nodes[i], r, constitutive)
				J[network.dimension*i,network.dimension*j+1],J[network.dimension*i+1,network.dimension*j+1],J[network.dimension*i+2,network.dimension*j+1]=find_jacobian_y_3d(network, network.interior_nodes[i], r, constitutive)
				J[network.dimension*i,network.dimension*j+2],J[network.dimension*i+1,network.dimension*j+2],J[network.dimension*i+2,network.dimension*j+2]=find_jacobian_z_3d(network, network.interior_nodes[i], r, constitutive)
			elif network.interior_nodes[i]==network.ridge_vertices[k][1] and network.ridge_vertices[k][0] in network.interior_nodes:
				r=network.ridge_vertices[k][0]
				j=network.interior_nodes.index(r)
				J[network.dimension*i,network.dimension*j],J[network.dimension*i+1,network.dimension*j],J[network.dimension*i+2,network.dimension*j]=find_jacobian_x_3d(network, network.interior_nodes[i], r, constitutive)
				J[network.dimension*i,network.dimension*j+1],J[network.dimension*i+1,network.dimension*j+1],J[network.dimension*i+2,network.dimension*j+1]=find_jacobian_y_3d(network, network.interior_nodes[i], r, constitutive)
				J[network.dimension*i,network.dimension*j+2],J[network.dimension*i+1,network.dimension*j+2],J[network.dimension*i+2,network.dimension*j+2]=find_jacobian_z_3d(network, network.interior_nodes[i], r, constitutive)
	return J


def iterative_newton_3d(network, constitutive,details):
	max_iter = 200
	epsilon = 1e-8
	for k in range(max_iter):
		F = write_vector_F_3d(network, constitutive)
		J = write_matrix_J_3d(network, constitutive)
		#print np.linalg.cond(J)
		#conversion in sparse
		J_sparse = csc_matrix(J)
		#start = time.time()
		diff = spsolve(J_sparse,-F)
		#print time.time()-start
		for i in range(len(network.interior_nodes)):
			j=network.interior_nodes[i]
			network.vertices[j,0] = network.vertices[j,0] + diff[network.dimension*i]
			network.vertices[j,1] = network.vertices[j,1] + diff[network.dimension*i+1]
			network.vertices[j,2] = network.vertices[j,2] + diff[network.dimension*i+2]
		if details == True:
			if k%5==0:
				print('Step: ', k, ', convergence: ',np.linalg.norm(diff))
		if np.linalg.norm(diff) < epsilon:
			print(['convergence!, nre iter:', k])
			break
		elif k== max_iter-1:
			raise ValueError 
#		else:
#			print(k, 'not converged')
	return network



### Global function of solving one step of the network tensile test depending on the dimension

def solve_force_balance(network, constitutive, details):
	if constitutive == 'linear':
		network = linear_scheme(network)
	else:
		if network.dimension == 2:
			network = iterative_newton(network, constitutive,details)
		elif network.dimension == 3:
			network = iterative_newton_3d(network, constitutive,details)
	return network

