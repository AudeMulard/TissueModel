import numpy as np
import matplotlib.pyplot as plt
np.set_printoptions(precision=0)
### Set boundary conditions on extoerior nodes

def new_bc(network, defo):
	new_positions = network.vertices
	for k in network.boundary_nodes_left:
	        new_positions[k,0] = new_positions[k,0] - defo
	for k in network.boundary_nodes_right:
	        new_positions[k,0] = new_positions[k,0] + defo
	for k in network.interior_nodes:
#		new_positions[k,0] = new_positions[k,0] * 1.05 # essential movement to not have a singular matrix at first step
		new_positions[k] = new_positions[k] * 1.05 # essential movement to not have a singular matrix at first step
	network.vertices=new_positions
	return network

### Calculate all the functions and jacobians

# Def of norm square
def length_square(x):
	return x[0]**2+x[1]**2


######################### LINEAR SCHEME #########################################
# Write the contribution of the node j to the equilibrium of the point i

def write_system(network):
	matrix = np.zeros((len(network.interior_nodes)*2,len(network.interior_nodes)*2))
	rest = np.zeros(len(network.interior_nodes)*2)
	Ef = network.Ef
	for i in range(len(network.interior_nodes)):
		for k in network.list_nodes_ridges[network.interior_nodes[i]]:
			matrix[2*i,2*i] += Ef
			matrix[2*i+1,2*i+1] +=Ef
			if network.interior_nodes[i]==network.ridge_vertices[k][0]:
				if network.ridge_vertices[k][1] in network.interior_nodes:
					r=network.ridge_vertices[k][1]
					j=network.interior_nodes.index(r)
					matrix[2*i,2*j]=-Ef
					matrix[2*i+1,2*j+1]=-Ef
					rest[2*i] += Ef*(network.vertices_ini[network.interior_nodes[i]][0]-network.vertices_ini[network.ridge_vertices[k][1]][0])
					rest[2*i+1] += Ef*(network.vertices_ini[network.interior_nodes[i]][1]-network.vertices_ini[network.ridge_vertices[k][1]][1])
				else:
					rest[2*i] = rest[2*i] + Ef*network.vertices[network.ridge_vertices[k][1]][0] + Ef*(network.vertices_ini[network.interior_nodes[i]][0]-network.vertices_ini[network.ridge_vertices[k][1]][0])
					rest[2*i+1] = rest[2*i+1]+ Ef*network.vertices[network.ridge_vertices[k][1]][1] + Ef*(network.vertices_ini[network.interior_nodes[i]][1]-network.vertices_ini[network.ridge_vertices[k][1]][1])
			elif network.interior_nodes[i]==network.ridge_vertices[k][1]:
				if network.ridge_vertices[k][0] in network.interior_nodes:
					r=network.ridge_vertices[k][0]
					j=network.interior_nodes.index(r)
					matrix[2*i,2*j]=-Ef
					matrix[2*i+1,2*j+1]=-Ef
					rest[2*i] += Ef*(network.vertices_ini[network.interior_nodes[i]][0]-network.vertices_ini[network.ridge_vertices[k][0]][0])
					rest[2*i+1] += Ef*(network.vertices_ini[network.interior_nodes[i]][1]-network.vertices_ini[network.ridge_vertices[k][0]][1])
				else:
					rest[2*i] = rest[2*i] + Ef*network.vertices[network.ridge_vertices[k][0]][0] + Ef*(network.vertices_ini[network.interior_nodes[i]][0]-network.vertices_ini[network.ridge_vertices[k][0]][0])
					rest[2*i+1] = rest[2*i+1]+ Ef*network.vertices[network.ridge_vertices[k][0]][1] + Ef*(network.vertices_ini[network.interior_nodes[i]][1]-network.vertices_ini[network.ridge_vertices[k][0]][1])
	return matrix, rest

def linear_scheme(network):
	matrix, rest = write_system(network)
	new_pos = np.linalg.solve(matrix,rest)
	print matrix
	for i in range(len(network.interior_nodes)):
		j=network.interior_nodes[i]
		network.vertices[j,0] = new_pos[2*i]
		network.vertices[j,1] = new_pos[2*i+1]
	return network

######################### NONLINEAR SCHEME #########################################
# calculate the whole jacobian at one point, differentiated by itself
# calculate the force at each point
def write_force(network, i, j, constitutive):
	if constitutive == 'exponential':
		lambda_f = length_square(network.vertices[i]-network.vertices[j])/length_square(network.vertices_ini[i]-network.vertices_ini[j])
		epsilon = 0.5*(lambda_f-1.)
		force=network.Ef*network.A/network.B*(np.exp(network.B*epsilon)-1.)
#		print i,j,lambda_f,force
		return force*(network.vertices[i]-network.vertices[j])
	elif constitutive == 'constant':
		return network.Ef*((network.vertices[i]-network.vertices[j])-(network.vertices_ini[i]-network.vertices_ini[j]))

def write_force_eq_point(network, i, constitutive):
	force_equation =0
	for j in network.list_nodes_ridges[network.interior_nodes[i]]:
		if network.interior_nodes[i]==network.ridge_vertices[j][0]:
			force_equation += write_force(network, network.interior_nodes[i], network.ridge_vertices[j][1], constitutive)
		if network.interior_nodes[i]==network.ridge_vertices[j][1]:
			force_equation += write_force(network, network.interior_nodes[i], network.ridge_vertices[j][0], constitutive)
	return force_equation

# Calculate jacobian for the force contributed by j on i, differentiated by x
def find_jacobian_x(network, i, j, constitutive):
	if constitutive == 'exponential':
		lambda_f = length_square(network.vertices[i]-network.vertices[j])/length_square(network.vertices_ini[i]-network.vertices_ini[j])
		jac = network.Ef*network.A/length_square(network.vertices_ini[i]-network.vertices_ini[j])*(network.vertices[i,0]-network.vertices[j,0])*np.exp(network.B*0.5*(lambda_f-1.))
		force = network.Ef*network.A/network.B*(np.exp(network.B*0.5*(lambda_f-1.))-1.)
		return jac*(network.vertices[i]-network.vertices[j]) + force*np.array([1,0])
#		return np.array([1,0])
	elif constitutive == 'constant':
		return network.Ef*np.array([1,0])

# Calculate jacobian for the force contributed by j on i, differentiated by y
def find_jacobian_y(network, i, j, constitutive):
	if constitutive == 'exponential':
		lambda_f = length_square(network.vertices[i]-network.vertices[j])/length_square(network.vertices_ini[i]-network.vertices_ini[j])
		jac = network.Ef*network.A/length_square(network.vertices_ini[i]-network.vertices_ini[j])*(network.vertices[i,1]-network.vertices[j,1])*np.exp(network.B*0.5*(lambda_f-1.))
		force = network.Ef*network.A/network.B*(np.exp(network.B*0.5*(lambda_f-1.))-1.)
		return jac*(network.vertices[i]-network.vertices[j]) + force*np.array([0,1])
#		return np.array([0,1])
	elif constitutive == 'constant':
		return network.Ef*np.array([0,1])

def write_jacobian_point(network, i, constitutive):
	jacobian_x=0
	jacobian_y=0
	for j in network.list_nodes_ridges[network.interior_nodes[i]]:
		if network.interior_nodes[i]==network.ridge_vertices[j][0]:
			jacobian_x += find_jacobian_x(network, network.interior_nodes[i], network.ridge_vertices[j][1], constitutive)
			jacobian_y += find_jacobian_y(network, network.interior_nodes[i], network.ridge_vertices[j][1], constitutive)
		if network.interior_nodes[i]==network.ridge_vertices[j][1]:
			jacobian_x += find_jacobian_x(network, network.interior_nodes[i], network.ridge_vertices[j][0], constitutive)
			jacobian_y += find_jacobian_y(network, network.interior_nodes[i], network.ridge_vertices[j][0], constitutive)
#		print network.interior_nodes[i], j, jacobian_x
	return np.array([jacobian_x, jacobian_y])

### Iterative Newton scheme to solve the nonlinear equations

# Calculate vector F

def write_vector_F(network, constitutive):
	F = np.zeros(len(network.interior_nodes)*2)
	for i in range(len(network.interior_nodes)):
		F[2*i] = np.array(write_force_eq_point(network,i, constitutive))[0]
		F[2*i+1] = np.array(write_force_eq_point(network,i, constitutive))[1]
	return F

# Calculate matrix J

def write_matrix_J(network, constitutive):
	J = np.zeros((len(network.interior_nodes)*2,len(network.interior_nodes)*2))
	for i in range(len(network.interior_nodes)):
		J[2*i,2*i]=write_jacobian_point(network, i, constitutive)[0,0]
		J[2*i,2*i+1]=write_jacobian_point(network, i, constitutive)[1,0]
		J[2*i+1, 2*i] = write_jacobian_point(network,i, constitutive)[0,1]
		J[2*i+1,2*i+1]=write_jacobian_point(network, i, constitutive)[1,1]
#		print i, write_jacobian_point(network, i, constitutive)
		for k in network.list_nodes_ridges[network.interior_nodes[i]]:
			if network.interior_nodes[i]==network.ridge_vertices[k][0] and network.ridge_vertices[k][1] in network.interior_nodes:
				r=network.ridge_vertices[k][1]
				j=network.interior_nodes.index(r)
				J[2*i,2*j]=-find_jacobian_x(network, network.interior_nodes[i], r, constitutive)[0]
				J[2*i,2*j+1]=-find_jacobian_y(network, network.interior_nodes[i], r, constitutive)[0]
				J[2*i+1,2*j]=-find_jacobian_x(network, network.interior_nodes[i], r, constitutive)[1]
				J[2*i+1,2*j+1]=-find_jacobian_y(network, network.interior_nodes[i], r, constitutive)[1]
			if network.interior_nodes[i]==network.ridge_vertices[k][1] and network.ridge_vertices[k][0] in network.interior_nodes:
				r=network.ridge_vertices[k][0]
				j=network.interior_nodes.index(r)
				J[2*i,2*j]=-find_jacobian_x(network, network.interior_nodes[i], r, constitutive)[0]
				J[2*i,2*j+1]=-find_jacobian_y(network, network.interior_nodes[i], r, constitutive)[0]
				J[2*i+1,2*j]=-find_jacobian_x(network, network.interior_nodes[i], r, constitutive)[1]
				J[2*i+1,2*j+1]=-find_jacobian_y(network, network.interior_nodes[i], r, constitutive)[1]
			#print i, k,J
	return J


def iterative_newton(network, constitutive):
	max_iter = 50
	epsilon = 1e-8
	for k in range(max_iter):
		F = write_vector_F(network, constitutive)
		J = write_matrix_J(network, constitutive)
		print J
		diff = np.linalg.solve(J,-F)
		for i in range(len(network.interior_nodes)):
			j=network.interior_nodes[i]
			network.vertices[j,0] = network.vertices[j,0] + diff[2*i]
			network.vertices[j,1] = network.vertices[j,1] + diff[2*i+1]
		if np.linalg.norm(diff) < epsilon:
			print('convergence!, nre iter:', k)
			break
		else:
			print(k, 'not converged')
	return network


### Global function of solving one step if the netowrk

def solve_force_balance(network, defo, constitutive, scheme):
	network = new_bc(network, defo)
	if scheme == 'nonlinear':
		network = iterative_newton(network, constitutive)
	if scheme == 'linear':
		network = linear_scheme(network)
	return network

