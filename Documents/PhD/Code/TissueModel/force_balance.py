import numpy as np
import matplotlib.pyplot as plt
#np.set_printoptions(precision=0)


### Set boundary conditions on exterior nodes

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

# Def of norm square
def length_square(network,x):	
	if network.dimension == 2:
		return x[0]**2+x[1]**2
	elif network.dimension == 3:
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

"""
######################### NONLINEAR SCHEME #########################################
# calculate the whole jacobian at one point, differentiated by itself
# calculate the force at each point
def write_force(network, i, j, constitutive):
	if constitutive == 'exponential':
		lambda_f = length_square(network,network.vertices[i]-network.vertices[j])/length_square(network,network.vertices_ini[i]-network.vertices_ini[j])
		epsilon = 0.5*(lambda_f-1.)
		force=network.Ef*network.A/network.B*(np.exp(network.B*epsilon)-1.)
		return force*(network.vertices[i]-network.vertices[j])
	elif constitutive == 'constant':
		return network.Ef*((network.vertices[i]-network.vertices[j])-(network.vertices_ini[i]-network.vertices_ini[j]))
	elif constitutive == 'linear':
		return network.Ef*(1-np.sqrt(length_square(network,network.vertices_ini[i]-network.vertices_ini[j]))/np.sqrt(length_square(network,network.vertices[i]-network.vertices[j])))*(network.vertices[i]-network.vertices[j])
	elif constitutive == 'linear2':
		return network.Ef*(np.sqrt(length_square(network,network.vertices[i]-network.vertices[j]))-np.sqrt(length_square(network,network.vertices_ini[i]-network.vertices_ini[j])))*(network.vertices[i]-network.vertices[j])/np.sqrt(length_square(network,network.vertices[i]-network.vertices[j]))

def write_force_eq_point(network, i, constitutive):
	force_equation =0
	for j in network.list_nodes_ridges[network.interior_nodes[i]]:
		if network.interior_nodes[i]==network.ridge_vertices[j][0]:
			force_equation += write_force(network, network.interior_nodes[i], network.ridge_vertices[j][1], constitutive)
		if network.interior_nodes[i]==network.ridge_vertices[j][1]:
			force_equation += write_force(network, network.interior_nodes[i], network.ridge_vertices[j][0], constitutive)
	return force_equation
"""

# Calculate jacobian for the force contributed by j on i, differentiated by x
def find_jacobian_x(network, i, j, constitutive):
	if constitutive == 'exponential':
		lambda_f = length_square(network,network.vertices[i]-network.vertices[j])/length_square(network,network.vertices_ini[i]-network.vertices_ini[j])
		jac = network.Ef*network.A/length_square(network,network.vertices_ini[i]-network.vertices_ini[j])*(network.vertices[i,0]-network.vertices[j,0])*np.exp(network.B*0.5*(lambda_f-1.))
		force = network.Ef*network.A/network.B*(np.exp(network.B*0.5*(lambda_f-1.))-1.)
		return jac*(network.vertices[i]-network.vertices[j]) + force*np.array([1,0])
	elif constitutive == 'constant':
		return network.Ef*np.array([1,0])
	elif constitutive == 'linear':
		length = np.sqrt(length_square(network,network.vertices[i]-network.vertices[j]))
		jac = network.Ef*np.sqrt(length_square(network,network.vertices_ini[i]-network.vertices_ini[j]))*(network.vertices[i][0]-network.vertices[j][0])/length**3
		force = network.Ef*(1-np.sqrt(length_square(network,network.vertices_ini[i]-network.vertices_ini[j]))/np.sqrt(length_square(network,network.vertices[i]-network.vertices[j])))
		return jac*(network.vertices[i]-network.vertices[j]) + force*np.array([1,0])
	elif constitutive == 'linear2':
		length = np.sqrt(length_square(network,network.vertices[i]-network.vertices[j]))
		jac = network.Ef*(network.vertices[i][0]-network.vertices[j][0])/length
		force = network.Ef*(1-np.sqrt(length_square(network,network.vertices_ini[i]-network.vertices_ini[j]))/np.sqrt(length_square(network,network.vertices[i]-network.vertices[j])))
		return jac*(network.vertices[i]-network.vertices[j])/np.sqrt(length_square(network,network.vertices[i]-network.vertices[j])) + force*np.array([1,0])

# Calculate jacobian for the force contributed by j on i, differentiated by y
def find_jacobian_y(network, i, j, constitutive):
	if constitutive == 'exponential':
		lambda_f = length_square(network,network.vertices[i]-network.vertices[j])/length_square(network,network.vertices_ini[i]-network.vertices_ini[j])
		jac = network.Ef*network.A/length_square(network,network.vertices_ini[i]-network.vertices_ini[j])*(network.vertices[i,1]-network.vertices[j,1])*np.exp(network.B*0.5*(lambda_f-1.))
		force = network.Ef*network.A/network.B*(np.exp(network.B*0.5*(lambda_f-1.))-1.)
		return jac*(network.vertices[i]-network.vertices[j]) + force*np.array([0,1])
#		return np.array([0,1])
	elif constitutive == 'constant':
		return network.Ef*np.array([0,1])
	elif constitutive == 'linear':
		length = np.sqrt(length_square(network,network.vertices[i]-network.vertices[j]))
		jac = network.Ef*np.sqrt(length_square(network,network.vertices_ini[i]-network.vertices_ini[j]))*(network.vertices[i][1]-network.vertices[j][1])/length**3
		force = network.Ef*(1-np.sqrt(length_square(network,network.vertices_ini[i]-network.vertices_ini[j]))/np.sqrt(length_square(network,network.vertices[i]-network.vertices[j])))
		return jac*(network.vertices[i]-network.vertices[j]) + force*np.array([0,1])
	elif constitutive == 'linear2':
		length = np.sqrt(length_square(network,network.vertices[i]-network.vertices[j]))
		jac = network.Ef*(network.vertices[i][1]-network.vertices[j][1])/length
		force = network.Ef*(1-np.sqrt(length_square(network,network.vertices_ini[i]-network.vertices_ini[j]))/np.sqrt(length_square(network,network.vertices[i]-network.vertices[j])))
		return jac*(network.vertices[i]-network.vertices[j])/np.sqrt(length_square(network,network.vertices[i]-network.vertices[j])) + force*np.array([0,1])

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
	return J


def iterative_newton(network, constitutive):
	max_iter = 100
	epsilon = 1e-8
	for k in range(max_iter):
		F = write_vector_F(network, constitutive)
		J = write_matrix_J(network, constitutive)
		diff = np.linalg.solve(J,-F)
		for i in range(len(network.interior_nodes)):
			j=network.interior_nodes[i]
			network.vertices[j,0] = network.vertices[j,0] + diff[2*i]
			network.vertices[j,1] = network.vertices[j,1] + diff[2*i+1]
		if np.linalg.norm(diff) < epsilon:
			print('convergence!, nre iter:', k)
			break
		elif k== max_iter-1:
			raise ValueError 
	return network



######################### NONLINEAR SCHEME 3D #########################################
# calculate the whole jacobian at one point, differentiated by itself
# calculate the force at each point
def write_force(network, i, j, constitutive):
	if constitutive == 'linear2':
		return network.Ef*(np.sqrt(length_square(network,network.vertices[i]-network.vertices[j]))-np.sqrt(length_square(network,network.vertices_ini[i]-network.vertices_ini[j])))*(network.vertices[i]-network.vertices[j])/np.sqrt(length_square(network,network.vertices[i]-network.vertices[j]))

def write_force_eq_point(network, i, constitutive):
	force_equation =0
	for j in network.list_nodes_ridges[network.interior_nodes[i]]:
		if network.interior_nodes[i]==network.ridge_vertices[j][0]:
			force_equation += write_force(network, network.interior_nodes[i], network.ridge_vertices[j][1], constitutive)
		if network.interior_nodes[i]==network.ridge_vertices[j][1]:
			force_equation += write_force(network, network.interior_nodes[i], network.ridge_vertices[j][0], constitutive)
	return force_equation

# Calculate jacobian for the force contributed by j on i, differentiated by x
def find_jacobian_x_3d(network, i, j, constitutive):
	if constitutive == 'linear2':
		length = np.sqrt(length_square(network,network.vertices[i]-network.vertices[j]))
		jac = network.Ef*(network.vertices[i][0]-network.vertices[j][0])/length
		force = network.Ef*(1-np.sqrt(length_square(network,network.vertices_ini[i]-network.vertices_ini[j]))/np.sqrt(length_square(network,network.vertices[i]-network.vertices[j])))
		return jac*(network.vertices[i]-network.vertices[j])/np.sqrt(length_square(network,network.vertices[i]-network.vertices[j])) + force*np.array([1,0,0])

# Calculate jacobian for the force contributed by j on i, differentiated by y
def find_jacobian_y_3d(network, i, j, constitutive):
	if constitutive == 'linear2':
		length = np.sqrt(length_square(network,network.vertices[i]-network.vertices[j]))
		jac = network.Ef*(network.vertices[i][1]-network.vertices[j][1])/length
		force = network.Ef*(1-np.sqrt(length_square(network,network.vertices_ini[i]-network.vertices_ini[j]))/np.sqrt(length_square(network,network.vertices[i]-network.vertices[j])))
		return jac*(network.vertices[i]-network.vertices[j])/np.sqrt(length_square(network,network.vertices[i]-network.vertices[j])) + force*np.array([0,1,0])

# Calculate jacobian for the force contributed by j on i, differentiated by z
def find_jacobian_z_3d(network, i, j, constitutive):
	if constitutive == 'linear2':
		length = np.sqrt(length_square(network,network.vertices[i]-network.vertices[j]))
		jac = network.Ef*(network.vertices[i][2]-network.vertices[j][2])/length
		force = network.Ef*(1-np.sqrt(length_square(network,network.vertices_ini[i]-network.vertices_ini[j]))/np.sqrt(length_square(network,network.vertices[i]-network.vertices[j])))
		return jac*(network.vertices[i]-network.vertices[j])/np.sqrt(length_square(network,network.vertices[i]-network.vertices[j])) + force*np.array([0,0,1])

def write_jacobian_point_3d(network, i, constitutive):
	jacobian_x=0
	jacobian_y=0
	jacobian_z=0
	for j in network.list_nodes_ridges[network.interior_nodes[i]]:
		if network.interior_nodes[i]==network.ridge_vertices[j][0]:
			jacobian_x += find_jacobian_x_3d(network, network.interior_nodes[i], network.ridge_vertices[j][1], constitutive)
			jacobian_y += find_jacobian_y_3d(network, network.interior_nodes[i], network.ridge_vertices[j][1], constitutive)
			jacobian_z += find_jacobian_z_3d(network, network.interior_nodes[i], network.ridge_vertices[j][1], constitutive)
		if network.interior_nodes[i]==network.ridge_vertices[j][1]:
			jacobian_x += find_jacobian_x_3d(network, network.interior_nodes[i], network.ridge_vertices[j][0], constitutive)
			jacobian_y += find_jacobian_y_3d(network, network.interior_nodes[i], network.ridge_vertices[j][0], constitutive)
			jacobian_z += find_jacobian_z_3d(network, network.interior_nodes[i], network.ridge_vertices[j][0], constitutive)
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
		for coord in coords:
			F[network.dimension*i+coord] = np.array(write_force_eq_point(network,i, constitutive))[coord]
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
		"""
		J[network.dimension*i,network.dimension*i]=write_jacobian_point_3d(network, i, constitutive)[0,0]
		J[network.dimension*i,network.dimension*i+1]=write_jacobian_point_3d(network, i, constitutive)[1,0]
		J[network.dimension*i,network.dimension*i+2]=write_jacobian_point_3d(network, i, constitutive)[2,0]
		J[network.dimension*i+1,network.dimension*i] = write_jacobian_point_3d(network,i, constitutive)[0,1]
		J[network.dimension*i+1,network.dimension*i+1]=write_jacobian_point_3d(network, i, constitutive)[1,1]
		J[network.dimension*i+1,network.dimension*i+2]=write_jacobian_point_3d(network, i, constitutive)[2,1]
		J[network.dimension*i+2,network.dimension*i] = write_jacobian_point_3d(network,i, constitutive)[0,2]
		J[network.dimension*i+2,network.dimension*i+1]=write_jacobian_point_3d(network, i, constitutive)[1,2]
		J[network.dimension*i+2,network.dimension*i+2]=write_jacobian_point_3d(network, i, constitutive)[2,2]
		"""
		for k in network.list_nodes_ridges[network.interior_nodes[i]]:
			if network.interior_nodes[i]==network.ridge_vertices[k][0] and network.ridge_vertices[k][1] in network.interior_nodes:
				r=network.ridge_vertices[k][1]
				j=network.interior_nodes.index(r)
				for coord1 in coords:
					J[network.dimension*i+coord1,network.dimension*j]=-find_jacobian_x_3d(network, network.interior_nodes[i], r, constitutive)[coord1]
					J[network.dimension*i+coord1,network.dimension*j+1]=-find_jacobian_y_3d(network, network.interior_nodes[i], r, constitutive)[coord1]
					J[network.dimension*i+coord1,network.dimension*j+2]=-find_jacobian_z_3d(network, network.interior_nodes[i], r, constitutive)[coord1]
			elif network.interior_nodes[i]==network.ridge_vertices[k][1] and network.ridge_vertices[k][0] in network.interior_nodes:
				r=network.ridge_vertices[k][0]
				j=network.interior_nodes.index(r)
				for coord1 in coords:
					J[network.dimension*i+coord1,network.dimension*j]=-find_jacobian_x_3d(network, network.interior_nodes[i], r, constitutive)[coord1]
					J[network.dimension*i+coord1,network.dimension*j+1]=-find_jacobian_y_3d(network, network.interior_nodes[i], r, constitutive)[coord1]
					J[network.dimension*i+coord1,network.dimension*j+2]=-find_jacobian_z_3d(network, network.interior_nodes[i], r, constitutive)[coord1]
				"""
				J[network.dimension*i,network.dimension*j]=-find_jacobian_x_3d(network, network.interior_nodes[i], r, constitutive)[0]
				J[network.dimension*i,network.dimension*j+1]=-find_jacobian_y_3d(network, network.interior_nodes[i], r, constitutive)[0]
				J[network.dimension*i,network.dimension*j+2]=-find_jacobian_z_3d(network, network.interior_nodes[i], r, constitutive)[0]
				J[network.dimension*i+1,network.dimension*j]=-find_jacobian_x_3d(network, network.interior_nodes[i], r, constitutive)[1]
				J[network.dimension*i+1,network.dimension*j+1]=-find_jacobian_y_3d(network, network.interior_nodes[i], r, constitutive)[1]
				J[network.dimension*i+1,network.dimension*j+2]=-find_jacobian_z_3d(network, network.interior_nodes[i], r, constitutive)[1]
				J[network.dimension*i+2,network.dimension*j]=-find_jacobian_x_3d(network, network.interior_nodes[i], r, constitutive)[2]
				J[network.dimension*i+2,network.dimension*j+1]=-find_jacobian_y_3d(network, network.interior_nodes[i], r, constitutive)[2]
				J[network.dimension*i+2,network.dimension*j+2]=-find_jacobian_z_3d(network, network.interior_nodes[i], r, constitutive)[2]
				
			elif network.interior_nodes[i]==network.ridge_vertices[k][1] and network.ridge_vertices[k][0] in network.interior_nodes:
				r=network.ridge_vertices[k][0]
				j=network.interior_nodes.index(r)
				for coord1 in coords:
					J[network.dimension*i+coord1,network.dimension*j]=-find_jacobian_x_3d(network, network.interior_nodes[i], r, constitutive)[coord1]
					J[network.dimension*i+coord1,network.dimension*j+1]=-find_jacobian_y_3d(network, network.interior_nodes[i], r, constitutive)[coord1]
					J[network.dimension*i+coord1,network.dimension*j+2]=-find_jacobian_z_3d(network, network.interior_nodes[i], r, constitutive)[coord1]
				
				J[network.dimension*i,network.dimension*j]=-find_jacobian_x_3d(network, network.interior_nodes[i], r, constitutive)[0]
				J[network.dimension*i,network.dimension*j+1]=-find_jacobian_y_3d(network, network.interior_nodes[i], r, constitutive)[0]
				J[network.dimension*i,network.dimension*j+2]=-find_jacobian_z_3d(network, network.interior_nodes[i], r, constitutive)[0]
				J[network.dimension*i+1,network.dimension*j]=-find_jacobian_x_3d(network, network.interior_nodes[i], r, constitutive)[1]
				J[network.dimension*i+1,network.dimension*j+1]=-find_jacobian_y_3d(network, network.interior_nodes[i], r, constitutive)[1]
				J[network.dimension*i+1,network.dimension*j+2]=-find_jacobian_z_3d(network, network.interior_nodes[i], r, constitutive)[1]
				J[network.dimension*i+2,network.dimension*j]=-find_jacobian_x_3d(network, network.interior_nodes[i], r, constitutive)[2]
				J[network.dimension*i+2,network.dimension*j+1]=-find_jacobian_y_3d(network, network.interior_nodes[i], r, constitutive)[2]
				J[network.dimension*i+2,network.dimension*j+2]=-find_jacobian_z_3d(network, network.interior_nodes[i], r, constitutive)[2]
				"""
	return J


def iterative_newton_3d(network, constitutive,step):
	max_iter = 200
	epsilon = 1e-10
	for k in range(max_iter):
		F = write_vector_F_3d(network, constitutive)
		J = write_matrix_J_3d(network, constitutive)
		diff = np.linalg.solve(J,-F)
		for i in range(len(network.interior_nodes)):
			j=network.interior_nodes[i]
			network.vertices[j,0] = network.vertices[j,0] + diff[network.dimension*i]
			network.vertices[j,1] = network.vertices[j,1] + diff[network.dimension*i+1]
			network.vertices[j,2] = network.vertices[j,2] + diff[network.dimension*i+2]
		if np.linalg.norm(diff) < epsilon:
			print(step,'convergence!, nre iter:', k)
			break
#		else:
#			print(k, 'not converged')
	return network



### Global function of solving one step if the netowrk

def solve_force_balance(network, defo, constitutive, scheme, side, step):
	if scheme == 'nonlinear':
		if network.dimension ==2:
			print(step)
			network = iterative_newton(network, constitutive)
		elif network.dimension == 3:
			network = iterative_newton_3d(network, constitutive, step)
	elif scheme == 'linear':
		network = linear_scheme(network)
	return network

