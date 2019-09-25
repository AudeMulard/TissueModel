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
	network.vertices=new_positions
	return network

### Calculate all the functions and jacobians

# Def of norm square
def length_square(network,x):
	return x[0]**2+x[1]**2


# Calculate jacobian for the force contributed by j on i, differentiated by x
def find_jacobian_x(network, i, j, constitutive):
	length = np.sqrt(length_square(network,network.vertices[i]-network.vertices[j]))
	jac = network.Ef*(network.vertices[i][0]-network.vertices[j][0])/length
	force = network.Ef*(1-np.sqrt(length_square(network,network.vertices_ini[i]-network.vertices_ini[j]))/np.sqrt(length_square(network,network.vertices[i]-network.vertices[j])))
	return jac*(network.vertices[i]-network.vertices[j])/np.sqrt(length_square(network,network.vertices[i]-network.vertices[j])) + force*np.array([1,0])

# Calculate jacobian for the force contributed by j on i, differentiated by y
def find_jacobian_y(network, i, j, constitutive):
	length = np.sqrt(length_square(network,network.vertices[i]-network.vertices[j]))
	jac = network.Ef*(network.vertices[i][1]-network.vertices[j][1])/length
	force = network.Ef*(1-np.sqrt(length_square(network,network.vertices_ini[i]-network.vertices_ini[j]))/np.sqrt(length_square(network,network.vertices[i]-network.vertices[j])))
	return jac*(network.vertices[i]-network.vertices[j])/np.sqrt(length_square(network,network.vertices[i]-network.vertices[j])) + force*np.array([0,1])

def write_jacobian_point(network, i, constitutive):
	jacobian_x=0
	jacobian_y=0
	if i < len(network.interior_nodes):
		for j in network.list_nodes_ridges[network.interior_nodes[i]]:
			if network.interior_nodes[i]==network.ridge_vertices[j][0]:
				jacobian_x += find_jacobian_x(network, network.interior_nodes[i], network.ridge_vertices[j][1], constitutive)
				jacobian_y += find_jacobian_y(network, network.interior_nodes[i], network.ridge_vertices[j][1], constitutive)
			if network.interior_nodes[i]==network.ridge_vertices[j][1]:
				jacobian_x += find_jacobian_x(network, network.interior_nodes[i], network.ridge_vertices[j][0], constitutive)
				jacobian_y += find_jacobian_y(network, network.interior_nodes[i], network.ridge_vertices[j][0], constitutive)
	elif i < 2*len(network.interior_nodes)+len(network.boundary_nodes_left):
		b = i-2*len(network.interior_nodes)
		for j in network.list_nodes_ridges[network.boundary_nodes_left[b]]:
			if network.boundary_nodes_left[b]==network.ridge_vertices[j][0]:
				jacobian_x += find_jacobian_x(network, network.boundary_nodes_left[b], network.ridge_vertices[j][1], constitutive)
				jacobian_y += find_jacobian_y(network, network.boundary_nodes_left[b], network.ridge_vertices[j][1], constitutive)
			if network.boundary_nodes_left[b]==network.ridge_vertices[j][1]:
				jacobian_x += find_jacobian_x(network, network.boundary_nodes_left[b], network.ridge_vertices[j][0], constitutive)
				jacobian_y += find_jacobian_y(network, network.boundary_nodes_left[b], network.ridge_vertices[j][0], constitutive)
	else:
		b = i-2*len(network.interior_nodes)-len(network.boundary_nodes_left)
		for j in network.list_nodes_ridges[network.boundary_nodes_right[b]]:
			if network.boundary_nodes_right[b]==network.ridge_vertices[j][0]:
				jacobian_x += find_jacobian_x(network, network.boundary_nodes_right[b], network.ridge_vertices[j][1], constitutive)
				jacobian_y += find_jacobian_y(network, network.boundary_nodes_right[b], network.ridge_vertices[j][1], constitutive)
			if network.boundary_nodes_right[b]==network.ridge_vertices[j][1]:
				jacobian_x += find_jacobian_x(network, network.boundary_nodes_right[b], network.ridge_vertices[j][0], constitutive)
				jacobian_y += find_jacobian_y(network, network.boundary_nodes_right[b], network.ridge_vertices[j][0], constitutive)
	return np.array([jacobian_x, jacobian_y])

### Iterative Newton scheme to solve the nonlinear equations
def write_force(network, i, j, constitutive):
	return network.Ef*(np.sqrt(length_square(network,network.vertices[i]-network.vertices[j]))-np.sqrt(length_square(network,network.vertices_ini[i]-network.vertices_ini[j])))*(network.vertices[i]-network.vertices[j])/np.sqrt(length_square(network,network.vertices[i]-network.vertices[j]))

def write_force_eq_point(network, i, constitutive):
	force_equation =0
	for j in network.list_nodes_ridges[network.interior_nodes[i]]:
		if network.interior_nodes[i]==network.ridge_vertices[j][0]:
			force_equation += write_force(network, network.interior_nodes[i], network.ridge_vertices[j][1], constitutive)
		if network.interior_nodes[i]==network.ridge_vertices[j][1]:
			force_equation += write_force(network, network.interior_nodes[i], network.ridge_vertices[j][0], constitutive)
	return force_equation

# Calculate vector F

def write_vector_F(network, constitutive):
	F = np.zeros(len(network.interior_nodes)*2+len(network.boundary_nodes_left)+len(network.boundary_nodes_right))
	for i in range(len(network.interior_nodes)):
		F[2*i] = np.array(write_force_eq_point(network,i, constitutive))[0]
		F[2*i+1] = np.array(write_force_eq_point(network,i, constitutive))[1]
	for i in range(len(network.boundary_nodes_left)):
		F[2*len(network.interior_nodes)+i]=np.array(write_force_eq_point(network,i, constitutive))[1]
	for i in range(len(network.boundary_nodes_right)):
		F[2*len(network.interior_nodes)+len(network.boundary_nodes_left)+i]=np.array(write_force_eq_point(network,i, constitutive))[1]		
	return F

# Calculate matrix J

def write_matrix_J(network, constitutive):
	J = np.zeros((len(network.interior_nodes)*2+len(network.boundary_nodes_left)+len(network.boundary_nodes_right),len(network.interior_nodes)*2+len(network.boundary_nodes_left)+len(network.boundary_nodes_right)))
	for i in range(len(network.interior_nodes)):
		J[2*i,2*i] = write_jacobian_point(network, i, constitutive)[0,0]
		J[2*i,2*i+1] = write_jacobian_point(network, i, constitutive)[1,0]
		J[2*i+1,2*i] = write_jacobian_point(network, i, constitutive)[0,1]
		J[2*i+1,2*i+1] = write_jacobian_point(network, i, constitutive)[1,1]
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
	for i in range(len(network.boundary_nodes_left)):
		b = 2*len(network.interior_nodes)+i
		J[b,b] = write_jacobian_point(network, b, constitutive)[0,0]
		J[b,b+1] = write_jacobian_point(network, b, constitutive)[1,0]
		for k in network.list_nodes_ridges[network.boundary_nodes_left[i]]:
			if network.boundary_nodes_left[i]==network.ridge_vertices[k][0] and network.ridge_vertices[k][1] in network.interior_nodes:
				r=network.ridge_vertices[k][1]
				j=network.interior_nodes.index(r)
				J[b,2*j]=-find_jacobian_x(network, network.interior_nodes[i], r, constitutive)[0]
				J[b,2*j+1]=-find_jacobian_y(network, network.interior_nodes[i], r, constitutive)[0]
			if network.boundary_nodes_left[i]==network.ridge_vertices[k][1] and network.ridge_vertices[k][0] in network.interior_nodes:
				r=network.ridge_vertices[k][0]
				j=network.interior_nodes.index(r)
				J[b,2*j]=-find_jacobian_x(network, network.interior_nodes[i], r, constitutive)[0]
				J[b,2*j+1]=-find_jacobian_y(network, network.interior_nodes[i], r, constitutive)[0]
	for i in range(len(network.boundary_nodes_right)):
		print b
		b = 2*len(network.interior_nodes)+len(network.boundary_nodes_left) +i
		J[b,b] = write_jacobian_point(network, b, constitutive)[0,0]
		J[b,b+1] = write_jacobian_point(network, b, constitutive)[1,0]
		for k in network.list_nodes_ridges[network.boundary_nodes_left[i]]:
			if network.boundary_nodes_left[i]==network.ridge_vertices[k][0] and network.ridge_vertices[k][1] in network.interior_nodes:
				r=network.ridge_vertices[k][1]
				j=network.interior_nodes.index(r)
				J[b,2*j]=-find_jacobian_x(network, network.interior_nodes[i], r, constitutive)[0]
				J[b,2*j+1]=-find_jacobian_y(network, network.interior_nodes[i], r, constitutive)[0]
			if network.boundary_nodes_left[i]==network.ridge_vertices[k][1] and network.ridge_vertices[k][0] in network.interior_nodes:
				r=network.ridge_vertices[k][0]
				j=network.interior_nodes.index(r)
				J[b,2*j]=-find_jacobian_x(network, network.interior_nodes[i], r, constitutive)[0]
				J[b,2*j+1]=-find_jacobian_y(network, network.interior_nodes[i], r, constitutive)[0]
	return J


def iterative_newton(network, constitutive):
	max_iter = 100
	epsilon = 1e-8
	for k in range(max_iter):
		F = write_vector_F(network, constitutive)
		J = write_matrix_J(network, constitutive)
		print J, len(J)
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

### Global function of solving one step if the netowrk

def solve_force_balance(network, defo, constitutive, scheme, side, step):
	print step
	network = iterative_newton(network, constitutive)
	return network

