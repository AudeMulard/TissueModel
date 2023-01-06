# -*- coding: mbcs -*-
from part import *
from material import *
from section import *
from assembly import *
from step import *
from interaction import *
from load import *
from mesh import *
from optimization import *
from job import *
from sketch import *
from visualization import *
from connectorBehavior import *
import csv,os
import numpy as np
from Plotting.information_network import load_network_info, sorted_ls, load_info_test
from Network_generation import *
import fnmatch

material_density=5 #lower can lead to instability: the higher, the higher the stable limit
bulk_viscosity = 1.2 #do not change

os.chdir('../Data/')# insert where the network geometry was recorded
os.chdir(sorted_ls('.')[-1])
os.chdir(sorted_ls('.')[-1])

filenames=sorted(fnmatch.filter(sorted_ls('.'), 'network_vertices_initial_*.csv'))
network = load_network_info(int(filenames[-1][-9:-4]))


ridges_to_delete = []
for i in range(len(network.ridge_vertices)):
	if network.ridge_vertices[i][0]==network.ridge_vertices[i][1]:
		ridges_to_delete.append(i)
ridges_to_delete=np.array(sorted(ridges_to_delete, reverse=True))
print(ridges_to_delete)
for k in ridges_to_delete:
	network.ridge_vertices=np.delete(network.ridge_vertices, int(k), axis=0)

tensile_test = load_info_test(int(filenames[-1][-9:-4]))
print(tensile_test.element_size)
space_discretization=tensile_test.space_discretization
traction_distance = tensile_test.traction_distance
iterations = abs(int(traction_distance / space_discretization))
myModel = mdb.models['Model-1']

vertices=network.vertices
ridge_vertices=network.ridge_vertices

step_time =25.#tried at 22, goes down on the end

test_number = int(np.random.rand()*10e8)
#### PART DEFINITION ####
def define_part():
	if int(network.dimension) == 2:
		for i in range(len(ridge_vertices)):
			ridge = ridge_vertices[i]
			partname = 'Part-' + str(i+1)
			myModel.Part(dimensionality=TWO_D_PLANAR, name=partname, type=
				DEFORMABLE_BODY)
			try:
				myModel.parts[partname].DatumPointByCoordinate(coords=(vertices[ridge[0]][0],vertices[ridge[0]][1], 0.0))
				myModel.parts[partname].DatumPointByCoordinate(coords=(vertices[ridge[1]][0],vertices[ridge[1]][1], 0.0))
			except IndexError:
				myModel.parts[partname].DatumPointByCoordinate(coords=(vertices[i][0],vertices[i][1], 0.0))
			myModel.parts[partname].WirePolyLine(mergeType=IMPRINT, meshable=
				ON, points=((myModel.parts[partname].datums[1], 
				myModel.parts[partname].datums[2]), ))
	elif int(network.dimension)==3:
		for i in range(len(ridge_vertices)):
			ridge = ridge_vertices[i]
			partname = 'Part-' + str(i+1)
			myModel.Part(dimensionality=THREE_D, name=partname, type=
				DEFORMABLE_BODY)
			try:
				myModel.parts[partname].DatumPointByCoordinate(coords=(vertices[ridge[0]][0],vertices[ridge[0]][1], vertices[ridge[0]][2]))
				myModel.parts[partname].DatumPointByCoordinate(coords=(vertices[ridge[1]][0],vertices[ridge[1]][1], vertices[ridge[1]][2]))
			except IndexError:
				myModel.parts[partname].DatumPointByCoordinate(coords=(vertices[i][0],vertices[i][1], vertices[i][2]))
			myModel.parts[partname].WirePolyLine(mergeType=IMPRINT, meshable=
				ON, points=((myModel.parts[partname].datums[1], 
				myModel.parts[partname].datums[2]), ))


#### MATERIAL AND SECTION DEFINITION ####

# Truss section
def define_material(network):
	myModel.Material(name='Material-2')
	myModel.materials['Material-2'].Elastic(table=((network.beam_young, network.beam_poisson), ))
	myModel.CircularProfile(name='Profile-1', r=network.beam_profile)
	myModel.BeamSection(consistentMassMatrix=False, integration=
		DURING_ANALYSIS, material='Material-2', name='Section-2', poissonRatio=0.0, 
		profile='Profile-1', temperatureVar=LINEAR)
	for i in range(len(ridge_vertices)):
		partname = 'Part-' + str(i+1)
		myModel.parts[partname].SectionAssignment(offset=0.0, 
			offsetField='', offsetType=MIDDLE_SURFACE, region=Region(
			edges=myModel.parts[partname].edges.getSequenceFromMask(mask=('[#1 ]', ), )),
			sectionName='Section-2', thicknessAssignment=
			FROM_SECTION)
	
	for i in range(len(ridge_vertices)):
		partname = 'Part-' + str(i+1)
		myModel.parts[partname].assignBeamSectionOrientation(method=
			N1_COSINES, n1=(0.0, 0.0, -1.0), region=Region(
			edges=myModel.parts[partname].edges.getSequenceFromMask(mask=('[#1 ]', ), )))
	mdb.models['Model-1'].rootAssembly.DatumCsysByDefault(CARTESIAN)
	mdb.models['Model-1'].materials['Material-2'].Density(table=((material_density, ), ))

#### ASSEMBLY ####
# Creation of instances


def assembly(network):
	list_node_label=[]
	for i in range(len(ridge_vertices)):
		partname = 'Part-' + str(i+1)
		instancename = 'Part-' + str(i+1) + '-1'
		myModel.rootAssembly.Instance(dependent=OFF, name=instancename, 
			part=myModel.parts[partname])
	for k in range(len(vertices)):
		ridge=network.list_nodes_ridges[k]
		instancename = 'Part-' + str(ridge[0]+1) + '-1' 
		try:
			coords = (vertices[k][0],vertices[k][1],vertices[k][2])
		except IndexError:
			coords = (vertices[k][0],vertices[k][1],0.0)
		list_node_label.append(myModel.rootAssembly.instances[instancename].vertices.findAt(coords).index)
	filename = 'node_label_%09d.csv' % test_number
	with open(filename,'w') as writeFile:
		writer = csv.writer(writeFile,delimiter=',')
		writer.writerow(list_node_label)
	return list_node_label


# Step Creation
def set_steps(network):
    mdb.models['Model-1'].ExplicitDynamicsStep(name='Step-1', previous='Initial', improvedDtMethod=ON,quadBulkViscosity=bulk_viscosity)
    mdb.models['Model-1'].steps['Step-1'].setValues(timePeriod=step_time, improvedDtMethod=ON)
    myModel.fieldOutputRequests['F-Output-1'].setValues(variables=('S', 'E', 'U', 'RF', 'CF', 'COORD'))
    myModel.fieldOutputRequests['F-Output-1'].setValues(
		numIntervals=iterations)

def set_boundary_conditions(network): 
    mdb.models['Model-1'].SmoothStepAmplitude(name='Amp-2', timeSpan=TOTAL, data=((
    0.0, 0.0), (step_time, 1.0)))
    for k in network.boundary_nodes_right:
        ridge=list_nodes_ridges[k]
        instancename = 'Part-' + str(ridge[0]+1) + '-1'
        try:
            coords = (vertices[k][0],vertices[k][1],vertices[k][2])
        except IndexError:
            coords = (vertices[k][0],vertices[k][1],0.0)
        vertice = myModel.rootAssembly.instances[instancename].vertices.findAt(coords)
        if vertice.index==0:
            myModel.DisplacementBC(amplitude='Amp-2', createStepName='Step-1', distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, name='BC-'+str(ridge[0]+1), region=Region(vertices=mdb.models['Model-1'].rootAssembly.instances[instancename].vertices[:1]), u1=traction_distance, u2=0.0,u3=0.0,ur3=UNSET)
        else:
            myModel.DisplacementBC(name='BC-'+str(ridge[0]+1), createStepName='Step-1', region=Region(vertices=mdb.models['Model-1'].rootAssembly.instances[instancename].vertices[1:]), u1=traction_distance,u2=0.0,u3=0.0,ur3=UNSET, amplitude='Amp-2', fixed=OFF, distributionType=UNIFORM, fieldName='', localCsys=None)
    for k in network.boundary_nodes_left:
        ridge=list_nodes_ridges[k]
        instancename = 'Part-' + str(ridge[0]+1) + '-1'
        try:
            coords = (vertices[k][0],vertices[k][1],vertices[k][2])
        except IndexError:
            coords = (vertices[k][0],vertices[k][1],0.0)
        vertice = myModel.rootAssembly.instances[instancename].vertices.findAt(coords)
        if vertice.index==0:
            myModel.PinnedBC(createStepName='Initial', localCsys=None, name='BC-'+str(ridge[0]+1), region=Region(vertices=mdb.models['Model-1'].rootAssembly.instances[instancename].vertices[:1]))
        else:
            myModel.PinnedBC(createStepName='Initial', localCsys=None, name='BC-'+str(ridge[0]+1), region=Region(vertices=mdb.models['Model-1'].rootAssembly.instances[instancename].vertices[1:]))


def define_mesh(mask):
	number_elements = []
	for i in range(len(ridge_vertices)):
		instancename = 'Part-' + str(i+1) + '-1'
		if int(network.dimension)==2:
			myModel.rootAssembly.setElementType(elemTypes=(ElemType(
				elemCode=B21, elemLibrary=EXPLICIT), ), regions=(
				myModel.rootAssembly.instances[instancename].edges.getSequenceFromMask(
				mask=('[#1 ]', ), ), ))
		elif int(network.dimension)==3:
				myModel.rootAssembly.setElementType(elemTypes=(ElemType(
				elemCode=B31, elemLibrary=EXPLICIT), ), regions=(
				myModel.rootAssembly.instances[instancename].edges.getSequenceFromMask(
				mask=('[#1 ]', ), ), ))
		myModel.rootAssembly.seedPartInstance(regions=(
			mdb.models['Model-1'].rootAssembly.instances[instancename], ), size=tensile_test.element_size)
		mdb.models['Model-1'].rootAssembly.generateMesh(regions=(
		mdb.models['Model-1'].rootAssembly.instances[instancename], ))
		number_elements.append(len(mdb.models['Model-1'].rootAssembly.instances[instancename].elements))
	filename = 'number_elements_%09d.csv' % test_number
	with open(filename,'w') as writeFile:
		writer = csv.writer(writeFile,delimiter=',')
		writer.writerow(number_elements)



list_nodes_ridges=[[] for i in range(len(vertices))]
for i in range(len(ridge_vertices)):
	list_nodes_ridges[ridge_vertices[i][0]].append(i)
	list_nodes_ridges[ridge_vertices[i][1]].append(i)


def create_connectors(network):
	connector_list=[]
	for k in range(len(list_nodes_ridges)):
		if int(network.dimension)==2: coords = (vertices[k][0],vertices[k][1],0.0)
		elif int(network.dimension)==3: coords = (vertices[k][0],vertices[k][1],vertices[k][2])
		list_ridge = list_nodes_ridges[k]
		if len(list_ridge) > 1:
			for i in range(len(list_ridge)-1):
				instancename1='Part-'+str(list_ridge[i]+1)+'-1'
				instancename2='Part-'+str(list_ridge[i+1]+1)+'-1'
				vertice1 = myModel.rootAssembly.instances[instancename1].vertices.findAt(coords)
				vertice2 = myModel.rootAssembly.instances[instancename2].vertices.findAt(coords)
				connector_list.append((vertice1, vertice2))
	myModel.rootAssembly.WirePolyLine(mergeType=IMPRINT, meshable=OFF
					, points=tuple(connector_list))
	mask = mdb.models['Model-1'].rootAssembly.edges.getMask()
	mdb.models['Model-1'].rootAssembly.Set(edges=
		mdb.models['Model-1'].rootAssembly.edges.getSequenceFromMask((mask[0], 
		), ), name='Wire-2-Set-1')
	mdb.models['Model-1'].MPCSection(name='ConnSect-2', mpcType=PIN_MPC, 
		userMode=DOF_MODE, userType=0)
	myModel.rootAssembly.SectionAssignment(region=
		myModel.rootAssembly.sets['Wire-2-Set-1'], sectionName=
		'ConnSect-2') # make two options? of connector coeff or non connector coeff?
	return mask

define_part()
define_material(network)
list_node_label=assembly(network)
set_steps(network)

mask = create_connectors(network)
network = network.sort_nodes()
set_boundary_conditions(network)



#### JOB #### 
job_name = 'Job-'+str(test_number)


def job():
	mdb.Job(atTime=None, contactPrint=OFF, description='', echoPrint=OFF, 
		explicitPrecision=SINGLE, getMemoryFromAnalysis=True, historyPrint=OFF, 
		memory=90, memoryUnits=PERCENTAGE, model='Model-1', modelPrint=OFF, 
		multiprocessingMode=MPI, name=job_name, nodalOutputPrecision=SINGLE, 
		numCpus=4, numDomains=4,numGPUs=0, queue=None, resultsFormat=ODB, scratch='', type=
		ANALYSIS, userSubroutine='', waitHours=0, waitMinutes=0)
	mdb.jobs[job_name].submit(consistencyChecking=OFF)
	mdb.jobs[job_name].waitForCompletion()


define_mesh(mask)
job()

def write_stress_report(odb,filename,network):
	picked_nodes =[]
	for node in network.boundary_nodes_right:
		part = network.list_nodes_ridges[node]
		instancename='PART-'+str(part[0]+1)+'-1'
		p = odb.rootAssembly.instances[instancename]
		for k in range(len(p.nodes)):
			if p.nodes[k].coordinates[0]>=1-10e-6:
				picked_nodes.append(p.nodes[k:k+1])
				break
	node_set_name='node_set_bc_10_' +str(len(network.vertices))
	odb.rootAssembly.NodeSet(name = node_set_name, nodes = picked_nodes)
	reports=session.xyDataListFromField(odb=odb, outputPosition=NODAL, variable=(('RF', NODAL, ((COMPONENT, 'RF1'), )), ),nodeSets=(node_set_name, ))
	x_data = sum(reports)
	x0 = session.xyDataObjects[x_data.name]
	session.writeXYReport(fileName=filename, xyData=(x0, ))

from odbAccess import *
o1 = session.openOdb(name=job_name+'.odb',readOnly=False)
session.viewports['Viewport: 1'].setValues(displayedObject=o1)
odb = session.odbs[job_name+'.odb']
nf = NumberFormat(numDigits=6, precision=0, format=SCIENTIFIC)
session.fieldReportOptions.setValues(reportFormat=COMMA_SEPARATED_VALUES, 
			numberFormat=nf)
for j in range(len(odb.steps)):
	stepname = 'Step-%d' % (j+1)
	k=0
	stress_data = 'stress_data_%02d_%09d.rpt' % (j+1,int(test_number))
	write_stress_report(odb,stress_data,network)
	for i in range(len(odb.steps[stepname].frames)):
		lastFrame=odb.steps[stepname].frames[i]
		name='network_vertices_%02d_%02d_%09d.csv' % (j+1,k,int(test_number))
		print(i)
		print(lastFrame.frameValue, (k)*0.1)
		session.writeFieldReport(fileName=name,append=OFF,sortItem='Node Label',
			odb=odb,step=0,frame=lastFrame,outputPosition=NODAL,variable=(('COORD', NODAL),))
		k+=1
