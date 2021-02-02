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


os.chdir('../Data/Study_networks/understanding/')
#os.chdir(sorted_ls('.')[-1])
#os.chdir(sorted_ls('.')[-1])
print(os.getcwd())


filenames=sorted(fnmatch.filter(sorted_ls('.'), 'network_vertices_initial_*.csv'))
network = load_network_info(706)
#network = load_network_info(523)
#network.ridge_vertices = np.delete(network.ridge_vertices, 248,axis=0)
#network.vertices = [[0.,0.],[1.0,0.]]
#network.ridge_vertices = [[0,1]]
network = network.create_ridge_node_list()
network = network.sort_nodes()
tensile_test = load_info_test(706)
space_discretization=tensile_test.space_discretization
traction_distance = tensile_test.traction_distance
iterations = abs(int(traction_distance / space_discretization))
myModel = mdb.models['Model-1']

vertices=network.vertices
ridge_vertices=network.ridge_vertices

test_number = int(np.random.rand()*10e8)

#### PART DEFINITION ####
def define_part():
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
    myModel.StaticStep(name='Step-1', previous='Initial',maxNumInc=500, minInc=1e-10, nlgeom=ON)
    #myModel.FieldOutputRequest(name='F-Output-3',createStepName='Step-1', variables=('COORD', 'S','E','SE'),numIntervals=
    #    iterations)
    myModel.fieldOutputRequests['F-Output-1'].setValues(variables=(
    'S', 'SEQUT', 'E', 'SE', 'U', 'RF', 'CF', 'COORD', 'MVF'))
    myModel.fieldOutputRequests['F-Output-1'].setValues(
    numIntervals=10)
    myModel.steps['Step-1'].setValues(stabilizationMethod=DISSIPATED_ENERGY_FRACTION, 
    continueDampingFactors=True, adaptiveDampingRatio=0.1)

def set_boundary_conditions(network): ## can be changed with the node label list
    list_vertices_right = []
    for k in network.boundary_nodes_right:
        ridge=list_nodes_ridges[k]
        #print(ridge)
        instancename = 'Part-' + str(ridge[0]+1) + '-1' 
        try:
            coords = (vertices[k][0],vertices[k][1],vertices[k][2])
        except IndexError:
            coords = (vertices[k][0],vertices[k][1],0.0)
        list_vertices_right.append(myModel.rootAssembly.instances[instancename].vertices.findAt(coords))
    list_vertices_left = []
    for k in network.boundary_nodes_left:
        ridge=list_nodes_ridges[k]
        instancename = 'Part-' + str(ridge[0]+1) + '-1' 
        try:
            coords = (vertices[k][0],vertices[k][1],vertices[k][2])
        except IndexError:
            coords = (vertices[k][0],vertices[k][1],0.0)
        list_vertices_left.append(myModel.rootAssembly.instances[instancename].vertices.findAt(coords))
    myModel.PinnedBC(createStepName='Initial', localCsys=None, name=
        'BC-1', region=Region(vertices=VertexArray(list_vertices_left)))
    myModel.DisplacementBC(amplitude=UNSET, createStepName='Step-1',
        distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, name=
        'BC-2', region=Region(vertices=VertexArray(list_vertices_right)), u1=traction_distance, u2=0.0,ur3=UNSET)


## to be deleted and adapted with network



def define_mesh(mask):
	number_elements = []
	for i in range(len(ridge_vertices)):
		instancename = 'Part-' + str(i+1) + '-1'
		myModel.rootAssembly.setElementType(elemTypes=(ElemType(
			elemCode=B21, elemLibrary=STANDARD), ), regions=(
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
		coords = (vertices[k][0],vertices[k][1],0.0)
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
	mdb.models['Model-1'].ConnectorSection(name='ConnSect-1', translationalType=
		JOIN)
	mdb.models['Model-1'].sections['ConnSect-1'].setValues(behaviorOptions=(
		ConnectorElasticity(table=((network.connector_coeff, ), ), independentComponents=(), 
		components=(6, )), ), rotationalType=ROTATION)
	mdb.models['Model-1'].sections['ConnSect-1'].behaviorOptions[0].ConnectorOptions(
		)
	myModel.rootAssembly.SectionAssignment(region=
		myModel.rootAssembly.sets['Wire-2-Set-1'], sectionName=
		'ConnSect-1')
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
		multiprocessingMode=DEFAULT, name=job_name, nodalOutputPrecision=SINGLE, 
		numCpus=1, numGPUs=0, queue=None, resultsFormat=ODB, scratch='', type=
		ANALYSIS, userSubroutine='', waitHours=0, waitMinutes=0)
	#mdb.jobs[job_name].submit(consistencyChecking=OFF)
	#mdb.jobs[job_name].waitForCompletion()


define_mesh(mask)
job()

def write_stress_report(odb,filename,network):
	picked_nodes =[]
	for node in network.boundary_nodes_right:
		part = network.list_nodes_ridges[node]
		instancename='PART-'+str(part[0]+1)+'-1'
		p = odb.rootAssembly.instances[instancename]
		picked_nodes.append(p.nodes[-1:])
	node_set_name='node_set_' +str(len(network.vertices))
	odb.rootAssembly.NodeSet(name = node_set_name, nodes = picked_nodes)
	#reports = session.xyDataListFromField(odb=odb, outputPosition=NODAL, variable=(('S', 
		#INTEGRATION_POINT, ((COMPONENT, 'S11'), )), ), nodeSets=(node_set_name, ))
	reports=session.xyDataListFromField(odb=odb, outputPosition=NODAL, variable=(('RF', NODAL, ((COMPONENT, 'RF1'), )), ),nodeSets=(node_set_name, ))
	"""
	reports_corrected=[]
	for report in reports:
		report = np.array(report)
		length = len(report)
		for k in range(1,length-1):
			if report[length-k][0] == report[length-(k+1)][0]:
				report= np.delete(report,length-k,axis=0)
		#print(report)
		reports_corrected.append(report)
	lengths = []
	for k in reports_corrected:
		lengths.append(len(k))
		if len(k)<max(lengths):
			print(k)
	strain , stress = [],[]
	import csv
	#print(max(lengths),min(lengths))
	with open(filename, 'w') as csvfile:
		spamwriter = csv.writer(csvfile)
		for i in range(max(lengths)):
			strain.append(np.array(reports_corrected[0])[i,0])
			sum=0
			for k in range(len(reports_corrected)):
				sum += np.array(reports_corrected[k])[i,1]
				if len(report[k]) == max(lengths): sum += np.array(report[k])[2*i,1]
				elif i < (max(lengths)-1)/2:
					if np.array(report[k])[i,0] == strain[-1]: sum += np.array(report[k])[2*i,1]
					else: sum += np.array(report[k])[2*i-1,1]
				else: sum += np.array(report[k])[2*i-1,1]
			spamwriter.writerow([strain[-1],sum])
			stress.append(sum)
	with open(filename, 'w') as writeFile:
		writeFile.write(strain,stress)
		writeFile.close()
	x_data = np.zeros(((max(lengths)-1)/2+1,2))
	for k in range(len(report)):
		print(list(report[k])[:((max(lengths)-1)):2])
		x_data += list(report[k])[:((max(lengths)-1)):2]
	for k in range(len(report))
	x_data = sum(report[::2])[::2]
	try:
		x_data = sum(report[::2])
	except XypError:
		x_data = sum(report[1::2])
	"""
	x_data = sum(reports)
	x0 = session.xyDataObjects[x_data.name]
	session.writeXYReport(fileName=filename, xyData=(x0, ))
"""
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
		if abs(lastFrame.frameValue- (k)*0.1) <10e-5: 
			session.writeFieldReport(fileName=name,append=OFF,sortItem='Node Label',
				odb=odb,step=0,frame=lastFrame,outputPosition=NODAL,variable=(('COORD', NODAL),))
			k+=1
"""