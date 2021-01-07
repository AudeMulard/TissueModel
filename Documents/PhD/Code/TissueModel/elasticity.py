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


os.chdir('../Data/')
os.chdir(sorted_ls('.')[-1])
os.chdir(sorted_ls('.')[-1])
print(os.getcwd())

#filenames=sorted(fnmatch.filter(os.listdir('.'), 'network_vertices_initial_*.csv'))
filenames=sorted(fnmatch.filter(os.listdir('.'), 'network_vertices_initial_*.csv'))
network = load_network_info(int(filenames[-1][-9:-4]))

#print(network.beam_profile)
tensile_test = load_info_test(int(filenames[-1][-9:-4]))
space_discretization=tensile_test.space_discretization
traction_distance = tensile_test.traction_distance
iterations = abs(int(traction_distance / space_discretization))
myModel = mdb.models['Model-1']

vertices=network.vertices
ridge_vertices=network.ridge_vertices

#partname = 'Part'
test_number = int(np.random.rand()*10e8)
#### PART DEFINITION ####
def define_part(network):
	for i in range(len(ridge_vertices)):
		ridge = ridge_vertices[i]
		partname = 'Part-' + str(i+1)
		myModel.Part(dimensionality=TWO_D_PLANAR, name=partname, type=
			DEFORMABLE_BODY)
		myModel.parts[partname].ReferencePoint(point=(0.0,0.0, 0.0))
		try:
			myModel.parts[partname].DatumPointByCoordinate(coords=(vertices[ridge[0]][0],vertices[ridge[0]][1], 0.0))
			myModel.parts[partname].DatumPointByCoordinate(coords=(vertices[ridge[1]][0],vertices[ridge[1]][1], 0.0))
		except IndexError:
			myModel.parts[partname].DatumPointByCoordinate(coords=(vertices[i][0],vertices[i][1], 0.0))
		myModel.parts[partname].WirePolyLine(mergeType=IMPRINT, meshable=
			ON, points=((myModel.parts[partname].datums[2], 
			myModel.parts[partname].datums[3]), ))
"""
	connector_list = []
	for i in range(len(ridge_vertices)):
		ridge = ridge_vertices[i]
		if ridge[0]+2<1000:
			num1 = 2
		else:
			num1 = 3
		if ridge[1]+2<1000:
			num2 = 2
		else:
			num2 = 3
		connector_list.append((myModel.parts[partname].datums[int(ridge[0]+num1)],myModel.parts[partname].datums[int(ridge[1]+num2)]))

	myModel.parts[partname].WirePolyLine(mergeType=IMPRINT, meshable=ON
					, points=tuple(connector_list))
"""


#### MATERIAL AND SECTION DEFINITION ####

# Truss section
def define_material(network):
	"""
	myModel.Material(name='Material-1')
	myModel.materials['Material-1'].Elastic(table=((network.truss_young, network.truss_poisson), ))

	myModel.TrussSection(area=network.truss_area, material='Material-1', name=
		'Section-1')

	for i in range(len(ridge_vertices)):
		partname = 'Part-' + str(i+1)
		myModel.parts[partname].SectionAssignment(offset=0.0, 
			offsetField='', offsetType=MIDDLE_SURFACE,  region=Region(
			edges=myModel.parts[partname].edges.getSequenceFromMask(mask=('[#1 ]', ), )), sectionName=
			'Section-1', thicknessAssignment=FROM_SECTION)
	
	mdb.models['Model-1'].parts[partname].Set(edges=
		mdb.models['Model-1'].parts[partname].edges.getSequenceFromMask((mask[0], 
		), ), name='Wire-2-Set-1')
	# Beam section
	"""
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
#instancename='Part_Instance'

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
    myModel.StaticStep(name='Step-2', previous='Step-1',maxNumInc=500, minInc=1e-10, nlgeom=ON)
    myModel.StaticStep(name='Step-3', previous='Step-2',maxNumInc=500, minInc=1e-10, nlgeom=ON)
    myModel.StaticStep(name='Step-4', previous='Step-3',maxNumInc=500, minInc=1e-10, nlgeom=ON)
    myModel.StaticStep(name='Step-5', previous='Step-4',maxNumInc=500, minInc=1e-10, nlgeom=ON)
    myModel.StaticStep(name='Step-6', previous='Step-5',maxNumInc=500, minInc=1e-10, nlgeom=ON)
    myModel.FieldOutputRequest(name='F-Output-3',createStepName='Step-1', variables=('COORD', 'S','E',),numIntervals=
        iterations)

def set_boundary_conditions(network): ## can be changed zith the node label list
    list_vertices_right = []
    list_node_label_right=[]
    list_node_label_left=[]
    for k in network.boundary_nodes_right:
        #k = list_node_label[node]#
        ridge=list_nodes_ridges[k]
        print(ridge)
        instancename = 'Part-' + str(ridge[0]+1) + '-1' 
        try:
            coords = (vertices[k][0],vertices[k][1],vertices[k][2])
        except IndexError:
            coords = (vertices[k][0],vertices[k][1],0.0)
        list_vertices_right.append(myModel.rootAssembly.instances[instancename].vertices.findAt(coords))
        list_node_label_right.append(myModel.rootAssembly.instances[instancename].vertices.findAt(coords).index)
    list_vertices_left = []
    for k in network.boundary_nodes_left:
        #k = list_node_label[node]
        ridge=list_nodes_ridges[k]
        instancename = 'Part-' + str(ridge[0]+1) + '-1' 
        try:
            coords = (vertices[k][0],vertices[k][1],vertices[k][2])
        except IndexError:
            coords = (vertices[k][0],vertices[k][1],0.0)
        list_vertices_left.append(myModel.rootAssembly.instances[instancename].vertices.findAt(coords))
        list_node_label_left.append(myModel.rootAssembly.instances[instancename].vertices.findAt(coords).index)
    myModel.PinnedBC(createStepName='Initial', localCsys=None, name=
        'BC-1', region=Region(vertices=VertexArray(list_vertices_left)))
    myModel.DisplacementBC(amplitude=UNSET, createStepName='Step-1',
        distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, name=
        'BC-2', region=Region(vertices=VertexArray(list_vertices_right)), u1=traction_distance, u2=0.0, u3=0.0, ur1=UNSET, ur2=UNSET, ur3=UNSET)
    myModel.boundaryConditions['BC-2'].setValuesInStep(stepName='Step-2', u1=0.0)
    myModel.boundaryConditions['BC-2'].setValuesInStep(stepName='Step-3', u1=traction_distance)
    myModel.boundaryConditions['BC-2'].setValuesInStep(stepName='Step-4', u1=0.0)
    myModel.boundaryConditions['BC-2'].setValuesInStep(stepName='Step-5', u1=traction_distance)
    myModel.boundaryConditions['BC-2'].setValuesInStep(stepName='Step-6', u1=0.0)
## to be deleted and adapted with network



def define_mesh(mask):
	number_elements = []
	for i in range(len(ridge_vertices)):
		instancename = 'Part-' + str(i+1) + '-1'
		myModel.rootAssembly.setElementType(elemTypes=(ElemType(
			elemCode=B21, elemLibrary=STANDARD), ), regions=(
			myModel.rootAssembly.instances[instancename].edges.getSequenceFromMask(
			mask=('[#1 ]', ), ), ))
		myModel.rootAssembly.seedPartInstance(deviationFactor=0.1, 
			minSizeFactor=0.1, regions=(
			mdb.models['Model-1'].rootAssembly.instances[instancename], ), size=tensile_test.element_size)
		mdb.models['Model-1'].rootAssembly.generateMesh(regions=(
		mdb.models['Model-1'].rootAssembly.instances[instancename], ))
		number_elements.append(len(mdb.models['Model-1'].rootAssembly.instances[instancename].elements))
	filename = 'number_elements_%09d.csv' % test_number
	with open(filename,'w') as writeFile:
		writer = csv.writer(writeFile,delimiter=',')
		writer.writerow(number_elements)


define_part(network)
#mask = myModel.parts[partname].edges.getMask()
define_material(network)
list_node_label=assembly(network)
set_steps(network)

list_nodes_ridges=[[] for i in range(len(vertices))]
for i in range(len(ridge_vertices)):
    list_nodes_ridges[ridge_vertices[i][0]].append(i)
    list_nodes_ridges[ridge_vertices[i][1]].append(i)

connector_list=[]
for k in range(len(list_nodes_ridges)):
	coords = (vertices[k][0],vertices[k][1],0.0)
	list_ridge = list_nodes_ridges[k]
	if len(list_ridge) > 1:
		for i in range(len(list_ridge)-1):
			instancename1='Part-'+str(list_ridge[i]+1)+'-1'
			#if len(list_ridge)> 2 and i == len(list_ridge)-1: i=-1
			#elif i == len(list_ridge)-1: 
			#    break
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
    ConnectorElasticity(table=((0.001, ), ), independentComponents=(), 
    components=(6, )), ), rotationalType=ROTATION)
mdb.models['Model-1'].sections['ConnSect-1'].behaviorOptions[0].ConnectorOptions(
    )
#myModel.rootAssembly.Set(edges=
#    myModel.rootAssembly.edges.getSequenceFromMask(('[#7 ]', ), )
#    , name='Wire-1-Set-1')
myModel.rootAssembly.SectionAssignment(region=
    myModel.rootAssembly.sets['Wire-2-Set-1'], sectionName=
    'ConnSect-1')


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
	mdb.jobs[job_name].submit(consistencyChecking=OFF)
	mdb.jobs[job_name].waitForCompletion()


define_mesh(mask)
job()

from odbAccess import *
o1 = session.openOdb(name=job_name+'.odb')
session.viewports['Viewport: 1'].setValues(displayedObject=o1)
odb = session.odbs[job_name+'.odb']
for j in range(len(odb.steps)):
	stepname = 'Step-%d' % (j+1)
	k=0
	for i in range(len(odb.steps[stepname].frames)):
		lastFrame=odb.steps[stepname].frames[i]
		name='network_vertices_%02d_%02d_%09d.csv' % (j+1,k,test_number)
		stress_data = 'stress_data_%02d_%02d_%09d.csv' % (j+1,k,test_number)
		print(i)
		print(lastFrame.frameValue, (k)*0.1)
		nf = NumberFormat(numDigits=6, precision=0, format=SCIENTIFIC)
		session.fieldReportOptions.setValues(reportFormat=COMMA_SEPARATED_VALUES, 
			numberFormat=nf)
		if abs(lastFrame.frameValue- (k)*0.1) <10e-5: 
			session.writeFieldReport(fileName=name,append=OFF,sortItem='Node Label',
				odb=odb,step=j,frame=lastFrame,outputPosition=NODAL,variable=(('COORD', NODAL),))
			k+=1
		#session.writeFieldReport(fileName=name,append=ON,sortItem='Node Label',
		#    odb=odb,step=0,frame=lastFrame,outputPosition=NODAL,variable=(('S', NODAL),))
			session.writeFieldReport(fileName=stress_data, append=OFF, sortItem='Node Label', 
				odb=odb, step=j, frame=lastFrame, outputPosition=NODAL, 
				variable=(('S',INTEGRATION_POINT,((COMPONENT, 'S11'), 
				)), 
				))
"""
while mdb.jobs[job_name].status != COMPLETED:
	try:
		define_mesh(mask)
		job()
		ouput_writing()
	except VisError:
		tensile_test.element_size = tensile_test.element_size/2
		define_mesh(mask)
		job()
		ouput_writing()
	tensile_test.element_size = tensile_test.element_size/2


#odb = session.odbs['C:/Temp/trans_p100_E8-8_SY200_N1-INIT-T0-T1.odb']
session.writeFieldReport(fileName='abaqus1.rpt', append=OFF,
    sortItem='Node Label', odb=odb, step=1, frame=odb.steps['Step-1'].frames[-1],
    outputPosition=INTEGRATION_POINT, variable=(('COORD', INTEGRATION_POINT), ))

mdb.jobs['Job-726'].waitForCompletion()

# Get coordinates at each step
from odbAccess import *
odb = openOdb(path='Job-726.odb')
for i in range(len(odb.steps['Step-1'].frames)):
    lastFrame=odb.steps['Step-1'].frames[i]
    with open('results1.txt','w') as writeFile:
        #writer = csv.writer(writeFile,delimiter=',')
        for i in range(len(lastFrame.fieldOutputs['COORD'].values)):
            stress_step +=lastFrame.fieldOutputs['COORD'].values[i].data[0]
            writeFile.write(str(lastFrame.fieldOutputs['COORD'].values[i].data[0]+','+str(lastFrame.fieldOutputs['COORD'].values[i].data[1])+ '\n'))


stress=[]
strain=[]
from odbAccess import *
odb = openOdb(path='Job-1.odb')
for i in range(len(odb.steps['Step-1'].frames)):
    lastFrame=odb.steps['Step-1'].frames[i]
	stress_step=0
	strain_step=0
	with open('results1.txt','w') as writeFile:
		#writer = csv.writer(writeFile,delimiter=',')
		for i in range(len(lastFrame.fieldOutputs['S'].values)):
			stress_step +=lastFrame.fieldOutputs['S'].values[i].data[0]
			strain_step +=lastFrame.fieldOutputs['E'].values[i].data[0]
			writeFile.write(str(lastFrame.fieldOutputs['S'].values[i].data[0])+ '\n')
	stress.append(stress_step/len(lastFrame.fieldOutputs['S'].values))
	strain.append(strain_step/len(lastFrame.fieldOutputs['S'].values))
"""