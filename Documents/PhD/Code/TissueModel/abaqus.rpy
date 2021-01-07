# -*- coding: mbcs -*-
#
# Abaqus/CAE Release 2018 replay file
# Internal Version: 2017_11_07-17.21.41 127140
# Run by am2548 on Thu Jan  7 08:52:52 2021
#

# from driverUtils import executeOnCaeGraphicsStartup
# executeOnCaeGraphicsStartup()
#: Executing "onCaeGraphicsStartup()" in the site directory ...
from abaqus import *
from abaqusConstants import *
session.Viewport(name='Viewport: 1', origin=(0.0, 0.0), width=136.345321655273, 
    height=155.590744018555)
session.viewports['Viewport: 1'].makeCurrent()
session.viewports['Viewport: 1'].maximize()
from caeModules import *
from driverUtils import executeOnCaeStartup
executeOnCaeStartup()
session.viewports['Viewport: 1'].partDisplay.geometryOptions.setValues(
    referenceRepresentation=ON)
execfile('new_solver_1.py', __main__.__dict__)
#: C:\Users\am2548\TissueModel\Documents\PhD\Code\Data_1\Study_networks\Jan-06-2021_0036
#: The section "ConnSect-1" has been assigned to 3139 wires or attachment lines.
session.viewports['Viewport: 1'].partDisplay.setValues(sectionAssignments=ON, 
    engineeringFeatures=ON)
session.viewports['Viewport: 1'].partDisplay.geometryOptions.setValues(
    referenceRepresentation=OFF)
a = mdb.models['Model-1'].rootAssembly
session.viewports['Viewport: 1'].setValues(displayedObject=a)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(
    optimizationTasks=OFF, geometricRestrictions=OFF, stopConditions=OFF)
mdb.jobs['Job-935729267'].submit(consistencyChecking=OFF)
#: The job input file "Job-935729267.inp" has been submitted for analysis.
#: Job Job-935729267: Analysis Input File Processor completed successfully.
session.viewports['Viewport: 1'].assemblyDisplay.setValues(step='Step-1')
session.viewports['Viewport: 1'].assemblyDisplay.setValues(
    adaptiveMeshConstraints=ON)
mdb.models['Model-1'].fieldOutputRequests['F-Output-1'].setValues(variables=(
    'S', 'SEQUT', 'E', 'SE', 'U', 'RF', 'CF', 'COORD', 'MVF'))
mdb.models['Model-1'].fieldOutputRequests['F-Output-1'].setValues(
    numIntervals=10)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(
    adaptiveMeshConstraints=OFF)
