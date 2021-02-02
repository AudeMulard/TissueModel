# -*- coding: mbcs -*-
#
# Abaqus/CAE Release 2018 replay file
# Internal Version: 2017_11_07-17.21.41 127140
# Run by am2548 on Fri Jan 29 19:07:10 2021
#

# from driverUtils import executeOnCaeGraphicsStartup
# executeOnCaeGraphicsStartup()
#: Executing "onCaeGraphicsStartup()" in the site directory ...
from abaqus import *
from abaqusConstants import *
session.Viewport(name='Viewport: 1', origin=(1.38281, 1.38426), width=203.55, 
    height=137.319)
session.viewports['Viewport: 1'].makeCurrent()
from driverUtils import executeOnCaeStartup
executeOnCaeStartup()
execfile('new_solver_2.py', __main__.__dict__)
#: C:\Users\am2548\TissueModel\Documents\PhD\Code\Data_1\Com_Vor_GN\Jan-29-2021_0005
#: The section "ConnSect-1" has been assigned to 4197 wires or attachment lines.
#: Model: C:/Users/am2548/TissueModel/Documents/PhD/Code/Data_1/Com_Vor_GN/Jan-29-2021_0005/Job-247131712.odb
#: Number of Assemblies:         1
#: Number of Assembly instances: 0
#: Number of Part instances:     2587
#: Number of Meshes:             2588
#: Number of Element Sets:       4
#: Number of Node Sets:          2
#: Number of Steps:              1
#: 0
#: (0.0, 0.0)
#: 1
#: (0.0399999991059303, 0.1)
#: 2
#: (0.0799999982118607, 0.2)
#: 3
#: (0.119999997317791, 0.3)
#: 4
#: (0.159999996423721, 0.4)
#: 5
#: (0.200000002980232, 0.5)
#: 6
#: (0.239999994635582, 0.6)
#: 7
#: (0.24432261288166, 0.7)
print 'RT script done'
#: RT script done
