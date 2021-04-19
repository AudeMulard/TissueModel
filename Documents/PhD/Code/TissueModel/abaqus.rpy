# -*- coding: mbcs -*-
#
# Abaqus/CAE Release 2018 replay file
# Internal Version: 2017_11_07-17.21.41 127140
# Run by am2548 on Mon Apr 19 12:42:59 2021
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
execfile('new_solver_4.py', __main__.__dict__)
#: C:\Users\am2548\TissueModel\Documents\PhD\Code\Data_1\default\Apr-19-2021_0025
#: The section "ConnSect-1" has been assigned to 1090 wires or attachment lines.
