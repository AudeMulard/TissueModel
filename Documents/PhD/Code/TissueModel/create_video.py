from abaqus import *
from abaqusConstants import *

job_name = 'Job-646906510'

from odbAccess import *
o1 = session.openOdb(name=job_name+'.odb',readOnly=False)
session.viewports['Viewport: 1'].setValues(displayedObject=o1)
session.viewports['Viewport: 1'].odbDisplay.basicOptions.setValues(renderBeamProfiles=ON,
    beamScaleFactor=1.0)
session.animationController.animationOptions.setValues(
    timeHistoryMode=TIME_BASED, timeIncrement=0.02, minTimeAutoCompute=True, 
    maxTimeAutoCompute=True)
session.animationController.setValues(animationType=TIME_HISTORY, viewports=(
    'Viewport: 1', ))
session.animationController.play(duration=UNLIMITED)

session.imageAnimationOptions.setValues(vpDecorations=ON, vpBackground=OFF, 
    compass=OFF,timeScale=1, frameRate=50)	

session.writeImageAnimation(
    fileName='C:/Users/am2548/TissueModel/Documents/PhD/Code/TissueModel/video_gm', 
    format=QUICKTIME, canvasObjects=(session.viewports['Viewport: 1'], ))