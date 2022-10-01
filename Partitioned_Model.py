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
from abaqus import getInput
from odbAccess import *
import numpy as np
import math
import regionToolset
from scipy.optimize import fsolve
import numpy as np
import os
session.journalOptions.setValues(replayGeometry=COORDINATE, recoverGeometry=COORDINATE)

setPath = r'D:/OneDrive - Indian Institute of Science/TRINA_TCPA/Week 23/Simulations2/'
os.chdir(setPath)


#Max Diameter of the helix
D = 0.75 #mm
#Diameter of the nylon wire
ds = 0.3 #mm
#Hegiht pf Helix
H = 150.0 #mm
#Length of Nylon
L = 500.0 #mm
#Number of Turns
N = 280.0

#If the avobe data is not present, explicitly provide a Pitch below
Pitch = H/N

#Number of turns of the coil
N_Rev = 2

Dmid = D-ds #Diameter of the helix
R = D/2.0
rs = ds/2.0
Rmid = Dmid/2.0
H_model = Pitch*N_Rev

#Number of Partitions (Equally Spaced
n_partition = 4
r_parts = []
for i in range(n_partition):
    r_parts = [rs*(n_partition-i)/(n_partition)]+r_parts
r_parts_temp = [0] + r_parts

#For unequal spacing, comment the 2 previous lines and use the 2 lines below
# r_parts = [rs/5, rs/2, rs/1.15, rs]; #ideally, the last part has a max radius of the wire
# n_partition = len(r_parts)
# r_parts_temp = [0] + r_parts



# ModelName = 'Model-1'
# jobName = ModelName
# outputDispName = ModelName

Temp_Part = 'TempPart'
Use_Part = 'Helix-1'


#Scaled Geometry
#Helix Shell
mdb.models['Model-1'].ConstrainedSketch(name='__profile__', sheetSize=1.0)
mdb.models['Model-1'].sketches['__profile__'].ConstructionLine(point1=(0.0, 
    -0.5), point2=(0.0, 0.5))
mdb.models['Model-1'].sketches['__profile__'].geometry.findAt((0.0, 0.0))
mdb.models['Model-1'].sketches['__profile__'].FixedConstraint(entity=
    mdb.models['Model-1'].sketches['__profile__'].geometry.findAt((0.0, 0.0), 
    ))
for i in range(n_partition):
    mdb.models['Model-1'].sketches['__profile__'].CircleByCenterPerimeter(center=(
        Rmid , 0.0), point1=((Rmid + r_parts[i]) , 0.0))
mdb.models['Model-1'].Part(dimensionality=THREE_D, name=Temp_Part, type=
    DEFORMABLE_BODY)
mdb.models['Model-1'].parts[Temp_Part].BaseShellRevolve(angle=360*N_Rev, 
    flipPitchDirection=OFF, flipRevolveDirection=OFF, moveSketchNormalToPath=
    OFF, pitch=Pitch, sketch=mdb.models['Model-1'].sketches['__profile__'])
del mdb.models['Model-1'].sketches['__profile__']

#Closing the shell
mdb.models['Model-1'].parts[Temp_Part].DatumPlaneByPrincipalPlane(offset=0.0, 
    principalPlane=XYPLANE)
mdb.models['Model-1'].parts[Temp_Part].DatumPlaneByRotation(angle= (N_Rev-math.floor(N_Rev))*360, axis=
    mdb.models['Model-1'].parts[Temp_Part].datums[1], plane=
    mdb.models['Model-1'].parts[Temp_Part].datums[2])
mdb.models['Model-1'].ConstrainedSketch(gridSpacing=0.08, name='__profile__', 
    sheetSize=3.3, transform=
    mdb.models['Model-1'].parts[Temp_Part].MakeSketchTransform(
    sketchPlane=mdb.models['Model-1'].parts[Temp_Part].datums[2], 
    sketchPlaneSide=SIDE1, 
    sketchUpEdge=mdb.models['Model-1'].parts[Temp_Part].datums[1], 
    sketchOrientation=RIGHT, origin=(0,0,0)))
mdb.models['Model-1'].parts[Temp_Part].projectReferencesOntoSketch(filter=
    COPLANAR_EDGES, sketch=mdb.models['Model-1'].sketches['__profile__'])
mdb.models['Model-1'].sketches['__profile__'].CircleByCenterPerimeter(center=(
    -Rmid , 0.0), point1=(-R , 0.0))
mdb.models['Model-1'].parts[Temp_Part].Shell(sketch=
    mdb.models['Model-1'].sketches['__profile__'], sketchOrientation=RIGHT, 
    sketchPlane=mdb.models['Model-1'].parts[Temp_Part].datums[2], 
    sketchPlaneSide=SIDE1, sketchUpEdge=
    mdb.models['Model-1'].parts[Temp_Part].datums[1])
del mdb.models['Model-1'].sketches['__profile__']

mdb.models['Model-1'].ConstrainedSketch(gridSpacing=0.08, name='__profile__', 
    sheetSize=3.3, transform=
    mdb.models['Model-1'].parts[Temp_Part].MakeSketchTransform(
    sketchPlane=mdb.models['Model-1'].parts[Temp_Part].datums[3], 
    sketchPlaneSide=SIDE1, 
    sketchUpEdge=mdb.models['Model-1'].parts[Temp_Part].datums[1], 
    sketchOrientation=RIGHT, origin=(0, -1*Pitch*N_Rev , 0)))
mdb.models['Model-1'].parts[Temp_Part].projectReferencesOntoSketch(filter=
    COPLANAR_EDGES, sketch=mdb.models['Model-1'].sketches['__profile__'])
mdb.models['Model-1'].sketches['__profile__'].CircleByCenterPerimeter(center=(
    -Rmid , 0.0), point1=(-R , 0.0))
mdb.models['Model-1'].parts[Temp_Part].Shell(sketch=
    mdb.models['Model-1'].sketches['__profile__'], sketchOrientation=RIGHT, 
    sketchPlane=mdb.models['Model-1'].parts[Temp_Part].datums[3], 
    sketchPlaneSide=SIDE1, sketchUpEdge=
    mdb.models['Model-1'].parts[Temp_Part].datums[1])
del mdb.models['Model-1'].sketches['__profile__']

#Create Solid
for i in range(n_partition):
    mdb.models['Model-1'].parts[Temp_Part].AddCells(faceList=(
        mdb.models['Model-1'].parts[Temp_Part].faces.findAt((  ((Rmid+r_parts_temp[i])+(Rmid+r_parts_temp[i+1])) /2, 0, 0.0), 
        ), ), flipped=False)
        



#Sets
set_list = []    
for i in range(n_partition):
    mdb.models['Model-1'].parts[Temp_Part].Set(name='Cell'+str(i+1),cells=mdb.models['Model-1'].parts[Temp_Part].cells.findAt(((  ((Rmid+r_parts_temp[i])+(Rmid+r_parts_temp[i+1])) /2, 0, 0.0),),))
    set_list.append(mdb.models['Model-1'].parts[Temp_Part].sets['Cell'+str(i+1)])
mdb.models['Model-1'].parts[Temp_Part].SetByBoolean(name='All_Cells', sets=tuple(set_list))

#Mesh
mdb.models['Model-1'].parts[Temp_Part].seedPart(deviationFactor=0.1, 
    minSizeFactor=0.1, size=ds/10)
    
for i in range(n_partition):
    mdb.models['Model-1'].parts['TempPart'].setMeshControls(elemShape=TET, regions=
        mdb.models['Model-1'].parts['TempPart'].cells.findAt(((
        ((Rmid+r_parts_temp[i])+(Rmid+r_parts_temp[i+1])) /2, 0, 0.0),),), technique=FREE)

mdb.models['Model-1'].parts[Temp_Part].generateMesh()

for i in range(n_partition):
    mdb.models['Model-1'].parts[Temp_Part].setElementType(elemTypes=(ElemType(
        elemCode=C3D20R, elemLibrary=STANDARD), ElemType(elemCode=C3D15, 
        elemLibrary=STANDARD), ElemType(elemCode=C3D10, elemLibrary=STANDARD, 
        secondOrderAccuracy=OFF, distortionControl=DEFAULT)), regions=(
        mdb.models['Model-1'].parts[Temp_Part].cells.findAt(((
        ((Rmid+r_parts_temp[i])+(Rmid+r_parts_temp[i+1])) /2, 0, 0.0),),), ))
        
mdb.models['Model-1'].rootAssembly.regenerate()





# #Create Scaled Part
# mdb.models['Model-1'].Part(compressFeatureList=ON, name= Use_Part, 
    # objectToCopy=mdb.models['Model-1'].parts[Temp_Part], scale=0.001)
    
# D = D*1e-3
# ds = ds*1e-3
# H = H*1e-3
# L = L*1e-3
# Pitch = H/N    
# r_parts = np.array(r_parts)*1e-3
# r_parts_temp = np.array(r_parts_temp)*1e-3
# Dmid = D-ds
# R = D/2.0
# rs = ds/2.0
# Rmid = Dmid/2.0
# H_model = Pitch*N_Rev

# #Sets
# set_list = []    
# for i in range(n_partition):
    # mdb.models['Model-1'].parts[Use_Part].Set(name='Cell'+str(i+1),cells=mdb.models['Model-1'].parts[Use_Part].cells.findAt(((  ((Rmid+r_parts_temp[i])+(Rmid+r_parts_temp[i+1])) /2, 0, 0.0),),))
    # set_list.append(mdb.models['Model-1'].parts[Use_Part].sets['Cell'+str(i+1)])
# mdb.models['Model-1'].parts[Use_Part].SetByBoolean(name='All_Cells', sets=tuple(set_list))

# #Mesh
# mdb.models['Model-1'].parts[Temp_Part].seedPart(deviationFactor=0.1, 
    # minSizeFactor=0.1, size=0.05e-3)
# mdb.models['Model-1'].parts[Use_Part].generateMesh()
# mdb.models['Model-1'].rootAssembly.regenerate()



