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
session.journalOptions.setValues(replayGeometry=COORDINATE, recoverGeometry=COORDINATE)

###
# This File Creates a Helix and sets dependant on the number of concentric parts needed.
# Sets named P1, P2, P3... represent the the different parts from the inner most to the outer most.
# All Units in SI
###

#Max Diameter (Measured experimentally
D = 0.75e-3

#Diameter of the nylon wire
ds = 0.3e-3

#Diameter of the helix
Dmid = D-ds

R = D/2.0
rs = ds/2.0
Rmid = Dmid/2.0


#Number of parts
N_parts = 3;

#Max Radius of parts
r_parts = [rs/5, rs/2, rs]; #ideally, the last part has a max radius of the wire

#Number of turns of the coil
N_Rev = 0.5 #Do not make this too large, it will take a very long time to solve


ModelName = 'TheMesh'
jobName = ModelName
outputDispName = ModelName
i = 0.33
w = 30e-3

D = 0.75e-3
ds = 0.3e-3
Dmid = D-ds

R = D/2.0
rs = ds/2.0
Rmid = Dmid/2.0


H = 150.0e-3
L = 500.0e-3
N = 280.0
#N = 17.0*H*100
Pitch = H/N

H_model = Pitch*N_Rev

h_conv = 50


L_Full = math.sqrt(H**2 + (N**2)*((math.pi)**2)*(D**2))
Vol_full = math.pi*(ds**2)*L_Full/4
#Vol_full = math.pi*(ds**2)*(math.pi*D*N)/4

mass_full = 0.09*1e-3*L_Full
rho = mass_full/Vol_full

R_ini = 14
Q = (i**2) * R_ini
q_model = Q/Vol_full

#L_ini = math.pi*D*N
A_ini = math.pi*(ds**2)/4

surf_i_model = i/A_ini

# Te = np.linspace(0, 200);
# T2 = Te+273.15

# rho_ini = R_ini*A_ini/L_Full
# #rho_ini = R_ini*A_ini/L_ini

# alpha_T = 0.003
# rho_T = rho_ini*(1+alpha_T*(T2-298))
# con_T = 1/rho_T

# # w = 30e-3
# w_model = w*9.81*H_model/H

# p_temp = (74.13*(w**2)-9.53*w+1.247)
# temp = (1/(1+np.exp(0.027*(Te-120))))
# E = (p_temp*(5.79e9*(1+temp)))

# poi = 0.4*np.ones_like(E)
# alpha22 = -1.3e-4*(1-np.cos(math.pi*(Te-23)/120));
# alpha1133 = 1.5e-4 + alpha22;

# arr_E = np.array((E,poi,T2)).T
# arr_alpha = np.array((alpha1133,alpha22,alpha1133,T2)).T
# #arr_alpha = np.array((alpha22,alpha22*1e-9,alpha22,T2)).T
# arr_con_T = np.array((con_T,T2)).T

# E_use = tuple([tuple(e) for e in arr_E])
# alpha_use = tuple([tuple(e) for e in arr_alpha])
# con_T_use = tuple([tuple(e) for e in arr_con_T])


mdb.Model(modelType=STANDARD_EXPLICIT, name=ModelName)


#Create Model

mdb.models[ModelName].ConstrainedSketch(name='__profile__', sheetSize=200.0)
mdb.models[ModelName].sketches['__profile__'].ConstructionLine(point1=(0.0, 
    -100.0), point2=(0.0, 100.0))
mdb.models[ModelName].sketches['__profile__'].FixedConstraint(entity=
    mdb.models[ModelName].sketches['__profile__'].geometry[2])
mdb.models[ModelName].sketches['__profile__'].CircleByCenterPerimeter(center=(
    Rmid*1e3, 0.0), point1=(R*1e3, 0.0))
mdb.models[ModelName].Part(dimensionality=THREE_D, name='Helix_Temp', type=
    DEFORMABLE_BODY)
mdb.models[ModelName].parts['Helix_Temp'].BaseSolidRevolve(angle=360*N_Rev, 
    flipPitchDirection=ON, flipRevolveDirection=OFF, moveSketchNormalToPath=ON, 
    pitch=Pitch*1e3, sketch=mdb.models[ModelName].sketches['__profile__'])
del mdb.models[ModelName].sketches['__profile__']



mdb.models[ModelName].Part(compressFeatureList=ON, name='Helix', 
    objectToCopy=mdb.models[ModelName].parts['Helix_Temp'], scale=0.001)

##
del mdb.models[ModelName].parts['Helix_Temp']
##


#Sets
mdb.models[ModelName].parts['Helix'].Set(name='cell_all',cells=mdb.models[ModelName].parts['Helix'].cells.findAt(((Rmid,0,0),),))

#mdb.models[ModelName].parts['Helix'].Set(elements=mdb.models[ModelName].parts['Helix'].sets['cell_all'].elements, name='All')

mdb.models[ModelName].parts['Helix'].Set(name='bottom_face',faces=mdb.models[ModelName].parts['Helix'].faces.findAt(((Rmid,0,0),),))

mdb.models[ModelName].parts['Helix'].Set(name='top_face',faces=mdb.models[ModelName].parts['Helix'].faces.findAt(((Rmid,H_model,0),),))

#mdb.models[ModelName].parts['Helix'].Set(elements=mdb.models[ModelName].parts['Helix'].sets['bottom_face'].elements, name='face1')

#Faces
#ConvectionSurface
mdb.models[ModelName].parts['Helix'].Surface(side1Faces=mdb.models[ModelName].parts['Helix'].faces.getByBoundingCylinder((0,-1*Pitch,0),(0,(Pitch*(N_Rev+1)),0),D+0.1e-3), name='All_Faces')
#Bottom Surface
mdb.models[ModelName].parts['Helix'].Surface(name='bottom_face', side1Faces= mdb.models[ModelName].parts['Helix'].faces.findAt(((Rmid, 0, 0), )))
#Top Surface
mdb.models[ModelName].parts['Helix'].Surface(name='top_face', side1Faces= mdb.models[ModelName].parts['Helix'].faces.findAt(((Rmid, H_model, 0), )))


#Mesh Model
mdb.models[ModelName].parts['Helix'].seedPart(deviationFactor=0.1, 
    minSizeFactor=0.1, size=D/40)
mdb.models[ModelName].parts['Helix'].setMeshControls(algorithm=MEDIAL_AXIS, regions = mdb.models[ModelName].parts['Helix'].cells.findAt(((Rmid,0,0),),))
mdb.models[ModelName].parts['Helix'].generateMesh()

#MeshElementType
mdb.models[ModelName].parts['Helix'].setElementType(elemTypes=(ElemType(
    elemCode=DC3D8E, elemLibrary=STANDARD), ElemType(elemCode=DC3D6E, elemLibrary=STANDARD), 
    ElemType(elemCode=DC3D4E, elemLibrary=STANDARD)), regions=(
    mdb.models[ModelName].parts['Helix'].cells.findAt(((Rmid,0,0), )), ))
    
mdb.models[ModelName].parts['Helix'].Set(nodes=mdb.models[ModelName].parts['Helix'].nodes.getByBoundingCylinder((0,Pitch*(N_Rev/2 - 0.1),0),(0,Pitch*(N_Rev/2+0.1),0),D+0.1e-3), name='Mid')
mdb.models[ModelName].parts['Helix'].Set(elements=mdb.models[ModelName].parts['Helix'].elements.getByBoundingCylinder((0,Pitch*(-0.1),0),(0,Pitch*(N_Rev+0.1),0),D+0.1e-3), name='Helix_Elements')



# t = np.linspace(0, N_Rev*2*np.pi,100*N_Rev);
r = Rmid;
pitch = Pitch;
h = pitch/(2*np.pi);

x = r*np.cos(t);
z = r*np.sin(t);
y = h*t;


#Finding the centroids of all elements
a = mdb.models[ModelName].parts['Helix']
elements = a.elements
centroid_list = []
for i, el in enumerate(elements):
    # print(el)
    region = regionToolset.Region(elements=elements[i:i+1])
    properties = a.getMassProperties(regions=region)
    centroid_list.append(list(properties['volumeCentroid']))
# print(centroid_list)




#Determining which part an element belongs to

r_parts = [0] + r_parts

for part in range(N_parts):

    lst = [];
    for elemid in range(1,len(centroid_list)+1):

        xyz = centroid_list[elemid-1]
        x1 = xyz[0]
        y1 = xyz[1]
        z1 = xyz[2]   
        m = 0;
        L = np.ones([11,4])*1000;
        for gue in np.arange(11): 
            def needed (t1,r,h,x1,y1,z1):
                y = x1*r*np.sin(t1)-z1*r*np.cos(t1)+h*h*t1-y1*h;
                return y        
            sol1 = fsolve(needed, gue, (r,h,x1,y1,z1))
            t1 = sol1;
            #if (x1*r*cos(t1)+y1*r*sin(t1)+h*h) > 0    
            x2 = r*np.cos(t1);
            z2 = r*np.sin(t1);
            y2 = h*t1;
            L[m,:] = [np.sqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2), x2, y2, z2];
            #plot3([x1 x2],[y1 y2],[z1 z2],':*','LineWidth',0.1);
            m = m+1;
            #end
        L = L[:,0].tolist()
        Lmin = min(L)
        print(Lmin)
        if Lmin > r_parts[part] and Lmin < r_parts[part+1]:
            lst.append(elemid)
        # else:
            # lst2.append(elemid)
        # k1 = L[:,0];
        # k = k1.tolist();
        # i = k.index(min(k));
        
    mdb.models[ModelName].parts['Helix'].SetFromElementLabels(name='P'+str(part+1), elementLabels= lst)
    # mdb.models[ModelName].parts['Helix'].SetFromElementLabels(name='P2', elementLabels= lst2)







#Orientation
# 1- radial out, 2- follows path, 3- vertically up

##Cylindrical Coordinate
mdb.models[ModelName].parts['Helix'].DatumCsysByThreePoints(coordSysType=
    CYLINDRICAL, name='Datum csys-1', origin=(0.0, 0.0, 0.0), point1=(0.0, 0.0, 
    1.0), point2=(1.0, 0.0, 1.0))
    
DatumId = mdb.models[ModelName].parts['Helix'].features['Datum csys-1'].id

mdb.models[ModelName].parts['Helix'].MaterialOrientation(
    additionalRotationField='', additionalRotationType=ROTATION_NONE, angle=0.0
    , axis=AXIS_3, fieldName='', localCsys=
    mdb.models[ModelName].parts['Helix'].datums[DatumId], orientationType=SYSTEM, 
    region=mdb.models[ModelName].parts['Helix'].sets['cell_all'], 
    stackDirection=STACK_ORIENTATION)



del mdb.models['Model-1']

##########################
# #####
    # mdb.models[ModelName].parts['Helix'].Surface(name='face', side1Faces=
        # mdb.models[ModelName].parts['Helix'].faces.findAt(((Rmid, 0, 0), )))


    # t = mdb.models[ModelName].parts['Helix'].sets['face1'].elements

    # print(t)
    # print(t[-1].label)
# ######

# mat1 = [];
# mat2 = [];

# for face in range(1,456):

    # lst1=list(range(1,32))
    # lst2=[]
    # for k in range(1,16):
        # temp1 = random.choice(lst1)
        # lst2.append(temp1)
        # lst1.remove(temp1)

    # nplst1 = ((face-1)*31)+np.array(lst1)  
    # for val in nplst1.tolist():
        # mat1.append(val)
    # nplst2 = ((face-1)*31)+np.array(lst2)    
    # for val in nplst2.tolist():
        # mat2.append(val)
    
# mdb.models['Model-1'].parts['Helix_2017_02'].SetFromElementLabels(name='Material1', elementLabels= mat1)
# mdb.models['Model-1'].parts['Helix_2017_02'].SetFromElementLabels(name='Material2', elementLabels= mat2)

# #########

# a = mdb.models['Model-1'].rootAssembly
# elements1 = element/selection of elements
# a.getMassProperties(regions=elements)["volumeCentroid"]

# ####
# import regionToolset
# a = mdb.models['Model-1'].rootAssembly
# elements = a.instances['Part-1-1'].elements
# centroid_list = []
# for i, el in enumerate(elements):
# print(el)
# region = regionToolset.Region(elements=elements[i:i+1])
# properties = a.getMassProperties(regions=region)
# centroid_list.append(list(properties['volumeCentroid']))
# print(centroid_list)



