# -*- coding: utf-8 -*-
"""
Created on Sun Jul 16 14:40:46 2023

@author: Kirk Ming Yeoh (e0546208@u.nus.edu)
"""

import os
import numpy as np
import math

### Update the following parameters
os.chdir('E:\\Kirk Ming Abaqus\\DFE2_Demo') # Directory where input files are, use double '\'
MacroInpName = 'Demo_3D_Macro.inp' # Name of macroscale input file
RVEInpName = 'Demo_3D_RVE.inp' # Name of RVE input file 
NewInpName = 'DFE2_3D.inp' # Name of new Direct FE2 input file


GP = [[-3**-0.5,-3**-0.5,-3**-0.5],[3**-0.5,-3**-0.5,-3**-0.5],[3**-0.5,3**-0.5,-3**-0.5],[-3**-0.5,3**-0.5,-3**-0.5],[-3**-0.5,-3**-0.5,3**-0.5],[3**-0.5,-3**-0.5,3**-0.5],[3**-0.5,3**-0.5,3**-0.5],[-3**-0.5,3**-0.5,3**-0.5]]

### Define functions
# Search inp array for a particular part's elements and nodes
def Search(inp,key1,out1,type1): #if dealine with nodes and coordinates, set type1 to 1 
    for i in range(len(inp)):
        if inp[i].count(key1)!=0:
            break    
    for temp_line in inp[i+1:]:
        if (temp_line == '') or (temp_line.count("*") != 0):
            break
        temp_line = (temp_line.replace(',',' ')).split()
        temp_line.pop(0) #Removes element number in connectivity and node number in nodal coordinates
        for i in range(len(temp_line)):
            if type1 == 1:
                temp_line[i] = float(temp_line[i])
            else:
                temp_line[i] = int(temp_line[i])-1 #All node/element labels are -1 for use in Python as indices
        out1.append(temp_line)   

# Removes corner nodes from an edge        
def TakeVertexOut(face):
    face.pop(0)
    face.pop(-1)
    return face   

# Calculates trilinear shape function values
def Trilin_Interpolation(tsi,eta,zeta):
    N1=float(0.125*(1-tsi)*(1-eta)*(1-zeta))
    N2=float(0.125*(1+tsi)*(1-eta)*(1-zeta))
    N3=float(0.125*(1+tsi)*(1+eta)*(1-zeta))
    N4=float(0.125*(1-tsi)*(1+eta)*(1-zeta))
    N5=float(0.125*(1-tsi)*(1-eta)*(1+zeta))
    N6=float(0.125*(1+tsi)*(1-eta)*(1+zeta))
    N7=float(0.125*(1+tsi)*(1+eta)*(1+zeta))
    N8=float(0.125*(1-tsi)*(1+eta)*(1+zeta))   
    return [N1,N2,N3,N4,N5,N6,N7,N8]

# Sorts nodes along an edge using their coordinates
def SortListofNodes1D(faceN,coordinate): # 0 for x; 1 for y; 2 for z; gives node numbers ready to be called with python (already -1)
    newlist = []
    oldlist = []
    for i in range(len(faceN)):
        oldlist.append(RVENodalCoord[faceN[i]][coordinate])
    
    orderedlist = sorted(oldlist)
    for j in range(len(orderedlist)):
        ind = oldlist.index(orderedlist[j])
        newlist.append(faceN[ind])
    
    return newlist

# Sorts nodes along a face using their coordinates
def SortListofNodes2D(faceN,coord1,coord2): # 0 for x; 1 for y; 2 for z; gives node numbers ready to be called with python (already -1)
    oldlistC1 = [] #first coordinate to be sorted along
    newlistN = []
    
    #Sort by first coodinate
    for i in range(len(faceN)):
        oldlistC1.append(RVENodalCoord[faceN[i]][coord1])
        
    newlistC1 = sorted(list(dict.fromkeys(oldlistC1)))
        
    for i in range(len(newlistC1)):
        C1 = newlistC1[i]
        sublistN = []
        sublistC2 = []
        for j in range(len(faceN)):
            C1N = RVENodalCoord[faceN[j]][coord1]
            
            if (C1N==C1):
                sublistN.append(faceN[j])
                sublistC2.append(RVENodalCoord[faceN[j]][coord2])
                
        newlistC2 = sorted(sublistC2)
        for j in range(len(sublistN)):
            Nindex = sublistC2.index(newlistC2[j])
            newlistN.append(sublistN[Nindex]) 
    
    return newlistN

# Removes nodes from a set that matches the given coordinates
def ExcludeNodes(faceN,coord1,coord2,coord3): # coord1, coord2 and coord3 are the coordinates to exclude, to be given as lists
    newlistN = []
    
    for i in range(len(faceN)):
        if RVENodalCoord[faceN[i]][0] not in coord1:
            if RVENodalCoord[faceN[i]][1] not in coord2:
                if RVENodalCoord[faceN[i]][2] not in coord3:
                    newlistN.append(faceN[i])
    
    return newlistN        


### Extracting information from Macro and RVE input files
# Macroscale input file
inp1 = []    
f1 = open(MacroInpName,'r')

while 1:
    line = f1.readline()
    if not line:
        break
    line = line.strip() #removes additional white spaces on left and right
    inp1.append(line)
    
f1.close()

# Removing 'generate' for easier processing
for i in reversed(range(len(inp1))):
    if (inp1[i].count('generate')!=0):
        Temp = (inp1[i+1].replace(',',' ')).split()
        k = 0 #term counter
        m = 0 #extra line counter
        for j in range(int(Temp[0]),int(Temp[1])+1,int(Temp[2])):
            if k==0:
                Temp2 = str(j)
                k = k+1
            elif k==16:
                inp1.insert(i+2+m,Temp2)
                m = m+1
                Temp2 = str(j)
                k = 1
            else:
                Temp2 = Temp2+', '+str(j)
                k = k+1
        inp1.insert(i+2+m,Temp2)
        inp1[i] = inp1[i][0:len(inp1[i])-10]
        del inp1[i+1]

# RVE input file
inp2 = []    
f2 = open(RVEInpName,'r')

while 1:
    line = f2.readline()
    if not line:
        break
    line = line.strip() #removes additional white spaces on left and right
    inp2.append(line)
    
f2.close()

# Removing 'generate' for easier processing
for i in reversed(range(len(inp2))):
    if (inp2[i].count('generate')!=0):
        Temp = (inp2[i+1].replace(',',' ')).split()
        k = 0 #term counter
        m = 0 #extra line counter
        for j in range(int(Temp[0]),int(Temp[1])+1,int(Temp[2])):
            if k==0:
                Temp2 = str(j)
                k = k+1
            elif k==16:
                inp2.insert(i+2+m,Temp2)
                m = m+1
                Temp2 = str(j)
                k = 1
            else:
                Temp2 = Temp2+', '+str(j)
                k = k+1
        inp2.insert(i+2+m,Temp2)
        inp2[i] = inp2[i][0:len(inp2[i])-10]
        del inp2[i+1]

# Extracting macroscale element info from old inp file
MacroNodalConnect,MacroNodalCoord = [],[]
Search(inp1,'*Element',MacroNodalConnect,0)
Search(inp1,'*Node',MacroNodalCoord,1)

StartConst = 'a'
for i in range(len(inp1)):
    if (inp1[i].count('*Instance,'))!=0:
        Line = inp1[i].split(', ')
        MacroInstName = Line[1][5:]
    if (inp1[i].count('*Material,'))!=0:
        inp1[i+2] = '1e-10,1e-10'
    if (inp1[i].count('** Constraint'))!=0 and (StartConst == 'a'):
        StartConst = i

# Extracting RVE info from old inp file
RVENodalCoord = []
Search(inp2,'*Node',RVENodalCoord,1)

Sections = []
Materials = []
StartEle = 'a'
for i in range(len(inp2)):
    if (inp2[i].count('*Element'))!=0 and (StartEle == 'a'):
        StartEle = i
    if (inp2[i].count('** Section'))!=0:
        Sections.append(i)
    if (inp2[i].count('*Material'))!=0:
        Materials.append(i)

RVEMats = open('RVEMats.dat','w')
for i in range(len(Materials)):
    for j in range(Materials[i]+1,len(inp2)):
        if (inp2[j].count('*Material'))!=0 or (inp2[j].count('**'))!=0:
            MatEnd = j
            break
        MatEnd = j+1 # If no further information provided beyond RVE materials
    for j in range(Materials[i],MatEnd):
        print>>RVEMats,inp2[j]
RVEMats.close()


### Processing the macroscale part information
# Sorting the nodal connectivity to match with DFE2 conventions
Nele = len(MacroNodalConnect)
NodalConnect = []
NodalCoordX = []
NodalCoordY = []
NodalCoordZ = []
for i in range(Nele):
    Nodes = []
    TempCoordX = []
    TempCoordY = []
    TempCoordZ = []
    Disp = [[] for x in range(8)]
    
    # Obtaining nodal information from this particular element
    for j in range(8):
        Nodes.append(MacroNodalConnect[i][j])
        TempCoordX.append(MacroNodalCoord[MacroNodalConnect[i][j]][0])
        TempCoordY.append(MacroNodalCoord[MacroNodalConnect[i][j]][1])
        TempCoordZ.append(MacroNodalCoord[MacroNodalConnect[i][j]][2])
    
    # Obtaining centroid of the element    
    X0 = sum(TempCoordX)/8
    Y0 = sum(TempCoordY)/8
    Z0 = sum(TempCoordZ)/8
    
    # Obtaining displacement vector of each node relative to centroid
    for j in range(8):
        Disp[j].append(MacroNodalCoord[MacroNodalConnect[i][j]][0]-X0)
        Disp[j].append(MacroNodalCoord[MacroNodalConnect[i][j]][1]-Y0)
        Disp[j].append(MacroNodalCoord[MacroNodalConnect[i][j]][2]-Z0)
        
    Order = [0 for x in range(8)]
    Set = [] # Nodes 1-4 with \zeta < 0
    for j in range(8):
        if Disp[j][2]<0:
            Set.append(j)            
    for j in range(4):
        if Disp[Set[j]][1]<0:
            if Disp[Set[j]][0]<0:
                Order[0] = Set[j]
            else:
                Order[1] = Set[j]
        elif Disp[Set[j]][1]>0:
            if Disp[Set[j]][0]<0:
                Order[3] = Set[j]
            else:
                Order[2] = Set[j]
                
    Set = [] # Nodes 5-8 with \zeta > 0
    for j in range(8):
        if Disp[j][2]>0:
            Set.append(j)            
    for j in range(4):
        if Disp[Set[j]][1]<0:
            if Disp[Set[j]][0]<0:
                Order[4] = Set[j]
            else:
                Order[5] = Set[j]
        elif Disp[Set[j]][1]>0:
            if Disp[Set[j]][0]<0:
                Order[7] = Set[j]
            else:
                Order[6] = Set[j]
                
    SortedConnect = []
    SortedCoordX = []
    SortedCoordY = []
    SortedCoordZ = []
    
    for j in range(8):
        SortedConnect.append(Nodes[Order[j]])
        SortedCoordX.append(MacroNodalCoord[Nodes[Order[j]]][0])
        SortedCoordY.append(MacroNodalCoord[Nodes[Order[j]]][1])
        SortedCoordZ.append(MacroNodalCoord[Nodes[Order[j]]][2])
        
    NodalConnect.append(SortedConnect)
    NodalCoordX.append(SortedCoordX)
    NodalCoordY.append(SortedCoordY)
    NodalCoordZ.append(SortedCoordZ)


### Processing the RVE part information
# Finding the smallest and largest nodal coordinate in all 3 directions
RVE_ListX = []
RVE_ListY = []
RVE_ListZ = []
for i in range(len(RVENodalCoord)):
    RVE_ListX.append(RVENodalCoord[i][0])
    RVE_ListY.append(RVENodalCoord[i][1])
    RVE_ListZ.append(RVENodalCoord[i][2])
xMin = min(RVE_ListX)
xMax = max(RVE_ListX)
yMin = min(RVE_ListY)
yMax = max(RVE_ListY)
zMin = min(RVE_ListZ)
zMax = max(RVE_ListZ)
del RVE_ListX
del RVE_ListY
del RVE_ListZ

# Sorting the RVE boundary nodes
FaceLNodes,FaceRNodes,FaceBaNodes,FaceFNodes,FaceBNodes,FaceTNodes = [],[],[],[],[],[]
EdgeLNodes, EdgeRNodes, EdgeBNodes, EdgeTNodes = [],[],[],[]
Edge1Nodes, Edge2Nodes, Edge3Nodes, Edge4Nodes = [],[],[],[]
for i in range(len(RVENodalCoord)):
    if RVENodalCoord[i][0] == xMin:
        FaceLNodes.append(i)
    if RVENodalCoord[i][0] == xMax:
        FaceRNodes.append(i)
    if RVENodalCoord[i][1] == yMin:
        FaceBaNodes.append(i)
    if RVENodalCoord[i][1] == yMax:
        FaceFNodes.append(i)
    if RVENodalCoord[i][2] == zMin:
        FaceBNodes.append(i)
    if RVENodalCoord[i][2] == zMax:
        FaceTNodes.append(i)
for i in range(len(FaceBaNodes)):
    if FaceBaNodes[i] in FaceLNodes:
        EdgeLNodes.append(FaceBaNodes[i])
        if FaceBaNodes[i] in FaceBNodes:
            V1 = FaceBaNodes[i]
        if FaceBaNodes[i] in FaceTNodes:
            V4 = FaceBaNodes[i]
    if FaceBaNodes[i] in FaceRNodes:
        EdgeRNodes.append(FaceBaNodes[i])
        if FaceBaNodes[i] in FaceBNodes:
            V2 = FaceBaNodes[i]
        if FaceBaNodes[i] in FaceTNodes:
            V3 = FaceBaNodes[i]
    if FaceBaNodes[i] in FaceBNodes:
        EdgeBNodes.append(FaceBaNodes[i])
    if FaceBaNodes[i] in FaceTNodes:
        EdgeTNodes.append(FaceBaNodes[i])
for i in range(len(FaceBNodes)):
    if FaceBNodes[i] in FaceLNodes:
        Edge1Nodes.append(FaceBNodes[i])
        if FaceBNodes[i] in FaceFNodes:
            V5 = FaceBNodes[i]
    if FaceBNodes[i] in FaceRNodes:
        Edge2Nodes.append(FaceBNodes[i])
        if FaceBNodes[i] in FaceFNodes:
            V6 = FaceBNodes[i]
for i in range(len(FaceTNodes)):
    if FaceTNodes[i] in FaceLNodes:
        Edge4Nodes.append(FaceTNodes[i])
        if FaceTNodes[i] in FaceFNodes:
            V8 = FaceTNodes[i]
    if FaceTNodes[i] in FaceRNodes:
        Edge3Nodes.append(FaceTNodes[i])
        if FaceTNodes[i] in FaceFNodes:
            V7 = FaceTNodes[i]
    
# Sorting RVE dimensions and offsets
B_RVE = xMax - xMin
H_RVE = yMax - yMin
T_RVE = zMax - zMin
OffsetX = (xMax + xMin)/2
OffsetY = (yMax + yMin)/2
OffsetZ = (zMax + zMin)/2
Offset = [OffsetX,OffsetY,OffsetZ]

# Adjusting RVE nodal coordinates to correspond to a part centered at the origin
for i in RVENodalCoord:
    for j in range(3):
        i[j] = i[j] - Offset[j]


### Generating the RVE placement in the macroscale mesh
RVEParts = open('RVEParts.dat','w')
Insts = open('Insts.dat','w')
SF = []
for i in range(Nele):
    C = np.array([
            [1,-1,-1,-1,1,1,1,-1],
            [1,1,-1,-1,-1,-1,1,1],
            [1,1,1,-1,1,-1,-1,-1],
            [1,-1,1,-1,-1,1,-1,1],
            [1,-1,-1,1,1,-1,-1,1],
            [1,1,-1,1,-1,1,-1,-1],
            [1,1,1,1,1,1,1,1],
            [1,-1,1,1,-1,-1,1,-1],])
    C_inv = np.linalg.inv(C)
    [a0,a1,a2,a3,a4,a5,a6,a7] = np.dot(C_inv,NodalCoordX[i])
    [b0,b1,b2,b3,b4,b5,b6,b7] = np.dot(C_inv,NodalCoordY[i])
    [d0,d1,d2,d3,d4,d5,d6,d7] = np.dot(C_inv,NodalCoordZ[i])
    
    for j in range(8):
        [tsi,eta,zeta] = GP[j]
        N1,N2,N3,N4,N5,N6,N7,N8 = Trilin_Interpolation(tsi,eta,zeta)
        
        RVE_X = N1*NodalCoordX[i][0] + N2*NodalCoordX[i][1] + N3*NodalCoordX[i][2] + N4*NodalCoordX[i][3] + N5*NodalCoordX[i][4] + N6*NodalCoordX[i][5] + N7*NodalCoordX[i][6] + N8*NodalCoordX[i][7]
        RVE_Y = N1*NodalCoordY[i][0] + N2*NodalCoordY[i][1] + N3*NodalCoordY[i][2] + N4*NodalCoordY[i][3] + N5*NodalCoordY[i][4] + N6*NodalCoordY[i][5] + N7*NodalCoordY[i][6] + N8*NodalCoordY[i][7]
        RVE_Z = N1*NodalCoordZ[i][0] + N2*NodalCoordZ[i][1] + N3*NodalCoordZ[i][2] + N4*NodalCoordZ[i][3] + N5*NodalCoordZ[i][4] + N6*NodalCoordZ[i][5] + N7*NodalCoordZ[i][6] + N8*NodalCoordZ[i][7]
        
        J = np.array([
                [a1+a4*eta+a5*zeta+a7*eta*zeta,b1+b4*eta+b5*zeta+b7*eta*zeta,d1+d4*eta+d5*zeta+d7*eta*zeta],
                [a2+a4*tsi+a6*zeta+a7*tsi*zeta,b2+b4*tsi+b6*zeta+b7*tsi*zeta,d2+d4*tsi+d6*zeta+d7*tsi*zeta],
                [a3+a5*tsi+a6*eta+a7*tsi*eta,b3+b5*tsi+b6*eta+b7*tsi*eta,d3+d5*tsi+d6*eta+d7*tsi*eta]])
       
        J_RVE = (abs(np.linalg.det(J)/(B_RVE*H_RVE*T_RVE)))**(1.0/3.0)
        
        if round(J_RVE,5) in SF: # If such an RVE was created previously, reuse it for the next instance
            N = SF.index(round(J_RVE,5))
        else: # If not, then we create another RVE with this current thickness scaling
            N = len(SF)
            SF.append(round(J_RVE,5))
            print>>RVEParts,'**'
            print>>RVEParts,'*Part, name=RVE-'+str(N+1)
            print>>RVEParts,'*Node'
            
            for k in range(len(RVENodalCoord)):
                print>>RVEParts,str(k+1)+', '+str(RVENodalCoord[k][0]*J_RVE)+', '+str(RVENodalCoord[k][1]*J_RVE)+', '+str(RVENodalCoord[k][2]*J_RVE)
                
            for k in range(StartEle,Sections[0]):
                print>>RVEParts,inp2[k]
                
            for k in range(len(Sections)):
                print>>RVEParts,inp2[Sections[k]]+'-'+str(N+1)
                print>>RVEParts,inp2[Sections[k]+1]
                print>>RVEParts,','
                
            print>>RVEParts,'*End Part'
            
        print>>Insts,'*Instance, name=Ele'+str(i+1)+'-RVE'+str(j+1)+', part=RVE-'+str(N+1)
        print>>Insts,str(RVE_X)+', '+str(RVE_Y)+', '+str(RVE_Z)
        print>>Insts,'*End Instance'
        print>>Insts,'**'

RVEParts.close()
Insts.close()


### Setting up the MPCs
Sets = open('Sets.dat','w')
Eqns = open('Eqns.dat','w')

# Pairing the nodes
EdgeLNodes = TakeVertexOut(SortListofNodes1D(EdgeLNodes,2))
EdgeRNodes = TakeVertexOut(SortListofNodes1D(EdgeRNodes,2))
    
EdgeBNodes = TakeVertexOut(SortListofNodes1D(EdgeBNodes,0))
EdgeTNodes = TakeVertexOut(SortListofNodes1D(EdgeTNodes,0))
    
Edge1Nodes = TakeVertexOut(SortListofNodes1D(Edge1Nodes,1))
Edge2Nodes = TakeVertexOut(SortListofNodes1D(Edge2Nodes,1))
Edge3Nodes = TakeVertexOut(SortListofNodes1D(Edge3Nodes,1))
Edge4Nodes = TakeVertexOut(SortListofNodes1D(Edge4Nodes,1))

FaceBaNodes = SortListofNodes2D(FaceBaNodes,0,2)
FaceFNodes = SortListofNodes2D(FaceFNodes,0,2)
PairingFacesBaF = []
for i in range(len(FaceBaNodes)):
    Temp = []
    Temp.append(FaceBaNodes[i])
    x1 = RVENodalCoord[FaceBaNodes[i]][0]
    z1 = RVENodalCoord[FaceBaNodes[i]][2]
    
    Dist = []
    for j in range(len(FaceFNodes)):
        x2 = RVENodalCoord[FaceFNodes[j]][0]
        z2 = RVENodalCoord[FaceFNodes[j]][2]
        
        Dist.append(math.sqrt(pow(x2-x1,2)+pow(z2-z1,2)))
        
    N = Dist.index(min(Dist))
    Temp.append(FaceFNodes[N])
    FaceFNodes.pop(N)
    PairingFacesBaF.append(Temp)    

FaceLNodes = SortListofNodes2D(ExcludeNodes(FaceLNodes,[],[RVENodalCoord[V3][1],RVENodalCoord[V6][1]],[RVENodalCoord[V3][2],RVENodalCoord[V6][2]]),1,2)
FaceRNodes = SortListofNodes2D(ExcludeNodes(FaceRNodes,[],[RVENodalCoord[V3][1],RVENodalCoord[V6][1]],[RVENodalCoord[V3][2],RVENodalCoord[V6][2]]),1,2)
PairingFacesLR = []
for i in range(len(FaceLNodes)):
    Temp = []
    Temp.append(FaceLNodes[i])
    y1 = RVENodalCoord[FaceLNodes[i]][1]
    z1 = RVENodalCoord[FaceLNodes[i]][2]
    
    Dist = []
    for j in range(len(FaceRNodes)):
        y2 = RVENodalCoord[FaceRNodes[j]][1]
        z2 = RVENodalCoord[FaceRNodes[j]][2]
        
        Dist.append(math.sqrt(pow(y2-y1,2)+pow(z2-z1,2)))
        
    N = Dist.index(min(Dist))
    Temp.append(FaceRNodes[N])
    FaceRNodes.pop(N)
    PairingFacesLR.append(Temp)

FaceBNodes = SortListofNodes2D(ExcludeNodes(FaceBNodes,[RVENodalCoord[V3][0],RVENodalCoord[V8][0]],[RVENodalCoord[V3][1],RVENodalCoord[V8][1]],[]),0,1)
FaceTNodes = SortListofNodes2D(ExcludeNodes(FaceTNodes,[RVENodalCoord[V3][0],RVENodalCoord[V8][0]],[RVENodalCoord[V3][1],RVENodalCoord[V8][1]],[]),0,1)
PairingFacesBT = []
for i in range(len(FaceBNodes)):
    Temp = []
    Temp.append(FaceBNodes[i])
    x1 = RVENodalCoord[FaceBNodes[i]][0]
    y1 = RVENodalCoord[FaceBNodes[i]][1]
    
    Dist = []
    for j in range(len(FaceTNodes)):
        x2 = RVENodalCoord[FaceTNodes[j]][0]
        y2 = RVENodalCoord[FaceTNodes[j]][1]
        
        Dist.append(math.sqrt(pow(x2-x1,2)+pow(y2-y1,2)))
        
    N = Dist.index(min(Dist))
    Temp.append(FaceTNodes[N])
    FaceTNodes.pop(N)
    PairingFacesBT.append(Temp)

# Calculating the coefficients and setting up the MPCs
for i in range(Nele):
    C = np.array([
            [1,-1,-1,-1,1,1,1,-1],
            [1,1,-1,-1,-1,-1,1,1],
            [1,1,1,-1,1,-1,-1,-1],
            [1,-1,1,-1,-1,1,-1,1],
            [1,-1,-1,1,1,-1,-1,1],
            [1,1,-1,1,-1,1,-1,-1],
            [1,1,1,1,1,1,1,1],
            [1,-1,1,1,-1,-1,1,-1],])
    C_inv = np.linalg.inv(C)
    [a0,a1,a2,a3,a4,a5,a6,a7] = np.dot(C_inv,NodalCoordX[i])
    [b0,b1,b2,b3,b4,b5,b6,b7] = np.dot(C_inv,NodalCoordY[i])
    [d0,d1,d2,d3,d4,d5,d6,d7] = np.dot(C_inv,NodalCoordZ[i])
    
    # Calling macroscale nodes into sets
    for j in range(8):
        print>>Sets,'*Nset, nset=Ele'+str(i+1)+'-N'+str(j+1)+', instance='+str(MacroInstName)
        print>>Sets,str(NodalConnect[i][j]+1)
        
    # At each macroscale GP
    for j in range(8): # 8 GPs
        [tsi,eta,zeta] = GP[j]
        
        J = np.array([
                [a1+a4*eta+a5*zeta+a7*eta*zeta,b1+b4*eta+b5*zeta+b7*eta*zeta,d1+d4*eta+d5*zeta+d7*eta*zeta],
                [a2+a4*tsi+a6*zeta+a7*tsi*zeta,b2+b4*tsi+b6*zeta+b7*tsi*zeta,d2+d4*tsi+d6*zeta+d7*tsi*zeta],
                [a3+a5*tsi+a6*eta+a7*tsi*eta,b3+b5*tsi+b6*eta+b7*tsi*eta,d3+d5*tsi+d6*eta+d7*tsi*eta]])
        J_inv = np.linalg.inv(J)
        J_RVE = (np.linalg.det(abs(J))/(B_RVE*H_RVE*T_RVE))**(1.0/3.0)
        dN1 = [-0.125*(1-eta)*(1-zeta),-0.125*(1-tsi)*(1-zeta),-0.125*(1-tsi)*(1-eta)]
        dN2 = [0.125*(1-eta)*(1-zeta),-0.125*(1+tsi)*(1-zeta),-0.125*(1+tsi)*(1-eta)]
        dN3 = [0.125*(1+eta)*(1-zeta),0.125*(1+tsi)*(1-zeta),-0.125*(1+tsi)*(1+eta)]
        dN4 = [-0.125*(1+eta)*(1-zeta),0.125*(1-tsi)*(1-zeta),-0.125*(1-tsi)*(1+eta)]
        dN5 = [-0.125*(1-eta)*(1+zeta),-0.125*(1-tsi)*(1+zeta),0.125*(1-tsi)*(1-eta)]
        dN6 = [0.125*(1-eta)*(1+zeta),-0.125*(1+tsi)*(1+zeta),0.125*(1+tsi)*(1-eta)]
        dN7 = [0.125*(1+eta)*(1+zeta),0.125*(1+tsi)*(1+zeta),0.125*(1+tsi)*(1+eta)]
        dN8 = [-0.125*(1+eta)*(1+zeta),0.125*(1-tsi)*(1+zeta),0.125*(1-tsi)*(1+eta)]
        N_NatDeriv = [dN1,dN2,dN3,dN4,dN5,dN6,dN7,dN8]
        N_GloDeriv = [[],[],[],[],[],[],[],[]]
        
        Shape_fn = Trilin_Interpolation(tsi,eta,zeta)
        
        # Calculating shape function gradients
        for k in range(8):
            N_GloDeriv[k] = np.dot(J_inv,np.transpose(np.array(N_NatDeriv[k])))
        
        # Calculating the scaled RVE dimensions
        dx = B_RVE*J_RVE
        dy = H_RVE*J_RVE
        dz = T_RVE*J_RVE
        
        # Calling sets and setting up the MPCs for the left and right boundaries of FaceBa
        for k in range(len(EdgeLNodes)):
            print>>Sets,'*Nset, nset=Ele'+str(i+1)+'-RVE'+str(j+1)+'-EdgeNodeL'+str(k+1)+', instance=Ele'+str(i+1)+'-RVE'+str(j+1)
            print>>Sets,str(EdgeLNodes[k]+1)
            print>>Sets,'*Nset, nset=Ele'+str(i+1)+'-RVE'+str(j+1)+'-EdgeNodeR'+str(k+1)+', instance=Ele'+str(i+1)+'-RVE'+str(j+1)
            print>>Sets,str(EdgeRNodes[k]+1)
            
            for dof in range(3):
                print>>Eqns,'** Constraint: Ele'+str(i+1)+'-RVE'+str(j+1)+'-EdgeLR'+str(k+1)+'-DOF'+str(dof+1)
                print>>Eqns,'*Equation'
                print>>Eqns,'10'
                print>>Eqns,'Ele'+str(i+1)+'-RVE'+str(j+1)+'-EdgeNodeR'+str(k+1)+', '+str(dof+1)+', -1.0'
                print>>Eqns,'Ele'+str(i+1)+'-RVE'+str(j+1)+'-EdgeNodeL'+str(k+1)+', '+str(dof+1)+', 1.0'
                for m in range(8):
                    print>>Eqns,'Ele'+str(i+1)+'-N'+str(m+1)+', '+str(dof+1)+', '+str(dx*N_GloDeriv[m][0])
                    
        # Calling sets and setting up the MPCs for top and bottom boundaries of FaceBa
        for k in range(len(EdgeBNodes)):
            print>>Sets,'*Nset, nset=Ele'+str(i+1)+'-RVE'+str(j+1)+'-EdgeNodeB'+str(k+1)+', instance=Ele'+str(i+1)+'-RVE'+str(j+1)
            print>>Sets,str(EdgeBNodes[k]+1)
            print>>Sets,'*Nset, nset=Ele'+str(i+1)+'-RVE'+str(j+1)+'-EdgeNodeT'+str(k+1)+', instance=Ele'+str(i+1)+'-RVE'+str(j+1)
            print>>Sets,str(EdgeTNodes[k]+1)
            
            for dof in range(3):
                print>>Eqns,'** Constraint: Ele'+str(i+1)+'-RVE'+str(j+1)+'-EdgeLR'+str(k+1)+'-DOF'+str(dof+1)
                print>>Eqns,'*Equation'
                print>>Eqns,'10'
                print>>Eqns,'Ele'+str(i+1)+'-RVE'+str(j+1)+'-EdgeNodeT'+str(k+1)+', '+str(dof+1)+', -1.0'
                print>>Eqns,'Ele'+str(i+1)+'-RVE'+str(j+1)+'-EdgeNodeB'+str(k+1)+', '+str(dof+1)+', 1.0'
                for m in range(8):
                    print>>Eqns,'Ele'+str(i+1)+'-N'+str(m+1)+', '+str(dof+1)+', '+str(dz*N_GloDeriv[m][2])
                    
        # Calling sets and setting up the MPCs for corners of FaceBa
        print>>Sets,'*Nset, nset=Ele'+str(i+1)+'-RVE'+str(j+1)+'-V1, instance=Ele'+str(i+1)+'-RVE'+str(j+1)
        print>>Sets,str(V1+1)
        print>>Sets,'*Nset, nset=Ele'+str(i+1)+'-RVE'+str(j+1)+'-V2, instance=Ele'+str(i+1)+'-RVE'+str(j+1)
        print>>Sets,str(V2+1)
        print>>Sets,'*Nset, nset=Ele'+str(i+1)+'-RVE'+str(j+1)+'-V3, instance=Ele'+str(i+1)+'-RVE'+str(j+1)
        print>>Sets,str(V3+1)
        print>>Sets,'*Nset, nset=Ele'+str(i+1)+'-RVE'+str(j+1)+'-V4, instance=Ele'+str(i+1)+'-RVE'+str(j+1)
        print>>Sets,str(V4+1)
        
        for dof in range(3):
            print>>Eqns,'** Constraint: Ele'+str(i+1)+'-RVE'+str(j+1)+'-V12-DOF'+str(dof+1)
            print>>Eqns,'*Equation'
            print>>Eqns,'10'
            print>>Eqns,'Ele'+str(i+1)+'-RVE'+str(j+1)+'-V2, '+str(dof+1)+', -1.0'
            print>>Eqns,'Ele'+str(i+1)+'-RVE'+str(j+1)+'-V1, '+str(dof+1)+', 1.0'
            for m in range(8):
                print>>Eqns,'Ele'+str(i+1)+'-N'+str(m+1)+', '+str(dof+1)+', '+str(dx*N_GloDeriv[m][0])
                
        for dof in range(3):
            print>>Eqns,'** Constraint: Ele'+str(i+1)+'-RVE'+str(j+1)+'-V23-DOF'+str(dof+1)
            print>>Eqns,'*Equation'
            print>>Eqns,'10'
            print>>Eqns,'Ele'+str(i+1)+'-RVE'+str(j+1)+'-V3, '+str(dof+1)+', -1.0'
            print>>Eqns,'Ele'+str(i+1)+'-RVE'+str(j+1)+'-V2, '+str(dof+1)+', 1.0'
            for m in range(8):
                print>>Eqns,'Ele'+str(i+1)+'-N'+str(m+1)+', '+str(dof+1)+', '+str(dz*N_GloDeriv[m][2])
                
        for dof in range(3):
            print>>Eqns,'** Constraint: Ele'+str(i+1)+'-RVE'+str(j+1)+'-V14-DOF'+str(dof+1)
            print>>Eqns,'*Equation'
            print>>Eqns,'10'
            print>>Eqns,'Ele'+str(i+1)+'-RVE'+str(j+1)+'-V4, '+str(dof+1)+', -1.0'
            print>>Eqns,'Ele'+str(i+1)+'-RVE'+str(j+1)+'-V1, '+str(dof+1)+', 1.0'
            for m in range(8):
                print>>Eqns,'Ele'+str(i+1)+'-N'+str(m+1)+', '+str(dof+1)+', '+str(dz*N_GloDeriv[m][2])
                
        for dof in range(3):
            print>>Eqns,'** Constraint: Ele'+str(i+1)+'-RVE'+str(j+1)+'-V1-DOF'+str(dof+1)
            print>>Eqns,'*Equation'
            print>>Eqns,'9'
            print>>Eqns,'Ele'+str(i+1)+'-RVE'+str(j+1)+'-V1, '+str(dof+1)+', -1.0'
            for m in range(8):
                print>>Eqns,'Ele'+str(i+1)+'-N'+str(m+1)+', '+str(dof+1)+', '+str(Shape_fn[m]-0.5*dx*N_GloDeriv[m][0]-0.5*dy*N_GloDeriv[m][1]-0.5*dz*N_GloDeriv[m][2])
            
        # Calling sets and setting up the MPCs for the front and back faces
        for k in range(len(PairingFacesBaF)):
            print>>Sets,'*Nset, nset=Ele'+str(i+1)+'-RVE'+str(j+1)+'-FaceNodeBa'+str(k+1)+', instance=Ele'+str(i+1)+'-RVE'+str(j+1)
            print>>Sets,str(PairingFacesBaF[k][0]+1)
            print>>Sets,'*Nset, nset=Ele'+str(i+1)+'-RVE'+str(j+1)+'-FaceNodeF'+str(k+1)+', instance=Ele'+str(i+1)+'-RVE'+str(j+1)
            print>>Sets,str(PairingFacesBaF[k][1]+1)
            
            for dof in range(3):
                print>>Eqns,'** Constraint: Ele'+str(i+1)+'-RVE'+str(j+1)+'-BaF'+str(k+1)+'-DOF'+str(dof+1)
                print>>Eqns,'*Equation'
                print>>Eqns,'10'
                print>>Eqns,'Ele'+str(i+1)+'-RVE'+str(j+1)+'-FaceNodeF'+str(k+1)+', '+str(dof+1)+', -1.0'
                print>>Eqns,'Ele'+str(i+1)+'-RVE'+str(j+1)+'-FaceNodeBa'+str(k+1)+', '+str(dof+1)+', 1.0'
                for m in range(8):
                    print>>Eqns,'Ele'+str(i+1)+'-N'+str(m+1)+', '+str(dof+1)+', '+str(dy*N_GloDeriv[m][1])
                    
        # Calling sets and setting up the MPCs for the edges parallel to y axis
        for k in range(len(Edge1Nodes)):
            print>>Sets,'*Nset, nset=Ele'+str(i+1)+'-RVE'+str(j+1)+'-Edge1Node'+str(k+1)+', instance=Ele'+str(i+1)+'-RVE'+str(j+1)
            print>>Sets,str(Edge1Nodes[k]+1)
            print>>Sets,'*Nset, nset=Ele'+str(i+1)+'-RVE'+str(j+1)+'-Edge2Node'+str(k+1)+', instance=Ele'+str(i+1)+'-RVE'+str(j+1)
            print>>Sets,str(Edge2Nodes[k]+1)
            print>>Sets,'*Nset, nset=Ele'+str(i+1)+'-RVE'+str(j+1)+'-Edge3Node'+str(k+1)+', instance=Ele'+str(i+1)+'-RVE'+str(j+1)
            print>>Sets,str(Edge3Nodes[k]+1)
            print>>Sets,'*Nset, nset=Ele'+str(i+1)+'-RVE'+str(j+1)+'-Edge4Node'+str(k+1)+', instance=Ele'+str(i+1)+'-RVE'+str(j+1)
            print>>Sets,str(Edge4Nodes[k]+1)
            
            for dof in range(3):
                print>>Eqns,'** Constraint: Ele'+str(i+1)+'-RVE'+str(j+1)+'-Edge12Node'+str(k+1)+'-DOF'+str(dof+1)
                print>>Eqns,'*Equation'
                print>>Eqns,'10'
                print>>Eqns,'Ele'+str(i+1)+'-RVE'+str(j+1)+'-Edge2Node'+str(k+1)+', '+str(dof+1)+', -1.0'
                print>>Eqns,'Ele'+str(i+1)+'-RVE'+str(j+1)+'-Edge1Node'+str(k+1)+', '+str(dof+1)+', 1.0'
                for m in range(8):
                    print>>Eqns,'Ele'+str(i+1)+'-N'+str(m+1)+', '+str(dof+1)+', '+str(dx*N_GloDeriv[m][0])
                    
            for dof in range(3):
                print>>Eqns,'** Constraint: Ele'+str(i+1)+'-RVE'+str(j+1)+'-Edge23Node'+str(k+1)+'-DOF'+str(dof+1)
                print>>Eqns,'*Equation'
                print>>Eqns,'10'
                print>>Eqns,'Ele'+str(i+1)+'-RVE'+str(j+1)+'-Edge3Node'+str(k+1)+', '+str(dof+1)+', -1.0'
                print>>Eqns,'Ele'+str(i+1)+'-RVE'+str(j+1)+'-Edge2Node'+str(k+1)+', '+str(dof+1)+', 1.0'
                for m in range(8):
                    print>>Eqns,'Ele'+str(i+1)+'-N'+str(m+1)+', '+str(dof+1)+', '+str(dz*N_GloDeriv[m][2])
                    
            for dof in range(3):
                print>>Eqns,'** Constraint: Ele'+str(i+1)+'-RVE'+str(j+1)+'-Edge14Node'+str(k+1)+'-DOF'+str(dof+1)
                print>>Eqns,'*Equation'
                print>>Eqns,'10'
                print>>Eqns,'Ele'+str(i+1)+'-RVE'+str(j+1)+'-Edge4Node'+str(k+1)+', '+str(dof+1)+', -1.0'
                print>>Eqns,'Ele'+str(i+1)+'-RVE'+str(j+1)+'-Edge1Node'+str(k+1)+', '+str(dof+1)+', 1.0'
                for m in range(8):
                    print>>Eqns,'Ele'+str(i+1)+'-N'+str(m+1)+', '+str(dof+1)+', '+str(dz*N_GloDeriv[m][2])
                    
        # Calling sets and setting up the MPCs for the left and right faces
        for k in range(len(PairingFacesLR)):
            print>>Sets,'*Nset, nset=Ele'+str(i+1)+'-RVE'+str(j+1)+'-FaceNodeL'+str(k+1)+', instance=Ele'+str(i+1)+'-RVE'+str(j+1)
            print>>Sets,str(PairingFacesLR[k][0]+1)
            print>>Sets,'*Nset, nset=Ele'+str(i+1)+'-RVE'+str(j+1)+'-FaceNodeR'+str(k+1)+', instance=Ele'+str(i+1)+'-RVE'+str(j+1)
            print>>Sets,str(PairingFacesLR[k][1]+1)
            
            for dof in range(3):
                print>>Eqns,'** Constraint: Ele'+str(i+1)+'-RVE'+str(j+1)+'-LR'+str(k+1)+'-DOF'+str(dof+1)
                print>>Eqns,'*Equation'
                print>>Eqns,'10'
                print>>Eqns,'Ele'+str(i+1)+'-RVE'+str(j+1)+'-FaceNodeR'+str(k+1)+', '+str(dof+1)+', -1.0'
                print>>Eqns,'Ele'+str(i+1)+'-RVE'+str(j+1)+'-FaceNodeL'+str(k+1)+', '+str(dof+1)+', 1.0'
                for m in range(8):
                    print>>Eqns,'Ele'+str(i+1)+'-N'+str(m+1)+', '+str(dof+1)+', '+str(dx*N_GloDeriv[m][0])
                    
        # Calling sets and setting up the MPCs for the top and bottom faces
        for k in range(len(PairingFacesBT)):
            print>>Sets,'*Nset, nset=Ele'+str(i+1)+'-RVE'+str(j+1)+'-FaceNodeB'+str(k+1)+', instance=Ele'+str(i+1)+'-RVE'+str(j+1)
            print>>Sets,str(PairingFacesBT[k][0]+1)
            print>>Sets,'*Nset, nset=Ele'+str(i+1)+'-RVE'+str(j+1)+'-FaceNodeT'+str(k+1)+', instance=Ele'+str(i+1)+'-RVE'+str(j+1)
            print>>Sets,str(PairingFacesBT[k][1]+1)
            
            for dof in range(3):
                print>>Eqns,'** Constraint: Ele'+str(i+1)+'-RVE'+str(j+1)+'-BT'+str(k+1)+'-DOF'+str(dof+1)
                print>>Eqns,'*Equation'
                print>>Eqns,'10'
                print>>Eqns,'Ele'+str(i+1)+'-RVE'+str(j+1)+'-FaceNodeT'+str(k+1)+', '+str(dof+1)+', -1.0'
                print>>Eqns,'Ele'+str(i+1)+'-RVE'+str(j+1)+'-FaceNodeB'+str(k+1)+', '+str(dof+1)+', 1.0'
                for m in range(8):
                    print>>Eqns,'Ele'+str(i+1)+'-N'+str(m+1)+', '+str(dof+1)+', '+str(dz*N_GloDeriv[m][2])
            
Sets.close()
Eqns.close()


### Writing the new Direct FE2 input file
RVEParts = open('RVEParts.dat','r')
Insts = open('Insts.dat','r')
Sets = open('Sets.dat','r')
Eqns = open('Eqns.dat','r')
RVEMats = open('RVEMats.dat','r')

RVEParts_lines = []
Insts_lines = []
Sets_lines = []
Eqns_lines = []
RVEMats_lines = []

while 1:
    line = RVEParts.readline()
    if not line:
        break
    line = line.strip()
    RVEParts_lines.append(line)

while 1:
    line = Insts.readline()
    if not line:
        break
    line = line.strip()
    Insts_lines.append(line)
    
while 1:
    line = Sets.readline()
    if not line:
        break
    line = line.strip()
    Sets_lines.append(line)
    
while 1:
    line = Eqns.readline()
    if not line:
        break
    line = line.strip()
    Eqns_lines.append(line)
    
while 1:
    line = RVEMats.readline()
    if not line:
        break
    line = line.strip()
    RVEMats_lines.append(line)

RVEParts.close()
Insts.close()
Sets.close()
Eqns.close()
RVEMats.close()           

f3 = open(NewInpName,'w')

Mark1 = inp1.index('*End Part')
for i in range(0,Mark1+1):
    print>>f3,inp1[i]

for i in range(len(RVEParts_lines)):
    print>>f3,RVEParts_lines[i] 
    
Mark2 = inp1.index('*End Instance')
for i in range(Mark1+1,Mark2+2):
    print>>f3,inp1[i]

for i in range(len(Insts_lines)):
    print>>f3,Insts_lines[i]

if StartConst == 'a':
    StartConst = inp1.index('*End Assembly')    
for i in range(Mark2+2,StartConst):
    print>>f3,inp1[i]

for i in range(len(Sets_lines)):
    print>>f3,Sets_lines[i]
    
for i in range(len(Eqns_lines)):
    print>>f3,Eqns_lines[i]                    
        
Mark3 = inp1.index('** ----------------------------------------------------------------')
for i in range(StartConst,Mark3):
    print>>f3,inp1[i]

for i in range(len(RVEMats_lines)):
    print>>f3,RVEMats_lines[i]
    
for i in range(Mark3,len(inp1)):
    print>>f3,inp1[i]
    
f3.close()

del RVEParts_lines
del Insts_lines
del Sets_lines
del Eqns_lines
del RVEMats_lines 
os.remove('RVEParts.dat')
os.remove('Insts.dat')
os.remove('Sets.dat')
os.remove('Eqns.dat')
os.remove('RVEMats.dat')


'''
Revision log

230714 Original  release

240916 Revision
Replaced 'remove' function with 'del' function in lines 147 and 182
'remove' function searches and deletes the first match, while 'del' function deletes the specific line as intended 

End of Revision
'''

















