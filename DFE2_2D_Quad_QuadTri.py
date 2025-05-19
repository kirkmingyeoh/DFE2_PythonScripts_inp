# -*- coding: utf-8 -*-
"""
Created on Fri Jul 14 15:07:06 2023

@author: Kirk Ming Yeoh (e0546208@u.nus.edu)
"""

import os
import numpy as np
import math

### Update the following parameters
os.chdir('E:\\Kirk Ming Abaqus\\DFE2_Demo') # directory where input files are, use double '\'
MacroInpName = 'Demo_2D_Macro.inp' # Name of macroscale input file
RVEInpName = 'Demo_2D_Micro.inp' # Name of RVE input file Demo_2D_RVE
NewInpName = 'DFE2_2D.inp' # Name of new Direct FE2 input file

'''
No further user input is required
'''

GP = [[-3**-0.5,-3**-0.5],[3**-0.5,-3**-0.5],[3**-0.5,3**-0.5],[-3**-0.5,3**-0.5]]

### Define functions
# Search inp list for a particular part's elements and nodes
def Search(inp,key1,out1,type1): # If dealing with nodes and coordinates, set type1 to 1 
    for n_inp_lines in range(len(inp)): # Read the input file line by line
        if inp[n_inp_lines].count(key1)!=0: # Stop at keyword key1
            break    
    for temp_line in inp[n_inp_lines+1:]: # Read from the line after the one containing keyword key1
        if (temp_line == '') or (temp_line.count("*") != 0): # Stop when the list of nodes or elements ends
            break
        temp_line = (temp_line.replace(',',' ')).split() # Split the line of nodal coordinates or element nodal connectivities into a list
        temp_line.pop(0) # Removes node number in nodal coordinates or element number in connectivity
        for n_list_terms in range(len(temp_line)):
            if type1 == 1:
                temp_line[n_list_terms] = float(temp_line[n_list_terms])
            else:
                temp_line[n_list_terms] = int(temp_line[n_list_terms])-1 # All node labels in connectivity list are -1 for use in Python as indices
        out1.append(temp_line) # Returns list of nodal cordinates or element nodal connectivity

# Removes corner nodes from an edge based on the smallest and largest coordinate
def TakeVertexOut(edge):
    del edge[0]
    del edge[-1]
    return edge # Returns a list of nodes on the edges without the corner node of the RVE

# Calculates bilinear shape function values
def Bilin_Interpolation(tsi,eta): # tsi and eta are the natural coordinates
    N1 = float(0.25*(1-tsi)*(1-eta))
    N2 = float(0.25*(1+tsi)*(1-eta))
    N3 = float(0.25*(1+tsi)*(1+eta))
    N4 = float(0.25*(1-tsi)*(1+eta))
    return [N1,N2,N3,N4] # Returns a list of shape function values

# Sorts nodes along an edge using their coordinates
def SortListofNodes1D(faceN,coordinate): # 0 for x; 1 for y; 2 for z; gives node numbers ready to be called with python (already -1)
    newlist = []
    oldlist = []
    for n_face_nodes in range(len(faceN)): # Obtain a list of coordinates for all the nodes
        oldlist.append(RVENodalCoord[faceN[n_face_nodes]][coordinate])
    
    orderedlist = sorted(oldlist) # Sort the nodal coordinates
    for n_oldlist_nodes in range(len(orderedlist)): # Sort the nodes based on the sorted list of coodinates
        ind = oldlist.index(orderedlist[n_oldlist_nodes])
        newlist.append(faceN[ind])
    
    return newlist # Returns a list of sorted nodes


### Extracting information from Macro and RVE input files
# Macroscale input file
inp1 = []    
f1 = open(MacroInpName,'r') # Open the macroscale input file as f1

while 1: # Read the macroscale input file line by line and store it
    line = f1.readline()
    if not line:
        break
    line = line.strip() # Removes additional white spaces on left and right
    inp1.append(line)
    
f1.close() 

# Removing 'generate' for easier processing
for n_inp1_lines in reversed(range(len(inp1))): # Read through the macroscale input file lines, done in reverse to avoid issues with line number when expanding the 'generate' keyword
    if (inp1[n_inp1_lines].count('generate')!=0): # Lines that compact node or element lists and contain the 'generate' keyword
        Temp = (inp1[n_inp1_lines+1].replace(',',' ')).split() # Split the key numbers in the compacted list containing the start, end and increment for the list
        n_terms = 0 # Term counter
        n_lines = 0 # Extra line counter
        for n_fulllist_terms in range(int(Temp[0]),int(Temp[1])+1,int(Temp[2])):
            if n_terms==0:
                Temp2 = str(n_fulllist_terms)
                n_terms = n_terms+1
            elif n_terms==16:
                inp1.insert(n_inp1_lines+2+n_lines,Temp2)
                n_lines = n_lines+1
                Temp2 = str(n_fulllist_terms)
                n_terms = 1
            else:
                Temp2 = Temp2+', '+str(n_fulllist_terms)
                n_terms = n_terms+1
        inp1.insert(n_inp1_lines+2+m,Temp2)
        inp1[n_inp1_lines] = inp1[n_inp1_lines][0:len(inp1[n_inp1_lines])-10]
        del inp1[n_inp1_lines+1]

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
    if (inp1[i].count('** Section'))!=0:
        if inp1[i+2] == ',':
            Thickness = 1.0
        else:
            Thickness = float(inp1[i+2].strip(','))
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
for i in range(Nele):
    Nodes = []
    X = []
    Y = []
    Ang = []
    TempConnect = []
    TempCoordX = []
    TempCoordY = []
    
    # Obtaining nodal information from this particular element
    for j in range(4):
        Nodes.append(MacroNodalConnect[i][j])
        X.append(MacroNodalCoord[MacroNodalConnect[i][j]][0])
        Y.append(MacroNodalCoord[MacroNodalConnect[i][j]][1])
    
    # Obtaining centroid of the element
    X0 = sum(X)/4
    Y0 = sum(Y)/4
    
    # Obtaining orientation of each node relative to centroid
    for j in range(4):
        X1 = X[j]-X0
        Y1 = Y[j]-Y0
        theta = math.atan(Y1/X1)
        if X1<0:
            theta = theta+math.pi
        if theta<0:
            theta = theta+2*(math.pi)
        Ang.append(theta*360/(2*(math.pi)))
    SAng = sorted(Ang)
    
    for j in range(4):
        Order1 = [2,3,0,1]
        N = Ang.index(SAng[Order1[j]])
        TempConnect.append(Nodes[N])
        TempCoordX.append(X[N])
        TempCoordY.append(Y[N])
        
    NodalConnect.append(TempConnect)
    NodalCoordX.append(TempCoordX)
    NodalCoordY.append(TempCoordY)


### Processing the RVE part information
# Finding the smallest and largest nodal coordinate in both directions
RVE_ListX = []
RVE_ListY = []
for i in range(len(RVENodalCoord)):
    RVE_ListX.append(RVENodalCoord[i][0])
    RVE_ListY.append(RVENodalCoord[i][1])
xMin = min(RVE_ListX)
xMax = max(RVE_ListX)
yMin = min(RVE_ListY)
yMax = max(RVE_ListY)
del RVE_ListX
del RVE_ListY

# Sorting the RVE boundary nodes
FaceLNodes,FaceRNodes,FaceBNodes,FaceTNodes = [],[],[],[]
for i in range(len(RVENodalCoord)):
    if RVENodalCoord[i][0] == xMin:
        FaceLNodes.append(i)
    if RVENodalCoord[i][0] == xMax:
        FaceRNodes.append(i)
    if RVENodalCoord[i][1] == yMin:
        FaceBNodes.append(i)
    if RVENodalCoord[i][1] == yMax:
        FaceTNodes.append(i)
for i in range(len(FaceLNodes)):
    if FaceLNodes[i] in FaceBNodes:
        V1 = FaceLNodes[i]
    if FaceLNodes[i] in FaceTNodes:
        V4 = FaceLNodes[i]
for i in range(len(FaceRNodes)):
    if FaceRNodes[i] in FaceBNodes:
        V2 = FaceRNodes[i]
    if FaceRNodes[i] in FaceTNodes:
        V3 = FaceRNodes[i]
    
# Sorting RVE dimensions and offsets
B_RVE = xMax - xMin
H_RVE = yMax - yMin
OffsetX = (xMax + xMin)/2
OffsetY = (yMax + yMin)/2
Offset = [OffsetX,OffsetY]

# Adjusting RVE nodal coordinates to correspond to a part centered at the origin
for i in RVENodalCoord:
    for j in range(2):
        i[j] = i[j] - Offset[j]
        

### Generating the RVE placement in the macroscale mesh
RVEParts = open('RVEParts.dat','w')
Insts = open('Insts.dat','w')
SF = []
for i in range(Nele):
    C = np.array([[1,-1,-1,1],[1,1,-1,-1],[1,1,1,1],[1,-1,1,-1]])
    C_inv = np.linalg.inv(C)
    [a0,a1,a2,a3] = np.dot(C_inv,NodalCoordX[i])
    [b0,b1,b2,b3] = np.dot(C_inv,NodalCoordY[i])
    
    for j in range(4):
        [tsi,eta] = GP[j]
        J = np.array([[a1+a3*eta,b1+b3*eta],[a2+a3*tsi,b2+b3*tsi]])
        J_RVE = abs(np.linalg.det(J)/(B_RVE*H_RVE)) # Scaling factor for thicknes
        
        [N1,N2,N3,N4] = Bilin_Interpolation(tsi,eta)
        RVE_X = N1*NodalCoordX[i][0] + N2*NodalCoordX[i][1] + N3*NodalCoordX[i][2] + N4*NodalCoordX[i][3]
        RVE_Y = N1*NodalCoordY[i][0] + N2*NodalCoordY[i][1] + N3*NodalCoordY[i][2] + N4*NodalCoordY[i][3]
        
        if round(J_RVE,5) in SF: # If such an RVE was created previously, reuse it for the next instance
            N = SF.index(round(J_RVE,5))
        else: # If not, then we create another RVE with this current thickness scaling
            N = len(SF)
            SF.append(round(J_RVE,5))
            print>>RVEParts,'**'
            print>>RVEParts,'*Part, name=RVE-'+str(N+1)
            print>>RVEParts,'*Node'
            
            for k in range(len(RVENodalCoord)):
                print>>RVEParts,str(k+1)+', '+str(RVENodalCoord[k][0])+', '+str(RVENodalCoord[k][1])
                
            for k in range(StartEle,Sections[0]):
                print>>RVEParts,inp2[k]
                
            for k in range(len(Sections)):
                print>>RVEParts,inp2[Sections[k]]+'-'+str(N+1)
                print>>RVEParts,inp2[Sections[k]+1]
                if inp2[Sections[k]+1].count('*Solid Section')!=0:
                    print>>RVEParts,str(Thickness*J_RVE)+','
                elif inp2[Sections[k]+1].count('*Cohesive Section')!=0:
                    Line = inp2[Sections[k]+2].split(',')
                    print>>RVEParts,str(Line[0])+','+str(Thickness*J_RVE)
                
            print>>RVEParts,'*End Part'
            
        print>>Insts,'*Instance, name=Ele'+str(i+1)+'-RVE'+str(j+1)+', part=RVE-'+str(N+1)
        print>>Insts,str(RVE_X)+', '+str(RVE_Y)+', 0.'
        print>>Insts,'*End Instance'
        print>>Insts,'**'
            
RVEParts.close()
Insts.close()

            
### Setting up the MPCs
Sets = open('Sets.dat','w')
Eqns = open('Eqns.dat','w')

# Pairing the nodes
FaceLNodes = TakeVertexOut(SortListofNodes1D(FaceLNodes,1))
FaceRNodes = TakeVertexOut(SortListofNodes1D(FaceRNodes,1))
PairingFacesLR = []
for i in range(len(FaceLNodes)):
    Temp = []
    Temp.append(FaceLNodes[i])
    Temp.append(FaceRNodes[i])
    PairingFacesLR.append(Temp)

FaceBNodes = TakeVertexOut(SortListofNodes1D(FaceBNodes,0))
FaceTNodes = TakeVertexOut(SortListofNodes1D(FaceTNodes,0))
PairingFacesBT = []
for i in range(len(FaceBNodes)):
    Temp = []
    Temp.append(FaceBNodes[i])
    Temp.append(FaceTNodes[i])
    PairingFacesBT.append(Temp)
    
# Calculating the coefficients and setting up the MPCs
for i in range(Nele):
    C = np.array([[1,-1,-1,1],[1,1,-1,-1],[1,1,1,1],[1,-1,1,-1]])
    C_inv = np.linalg.inv(C)
    [a0,a1,a2,a3] = np.dot(C_inv,NodalCoordX[i])
    [b0,b1,b2,b3] = np.dot(C_inv,NodalCoordY[i])
    
    # Calling macroscale nodes into sets
    for j in range(4):
        print>>Sets,'*Nset, nset=Ele'+str(i+1)+'-N'+str(j+1)+', instance='+str(MacroInstName)
        print>>Sets,str(NodalConnect[i][j]+1)
    
    # At each macroscale GP
    for j in range(4): # 4 GPs
        [tsi,eta] = GP[j]
        
        J = np.array([[a1+a3*eta,b1+b3*eta],[a2+a3*tsi,b2+b3*tsi]])
        J_inv = np.linalg.inv(J)
        dN1 = [-0.25*(1-eta),-0.25*(1-tsi)]
        dN2 = [0.25*(1-eta),-0.25*(1+tsi)]
        dN3 = [0.25*(1+eta),0.25*(1+tsi)]
        dN4 = [-0.25*(1+eta),0.25*(1-tsi)]
        N_NatDeriv = [dN1,dN2,dN3,dN4]
        N_GloDeriv = [[],[],[],[]]
        
        Shape_fn = Bilin_Interpolation(tsi,eta)
        
        # Calculating shape function gradients
        for k in range(4):
            N_GloDeriv[k] = np.dot(J_inv,np.transpose(np.array(N_NatDeriv[k])))
        
        # Calling sets and setting up the MPCs for left and right faces
        for k in range(len(PairingFacesLR)):
            print>>Sets,'*Nset, nset=Ele'+str(i+1)+'-RVE'+str(j+1)+'-FaceNodeL'+str(k+1)+', instance=Ele'+str(i+1)+'-RVE'+str(j+1)
            print>>Sets,str(PairingFacesLR[k][0]+1)
            print>>Sets,'*Nset, nset=Ele'+str(i+1)+'-RVE'+str(j+1)+'-FaceNodeR'+str(k+1)+', instance=Ele'+str(i+1)+'-RVE'+str(j+1)
            print>>Sets,str(PairingFacesLR[k][1]+1)
            
            for dof in range(2):
                print>>Eqns,'** Constraint: Ele'+str(i+1)+'-RVE'+str(j+1)+'-LR'+str(k+1)+'-DOF'+str(dof+1)
                print>>Eqns,'*Equation'
                print>>Eqns,'6'
                print>>Eqns,'Ele'+str(i+1)+'-RVE'+str(j+1)+'-FaceNodeR'+str(k+1)+', '+str(dof+1)+', -1.0'
                print>>Eqns,'Ele'+str(i+1)+'-RVE'+str(j+1)+'-FaceNodeL'+str(k+1)+', '+str(dof+1)+', 1.0'
                for m in range(4):
                    print>>Eqns,'Ele'+str(i+1)+'-N'+str(m+1)+', '+str(dof+1)+', '+str(B_RVE*N_GloDeriv[m][0])
                    
        # Calling sets and setting up the MPCs for top and bottom faces
        for k in range(len(PairingFacesBT)):
            print>>Sets,'*Nset, nset=Ele'+str(i+1)+'-RVE'+str(j+1)+'-FaceNodeB'+str(k+1)+', instance=Ele'+str(i+1)+'-RVE'+str(j+1)
            print>>Sets,str(PairingFacesBT[k][0]+1)
            print>>Sets,'*Nset, nset=Ele'+str(i+1)+'-RVE'+str(j+1)+'-FaceNodeT'+str(k+1)+', instance=Ele'+str(i+1)+'-RVE'+str(j+1)
            print>>Sets,str(PairingFacesBT[k][1]+1)
            
            for dof in range(2):
                print>>Eqns,'** Constraint: Ele'+str(i+1)+'-RVE'+str(j+1)+'-BT'+str(k+1)+'-DOF'+str(dof+1)
                print>>Eqns,'*Equation'
                print>>Eqns,'6'
                print>>Eqns,'Ele'+str(i+1)+'-RVE'+str(j+1)+'-FaceNodeT'+str(k+1)+', '+str(dof+1)+', -1.0'
                print>>Eqns,'Ele'+str(i+1)+'-RVE'+str(j+1)+'-FaceNodeB'+str(k+1)+', '+str(dof+1)+', 1.0'
                for m in range(4):
                    print>>Eqns,'Ele'+str(i+1)+'-N'+str(m+1)+', '+str(dof+1)+', '+str(H_RVE*N_GloDeriv[m][1])
                    
        # Calling sets and setting up the MPCs for vertices
        print>>Sets,'*Nset, nset=Ele'+str(i+1)+'-RVE'+str(j+1)+'-V1, instance=Ele'+str(i+1)+'-RVE'+str(j+1)
        print>>Sets,str(V1+1)
        print>>Sets,'*Nset, nset=Ele'+str(i+1)+'-RVE'+str(j+1)+'-V2, instance=Ele'+str(i+1)+'-RVE'+str(j+1)
        print>>Sets,str(V2+1)
        print>>Sets,'*Nset, nset=Ele'+str(i+1)+'-RVE'+str(j+1)+'-V3, instance=Ele'+str(i+1)+'-RVE'+str(j+1)
        print>>Sets,str(V3+1)
        print>>Sets,'*Nset, nset=Ele'+str(i+1)+'-RVE'+str(j+1)+'-V4, instance=Ele'+str(i+1)+'-RVE'+str(j+1)
        print>>Sets,str(V4+1)
        
        for dof in range(2):
            print>>Eqns,'** Constraint: Ele'+str(i+1)+'-RVE'+str(j+1)+'-V23-DOF'+str(dof+1)
            print>>Eqns,'*Equation'
            print>>Eqns,'6'
            print>>Eqns,'Ele'+str(i+1)+'-RVE'+str(j+1)+'-V3, '+str(dof+1)+', -1.0'
            print>>Eqns,'Ele'+str(i+1)+'-RVE'+str(j+1)+'-V2, '+str(dof+1)+', 1.0'
            for m in range(4):
                print>>Eqns,'Ele'+str(i+1)+'-N'+str(m+1)+', '+str(dof+1)+', '+str(H_RVE*N_GloDeriv[m][1])
                
        for dof in range(2):
            print>>Eqns,'** Constraint: Ele'+str(i+1)+'-RVE'+str(j+1)+'-V14-DOF'+str(dof+1)
            print>>Eqns,'*Equation'
            print>>Eqns,'6'
            print>>Eqns,'Ele'+str(i+1)+'-RVE'+str(j+1)+'-V4, '+str(dof+1)+', -1.0'
            print>>Eqns,'Ele'+str(i+1)+'-RVE'+str(j+1)+'-V1, '+str(dof+1)+', 1.0'
            for m in range(4):
                print>>Eqns,'Ele'+str(i+1)+'-N'+str(m+1)+', '+str(dof+1)+', '+str(H_RVE*N_GloDeriv[m][1])
                
        for dof in range(2):
            print>>Eqns,'** Constraint: Ele'+str(i+1)+'-RVE'+str(j+1)+'-V12-DOF'+str(dof+1)
            print>>Eqns,'*Equation'
            print>>Eqns,'6'
            print>>Eqns,'Ele'+str(i+1)+'-RVE'+str(j+1)+'-V2, '+str(dof+1)+', -1.0'
            print>>Eqns,'Ele'+str(i+1)+'-RVE'+str(j+1)+'-V1, '+str(dof+1)+', 1.0'
            for m in range(4):
                print>>Eqns,'Ele'+str(i+1)+'-N'+str(m+1)+', '+str(dof+1)+', '+str(B_RVE*N_GloDeriv[m][0])
                
        for dof in range(2):
            print>>Eqns,'** Constraint: Ele'+str(i+1)+'-RVE'+str(j+1)+'-V1-DOF'+str(dof+1)
            print>>Eqns,'*Equation'
            print>>Eqns,'5'
            print>>Eqns,'Ele'+str(i+1)+'-RVE'+str(j+1)+'-V1, '+str(dof+1)+', -1.0'
            for m in range(4):
                print>>Eqns,'Ele'+str(i+1)+'-N'+str(m+1)+', '+str(dof+1)+', '+str(Shape_fn[m]-0.5*B_RVE*N_GloDeriv[m][0]-0.5*H_RVE*N_GloDeriv[m][1])

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
Replaced 'remove' function with 'del' function in lines 102 and 137
'remove' function searches and deletes the first match, while 'del' function deletes the specific line as intended 

Added '+1' to line 180 of original code (currently line 185)
When reading RVE material definitions, current code will miss the last line of the entire RVE input file

End of 240916 Revision

250519 Revision (WIP)
Replaced variables with explanatory names
Added additional comments for clarity

End of 250519 Revision

Proposed Revisions (yet to be implemented)
Generalise RVE Sections portion to account for possible material orientations
Account for RVE sets and surfaces from the microscale input file, found from Part?
Account for RVE level interactions and constraints, found from Instance, including the rearrangement? adopt from code for Yuhao and ZB
Account for BCs applied at Initial step when writing RVE materials, which appears before the ------- line? adopt from code for Yuhao and ZB
Allow for (partial) reduced integration by changing the list of GP? adopt from code for ZB
Account for multiple macroscale parts, multiple types of RVE and different RVE orientations? adopt from code for Yuhao
Relabel J as J^T for better clarity
Account for X1=0 cases when sorting for the macroscale element connectivity
'''     
        
            
        
            
            


            





























