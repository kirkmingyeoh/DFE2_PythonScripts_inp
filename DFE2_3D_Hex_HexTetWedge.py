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

'''
No further user input is required
'''

GP = [[-3**-0.5,-3**-0.5,-3**-0.5],[3**-0.5,-3**-0.5,-3**-0.5],[3**-0.5,3**-0.5,-3**-0.5],[-3**-0.5,3**-0.5,-3**-0.5],[-3**-0.5,-3**-0.5,3**-0.5],[3**-0.5,-3**-0.5,3**-0.5],[3**-0.5,3**-0.5,3**-0.5],[-3**-0.5,3**-0.5,3**-0.5]]

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

# Calculates trilinear shape function values
def Trilin_Interpolation(tsi,eta,zeta): # tsi, eta and zeta are the natural coordinates
    N1 = float(0.125*(1-tsi)*(1-eta)*(1-zeta))
    N2 = float(0.125*(1+tsi)*(1-eta)*(1-zeta))
    N3 = float(0.125*(1+tsi)*(1+eta)*(1-zeta))
    N4 = float(0.125*(1-tsi)*(1+eta)*(1-zeta))
    N5 = float(0.125*(1-tsi)*(1-eta)*(1+zeta))
    N6 = float(0.125*(1+tsi)*(1-eta)*(1+zeta))
    N7 = float(0.125*(1+tsi)*(1+eta)*(1+zeta))
    N8 = float(0.125*(1-tsi)*(1+eta)*(1+zeta))   
    return [N1,N2,N3,N4,N5,N6,N7,N8] # Returns a list of shape function values

# Sorts nodes along an edge using their coordinates
def SortListofNodes1D(faceN,coordinate): # Coordinates: 0 for x; 1 for y; 2 for z
    newlist = []
    oldlist = []
    for n_face_nodes in range(len(faceN)): # Obtain a list of coordinates for all the nodes
        oldlist.append(RVENodalCoord[faceN[n_face_nodes]][coordinate])
    
    orderedlist = sorted(oldlist) # Sort the nodal coordinates
    for n_oldlist_nodes in range(len(orderedlist)): # Sort the nodes based on the sorted list of coodinates
        ind = oldlist.index(orderedlist[n_oldlist_nodes])
        newlist.append(faceN[ind])
    
    return newlist # Returns a list of sorted nodes in terms of node numbers ready to be called with python (already -1)

# Sorts nodes on a face using their coordinates
def SortListofNodes2D(faceN,coord1,coord2): # 0 for x; 1 for y; 2 for z; gives node numbers ready to be called with python (already -1)
    oldlistC1 = [] 
    newlistN = []
    
    # Sort by first coordinate
    for n_face_nodes in range(len(faceN)): # Obtain a list of coordinates along the first direction for all the nodes
        oldlistC1.append(RVENodalCoord[faceN[n_face_nodes]][coord1])
        
    newlistC1 = sorted(list(dict.fromkeys(oldlistC1))) # Sort the nodal coordinates along the first direction, removing repeats

    # Sort by second coordinate for each new unique first coordinate
    for n_newlistC1_coords in range(len(newlistC1)): 
        C1 = newlistC1[n_newlistC1_coords]
        sublistN = []
        sublistC2 = []
        for n2_face_nodes in range(len(faceN)): # Obtain the node labels and second coordinates for nodes with the first coordinate equal to C1 
            C1N = RVENodalCoord[faceN[n2_face_nodes]][coord1]
            
            if (C1N==C1):
                sublistN.append(faceN[n2_face_nodes])
                sublistC2.append(RVENodalCoord[faceN[n2_face_nodes]][coord2])
                
        newlistC2 = sorted(sublistC2) # Sort the list of second coordinates
        for n_newlistC1_nodes in range(len(sublistN)):
            Nindex = sublistC2.index(newlistC2[n_newlistC1_nodes])
            newlistN.append(sublistN[Nindex]) 
    
    return newlistN # Returns a list of sorted nodes

# Removes nodes from a set that matches a given set of coordinates
def ExcludeNodes(faceN,coord1,coord2,coord3): # coord1, coord2 and coord3 are the coordinates to exclude, to be given as lists
    newlistN = []
    
    for n_face_nodes in range(len(faceN)):
        if RVENodalCoord[faceN[n_face_nodes]][0] not in coord1: # If the node does not have an x coordinate matching any in the list
            if RVENodalCoord[faceN[n_face_nodes]][1] not in coord2: # If the node does not have a y coordinate matching any in the list
                if RVENodalCoord[faceN[n_face_nodes]][2] not in coord3: # If the node does not have a z coordinate matching any in the list
                    newlistN.append(faceN[n_face_nodes])
    
    return newlistN # Returns a list of nodes that do not have coordinates matching the given set        


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
            if n_terms==0: # Start of the list
                Temp2 = str(n_fulllist_terms)
                n_terms = n_terms+1
            elif n_terms==16: # 16th term of the list, where a new line is required
                inp1.insert(n_inp1_lines+2+n_lines,Temp2) # Insert the current filled line into the macroscale input file lines
                n_lines = n_lines+1
                Temp2 = str(n_fulllist_terms) # Start a new line with the next term
                n_terms = 1
            else: # All other terms
                Temp2 = Temp2+', '+str(n_fulllist_terms)
                n_terms = n_terms+1
        inp1.insert(n_inp1_lines+2+n_lines,Temp2) # Insert the final line into the macroscale input file lines
        inp1[n_inp1_lines] = inp1[n_inp1_lines][0:len(inp1[n_inp1_lines])-10] # Remove the 'generate' keyword from the opening line of the set or surface list
        del inp1[n_inp1_lines+1] # Remove the original compacted list

# RVE input file
inp2 = []    
f2 = open(RVEInpName,'r') # Open the RVE input file as f2

while 1: # Read the RVE input file line by line and store it
    line = f2.readline()
    if not line:
        break
    line = line.strip() # Removes additional white spaces on left and right
    inp2.append(line)
    
f2.close()

# Removing 'generate' for easier processing
for n_inp2_lines in reversed(range(len(inp2))): # Read through the RVE input file lines, done in reverse to avoid issues with line number when expanding the 'generate' keyword
    if (inp2[n_inp2_lines].count('generate')!=0): # Lines that compact node or element lists and contain the 'generate' keyword
        Temp = (inp2[n_inp2_lines+1].replace(',',' ')).split() # Split the key numbers in the compacted list containing the start, end and increment for the list
        n_terms = 0 # Term counter
        n_lines = 0 # Extra line counter
        for n_fulllist_terms in range(int(Temp[0]),int(Temp[1])+1,int(Temp[2])):
            if n_terms==0: # Start of the list
                Temp2 = str(n_fulllist_terms)
                n_terms = n_terms+1
            elif n_terms==16: # 16th term of the list, where a new line is required
                inp2.insert(n_inp2_lines+2+n_lines,Temp2) # Insert the current filled line into the RVE input file lines
                n_lines = n_lines+1
                Temp2 = str(n_fulllist_terms) # Start a new line with the next term
                n_terms = 1
            else: # All other terms
                Temp2 = Temp2+', '+str(n_fulllist_terms)
                n_terms = n_terms+1
        inp2.insert(n_inp2_lines+2+n_lines,Temp2) # Insert the final line into the RVE input file lines
        inp2[n_inp2_lines] = inp2[n_inp2_lines][0:len(inp2[n_inp2_lines])-10] # Remove the 'generate' keyword from the opening line of the set or surface list
        del inp2[n_inp2_lines+1] # Remove the original compacted list

# Extracting macroscale element info from old inp file
MacroNodalConnect,MacroNodalCoord = [],[]
Search(inp1,'*Element',MacroNodalConnect,0) # Search for macroscale elements' nodal connectivity and store it
Search(inp1,'*Node',MacroNodalCoord,1) # Search for macroscale nodal coordinates and store it

StartConst = 'a' # Marker to indicate start of Constraints, if any
for n_inp1_lines in range(len(inp1)): # Search through all lines in the macroscale input file
    if (inp1[n_inp1_lines].count('*Instance,'))!=0: # Find the line containing the macroscale Instance, using the keyword '*Instance,'
        Line = inp1[n_inp1_lines].split(', ') 
        MacroInstName = Line[1][5:] # Extract the macroscale Instance name
    if (inp1[n_inp1_lines].count('*Material,'))!=0: # Find the line containing the macroscale Material definition, using the keyword '*Material,'
        inp1[n_inp1_lines+2] = '1e-10,1e-10' # Replace the macroscale Material with null definitions
    if (inp1[n_inp1_lines].count('** Constraint'))!=0 and (StartConst == 'a'): # Find the line defining the start of macroscale Constraints if it has not been found, using the keyword '** Constraint'
        StartConst = n_inp1_lines # Mark the starting line for the macroscale Constraints

# Extracting RVE info from old inp file
RVENodalCoord = []
Search(inp2,'*Node',RVENodalCoord,1) # Search for RVE nodal coordinates and store it

Sections = [] # List to store first line number of each RVE Sections
Materials = [] # List to store first line number of each RVE Material
StartEle = 'a' # Marker to indicate start of RVE Elements
for n_inp2_lines in range(len(inp2)): # Search through all lines in the RVE input file
    if (inp2[n_inp2_lines].count('*Element'))!=0 and (StartEle == 'a'): # Find the line defining the start of RVE elements if it has not been found, using the keyword '*Element'
        StartEle = n_inp2_lines # Mark the starting line for the RVE elements
    if (inp2[n_inp2_lines].count('** Section'))!=0: # Find the lines defining the start of each RVE Section, using the keyword '** Section'
        Sections.append(n_inp2_lines) # Store the first line number for each RVE Section
    if (inp2[n_inp2_lines].count('*Material'))!=0: # Find the lines defining the start of each RVE Material, using the keyword '*Material'
        Materials.append(n_inp2_lines) # Store the first line number for each RVE Material

RVEMats = open('RVEMats.dat','w') # Open a temporary file to store information on RVE Materials
for n_inp2_lines in range(len(Materials)): # Looping through each RVE Material
    for n2_inp2_lines in range(Materials[n_inp2_lines]+1,len(inp2)): # Search for the end of each RVE Material definition, starting from the line after the first line
        if (inp2[n2_inp2_lines].count('*Material'))!=0 or (inp2[n2_inp2_lines].count('**'))!=0: # Start of next RVE Material or other definition, marked with the keywords '*Material' or '**' respectively 
            MatEnd = n2_inp2_lines # Mark the end of the current RVE Material
            break
        MatEnd = n2_inp2_lines+1 # Current RVE Material ends with the RVE input file if no further information provided beyond RVE Materials, such as Step
    for n2_inp2_lines in range(Materials[n_inp2_lines],MatEnd):
        print>>RVEMats,inp2[n2_inp2_lines] # Print all RVE Materials to the temporary file
RVEMats.close()


### Processing the macroscale part information
# Sorting the nodal connectivity to match with DFE2 conventions
N_macro_eles = len(MacroNodalConnect)
NodalConnect = []
NodalCoordX = []
NodalCoordY = []
NodalCoordZ = []
for i in range(N_macro_eles):
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
for i in range(N_macro_eles):
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
for i in range(N_macro_eles):
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
        J_RVE = (abs(np.linalg.det(J)/(B_RVE*H_RVE*T_RVE)))**(1.0/3.0)
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

Proposed Revisions (yet to be implemented)
Generalise RVE Sections portion to account for possible material orientations
Account for RVE sets and surfaces from the microscale input file, found from Part?
Account for RVE level interactions and constraints, found from Instance, including the rearrangement? adopt from code for Yuhao and ZB
Account for BCs applied at Initial step when writing RVE materials, which appears before the ------- line? adopt from code for Yuhao and ZB
Allow for (partial) reduced integration by changing the list of GP? adopt from code for ZB
Account for multiple macroscale parts, multiple types of RVE and different RVE orientations? adopt from code for Yuhao
Relabel J as J^T for better clarity
'''

















