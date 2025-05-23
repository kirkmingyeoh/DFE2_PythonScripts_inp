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
for n_inp2_lines in range(len(Materials)): # Loop through each RVE Material
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
N_macro_eles = len(MacroNodalConnect) # Total number of macroscale elements
NodalConnect = [] # List to store sorted macroscale element nodal connectivities
NodalCoordX = [] # List to store nodal x coordinates based on sorted macroscale element nodal connectivities
NodalCoordY = [] # List to store nodal y coordinates based on sorted macroscale element nodal connectivities
NodalCoordZ = [] # List to store nodal z coordinates based on sorted macroscale element nodal connectivities
for n_macro_eles in range(N_macro_eles): # Loop through all macroscale elements
    Nodes = [] # List to store nodes labels for all nodes of the current macroscale element
    TempCoordX = [] # List to store x coordinates for all nodes of the current macroscale element
    TempCoordY = [] # List to store y coordinates for all nodes of the current macroscale element
    TempCoordZ = [] # List to store z coordinates for all nodes of the current macroscale element
    Disp = [[] for x in range(len(MacroNodalConnect[n_macro_eles]))] # List to store displacement vector of all nodes wrt the centroid of the current macroscale element
    
    # Obtaining nodal information from this particular element
    for n_macroele_nodes in range(len(MacroNodalConnect[n_macro_eles])): # Loop through all nodes of the current macroscale element
        Nodes.append(MacroNodalConnect[n_macro_eles][n_macroele_nodes]) # Store the node number 
        TempCoordX.append(MacroNodalCoord[MacroNodalConnect[n_macro_eles][n_macroele_nodes]][0]) # Store the x coordinate
        TempCoordY.append(MacroNodalCoord[MacroNodalConnect[n_macro_eles][n_macroele_nodes]][1]) # Store the y coordinate
        TempCoordZ.append(MacroNodalCoord[MacroNodalConnect[n_macro_eles][n_macroele_nodes]][2]) # Store the z coordinate
    
    # Obtaining centroid of the element    
    X0 = sum(TempCoordX)/len(MacroNodalConnect[n_macro_eles]) # x coordinate of centroid
    Y0 = sum(TempCoordY)/len(MacroNodalConnect[n_macro_eles]) # y coordinate of centroid
    Z0 = sum(TempCoordZ)/len(MacroNodalConnect[n_macro_eles]) # z coordinate of centroid
    
    # Obtaining displacement vector of each node relative to centroid
    for n_macroele_nodes in range(len(MacroNodalConnect[n_macro_eles])): # Loop through all nodes of the current macroscale element
        Disp[n_macroele_nodes].append(MacroNodalCoord[MacroNodalConnect[n_macro_eles][n_macroele_nodes]][0]-X0) # x component of displacement vector wrt centroid
        Disp[n_macroele_nodes].append(MacroNodalCoord[MacroNodalConnect[n_macro_eles][n_macroele_nodes]][1]-Y0) # y component of displacement vector wrt centroid
        Disp[n_macroele_nodes].append(MacroNodalCoord[MacroNodalConnect[n_macro_eles][n_macroele_nodes]][2]-Z0) # z component of displacement vector wrt centroid
        
    Order = [0 for x in range(len(MacroNodalConnect[n_macro_eles]))] # List to store the sorted nodes
    # The 8 nodes are first split into two groups of 4, based on whether their \zeta natural coordinates are positive or negative
    # This works on the assumption that the macroscale element is not too distorted relative to the master cubic element 
    # After sorting based on their \zeta natural coordinate, the groups of nodes are each sorted based on their \tsi and \eta natural coordinates
    # When sorting in plane within each group, the first node has natural coordinates [-1,-1] and the order goes in a counter-clockwise fashion
    
    Set = [] # List to store nodes 1-4 with \zeta < 0
    for n_macroele_nodes in range(len(MacroNodalConnect[n_macro_eles])): # Loop through all nodes of the current macroscale element
        if Disp[n_macroele_nodes][2]<0: # If the z component of the displacement vector is negative
            Set.append(n_macroele_nodes) # Add the node to the set
    for n_macroele_nodes in range(len(MacroNodalConnect[n_macro_eles])/2): # Loop across the 4 chosen nodes with \zeta<0
        # Sort the 4 nodes on the face based on their x and y components of the displacement vector
        if Disp[Set[n_macroele_nodes]][1]<0: 
            if Disp[Set[n_macroele_nodes]][0]<0:
                Order[0] = Set[n_macroele_nodes] # The first node has x and y components which are both negative
            else:
                Order[1] = Set[n_macroele_nodes] # The second node has a positive x component and a negative y component
        elif Disp[Set[n_macroele_nodes]][1]>0:
            if Disp[Set[n_macroele_nodes]][0]<0:
                Order[3] = Set[n_macroele_nodes] # The fourth node has a negative x component and a positive y component 
            else:
                Order[2] = Set[n_macroele_nodes] # The third node has x and y components which are both positive 
                
    Set = [] # List to store nodes 5-8 with \zeta > 0
    for n_macroele_nodes in range(len(MacroNodalConnect[n_macro_eles])): # Loop through all nodes of the current macroscale element
        if Disp[n_macroele_nodes][2]>0: # If the z component of the displacement vector is positive
            Set.append(n_macroele_nodes) # Add the node to the set        
    for n_macroele_nodes in range(len(MacroNodalConnect[n_macro_eles])/2): # Loop across the 4 chosen nodes with \zeta<0
        # Sort the 4 nodes on the face based on their x and y components of the displacement vector
        if Disp[Set[n_macroele_nodes]][1]<0:
            if Disp[Set[n_macroele_nodes]][0]<0:
                Order[4] = Set[n_macroele_nodes] # The first node has x and y components which are both negative
            else:
                Order[5] = Set[n_macroele_nodes] # The second node has a positive x component and a negative y component
        elif Disp[Set[n_macroele_nodes]][1]>0:
            if Disp[Set[n_macroele_nodes]][0]<0:
                Order[7] = Set[n_macroele_nodes] # The fourth node has a negative x component and a positive y component 
            else:
                Order[6] = Set[n_macroele_nodes] # The third node has x and y components which are both positive
                
    SortedConnect = [] # List to store sorted node labels for all nodes of the current macroscale element
    SortedCoordX = [] # List to store sorted x coordinates for all nodes of the current macroscale element
    SortedCoordY = [] # List to store sorted y coordinates for all nodes of the current macroscale element
    SortedCoordZ = [] # List to store sorted z coordinates for all nodes of the current macroscale element
    
    for n_macroele_nodes in range(len(MacroNodalConnect[n_macro_eles])): # Loop through all nodes of the current macroscale element
        SortedConnect.append(Nodes[Order[n_macroele_nodes]]) # Extract the corresponding node label
        SortedCoordX.append(MacroNodalCoord[Nodes[Order[n_macroele_nodes]]][0]) # Extract the corresponding node x coordinate
        SortedCoordY.append(MacroNodalCoord[Nodes[Order[n_macroele_nodes]]][1]) # Extract the corresponding node y coordinate
        SortedCoordZ.append(MacroNodalCoord[Nodes[Order[n_macroele_nodes]]][2]) # Extract the corresponding node z coordinate
        
    NodalConnect.append(SortedConnect) # Add the sorted node for the current macroscale element to the main list
    NodalCoordX.append(SortedCoordX) # Add the sorted x coordinates for the nodes of the current macroscale element to the main list
    NodalCoordY.append(SortedCoordY) # Add the sorted y coordinates for the nodes of the current macroscale element to the main list
    NodalCoordZ.append(SortedCoordZ) # Add the sorted z coordinates for the nodes of the current macroscale element to the main list


### Processing the RVE part information
# Finding the smallest and largest nodal coordinate in all 3 directions
# Assumes a cuboidal RVE aligned along the global x, y and z directions
RVE_ListX = [] # Temporary list to store x coordinates of all RVE nodes
RVE_ListY = [] # Temporary list to store y coordinates of all RVE nodes
RVE_ListZ = [] # Temporary list to store z coordinates of all RVE nodes
for n_RVE_nodes in range(len(RVENodalCoord)): # Loop through all RVE nodes
    RVE_ListX.append(RVENodalCoord[n_RVE_nodes][0]) # Store the x coordinates
    RVE_ListY.append(RVENodalCoord[n_RVE_nodes][1]) # Store the y coordinates
    RVE_ListZ.append(RVENodalCoord[n_RVE_nodes][2]) # Store the z coordinates
xMin = min(RVE_ListX) # Smallest x coordinate
xMax = max(RVE_ListX) # Largest x coordinate
yMin = min(RVE_ListY) # Smallest y coordinate
yMax = max(RVE_ListY) # Largest y coordinate
zMin = min(RVE_ListZ) # Smallest z coordinate
zMax = max(RVE_ListZ) # Largest z coordinate
del RVE_ListX # Remove the temporary list
del RVE_ListY # Remove the temporary list
del RVE_ListZ # Remove the temporary list

# Sorting the RVE boundary nodes
FaceLNodes,FaceRNodes,FaceBaNodes,FaceFNodes,FaceBNodes,FaceTNodes = [],[],[],[],[],[] # List to store the nodes on the RVE boundary faces
EdgeLNodes, EdgeRNodes, EdgeBNodes, EdgeTNodes = [],[],[],[] # List to store the nodes on the RVE boundary edges attached to the back face
Edge1Nodes, Edge2Nodes, Edge3Nodes, Edge4Nodes = [],[],[],[] # List to store the nodes on the RVE boundary edges
for n_RVE_nodes in range(len(RVENodalCoord)): # Loop through all RVE nodes
    if RVENodalCoord[n_RVE_nodes][0] == xMin: # If the nodal x coordinate matches the RVE smallest x coordinate
        FaceLNodes.append(n_RVE_nodes) # Store the node the left face list
    if RVENodalCoord[n_RVE_nodes][0] == xMax:  # If the nodal x coordinate matches the RVE largest x coordinate
        FaceRNodes.append(n_RVE_nodes) # Store the node the right face list
    if RVENodalCoord[n_RVE_nodes][1] == yMin: # If the nodal y coordinate matches the RVE smallest y coordinate
        FaceBaNodes.append(n_RVE_nodes) # Store the node the back face list
    if RVENodalCoord[n_RVE_nodes][1] == yMax:  # If the nodal y coordinate matches the RVE largest y coordinate
        FaceFNodes.append(n_RVE_nodes) # Store the node the front face list
    if RVENodalCoord[n_RVE_nodes][2] == zMin: # If the nodal z coordinate matches the RVE smallest z coordinate
        FaceBNodes.append(n_RVE_nodes) # Store the node the bottom face list
    if RVENodalCoord[n_RVE_nodes][2] == zMax:  # If the nodal z coordinate matches the RVE largest z coordinate
        FaceTNodes.append(n_RVE_nodes) # Store the node the top face list
for n_FaceBa_nodes in range(len(FaceBaNodes)): # Loop through all nodes in the back face list
    if FaceBaNodes[n_FaceBa_nodes] in FaceLNodes: # If the node is also in the left face list
        EdgeLNodes.append(FaceBaNodes[n_FaceBa_nodes]) # Store the node in the left edge list
        if FaceBaNodes[n_FaceBa_nodes] in FaceBNodes: # If the node is also in the bottom face list
            V1 = FaceBaNodes[n_FaceBa_nodes] # Store the node as vertex V1
        if FaceBaNodes[n_FaceBa_nodes] in FaceTNodes: # If the node is also in the top face list
            V4 = FaceBaNodes[n_FaceBa_nodes] # Store the node as vertex V4
    if FaceBaNodes[n_FaceBa_nodes] in FaceRNodes: # If the node is also in the right face list
        EdgeRNodes.append(FaceBaNodes[n_FaceBa_nodes]) # Store the node in the right edge list
        if FaceBaNodes[n_FaceBa_nodes] in FaceBNodes: # If the node is also in the bottom face list 
            V2 = FaceBaNodes[n_FaceBa_nodes] # Store the node as vertex V2
        if FaceBaNodes[n_FaceBa_nodes] in FaceTNodes: # If the node is also in the top face list
            V3 = FaceBaNodes[n_FaceBa_nodes] # Store the node as vertex V3
    if FaceBaNodes[n_FaceBa_nodes] in FaceBNodes: # If the node is also in the bottom face list
        EdgeBNodes.append(FaceBaNodes[n_FaceBa_nodes]) # Store the node in the bottom edge list
    if FaceBaNodes[n_FaceBa_nodes] in FaceTNodes: # If the node is also in the top face list
        EdgeTNodes.append(FaceBaNodes[n_FaceBa_nodes]) # Store the node in the top edge list
for n_FaceB_nodes in range(len(FaceBNodes)): # If the node is also in the bottom face list
    if FaceBNodes[n_FaceB_nodes] in FaceLNodes: # If the node is also in the left face list
        Edge1Nodes.append(FaceBNodes[n_FaceB_nodes]) # Store the node in the 1 edge list
        if FaceBNodes[n_FaceB_nodes] in FaceFNodes: # If the node is also in the front face list
            V5 = FaceBNodes[n_FaceB_nodes] # Store the node as vertex V5
    if FaceBNodes[n_FaceB_nodes] in FaceRNodes: # If the node is also in the right face list
        Edge2Nodes.append(FaceBNodes[n_FaceB_nodes]) # Store the node in the 2 edge list
        if FaceBNodes[n_FaceB_nodes] in FaceFNodes: # If the node is also in the front face list
            V6 = FaceBNodes[n_FaceB_nodes] # Store the node as vertex V6
for n_FaceT_nodes in range(len(FaceTNodes)): # If the node is also in the top face list
    if FaceTNodes[n_FaceT_nodes] in FaceLNodes: # If the node is also in the left face list
        Edge4Nodes.append(FaceTNodes[n_FaceT_nodes]) # Store the node in the 4 edge list
        if FaceTNodes[n_FaceT_nodes] in FaceFNodes: # If the node is also in the front face list
            V8 = FaceTNodes[n_FaceT_nodes] # Store the node as vertex V8
    if FaceTNodes[n_FaceT_nodes] in FaceRNodes: # If the node is also in the right face list
        Edge3Nodes.append(FaceTNodes[n_FaceT_nodes]) # Store the node in the 3 edge list
        if FaceTNodes[n_FaceT_nodes] in FaceFNodes: # If the node is also in the front face list
            V7 = FaceTNodes[n_FaceT_nodes] # Store the node as vertex V7
    
# Sorting RVE dimensions and offsets
B_RVE = xMax - xMin # RVE dimension along x direction
H_RVE = yMax - yMin # RVE dimension along y direction
T_RVE = zMax - zMin # RVE dimension along z direction
OffsetX = (xMax + xMin)/2.0 # RVE centroid x coordinate
OffsetY = (yMax + yMin)/2.0 # RVE centroid y coordinate
OffsetZ = (zMax + zMin)/2.0 # RVE centroid z coordinate
Offset = [OffsetX,OffsetY,OffsetZ] # RVE centroid coordinates

# Adjusting RVE nodal coordinates to correspond to a part centered at the origin
for n_RVE_node in RVENodalCoord: # Loop through all RVE nodes
    for n_nodal_coord in range(3): # Loop through all coordinates of each node
        n_RVE_node[n_nodal_coord] = n_RVE_node[n_nodal_coord] - Offset[n_nodal_coord] # Subtracting the offset from all RVE nodal coordinates


### Generating the RVE placement in the macroscale mesh
RVEParts = open('RVEParts.dat','w') # Open a temporary file to store information on RVE Parts
Insts = open('Insts.dat','w') # Open a temporary file to store information on RVE Instances
SF = [] # List to store required RVE volume scaling factors
for n_macro_eles in range(N_macro_eles): # Loop through all macroscale elements
    # Calculate the mapping function between natural and global coordinates
    # x = a0 + a1*tsi + a2*eta + a3*zeta + a4*tsi*eta + a5*tsi*zeta + a6*eta*zeta + a7*tsi*eta*zeta
    # y = b0 + b1*tsi + b2*eta + b3*zeta + b4*tsi*eta + b5*tsi*zeta + b6*eta*zeta + b7*tsi*eta*zeta
    # z = d0 + d1*tsi + d2*eta + d3*zeta + d4*tsi*eta + d5*tsi*zeta + d6*eta*zeta + d7*tsi*eta*zeta
    # C is matrix where [C]*{a0;a1;a2;a3;a4;a5;a6;a7} = {x1;x2;x3;x4;x5;x6;x7;x8}, repeated for b-y and d-z, based on the equations above
    # 1,2,3,4,5,6,7,8 for x,y,z refer to the 8 macroscale nodes
    C = np.array([
            [1,-1,-1,-1,1,1,1,-1],
            [1,1,-1,-1,-1,-1,1,1],
            [1,1,1,-1,1,-1,-1,-1],
            [1,-1,1,-1,-1,1,-1,1],
            [1,-1,-1,1,1,-1,-1,1],
            [1,1,-1,1,-1,1,-1,-1],
            [1,1,1,1,1,1,1,1],
            [1,-1,1,1,-1,-1,1,-1],]) # Matrix of natural coordinates
    C_inv = np.linalg.inv(C) # Inverse of matrix C
    [a0,a1,a2,a3,a4,a5,a6,a7] = np.dot(C_inv,NodalCoordX[n_macro_eles]) # Coefficient a
    [b0,b1,b2,b3,b4,b5,b6,b7] = np.dot(C_inv,NodalCoordY[n_macro_eles]) # Coefficient b
    [d0,d1,d2,d3,d4,d5,d6,d7] = np.dot(C_inv,NodalCoordZ[n_macro_eles]) # Coefficient d
    
    for n_macroele_GPs in range(len(GP)): # Loop through all integration points of the macroscale element
        [tsi,eta,zeta] = GP[n_macroele_GPs] # Natural coordinates of the current integration point
        N1,N2,N3,N4,N5,N6,N7,N8 = Trilin_Interpolation(tsi,eta,zeta) # Shape function values corresponding to the current integration point
        
        RVE_X = N1*NodalCoordX[n_macro_eles][0] + N2*NodalCoordX[n_macro_eles][1] + N3*NodalCoordX[n_macro_eles][2] + N4*NodalCoordX[n_macro_eles][3] + N5*NodalCoordX[n_macro_eles][4] + N6*NodalCoordX[n_macro_eles][5] + N7*NodalCoordX[n_macro_eles][6] + N8*NodalCoordX[n_macro_eles][7] # x coordinate of the current integration point 
        RVE_Y = N1*NodalCoordY[n_macro_eles][0] + N2*NodalCoordY[n_macro_eles][1] + N3*NodalCoordY[n_macro_eles][2] + N4*NodalCoordY[n_macro_eles][3] + N5*NodalCoordY[n_macro_eles][4] + N6*NodalCoordY[n_macro_eles][5] + N7*NodalCoordY[n_macro_eles][6] + N8*NodalCoordY[n_macro_eles][7] # y coordinate of the current integration point
        RVE_Z = N1*NodalCoordZ[n_macro_eles][0] + N2*NodalCoordZ[n_macro_eles][1] + N3*NodalCoordZ[n_macro_eles][2] + N4*NodalCoordZ[n_macro_eles][3] + N5*NodalCoordZ[n_macro_eles][4] + N6*NodalCoordZ[n_macro_eles][5] + N7*NodalCoordZ[n_macro_eles][6] + N8*NodalCoordZ[n_macro_eles][7] # z coordinate of the current integration point
        
        J = np.array([
                [a1+a4*eta+a5*zeta+a7*eta*zeta,b1+b4*eta+b5*zeta+b7*eta*zeta,d1+d4*eta+d5*zeta+d7*eta*zeta],
                [a2+a4*tsi+a6*zeta+a7*tsi*zeta,b2+b4*tsi+b6*zeta+b7*tsi*zeta,d2+d4*tsi+d6*zeta+d7*tsi*zeta],
                [a3+a5*tsi+a6*eta+a7*tsi*eta,b3+b5*tsi+b6*eta+b7*tsi*eta,d3+d5*tsi+d6*eta+d7*tsi*eta]]) # Jacobian matrix of the current integration point
       
        J_RVE = (abs(np.linalg.det(J)/(B_RVE*H_RVE*T_RVE)))**(1.0/3.0) # Scaling factor for RVE volume at the current integration point
        
        if round(J_RVE,5) in SF: # Check if this volume scaling factor has been used before, reuse the same RVE Part if so == 
            Ind_SF = SF.index(round(J_RVE,5)) # Index of the first instance of this volume scaling factor
        else: # Create a new RVE with this current volume scaling factor if it has not been used before
            Ind_SF = len(SF) # Index the new volume scaling factor as the next term in the list
            SF.append(round(J_RVE,5)) # Add the new volume scaling factor to the list

            # Print headers of the new Part
            print>>RVEParts,'**'
            print>>RVEParts,'*Part, name=RVE-'+str(Ind_SF+1) # Name the part with the new index
            print>>RVEParts,'*Node'

            # Print the nodal coordinates of the new Part
            for n_RVE_nodes in range(len(RVENodalCoord)): # Loop through all RVE nodes
                print>>RVEParts,str(n_RVE_nodes+1)+', '+str(RVENodalCoord[n_RVE_nodes][0]*J_RVE)+', '+str(RVENodalCoord[n_RVE_nodes][1]*J_RVE)+', '+str(RVENodalCoord[n_RVE_nodes][2]*J_RVE) # RVE nodal coordinate values are scaled based on the volume scaling factor required

            # Print the elements of the new Part
            for n_inp2_lines in range(StartEle,Sections[0]):  # Loop through the RVE input file lines from the start of elements till the first Section
                print>>RVEParts,inp2[n_inp2_lines]

            # Print the Sections of the new Part
            for n_RVE_sections in range(len(Sections)):  # Loop through the RVE Sections
                print>>RVEParts,inp2[Sections[n_RVE_sections]]+'-'+str(Ind_SF+1) # Print the Section name, with the SF index added
                print>>RVEParts,inp2[Sections[n_RVE_sections]+1] # Print the next line of the Section
                print>>RVEParts,','
                
            # Print the end of the new Part
            print>>RVEParts,'*End Part'

        # Print the Instance of the new Part
        print>>Insts,'*Instance, name=Ele'+str(n_macro_eles+1)+'-RVE'+str(n_macroele_GPs+1)+', part=RVE-'+str(Ind_SF+1)  # Name of the Instance and the Part it refers to
        print>>Insts,str(RVE_X)+', '+str(RVE_Y)+', '+str(RVE_Z) # Translation required for the Instance, from the origin to the macroscale integration point coordinates
        print>>Insts,'*End Instance'
        print>>Insts,'**'

RVEParts.close()
Insts.close()


### Setting up the MPCs
Sets = open('Sets.dat','w') # Open a temporary file to store information on DFE2 Sets
Eqns = open('Eqns.dat','w') # Open a temporary file to store information on DFE2 MPCs

# Pairing the nodes
# Assume the RVE mesh is perfectly periodic
# Left and right edges
EdgeLNodes = TakeVertexOut(SortListofNodes1D(EdgeLNodes,2)) # Sort the left edge nodes based on their z coordinates and remove the two corner nodes
EdgeRNodes = TakeVertexOut(SortListofNodes1D(EdgeRNodes,2)) # Sort the right edge nodes based on their z coordinates and remove the two corner nodes

# Bottom and top edges
EdgeBNodes = TakeVertexOut(SortListofNodes1D(EdgeBNodes,0)) # Sort the bottom edge nodes based on their x coordinates and remove the two corner nodes
EdgeTNodes = TakeVertexOut(SortListofNodes1D(EdgeTNodes,0)) # Sort the top edge nodes based on their x coordinates and remove the two corner nodes

# 1, 2, 3, and 4 edges
Edge1Nodes = TakeVertexOut(SortListofNodes1D(Edge1Nodes,1)) # Sort the 1 edge nodes based on their y coordinates and remove the two corner nodes
Edge2Nodes = TakeVertexOut(SortListofNodes1D(Edge2Nodes,1)) # Sort the 2 edge nodes based on their y coordinates and remove the two corner nodes
Edge3Nodes = TakeVertexOut(SortListofNodes1D(Edge3Nodes,1)) # Sort the 3 edge nodes based on their y coordinates and remove the two corner nodes
Edge4Nodes = TakeVertexOut(SortListofNodes1D(Edge4Nodes,1)) # Sort the 4 edge nodes based on their y coordinates and remove the two corner nodes

# Back and front faces
FaceBaNodes = SortListofNodes2D(FaceBaNodes,0,2) # Sort the back face nodes based on their x and z coordinates
FaceFNodes = SortListofNodes2D(FaceFNodes,0,2) # Sort the front face nodes based on their x and z coordinates
PairingFacesBaF = [] # List to store the paired back and front face nodes
for n_FaceBa_nodes in range(len(FaceBaNodes)): # Loop through the back face nodes
    Temp = [] # Temporary list to store each pair of back and front face nodes
    Temp.append(FaceBaNodes[n_FaceBa_nodes]) # Back face node number of the pair
    x1 = RVENodalCoord[FaceBaNodes[n_FaceBa_nodes]][0] # x coordinate of the back face node of the pair
    z1 = RVENodalCoord[FaceBaNodes[n_FaceBa_nodes]][2] # z coordinate of the back face node of the pair

    # Find the front face node for pairing based on the shortest in-plane distance with the current back face node
    Dist = [] # Temporary list to store in-plane distance between each front face node with the current back face node
    for n_FaceF_nodes in range(len(FaceFNodes)): # Loop through the front face nodes
        x2 = RVENodalCoord[FaceFNodes[n_FaceF_nodes]][0] # x coordinate of the front face node
        z2 = RVENodalCoord[FaceFNodes[n_FaceF_nodes]][2] # z coordinate of the front face node
        
        Dist.append(math.sqrt(pow(x2-x1,2)+pow(z2-z1,2))) # In-plane distance between current front face node with the current back face node
        
    Ind_dist = Dist.index(min(Dist)) # Index the front face node with the smallest in-plane distance from the current back face node as the node to be paired
    Temp.append(FaceFNodes[Ind_dist]) # Front face node number of the pair
    FaceFNodes.pop(Ind_dist) # Remove the paired front face node from the list to reduce pairing time for subsequent back face nodes
    PairingFacesBaF.append(Temp) # Store the pair of node numbers into the back-front list

FaceLNodes = SortListofNodes2D(ExcludeNodes(FaceLNodes,[],[RVENodalCoord[V3][1],RVENodalCoord[V6][1]],[RVENodalCoord[V3][2],RVENodalCoord[V6][2]]),1,2) # Sort the left face nodes based on their y and z coordinates, excluding nodes on the edges using the y and z coordinates of V3 and V6
FaceRNodes = SortListofNodes2D(ExcludeNodes(FaceRNodes,[],[RVENodalCoord[V3][1],RVENodalCoord[V6][1]],[RVENodalCoord[V3][2],RVENodalCoord[V6][2]]),1,2) # Sort the right face nodes based on their y and z coordinates, excluding nodes on the edges using the y and z coordinates of V3 and V6
PairingFacesLR = [] # List to store the paired left and right face nodes
for n_FaceL_nodes in range(len(FaceLNodes)): # Loop through the left face nodes
    Temp = [] # Temporary list to store each pair of left and right face nodes
    Temp.append(FaceLNodes[n_FaceL_nodes]) # Left face node number of the pair
    y1 = RVENodalCoord[FaceLNodes[n_FaceL_nodes]][1] # y coordinate of the left face node of the pair
    z1 = RVENodalCoord[FaceLNodes[n_FaceL_nodes]][2] # z coordinate of the left face node of the pair

    # Find the right face node for pairing based on the shortest in-plane distance with the current left face node
    Dist = [] # Temporary list to store in-plane distance between each right face node with the current left face node
    for n_FaceR_nodes in range(len(FaceRNodes)): # Loop through the right face nodes
        y2 = RVENodalCoord[FaceRNodes[j]][1] # y coordinate of the right face node
        z2 = RVENodalCoord[FaceRNodes[j]][2] # z coordinate of the right face node
        
        Dist.append(math.sqrt(pow(y2-y1,2)+pow(z2-z1,2))) # In-plane distance between current right face node with the current left face node
        
    Ind_dist = Dist.index(min(Dist)) # Index the right face node with the smallest in-plane distance from the current left face node as the node to be paired
    Temp.append(FaceRNodes[Ind_dist]) # Right face node number of the pair
    FaceRNodes.pop(Ind_dist) # Remove the paired right face node from the list to reduce pairing time for subsequent left face nodes
    PairingFacesLR.append(Temp) # Store the pair of node numbers into the left-right list

FaceBNodes = SortListofNodes2D(ExcludeNodes(FaceBNodes,[RVENodalCoord[V3][0],RVENodalCoord[V8][0]],[RVENodalCoord[V3][1],RVENodalCoord[V8][1]],[]),0,1) # Sort the bottom face nodes based on their x and y coordinates, excluding nodes on the edges using the x and y coordinates of V3 and V8
FaceTNodes = SortListofNodes2D(ExcludeNodes(FaceTNodes,[RVENodalCoord[V3][0],RVENodalCoord[V8][0]],[RVENodalCoord[V3][1],RVENodalCoord[V8][1]],[]),0,1) # Sort the top face nodes based on their x and y coordinates, excluding nodes on the edges using the x and y coordinates of V3 and V8
PairingFacesBT = [] # List to store the paired bottom and top face nodes
for n_FaceB_nodes in range(len(FaceBNodes)): # Loop through the bottom face nodes
    Temp = [] # Temporary list to store each pair of bottom and top face nodes
    Temp.append(FaceBNodes[n_FaceB_nodes]) # Bottom face node number of the pair
    x1 = RVENodalCoord[FaceBNodes[n_FaceB_nodes]][0] # x coordinate of the bottom face node of the pair
    y1 = RVENodalCoord[FaceBNodes[n_FaceB_nodes]][1] # y coordinate of the bottom face node of the pair
    
    # Find the top face node for pairing based on the shortest in-plane distance with the current bottom face node
    Dist = [] # Temporary list to store in-plane distance between each top face node with the current bottom face node
    for n_FaceT_nodes in range(len(FaceTNodes)): # Loop through the top face nodes
        x2 = RVENodalCoord[FaceTNodes[n_FaceT_nodes]][0] # x coordinate of the top face node of the pair
        y2 = RVENodalCoord[FaceTNodes[n_FaceT_nodes]][1] # y coordinate of the top face node of the pair
        
        Dist.append(math.sqrt(pow(x2-x1,2)+pow(y2-y1,2))) # In-plane distance between current top face node with the current bottom face node
        
    Ind_dist = Dist.index(min(Dist)) # Index the top face node with the smallest in-plane distance from the current bottom face node as the node to be paired
    Temp.append(FaceTNodes[Ind_dist]) # Top face node number of the pair
    FaceTNodes.pop(Ind_dist) # Remove the paired top face node from the list to reduce pairing time for subsequent bottom face nodes
    PairingFacesBT.append(Temp) # Store the pair of node numbers into the bottom-top list

# Calculating the coefficients and setting up the MPCs
# For each macroscale element
for n_macro_eles in range(N_macro_eles): # Loop through all macroscale elements
    # Same mapping function between natural and global coordinates
    C = np.array([
            [1,-1,-1,-1,1,1,1,-1],
            [1,1,-1,-1,-1,-1,1,1],
            [1,1,1,-1,1,-1,-1,-1],
            [1,-1,1,-1,-1,1,-1,1],
            [1,-1,-1,1,1,-1,-1,1],
            [1,1,-1,1,-1,1,-1,-1],
            [1,1,1,1,1,1,1,1],
            [1,-1,1,1,-1,-1,1,-1],]) # Same matrix of natural coordinates based on the mapping function between natural and global coordinates
    C_inv = np.linalg.inv(C) # Inverse of matrix C
    [a0,a1,a2,a3,a4,a5,a6,a7] = np.dot(C_inv,NodalCoordX[n_macro_eles]) # Coefficient a
    [b0,b1,b2,b3,b4,b5,b6,b7] = np.dot(C_inv,NodalCoordY[n_macro_eles]) # Coefficient b
    [d0,d1,d2,d3,d4,d5,d6,d7] = np.dot(C_inv,NodalCoordZ[n_macro_eles]) # Coefficient c
    
    # Call macroscale nodes into sets
    for n_macroele_nodes in range(8): # Loop through all nodes of the macroscale element
        print>>Sets,'*Nset, nset=Ele'+str(n_macro_eles+1)+'-N'+str(n_macroele_nodes+1)+', instance='+str(MacroInstName) # Create a Set for the macroscale node
        print>>Sets,str(NodalConnect[n_macro_eles][n_macroele_nodes]+1) # Node number of the macroscale node
        
    # For each macroscale integration point
    for n_macroele_GPs in range(len(GP)): # Loop through all integration points of the macroscale element
        [tsi,eta,zeta] = GP[n_macroele_GPs] # Natural coordinates of the current integration point
        J = np.array([
                [a1+a4*eta+a5*zeta+a7*eta*zeta,b1+b4*eta+b5*zeta+b7*eta*zeta,d1+d4*eta+d5*zeta+d7*eta*zeta],
                [a2+a4*tsi+a6*zeta+a7*tsi*zeta,b2+b4*tsi+b6*zeta+b7*tsi*zeta,d2+d4*tsi+d6*zeta+d7*tsi*zeta],
                [a3+a5*tsi+a6*eta+a7*tsi*eta,b3+b5*tsi+b6*eta+b7*tsi*eta,d3+d5*tsi+d6*eta+d7*tsi*eta]]) # Jacobian matrix of the current integration point
        J_inv = np.linalg.inv(J) # Inverse of the Jacobian matrix of the current integration point
        J_RVE = (abs(np.linalg.det(J)/(B_RVE*H_RVE*T_RVE)))**(1.0/3.0) # Scaling factor for RVE volume at the current integration point

        # Shape function values at the current macroscale integration point
        Shape_fn = Trilin_Interpolation(tsi,eta,zeta)
        
        # Expressions for the shape function gradients
        # Derived by differentiating the shape functions wrt the natural coordinates
        dN1 = [-0.125*(1-eta)*(1-zeta),-0.125*(1-tsi)*(1-zeta),-0.125*(1-tsi)*(1-eta)] # Gradients of the shape function wrt tsi, eta and zeta for macroscale node 1
        dN2 = [0.125*(1-eta)*(1-zeta),-0.125*(1+tsi)*(1-zeta),-0.125*(1+tsi)*(1-eta)] # Gradients of the shape function wrt tsi, eta and zeta for macroscale node 2
        dN3 = [0.125*(1+eta)*(1-zeta),0.125*(1+tsi)*(1-zeta),-0.125*(1+tsi)*(1+eta)] # Gradients of the shape function wrt tsi, eta and zeta for macroscale node 3
        dN4 = [-0.125*(1+eta)*(1-zeta),0.125*(1-tsi)*(1-zeta),-0.125*(1-tsi)*(1+eta)] # Gradients of the shape function wrt tsi, eta and zeta for macroscale node 4
        dN5 = [-0.125*(1-eta)*(1+zeta),-0.125*(1-tsi)*(1+zeta),0.125*(1-tsi)*(1-eta)] # Gradients of the shape function wrt tsi, eta and zeta for macroscale node 5
        dN6 = [0.125*(1-eta)*(1+zeta),-0.125*(1+tsi)*(1+zeta),0.125*(1+tsi)*(1-eta)] # Gradients of the shape function wrt tsi, eta and zeta for macroscale node 6
        dN7 = [0.125*(1+eta)*(1+zeta),0.125*(1+tsi)*(1+zeta),0.125*(1+tsi)*(1+eta)] # Gradients of the shape function wrt tsi, eta and zeta for macroscale node 7
        dN8 = [-0.125*(1+eta)*(1+zeta),0.125*(1-tsi)*(1+zeta),0.125*(1-tsi)*(1+eta)] # Gradients of the shape function wrt tsi, eta and zeta for macroscale node 8
        N_NatDeriv = [dN1,dN2,dN3,dN4,dN5,dN6,dN7,dN8] # List of shape function gradients wrt tsi, eta and zeta
        N_GloDeriv = [[],[],[],[],[],[],[],[]] # List to store shape function gradients wrt to x, y and z

        # Calculate shape function gradients along x, y and z directions
         # Obtained by multiplying the inverse of the Jacobian matrix with the shape function gradients wrt to tsi, eta and zeta
        for n_macroele_nodes in range(8): # Loop through all nodes of the macroscale element
            N_GloDeriv[n_macroele_nodes] = np.dot(J_inv,np.transpose(np.array(N_NatDeriv[n_macroele_nodes]))) # Matrix multiplication between the inverse of the Jacobian matrix and shape function gradients wrt tsi, eta and zeta
        
        # Calculate the scaled RVE dimensions
        dx = B_RVE*J_RVE # Scaled dimension along x
        dy = H_RVE*J_RVE # Scaled dimension along y
        dz = T_RVE*J_RVE # Scaled dimension along z
        
        # Call Sets and set up the MPCs for the left and right edges of FaceBa
        for n_EdgeLR_nodepairs in range(len(EdgeLNodes)): # Loop through all left-right edge node pairs
            print>>Sets,'*Nset, nset=Ele'+str(n_macro_eles+1)+'-RVE'+str(n_macroele_GPs+1)+'-EdgeNodeL'+str(n_EdgeLR_nodepairs+1)+', instance=Ele'+str(n_macro_eles+1)+'-RVE'+str(n_macroele_GPs+1) # Create a Set for the left edge node of the pair
            print>>Sets,str(EdgeLNodes[n_FaceLR_nodepairs]+1) # Left face node number of the pair
            print>>Sets,'*Nset, nset=Ele'+str(n_macro_eles+1)+'-RVE'+str(n_macroele_GPs+1)+'-EdgeNodeR'+str(n_EdgeLR_nodepairs+1)+', instance=Ele'+str(n_macro_eles+1)+'-RVE'+str(n_macroele_GPs+1) # Create a Set for the right edge node of the pair
            print>>Sets,str(EdgeRNodes[n_FaceLR_nodepairs]+1) # Right edge node number of the pair
            
            for n_RVEnode_dofs in range(3): # Loop through all DOFs of the nodes
                print>>Eqns,'** Constraint: Ele'+str(n_macro_eles+1)+'-RVE'+str(n_macroele_GPs+1)+'-EdgeLR'+str(n_EdgeLR_nodepairs+1)+'-DOF'+str(n_RVEnode_dofs+1) # Create an Equation type Constraint for the DOF
                print>>Eqns,'*Equation'
                print>>Eqns,'10' # Number of terms in the Constraint
                print>>Eqns,'Ele'+str(n_macro_eles+1)+'-RVE'+str(n_macroele_GPs+1)+'-EdgeNodeR'+str(k+1)+', '+str(n_RVEnode_dofs+1)+', -1.0' # Right edge node DOF, also the first DOF which will be removed
                print>>Eqns,'Ele'+str(n_macro_eles+1)+'-RVE'+str(n_macroele_GPs+1)+'-EdgeNodeL'+str(k+1)+', '+str(n_RVEnode_dofs+1)+', 1.0' # Left edge node DOF
                for n_macroele_nodes in range(8): # Loop through all macroscale nodes
                    # Coefficient of the macroscale node term obtained by multiplying RVE dimension and macroscale shape function gradient along x direction
                    print>>Eqns,'Ele'+str(n_macro_eles+1)+'-N'+str(n_macroele_nodes+1)+', '+str(n_RVEnode_dofs+1)+', '+str(dx*N_GloDeriv[n_macroele_nodes][0]) # Macroscale node DOF
                    
        # Call Sets and set up the MPCs for bottom and top edges of FaceBa
        for n_FaceBT_nodepairs in range(len(EdgeBNodes)): # Loop through all bottom-top edge node pairs
            print>>Sets,'*Nset, nset=Ele'+str(n_macro_eles+1)+'-RVE'+str(n_macroele_GPs+1)+'-EdgeNodeB'+str(n_FaceBT_nodepairs+1)+', instance=Ele'+str(n_macro_eles+1)+'-RVE'+str(n_macroele_GPs+1) # Create a Set for the bottom edge node of the pair
            print>>Sets,str(EdgeBNodes[n_FaceBT_nodepairs]+1) # Bottom face node number of the pair
            print>>Sets,'*Nset, nset=Ele'+str(n_macro_eles+1)+'-RVE'+str(n_macroele_GPs+1)+'-EdgeNodeT'+str(n_FaceBT_nodepairs+1)+', instance=Ele'+str(n_macro_eles+1)+'-RVE'+str(n_macroele_GPs+1) # Create a Set for the top edge node of the pair
            print>>Sets,str(EdgeTNodes[n_FaceBT_nodepairs]+1) # Top face node number of the pair
            
            for n_RVEnode_dofs in range(3):
                print>>Eqns,'** Constraint: Ele'+str(n_macro_eles+1)+'-RVE'+str(n_macroele_GPs+1)+'-EdgeLR'+str(n_FaceBT_nodepairs+1)+'-DOF'+str(n_RVEnode_dofs+1)
                print>>Eqns,'*Equation'
                print>>Eqns,'10' # Number of terms in the Constraint
                print>>Eqns,'Ele'+str(n_macro_eles+1)+'-RVE'+str(n_macroele_GPs+1)+'-EdgeNodeT'+str(n_FaceBT_nodepairs+1)+', '+str(n_RVEnode_dofs+1)+', -1.0' # Top edge node DOF, also the first DOF which will be removed
                print>>Eqns,'Ele'+str(n_macro_eles+1)+'-RVE'+str(n_macroele_GPs+1)+'-EdgeNodeB'+str(n_FaceBT_nodepairs+1)+', '+str(n_RVEnode_dofs+1)+', 1.0' # Bottom edge node DOF
                for n_macroele_nodes in range(8): # Loop through all macroscale nodes
                    # Coefficient of the macroscale node term obtained by multiplying RVE dimension and macroscale shape function gradient along z direction
                    print>>Eqns,'Ele'+str(n_macro_eles+1)+'-N'+str(n_macroele_nodes+1)+', '+str(n_RVEnode_dofs+1)+', '+str(dz*N_GloDeriv[n_macroele_nodes][2]) # Macroscale node DOF
                    
        # Call Sets and set up the MPCs for vertices of FaceBa
        print>>Sets,'*Nset, nset=Ele'+str(n_macro_eles+1)+'-RVE'+str(n_macroele_GPs+1)+'-V1, instance=Ele'+str(n_macro_eles+1)+'-RVE'+str(n_macroele_GPs+1) # Create a Set for the node V1
        print>>Sets,str(V1+1) # Node number for V1
        print>>Sets,'*Nset, nset=Ele'+str(n_macro_eles+1)+'-RVE'+str(n_macroele_GPs+1)+'-V2, instance=Ele'+str(n_macro_eles+1)+'-RVE'+str(n_macroele_GPs+1) # Create a Set for the node V2
        print>>Sets,str(V2+1) # Node number for V2
        print>>Sets,'*Nset, nset=Ele'+str(n_macro_eles+1)+'-RVE'+str(n_macroele_GPs+1)+'-V3, instance=Ele'+str(n_macro_eles+1)+'-RVE'+str(n_macroele_GPs+1) # Create a Set for the node V3
        print>>Sets,str(V3+1) # Node number for V3
        print>>Sets,'*Nset, nset=Ele'+str(n_macro_eles+1)+'-RVE'+str(n_macroele_GPs+1)+'-V4, instance=Ele'+str(n_macro_eles+1)+'-RVE'+str(n_macroele_GPs+1) # Create a Set for the node V4
        print>>Sets,str(V4+1) # Node number for V4

        # Nodes V1 and V2
        for n_RVEnode_dofs in range(3): # Loop through all DOFs of the nodes
            print>>Eqns,'** Constraint: Ele'+str(n_macro_eles+1)+'-RVE'+str(n_macroele_GPs+1)+'-V12-DOF'+str(n_RVEnode_dofs+1) # Create an Equation type Constraint for the DOF
            print>>Eqns,'*Equation'
            print>>Eqns,'10' # Number of terms in the Constraint
            print>>Eqns,'Ele'+str(n_macro_eles+1)+'-RVE'+str(n_macroele_GPs+1)+'-V2, '+str(n_RVEnode_dofs+1)+', -1.0' # V2 node DOF, also the first DOF which will be removed
            print>>Eqns,'Ele'+str(n_macro_eles+1)+'-RVE'+str(n_macroele_GPs+1)+'-V1, '+str(n_RVEnode_dofs+1)+', 1.0' # V1 node DOF
            for n_macroele_nodes in range(8): # Loop through all macroscale nodes
                # Coefficient of the macroscale node term obtained by multiplying RVE dimension and macroscale shape function gradient along x direction
                print>>Eqns,'Ele'+str(n_macro_eles+1)+'-N'+str(n_macroele_nodes+1)+', '+str(n_RVEnode_dofs+1)+', '+str(dx*N_GloDeriv[n_macroele_nodes][0]) # Macroscale node DOF

        # Nodes V2 and V3
        for n_RVEnode_dofs in range(3): # Loop through all DOFs of the nodes
            print>>Eqns,'** Constraint: Ele'+str(n_macro_eles+1)+'-RVE'+str(n_macroele_GPs+1)+'-V23-DOF'+str(n_RVEnode_dofs+1) # Create an Equation type Constraint for the DOF
            print>>Eqns,'*Equation'
            print>>Eqns,'10' # Number of terms in the Constraint
            print>>Eqns,'Ele'+str(n_macro_eles+1)+'-RVE'+str(n_macroele_GPs+1)+'-V3, '+str(n_RVEnode_dofs+1)+', -1.0' # V3 node DOF, also the first DOF which will be removed
            print>>Eqns,'Ele'+str(n_macro_eles+1)+'-RVE'+str(n_macroele_GPs+1)+'-V2, '+str(n_RVEnode_dofs+1)+', 1.0' # V2 node DOF
            for n_macroele_nodes in range(8): # Loop through all macroscale nodes
                # Coefficient of the macroscale node term obtained by multiplying RVE dimension and macroscale shape function gradient along z direction
                print>>Eqns,'Ele'+str(n_macro_eles+1)+'-N'+str(n_macroele_nodes+1)+', '+str(n_RVEnode_dofs+1)+', '+str(dz*N_GloDeriv[n_macroele_nodes][2])
                
        # Nodes V1 and V4
        for n_RVEnode_dofs in range(3): # Loop through all DOFs of the nodes
            print>>Eqns,'** Constraint: Ele'+str(n_macro_eles+1)+'-RVE'+str(n_macroele_GPs+1)+'-V14-DOF'+str(n_RVEnode_dofs+1) # Create an Equation type Constraint for the DOF
            print>>Eqns,'*Equation'
            print>>Eqns,'10' # Number of terms in the Constraint
            print>>Eqns,'Ele'+str(n_macro_eles+1)+'-RVE'+str(n_macroele_GPs+1)+'-V4, '+str(n_RVEnode_dofs+1)+', -1.0' # V4 node DOF, also the first DOF which will be removed
            print>>Eqns,'Ele'+str(n_macro_eles+1)+'-RVE'+str(n_macroele_GPs+1)+'-V1, '+str(n_RVEnode_dofs+1)+', 1.0' # V1 node DOF
            for n_macroele_nodes in range(8): # Loop through all macroscale nodes
                # Coefficient of the macroscale node term obtained by multiplying RVE dimension and macroscale shape function gradient along z direction
                print>>Eqns,'Ele'+str(n_macro_eles+1)+'-N'+str(n_macroele_nodes+1)+', '+str(n_RVEnode_dofs+1)+', '+str(dz*N_GloDeriv[n_macroele_nodes][2])
                
        # Set the rigid body constraint for V1
        for n_RVEnode_dofs in range(3): # Loop through all DOFs of node V1
            print>>Eqns,'** Constraint: Ele'+str(n_macro_eles+1)+'-RVE'+str(n_macroele_GPs+1)+'-V1-DOF'+str(n_RVEnode_dofs+1) # Create an Equation type Constraint for the DOF
            print>>Eqns,'*Equation'
            print>>Eqns,'9'
            print>>Eqns,'Ele'+str(n_macro_eles+1)+'-RVE'+str(n_macroele_GPs+1)+'-V1, '+str(n_RVEnode_dofs+1)+', -1.0'
            for m in range(8):
                print>>Eqns,'Ele'+str(n_macro_eles+1)+'-N'+str(m+1)+', '+str(n_RVEnode_dofs+1)+', '+str(Shape_fn[m]-0.5*dx*N_GloDeriv[m][0]-0.5*dy*N_GloDeriv[m][1]-0.5*dz*N_GloDeriv[m][2])
            
        # Calling sets and setting up the MPCs for the front and back faces
        for k in range(len(PairingFacesBaF)):
            print>>Sets,'*Nset, nset=Ele'+str(n_macro_eles+1)+'-RVE'+str(n_macroele_GPs+1)+'-FaceNodeBa'+str(k+1)+', instance=Ele'+str(n_macro_eles+1)+'-RVE'+str(n_macroele_GPs+1)
            print>>Sets,str(PairingFacesBaF[k][0]+1)
            print>>Sets,'*Nset, nset=Ele'+str(n_macro_eles+1)+'-RVE'+str(n_macroele_GPs+1)+'-FaceNodeF'+str(k+1)+', instance=Ele'+str(n_macro_eles+1)+'-RVE'+str(n_macroele_GPs+1)
            print>>Sets,str(PairingFacesBaF[k][1]+1)
            
            for dof in range(3):
                print>>Eqns,'** Constraint: Ele'+str(n_macro_eles+1)+'-RVE'+str(n_macroele_GPs+1)+'-BaF'+str(k+1)+'-DOF'+str(dof+1)
                print>>Eqns,'*Equation'
                print>>Eqns,'10'
                print>>Eqns,'Ele'+str(n_macro_eles+1)+'-RVE'+str(n_macroele_GPs+1)+'-FaceNodeF'+str(k+1)+', '+str(dof+1)+', -1.0'
                print>>Eqns,'Ele'+str(n_macro_eles+1)+'-RVE'+str(n_macroele_GPs+1)+'-FaceNodeBa'+str(k+1)+', '+str(dof+1)+', 1.0'
                for m in range(8):
                    print>>Eqns,'Ele'+str(n_macro_eles+1)+'-N'+str(m+1)+', '+str(dof+1)+', '+str(dy*N_GloDeriv[m][1])
                    
        # Calling sets and setting up the MPCs for the edges parallel to y axis
        for k in range(len(Edge1Nodes)):
            print>>Sets,'*Nset, nset=Ele'+str(n_macro_eles+1)+'-RVE'+str(n_macroele_GPs+1)+'-Edge1Node'+str(k+1)+', instance=Ele'+str(n_macro_eles+1)+'-RVE'+str(n_macroele_GPs+1)
            print>>Sets,str(Edge1Nodes[k]+1)
            print>>Sets,'*Nset, nset=Ele'+str(n_macro_eles+1)+'-RVE'+str(n_macroele_GPs+1)+'-Edge2Node'+str(k+1)+', instance=Ele'+str(n_macro_eles+1)+'-RVE'+str(n_macroele_GPs+1)
            print>>Sets,str(Edge2Nodes[k]+1)
            print>>Sets,'*Nset, nset=Ele'+str(n_macro_eles+1)+'-RVE'+str(n_macroele_GPs+1)+'-Edge3Node'+str(k+1)+', instance=Ele'+str(n_macro_eles+1)+'-RVE'+str(n_macroele_GPs+1)
            print>>Sets,str(Edge3Nodes[k]+1)
            print>>Sets,'*Nset, nset=Ele'+str(n_macro_eles+1)+'-RVE'+str(n_macroele_GPs+1)+'-Edge4Node'+str(k+1)+', instance=Ele'+str(n_macro_eles+1)+'-RVE'+str(n_macroele_GPs+1)
            print>>Sets,str(Edge4Nodes[k]+1)
            
            for dof in range(3):
                print>>Eqns,'** Constraint: Ele'+str(n_macro_eles+1)+'-RVE'+str(n_macroele_GPs+1)+'-Edge12Node'+str(k+1)+'-DOF'+str(dof+1)
                print>>Eqns,'*Equation'
                print>>Eqns,'10'
                print>>Eqns,'Ele'+str(n_macro_eles+1)+'-RVE'+str(n_macroele_GPs+1)+'-Edge2Node'+str(k+1)+', '+str(dof+1)+', -1.0'
                print>>Eqns,'Ele'+str(n_macro_eles+1)+'-RVE'+str(n_macroele_GPs+1)+'-Edge1Node'+str(k+1)+', '+str(dof+1)+', 1.0'
                for m in range(8):
                    print>>Eqns,'Ele'+str(n_macro_eles+1)+'-N'+str(m+1)+', '+str(dof+1)+', '+str(dx*N_GloDeriv[m][0])
                    
            for dof in range(3):
                print>>Eqns,'** Constraint: Ele'+str(n_macro_eles+1)+'-RVE'+str(n_macroele_GPs+1)+'-Edge23Node'+str(k+1)+'-DOF'+str(dof+1)
                print>>Eqns,'*Equation'
                print>>Eqns,'10'
                print>>Eqns,'Ele'+str(n_macro_eles+1)+'-RVE'+str(n_macroele_GPs+1)+'-Edge3Node'+str(k+1)+', '+str(dof+1)+', -1.0'
                print>>Eqns,'Ele'+str(n_macro_eles+1)+'-RVE'+str(n_macroele_GPs+1)+'-Edge2Node'+str(k+1)+', '+str(dof+1)+', 1.0'
                for m in range(8):
                    print>>Eqns,'Ele'+str(n_macro_eles+1)+'-N'+str(m+1)+', '+str(dof+1)+', '+str(dz*N_GloDeriv[m][2])
                    
            for dof in range(3):
                print>>Eqns,'** Constraint: Ele'+str(n_macro_eles+1)+'-RVE'+str(n_macroele_GPs+1)+'-Edge14Node'+str(k+1)+'-DOF'+str(dof+1)
                print>>Eqns,'*Equation'
                print>>Eqns,'10'
                print>>Eqns,'Ele'+str(n_macro_eles+1)+'-RVE'+str(n_macroele_GPs+1)+'-Edge4Node'+str(k+1)+', '+str(dof+1)+', -1.0'
                print>>Eqns,'Ele'+str(n_macro_eles+1)+'-RVE'+str(n_macroele_GPs+1)+'-Edge1Node'+str(k+1)+', '+str(dof+1)+', 1.0'
                for m in range(8):
                    print>>Eqns,'Ele'+str(n_macro_eles+1)+'-N'+str(m+1)+', '+str(dof+1)+', '+str(dz*N_GloDeriv[m][2])
                    
        # Calling sets and setting up the MPCs for the left and right faces
        for k in range(len(PairingFacesLR)):
            print>>Sets,'*Nset, nset=Ele'+str(n_macro_eles+1)+'-RVE'+str(n_macroele_GPs+1)+'-FaceNodeL'+str(k+1)+', instance=Ele'+str(n_macro_eles+1)+'-RVE'+str(n_macroele_GPs+1)
            print>>Sets,str(PairingFacesLR[k][0]+1)
            print>>Sets,'*Nset, nset=Ele'+str(n_macro_eles+1)+'-RVE'+str(n_macroele_GPs+1)+'-FaceNodeR'+str(k+1)+', instance=Ele'+str(n_macro_eles+1)+'-RVE'+str(n_macroele_GPs+1)
            print>>Sets,str(PairingFacesLR[k][1]+1)
            
            for dof in range(3):
                print>>Eqns,'** Constraint: Ele'+str(n_macro_eles+1)+'-RVE'+str(n_macroele_GPs+1)+'-LR'+str(k+1)+'-DOF'+str(dof+1)
                print>>Eqns,'*Equation'
                print>>Eqns,'10'
                print>>Eqns,'Ele'+str(n_macro_eles+1)+'-RVE'+str(n_macroele_GPs+1)+'-FaceNodeR'+str(k+1)+', '+str(dof+1)+', -1.0'
                print>>Eqns,'Ele'+str(n_macro_eles+1)+'-RVE'+str(n_macroele_GPs+1)+'-FaceNodeL'+str(k+1)+', '+str(dof+1)+', 1.0'
                for m in range(8):
                    print>>Eqns,'Ele'+str(n_macro_eles+1)+'-N'+str(m+1)+', '+str(dof+1)+', '+str(dx*N_GloDeriv[m][0])
                    
        # Calling sets and setting up the MPCs for the top and bottom faces
        for k in range(len(PairingFacesBT)):
            print>>Sets,'*Nset, nset=Ele'+str(n_macro_eles+1)+'-RVE'+str(n_macroele_GPs+1)+'-FaceNodeB'+str(k+1)+', instance=Ele'+str(n_macro_eles+1)+'-RVE'+str(n_macroele_GPs+1)
            print>>Sets,str(PairingFacesBT[k][0]+1)
            print>>Sets,'*Nset, nset=Ele'+str(n_macro_eles+1)+'-RVE'+str(n_macroele_GPs+1)+'-FaceNodeT'+str(k+1)+', instance=Ele'+str(n_macro_eles+1)+'-RVE'+str(n_macroele_GPs+1)
            print>>Sets,str(PairingFacesBT[k][1]+1)
            
            for dof in range(3):
                print>>Eqns,'** Constraint: Ele'+str(n_macro_eles+1)+'-RVE'+str(n_macroele_GPs+1)+'-BT'+str(k+1)+'-DOF'+str(dof+1)
                print>>Eqns,'*Equation'
                print>>Eqns,'10'
                print>>Eqns,'Ele'+str(n_macro_eles+1)+'-RVE'+str(n_macroele_GPs+1)+'-FaceNodeT'+str(k+1)+', '+str(dof+1)+', -1.0'
                print>>Eqns,'Ele'+str(n_macro_eles+1)+'-RVE'+str(n_macroele_GPs+1)+'-FaceNodeB'+str(k+1)+', '+str(dof+1)+', 1.0'
                for m in range(8):
                    print>>Eqns,'Ele'+str(n_macro_eles+1)+'-N'+str(m+1)+', '+str(dof+1)+', '+str(dz*N_GloDeriv[m][2])
            
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
Replaced 'remove' function with 'del' function
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

















