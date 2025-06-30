# -*- coding: utf-8 -*-
"""
Created on Mon Jun 23 13:51:30 2025

@author: Kirk Ming Yeoh (e0546208@u.nus.edu)

This script contains all the modules and functions required for running the Direct FE2 set up scripts
"""

### Import required libraries
import math
import numpy as np
import os

### Data extraction modules
# Read information from user-defined input files and remove 'generate' keyword for easier processing
def ReadInp(Inp,InpLines):
    f_inp = open(Inp,'r')# Open the macroscale input file as f_inp
    
    while 1: # Read the macroscale input file line by line and store it
        line = f_inp.readline()
        if not line:
            break
        line = line.strip()
        InpLines.append(line)
        
    f_inp.close()
    
    for n_inp_lines in reversed(range(len(InpLines))): # Read through the macroscale input file lines, done in reverse to avoid issues with line number when expanding the 'generate' keyword
        if (InpLines[n_inp_lines].count('generate')!=0): # Lines that compact node or element lists and contain the 'generate' keyword
            Temp = (InpLines[n_inp_lines+1].replace(',',' ')).split() # Split the key numbers in the compacted list containing the start, end and increment for the list
            n_terms = 0 # Term counter
            n_lines = 0 # Extra line counter
            for n_fulllist_terms in range(int(Temp[0]),int(Temp[1])+1,int(Temp[2])):
                if n_terms == 0: # Start of the list
                    Temp2 = str(n_fulllist_terms)
                    n_terms = n_terms+1
                elif n_terms == 16: # 16th term of the list, where a new line is required
                    InpLines.insert(n_inp_lines+2+n_lines,Temp2) # Insert the current filled line into the macroscale input file lines
                    n_lines = n_lines+1
                    Temp2 = str(n_fulllist_terms) # Start a new line with the next term
                    n_terms = 1
                else: # All other terms
                    Temp2 = Temp2+', '+str(n_fulllist_terms)
                    n_terms = n_terms+1
            InpLines.insert(n_inp_lines+2+n_lines,Temp2) # Insert the final line into the macroscale input file lines
            InpLines[n_inp_lines] = InpLines[n_inp_lines][0:len(InpLines[n_inp_lines])-10] # Remove the 'generate' keyword from the opening line of the set or surface list
            del InpLines[n_inp_lines+1] # Remove the original compacted list

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

# Extract key macroscale information
def MacroInfo(inp1,RVE_Dim):
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
        if RVE_Dim==2:
            if (inp1[n_inp1_lines].count('** Section'))!=0: # Find the line containing the macroscale Section thickness, using the keyword '** Section'
                if inp1[n_inp1_lines+2] == ',': 
                    Thickness = 1.0 # Set the thickness to 1.0 if it is not defined
                else:
                    Thickness = float(inp1[n_inp1_lines+2].strip(',')) # Extract the thickness value if it is defined
            
    if RVE_Dim == 2:
        return MacroNodalConnect,MacroNodalCoord,MacroInstName,StartConst,Thickness
    elif RVE_Dim == 3:
        return MacroNodalConnect,MacroNodalCoord,MacroInstName,StartConst
        
# Sort the macroscale nodes
def SortMacroNodes(MacroNodalConnect,MacroNodalCoord,Macro_Dim,Sort):
    NodalConnect,NodalCoordX,NodalCoordY = [],[],[]
    if Macro_Dim == 3:
        NodalCoordZ = []
    
    for n_macro_eles in range(len(MacroNodalConnect)):
        if Sort == '': # Use Abaqus' default nodal connectivity
            TempConnect,TempCoordX,TempCoordY = [],[],[]
            if Macro_Dim == 3:
                TempCoordZ = []
            for n_macroele_nodes in range(len(MacroNodalConnect[n_macro_eles])):
                TempConnect.append(MacroNodalConnect[n_macro_eles][n_macroele_nodes])
                TempCoordX.append(MacroNodalCoord[MacroNodalConnect[n_macro_eles][n_macroele_nodes]][0])
                TempCoordY.append(MacroNodalCoord[MacroNodalConnect[n_macro_eles][n_macroele_nodes]][1])
                if Macro_Dim == 3:
                    TempCoordZ.append(MacroNodalCoord[MacroNodalConnect[n_macro_eles][n_macroele_nodes]][2])
            NodalConnect.append(TempConnect)
            NodalCoordX.append(TempCoordX)
            NodalCoordY.append(TempCoordY)
            if Macro_Dim == 3:
                NodalCoordZ.append(TempCoordZ)
            
#        elif Sort == 'global':
#            sss
#        elif Sort == 'user':
#            xxx
    
    if Macro_Dim == 2:
        return NodalConnect,NodalCoordX,NodalCoordY
    elif Macro_Dim == 3:
        return NodalConnect,NodalCoordX,NodalCoordY,NodalCoordZ
        
# Extract key RVE information
def RVEInfo(inp2,RVE_Dim):
    RVENodalConnect,RVENodalCoord = [],[]
    Search(inp2,'*Element',RVENodalConnect,0) # Search for RVE elements' nodal connectivity and store it
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
    
    RVEMats = open('RVEMats.dat','w') # Create a temporary file to store information on RVE Materials
    for n_inp2_lines in range(len(Materials)): # Loop through each RVE Material
        for n2_inp2_lines in range(Materials[n_inp2_lines]+1,len(inp2)): # Search for the end of each RVE Material definition, starting from the line after the first line
            if (inp2[n2_inp2_lines].count('*Material'))!=0 or (inp2[n2_inp2_lines].count('**'))!=0: # Start of next RVE Material or other definition, marked with the keywords '*Material' or '**' respectively 
                MatEnd = n2_inp2_lines # Mark the end of the current RVE Material
                break
            MatEnd = n2_inp2_lines+1 # Current RVE Material ends with the RVE input file if no further information provided beyond RVE Materials, such as Step
        for n2_inp2_lines in range(Materials[n_inp2_lines],MatEnd):
            RVEMats.write(inp2[n2_inp2_lines]+'\n')
    RVEMats.close()
    
    # Finding the smallest and largest nodal coordinate in all 3 directions
    # Assumes a cuboidal RVE aligned along the global x, y and z directions
    RVE_ListX = [] # Temporary list to store x coordinates of all RVE nodes
    RVE_ListY = [] # Temporary list to store y coordinates of all RVE nodes
    if RVE_Dim == 3:
        RVE_ListZ = [] # Temporary list to store z coordinates of all RVE nodes
    for n_RVE_nodes in range(len(RVENodalCoord)): # Loop through all RVE nodes
        RVE_ListX.append(RVENodalCoord[n_RVE_nodes][0]) # Store the x coordinates
        RVE_ListY.append(RVENodalCoord[n_RVE_nodes][1]) # Store the y coordinates
        if RVE_Dim == 3:
            RVE_ListZ.append(RVENodalCoord[n_RVE_nodes][2]) # Store the z coordinates
    xMin = min(RVE_ListX) # Smallest x coordinate
    xMax = max(RVE_ListX) # Largest x coordinate
    yMin = min(RVE_ListY) # Smallest y coordinate
    yMax = max(RVE_ListY) # Largest y coordinate
    if RVE_Dim == 3:
        zMin = min(RVE_ListZ) # Smallest z coordinate
        zMax = max(RVE_ListZ) # Largest z coordinate
    del RVE_ListX # Remove the temporary list
    del RVE_ListY # Remove the temporary list
    if RVE_Dim == 3:
        del RVE_ListZ # Remove the temporary list
    
    # Sorting RVE dimensions and offsets
    B_RVE = xMax - xMin # RVE dimension along x direction
    H_RVE = yMax - yMin # RVE dimension along y direction
    if RVE_Dim == 3:
        T_RVE = zMax - zMin # RVE dimension along z direction
    OffsetX = (xMax + xMin)/2.0 # RVE centroid x coordinate
    OffsetY = (yMax + yMin)/2.0 # RVE centroid y coordinate
    if RVE_Dim == 3:
        OffsetZ = (zMax + zMin)/2.0 # RVE centroid z coordinate
    if RVE_Dim == 2:
        Offset = [OffsetX,OffsetY] # RVE centroid coordinates
    elif RVE_Dim == 3:
        Offset = [OffsetX,OffsetY,OffsetZ] # RVE centroid coordinates
        
    # Adjusting RVE nodal coordinates to correspond to a part centered at the origin
    for n_RVE_node in RVENodalCoord: # Loop through all RVE nodes
        for n_nodal_coord in range(RVE_Dim): # Loop through all coordinates of each node
            n_RVE_node[n_nodal_coord] = n_RVE_node[n_nodal_coord] - Offset[n_nodal_coord] # Subtracting the offset from all RVE nodal coordinates
    
    if RVE_Dim == 2:
        MaxCoord = [xMin-Offset[0],xMax-Offset[0],yMin-Offset[1],yMax-Offset[1]] # List to store adjusted maximum and minimum coordinates
        return RVENodalConnect,RVENodalCoord,Sections,Materials,StartEle,B_RVE,H_RVE,Offset,MaxCoord
    elif RVE_Dim == 3:
        MaxCoord = [xMin-Offset[0],xMax-Offset[0],yMin-Offset[1],yMax-Offset[1],zMin-Offset[2],zMax-Offset[2]] # List to store adjusted maximum and minimum coordinates
        return RVENodalConnect,RVENodalCoord,Sections,Materials,StartEle,B_RVE,H_RVE,T_RVE,Offset,MaxCoord
        
# isclose equivalent comparison for floating point numbers
# Some versions of Abaqus Python IDE have older numpy modules which do not have the standard isclose()
def FPisclose(FP1,FP2,tolerance):
    if (abs(FP1-FP2) <= tolerance):
        return 1
    else:
        return 0

# Sort RVE boundary nodes                    
def RVENodes(RVENodalCoord,RVE_Dim,MaxCoord,Tol,Offset):
    if RVE_Dim == 2:
        FaceLNodes,FaceRNodes,FaceBNodes,FaceTNodes = [],[],[],[] # Lists to store the nodes on the RVE boundaries
        for n_RVE_nodes in range(len(RVENodalCoord)): # Loop through all RVE nodes
            if FPisclose(RVENodalCoord[n_RVE_nodes][0],MaxCoord[0],Tol): # If the nodal x coordinate matches the RVE smallest x coordinate
                FaceLNodes.append(n_RVE_nodes) # Store the node the left face list
            if FPisclose(RVENodalCoord[n_RVE_nodes][0],MaxCoord[1],Tol): # If the nodal x coordinate matches the RVE largest x coordinate
                FaceRNodes.append(n_RVE_nodes) # Store the node the right face list
            if FPisclose(RVENodalCoord[n_RVE_nodes][1],MaxCoord[2],Tol): # If the nodal y coordinate matches the RVE smallest y coordinate
                FaceBNodes.append(n_RVE_nodes) # Store the node the bottom face list
            if FPisclose(RVENodalCoord[n_RVE_nodes][1],MaxCoord[3],Tol): # If the nodal y coordinate matches the RVE largest y coordinate
                FaceTNodes.append(n_RVE_nodes) # Store the node the front face list
        for n_FaceL_nodes in range(len(FaceLNodes)): # Loop through all nodes in the left face list
            if FaceLNodes[n_FaceL_nodes] in FaceBNodes: # If the node is also in the bottom face list
                V1 = FaceLNodes[n_FaceL_nodes] # Store the node as vertex V1
            if FaceLNodes[n_FaceL_nodes] in FaceTNodes: # If the node is also in the top face list
                V4 = FaceLNodes[n_FaceL_nodes] # Store the node as vertex V4
        for n_FaceR_nodes in range(len(FaceRNodes)):
            if FaceRNodes[n_FaceR_nodes] in FaceBNodes: # If the node is also in the bottom face list
                V2 = FaceRNodes[n_FaceR_nodes] # Store the node as vertex V2
            if FaceRNodes[n_FaceR_nodes] in FaceTNodes: # If the node is also in the top face list
                V3 = FaceRNodes[n_FaceR_nodes] # Store the node as vertex V3
        return FaceLNodes,FaceRNodes,FaceBNodes,FaceTNodes,V1,V2,V3,V4 # Returns the lists of boundary nodes
    elif RVE_Dim == 3:
        FaceLNodes,FaceRNodes,FaceBaNodes,FaceFNodes,FaceBNodes,FaceTNodes = [],[],[],[],[],[] # Lists to store the nodes on the RVE boundary faces
        EdgeLBaNodes,EdgeLFNodes,EdgeLBNodes,EdgeLTNodes = [],[],[],[] # Lists to store nodes on the edges of the left face
        EdgeRBaNodes,EdgeRFNodes,EdgeRBNodes,EdgeRTNodes = [],[],[],[] # List to store nodes on the edges of the right face
        EdgeBaBNodes,EdgeBaTNodes,EdgeFBNodes,EdgeFTNodes = [],[],[],[] # List to store nodes on the remaining edges of the front and back faces
        for n_RVE_nodes in range(len(RVENodalCoord)): # Loop through all RVE nodes
            if FPisclose(RVENodalCoord[n_RVE_nodes][0],MaxCoord[0],Tol): # If the nodal x coordinate matches the RVE smallest x coordinate
                FaceLNodes.append(n_RVE_nodes) # Store the node the left face list
            if FPisclose(RVENodalCoord[n_RVE_nodes][0],MaxCoord[1],Tol):  # If the nodal x coordinate matches the RVE largest x coordinate
                FaceRNodes.append(n_RVE_nodes) # Store the node the right face list
            if FPisclose(RVENodalCoord[n_RVE_nodes][1],MaxCoord[2],Tol): # If the nodal y coordinate matches the RVE smallest y coordinate
                FaceBaNodes.append(n_RVE_nodes) # Store the node the back face list
            if FPisclose(RVENodalCoord[n_RVE_nodes][1],MaxCoord[3],Tol):  # If the nodal y coordinate matches the RVE largest y coordinate
                FaceFNodes.append(n_RVE_nodes) # Store the node the front face list
            if FPisclose(RVENodalCoord[n_RVE_nodes][2],MaxCoord[4],Tol): # If the nodal z coordinate matches the RVE smallest z coordinate
                FaceBNodes.append(n_RVE_nodes) # Store the node the bottom face list
            if FPisclose(RVENodalCoord[n_RVE_nodes][2],MaxCoord[5],Tol):  # If the nodal z coordinate matches the RVE largest z coordinate
                FaceTNodes.append(n_RVE_nodes) # Store the node the top face list
        for n_FaceL_nodes in range(len(FaceLNodes)): # Loop through all nodes in the left face list
            if FaceLNodes[n_FaceL_nodes] in FaceBaNodes: # If the node is also in the back face list
                EdgeLBaNodes.append(FaceLNodes[n_FaceL_nodes]) # Store the node in the left-back edge list
                if FaceLNodes[n_FaceL_nodes] in FaceBNodes: # If the node is also in the bottom face list
                    V1 = FaceLNodes[n_FaceL_nodes] # Store the node as vertex V1
                if FaceLNodes[n_FaceL_nodes] in FaceTNodes: # If the node is also in the top face list
                    V5 = FaceLNodes[n_FaceL_nodes] # Store the node as vertex V5
            if FaceLNodes[n_FaceL_nodes] in FaceFNodes: # If the node is also in the front face list
                EdgeLFNodes.append(FaceLNodes[n_FaceL_nodes]) # Store the node in the left-front edge list
                if FaceLNodes[n_FaceL_nodes] in FaceBNodes: # If the node is also in the bottom face list
                    V4 = FaceLNodes[n_FaceL_nodes] # Store the node as vertex V4
                if FaceLNodes[n_FaceL_nodes] in FaceTNodes: # If the node is also in the top face list
                    V8 = FaceLNodes[n_FaceL_nodes] # Store the node as vertex V8
            if FaceLNodes[n_FaceL_nodes] in FaceBNodes: # If the node is also in the bottom face list
                EdgeLBNodes.append(FaceLNodes[n_FaceL_nodes]) # Store the node in the left-bottom edge list
            if FaceLNodes[n_FaceL_nodes] in FaceTNodes: # If the node is also in the top face list
                EdgeLTNodes.append(FaceLNodes[n_FaceL_nodes]) # Store the node in the left-top edge list
        for n_FaceR_nodes in range(len(FaceRNodes)): # Loop through all nodes in the right face list
            if FaceRNodes[n_FaceR_nodes] in FaceBaNodes: # If the node is also in the back face list
                EdgeRBaNodes.append(FaceRNodes[n_FaceR_nodes]) # Store the node in the right-back edge list
                if FaceRNodes[n_FaceR_nodes] in FaceBNodes: # If the node is also in the bottom face list
                    V2 = FaceRNodes[n_FaceR_nodes] # Store the node as vertex V2
                if FaceRNodes[n_FaceR_nodes] in FaceTNodes: # If the node is also in the top face list
                    V6 = FaceRNodes[n_FaceR_nodes] # Store the node as vertex V6
            if FaceRNodes[n_FaceR_nodes] in FaceFNodes: # If the node is also in the front face list
                EdgeRFNodes.append(FaceRNodes[n_FaceR_nodes]) # Store the node in the right-front edge list
                if FaceRNodes[n_FaceR_nodes] in FaceBNodes: # If the node is also in the bottom face list
                    V3 = FaceRNodes[n_FaceR_nodes] # Store the node as vertex V3
                if FaceRNodes[n_FaceR_nodes] in FaceTNodes: # If the node is also in the top face list
                    V7 = FaceRNodes[n_FaceR_nodes] # Store the node as vertex V7
            if FaceRNodes[n_FaceR_nodes] in FaceBNodes: # If the node is also in the bottom face list
                EdgeRBNodes.append(FaceRNodes[n_FaceR_nodes]) # Store the node in the right-bottom edge list
            if FaceRNodes[n_FaceR_nodes] in FaceTNodes: # If the node is also in the top face list
                EdgeRTNodes.append(FaceRNodes[n_FaceR_nodes]) # Store the node in the right-top edge list
        for n_FaceBa_nodes in range(len(FaceBaNodes)): # Loop through all nodes in the back face list
            if FaceBaNodes[n_FaceBa_nodes] in FaceBNodes: # If the node is also in the bottom face list
                EdgeBaBNodes.append(FaceBaNodes[n_FaceBa_nodes]) # Store the node in the back-bottom edge list
            if FaceBaNodes[n_FaceBa_nodes] in FaceTNodes: # If the node is also in the top face list
                EdgeBaTNodes.append(FaceBaNodes[n_FaceBa_nodes]) # Store the node in the back-top edge list
        for n_FaceF_nodes in range(len(FaceFNodes)): # Loop through all nodes in the front face list
            if FaceFNodes[n_FaceF_nodes] in FaceBNodes: # If the node is also in the bottom face list
                EdgeFBNodes.append(FaceFNodes[n_FaceF_nodes]) # Store the node in the front-bottom edge list
            if FaceFNodes[n_FaceF_nodes] in FaceTNodes: # If the node is also in the top face list
                EdgeFTNodes.append(FaceFNodes[n_FaceF_nodes]) # Store the node in the front-top edge list
        return FaceLNodes,FaceRNodes,FaceBaNodes,FaceFNodes,FaceBNodes,FaceTNodes,EdgeLBaNodes,EdgeLFNodes,EdgeLBNodes,EdgeLTNodes,EdgeRBaNodes,EdgeRFNodes,EdgeRBNodes,EdgeRTNodes,EdgeBaBNodes,EdgeBaTNodes,EdgeFBNodes,EdgeFTNodes,V1,V2,V3,V4,V5,V6,V7,V8 # Returns the lists of boundary nodes
        

### Data pre-processing modules
# Sorts nodes along an edge using their coordinates
def SortListofNodes1D(faceN,coordinate,RVENodalCoord): # Coordinates: 0 for x; 1 for y; 2 for z
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
def SortListofNodes2D(faceN,coord1,coord2,RVENodalCoord): # 0 for x; 1 for y; 2 for z; gives node numbers ready to be called with python (already -1)
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
def ExcludeNodes(faceN,coord1,coord2,coord3,RVENodalCoord,RVE_Dim): # coord1, coord2 and coord3 are the coordinates to exclude, to be given as lists
    newlistN = []
    
    for n_face_nodes in range(len(faceN)):
        if RVENodalCoord[faceN[n_face_nodes]][0] not in coord1: # If the node does not have an x coordinate matching any in the list
            if RVENodalCoord[faceN[n_face_nodes]][1] not in coord2: # If the node does not have a y coordinate matching any in the list
                if RVE_Dim == 3:
                    if RVENodalCoord[faceN[n_face_nodes]][2] not in coord3: # If the node does not have a z coordinate matching any in the list
                        newlistN.append(faceN[n_face_nodes])
                else:
                    newlistN.append(faceN[n_face_nodes])
    
    return newlistN # Returns a list of nodes that do not have coordinates matching the given set 

# Pairs nodes on parallel faces
def PairFaceNodes(Face1Nodes,Face2Nodes,Coord1,Coord2,RVENodalCoord):
    PairingFaces = [] # List to store the paired face nodes
    for n_Face1_nodes in range(len(Face1Nodes)): # Loop through the first face nodes
        Temp = [] # Temporary list to store each pair of face nodes
        Temp.append(Face1Nodes[n_Face1_nodes]) # First face node number of the pair
        C11 = RVENodalCoord[Face1Nodes[n_Face1_nodes]][Coord1] # First coordinate of the first face node of the pair
        C12 = RVENodalCoord[Face1Nodes[n_Face1_nodes]][Coord2] # Second coordinate of the first face node of the pair
        
        # Find the second face node for pairing based on the shortest in-plane distance with the current first face node
        Dist = [] # Temporary list to store in-plane distance between each second face node with the current first face node
        for n_Face2_nodes in range(len(Face2Nodes)): # Loop through the second face nodes
            C21 = RVENodalCoord[Face2Nodes[n_Face2_nodes]][Coord1] # First coordinate of the second face node
            C22 = RVENodalCoord[Face2Nodes[n_Face2_nodes]][Coord2] # Second coordinate of the second face node
            Dist.append(math.sqrt(pow(C21-C11,2)+pow(C22-C12,2))) # In-plane distance between current second face node with the current first face node
            
        Ind_dist = Dist.index(min(Dist)) # Index the second face node with the smallest in-plane distance from the current first face node as the node to be paired
        Temp.append(Face2Nodes[Ind_dist]) # Second face node number of the pair
        del Face2Nodes[Ind_dist] # Remove the paired second face node from the list to reduce pairing time for subsequent first face nodes
        PairingFaces.append(Temp) # Store the pair of node numbers into the list
        
    return PairingFaces # Returns a list of paired nodes
    

### Direct FE2 set up modules
# Calculate the mapping function between natural and global coordinates
def NatGloCoord(Macro_Dim,MacroEle_Type,NodalCoordX,NodalCoordY,NodalCoordZ):
    if Macro_Dim == 2:
        if MacroEle_Type == 'Quad-Lin':
            # x = a0 + a1*tsi + a2*eta + a3*tsi*eta
            # y = b0 + b1*tsi + b2*eta + b3*tsi*eta
            # C is matrix where [C]*{Coeff_a} = {x1;x2;x3;x4}, repeated for b-y, based on the equations above
            # 1,2,3,4 for x and y refer to the 4 macroscale nodes
            C = np.array([
                    [1,-1,-1,1],
                    [1,1,-1,-1],
                    [1,1,1,1],
                    [1,-1,1,-1]]) # Matrix of natural coordinates
        C_inv = np.linalg.inv(C) # Inverse of matrix C
        Coeff_a = np.dot(C_inv,NodalCoordX) # Coefficient a
        Coeff_b = np.dot(C_inv,NodalCoordY) # Coefficient b
        return Coeff_a,Coeff_b
    
    elif Macro_Dim == 3:
        if MacroEle_Type == 'Hex-Lin':
            # x = a0 + a1*tsi + a2*eta + a3*zeta + a4*tsi*eta + a5*tsi*zeta + a6*eta*zeta + a7*tsi*eta*zeta
            # y = b0 + b1*tsi + b2*eta + b3*zeta + b4*tsi*eta + b5*tsi*zeta + b6*eta*zeta + b7*tsi*eta*zeta
            # z = d0 + d1*tsi + d2*eta + d3*zeta + d4*tsi*eta + d5*tsi*zeta + d6*eta*zeta + d7*tsi*eta*zeta
            # C is matrix where [C]*{Coeff_a} = {x1;x2;x3;x4;x5;x6;x7;x8}, repeated for b-y and d-z, based on the equations above
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
        Coeff_a = np.dot(C_inv,NodalCoordX) # Coefficient a
        Coeff_b = np.dot(C_inv,NodalCoordY) # Coefficient b
        Coeff_d = np.dot(C_inv,NodalCoordZ) # Coefficient d
        return Coeff_a,Coeff_b,Coeff_d

# Calculate parameters related to the Jacobian of the current macroscale integration point    
def MacroEleJ(MacroEle_Type,Coeff_a,Coeff_b,Coeff_d,tsi,eta,zeta):
    if MacroEle_Type == 'Quad-Lin':
        [a0,a1,a2,a3] = Coeff_a
        [b0,b1,b2,b3] = Coeff_b
        # Jacobian matrix of the current integration point
        J = np.array([
                [a1+a3*eta,a2+a3*tsi],
                [b1+b3*eta,b2+b3*tsi]])
        
    if MacroEle_Type == 'Hex-Lin':
        [a0,a1,a2,a3,a4,a5,a6,a7] = Coeff_a
        [b0,b1,b2,b3,b4,b5,b6,b7] = Coeff_b
        [d0,d1,d2,d3,d4,d5,d6,d7] = Coeff_d
        # Jacobian matrix of the current integration point
        J = np.array([
                [a1+a4*eta+a5*zeta+a7*eta*zeta,a2+a4*tsi+a6*zeta+a7*tsi*zeta,a3+a5*tsi+a6*eta+a7*tsi*eta],
                [b1+b4*eta+b5*zeta+b7*eta*zeta,b2+b4*tsi+b6*zeta+b7*tsi*zeta,b3+b5*tsi+b6*eta+b7*tsi*eta],
                [d1+d4*eta+d5*zeta+d7*eta*zeta,d2+d4*tsi+d6*zeta+d7*tsi*zeta,d3+d5*tsi+d6*eta+d7*tsi*eta]])
    
    J_T_inv = np.linalg.inv(np.transpose(J)) # Inverse of transpose of Jacobian matrix
    J_det = abs(np.linalg.det(J)) # Determinant of Jacobian matrix, for RVE volume scaling
    return J_T_inv,J_det

# Calculate the shape function values at the current macroscale integration point
def ShapeFn(MacroEle_Type,tsi,eta,zeta):
    if MacroEle_Type == 'Quad-Lin':
        N1 = float(0.25*(1-tsi)*(1-eta))
        N2 = float(0.25*(1+tsi)*(1-eta))
        N3 = float(0.25*(1+tsi)*(1+eta))
        N4 = float(0.25*(1-tsi)*(1+eta))
        return [N1,N2,N3,N4]
            
    elif MacroEle_Type == 'Hex-Lin':
        N1=float(0.125*(1-tsi)*(1-eta)*(1-zeta))
        N2=float(0.125*(1+tsi)*(1-eta)*(1-zeta))
        N3=float(0.125*(1+tsi)*(1+eta)*(1-zeta))
        N4=float(0.125*(1-tsi)*(1+eta)*(1-zeta))
        N5=float(0.125*(1-tsi)*(1-eta)*(1+zeta))
        N6=float(0.125*(1+tsi)*(1-eta)*(1+zeta))
        N7=float(0.125*(1+tsi)*(1+eta)*(1+zeta))
        N8=float(0.125*(1-tsi)*(1+eta)*(1+zeta))   
        return [N1,N2,N3,N4,N5,N6,N7,N8]

# Generate RVE Part and Instance
def RVEPartInst(J_RVE,SF,RVEParts,RVENodalCoord,RVE_Dim,StartEle,Sections,inp2,Thickness,Insts,n_macro_eles,n_macroele_GPs,MacroEle_N,NodalCoordX,NodalCoordY,NodalCoordZ):
    if round(J_RVE,5) in SF: # Check if this volume scaling factor has been used before, reuse the same RVE Part if so 
        Ind_SF = SF.index(round(J_RVE,5)) # Index of the first instance of this volume scaling factor
    else: # Create a new RVE with this current volume scaling factor if it has not been used before
        Ind_SF = len(SF) # Index the new volume scaling factor as the next term in the list
        SF.append(round(J_RVE,5)) # Add the new volume scaling factor to the list
        
        # Write headers of the new Part
        RVEParts.write('**'+'\n')
        RVEParts.write('*Part, name=RVE-'+str(Ind_SF+1)+'\n') # Name the part with the new index
        RVEParts.write('*Node'+'\n')
        
        # Write the nodal coordinates of the new Part
        for n_RVE_nodes in range(len(RVENodalCoord)): # Loop through all RVE nodes
            if RVE_Dim == 2:
                RVEParts.write(str(n_RVE_nodes+1)+', '+str(RVENodalCoord[n_RVE_nodes][0])+', '+str(RVENodalCoord[n_RVE_nodes][1])+'\n')
            elif RVE_Dim == 3:
                RVEParts.write(str(n_RVE_nodes+1)+', '+str(RVENodalCoord[n_RVE_nodes][0]*J_RVE)+', '+str(RVENodalCoord[n_RVE_nodes][1]*J_RVE)+', '+str(RVENodalCoord[n_RVE_nodes][2]*J_RVE)+'\n')
        
        # Write the elements and sets of the new Part
        for n_inp2_lines in range(StartEle,Sections[0]): # Loop through the RVE input file lines from the start of elements till the first Section
            RVEParts.write(inp2[n_inp2_lines]+'\n')
            
        # Write the Sections of the new Part
        for n_RVE_sections in range(len(Sections)):  # Loop through the RVE Sections
            RVEParts.write(inp2[Sections[n_RVE_sections]]+'-'+str(Ind_SF+1)+'\n') # Write the Section name, with the SF index added
            RVEParts.write(inp2[Sections[n_RVE_sections]+1]+'\n') # Write the next line of the Section
            if RVE_Dim == 2:
                if inp2[Sections[n_RVE_sections]+1].count('*Solid Section')!=0: # If the Section is a Solid Section
                    RVEParts.write(str(Thickness*J_RVE)+','+'\n') # Apply the SF to the thickness
                elif inp2[Sections[n_RVE_sections]+1].count('*Cohesive Section')!=0: # If the Section is a Cohesive Section
                    Line = inp2[Sections[n_RVE_sections]+2].split(',')
                    RVEParts.write(str(Line[0])+','+str(Thickness*J_RVE)+'\n') # Apply the SF to the thickness
            elif RVE_Dim == 3:
                RVEParts.write(','+'\n')
        
        # Write the end of the new Part
        RVEParts.write('*End Part'+'\n')
        
    # Write the Instance of the new Part
    Insts.write('*Instance, name=Ele'+str(n_macro_eles+1)+'-RVE'+str(n_macroele_GPs+1)+', part=RVE-'+str(Ind_SF+1)+'\n') # Name of the Instance and the Part it refers to
    # Determining the location of the Instance
    RVE_X = 0.0
    RVE_Y = 0.0
    RVE_Z = 0.0
    for n_macroele_nodes in range(len(MacroEle_N)):
        RVE_X += MacroEle_N[n_macroele_nodes]*NodalCoordX[n_macroele_nodes]
        RVE_Y += MacroEle_N[n_macroele_nodes]*NodalCoordY[n_macroele_nodes]
        if RVE_Dim == 3:
            RVE_Z += MacroEle_N[n_macroele_nodes]*NodalCoordZ[n_macroele_nodes]
    Insts.write(str(RVE_X)+', '+str(RVE_Y)+', '+str(RVE_Z)+'\n')  # Translation required for the Instance, from the origin to the macroscale integration point coordinates
    Insts.write('*End Instance'+'\n')
    Insts.write('**'+'\n')

# Call macroscale nodes into Sets
def MacroNodeSet(MacroNodalConnect,Sets,n_macro_eles,MacroInstName,NodalConnect):
    for n_macroele_nodes in range(len(MacroNodalConnect[n_macro_eles])): # Loop through all nodes of the macroscale element
        Sets.write('*Nset, nset=Ele'+str(n_macro_eles+1)+'-N'+str(n_macroele_nodes+1)+', instance='+str(MacroInstName)+'\n') # Create a Set for the macroscale node
        Sets.write(str(NodalConnect[n_macro_eles][n_macroele_nodes]+1)+'\n') # Node number of the macroscale node 

# Calculate shape function gradients along global coordinates
def ShapeFnDeriv(MacroEle_Type,tsi,eta,zeta,J_T_inv):
    if MacroEle_Type == 'Quad-Lin':
        # Expressions for the shape function gradients
        # Derived by differentiating the shape functions wrt the natural coordinates
        dN1 = [-0.25*(1-eta),-0.25*(1-tsi)] # Gradients of the shape function wrt tsi and eta for macroscale node 1
        dN2 = [0.25*(1-eta),-0.25*(1+tsi)] # Gradients of the shape function wrt tsi and eta for macroscale node 2
        dN3 = [0.25*(1+eta),0.25*(1+tsi)] # Gradients of the shape function wrt tsi and eta for macroscale node 3
        dN4 = [-0.25*(1+eta),0.25*(1-tsi)] # Gradients of the shape function wrt tsi and eta for macroscale node 4
        N_NatDeriv = [dN1,dN2,dN3,dN4] # List of shape function gradients wrt tsi and eta
        N_GloDeriv = [[],[],[],[]] # List to store shape function gradients wrt to x and y           
        
    elif MacroEle_Type == 'Hex-Lin':
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

    # Calculate shape function gradients along the global coordinates            
    # Obtained by multiplying the inverse of the Jacobian matrix with the shape function gradients wrt to the natural coordinates
    for n_macroele_nodes in range(len(N_NatDeriv)): # Loop through the shape function gradients of all nodes of the macroscale element
        N_GloDeriv[n_macroele_nodes] = np.dot(J_T_inv,np.transpose(np.array(N_NatDeriv[n_macroele_nodes]))) # Matrix multiplication between the inverse of the Jacobian matrix and shape function gradients wrt the natural coordinates
    return N_GloDeriv                    

# Call Sets and set up MPCs
def RVENodeSetMPC(RVE_Dim,J_RVE,Rot,Ori_Dist_Vect,NodeGroupList,Sets,n_macro_eles,n_macroele_GPs,NodeGroup_Label,MPC_Type,Eqns,N_GloDeriv,MacroEle_N):
    # Remove any volume scaling for 2D cases
    if RVE_Dim == 2:
        J_RVE = 1.0     
    
    # Call Sets and set up the MPCs for the node pair
    for n_nodegroups in range(len(NodeGroupList)):
        # Call Sets for the master node of the pair
        Sets.write('*Nset, nset=Ele'+str(n_macro_eles+1)+'-RVE'+str(n_macroele_GPs+1)+'-'+str(NodeGroup_Label[0])+str(NodeGroup_Label[1])+'Node'+str(n_nodegroups+1)+', instance=Ele'+str(n_macro_eles+1)+'-RVE'+str(n_macroele_GPs+1)+'\n')
        Sets.write(str(NodeGroupList[n_nodegroups][0]+1)+'\n')
        
        # Call Sets for the slave node of each pair and set up the MPC for the node pair
        for n_nodepair_terms in range(len(NodeGroupList[n_nodegroups])-1):
            Sets.write('*Nset, nset=Ele'+str(n_macro_eles+1)+'-RVE'+str(n_macroele_GPs+1)+'-'+str(NodeGroup_Label[0])+str(NodeGroup_Label[n_nodepair_terms+2])+'Node'+str(n_nodegroups+1)+', instance=Ele'+str(n_macro_eles+1)+'-RVE'+str(n_macroele_GPs+1)+'\n')
            Sets.write(str(NodeGroupList[n_nodegroups][n_nodepair_terms+1]+1)+'\n')
            
            Dist_Vect = NodePairDist(Rot,Ori_Dist_Vect[n_nodepair_terms],J_RVE)
            NodePair_Label = [NodeGroup_Label[0],NodeGroup_Label[1],NodeGroup_Label[n_nodepair_terms+2]]
            MPC(MPC_Type,RVE_Dim,Eqns,n_macro_eles,n_macroele_GPs,NodePair_Label,n_nodegroups,Dist_Vect,N_GloDeriv,MacroEle_N)
            
        # Set up MPC for rigid body motion, applied on corner V1
        if NodeGroup_Label[0] == 'V':
            Dist_Vect = NodePairDist(Rot,Ori_Dist_Vect[-1],J_RVE)
            NodePair_Label = [NodeGroup_Label[0],NodeGroup_Label[1]]
            MPC(MPC_Type+1,RVE_Dim,Eqns,n_macro_eles,n_macroele_GPs,NodePair_Label,n_nodegroups,Dist_Vect,N_GloDeriv,MacroEle_N)

# Calculate the rotated (if any) and scaled (if any) distances for the node pair
def NodePairDist(Rot,Ori_Dist_Vect,J_RVE):
    if Rot == '':
        Dist_Vect = np.array(Ori_Dist_Vect)*J_RVE # Scaled distance vector for the node pair
    else:
        Dist_Vect = np.dot(Rot,np.array(Ori_Dist_Vect))*J_RVE # Rotated and scaled distance vector for the node pair  
    return Dist_Vect    

# Call and write Sets for RVE nodes pairs
def RVENodeSet(NodeGroupList,Sets,n_macro_eles,n_macroele_GPs,NodeGroup_Label):
    for n_nodegroups in range(len(NodeGroupList)):
        # Call Sets for the master node of the group
        Sets.write('*Nset, nset=Ele'+str(n_macro_eles+1)+'-RVE'+str(n_macroele_GPs+1)+'-'+str(NodeGroup_Label[0])+str(NodeGroup_Label[1])+'Node'+str(n_nodegroups+1)+', instance=Ele'+str(n_macro_eles+1)+'-RVE'+str(n_macroele_GPs+1)+'\n')
        Sets.write(str(NodeGroupList[n_nodegroups][0]+1)+'\n')
        
        # Call Sets for the slave node of each pair
        for n_nodepair_terms in range(len(NodeGroupList[n_nodegroups])-1):
            Sets.write('*Nset, nset=Ele'+str(n_macro_eles+1)+'-RVE'+str(n_macroele_GPs+1)+'-'+str(NodeGroup_Label[0])+str(NodeGroup_Label[n_nodepair_terms+2])+'Node'+str(n_nodegroups+1)+', instance=Ele'+str(n_macro_eles+1)+'-RVE'+str(n_macroele_GPs+1)+'\n')
            Sets.write(str(NodeGroupList[n_nodegroups][n_nodepair_terms+1]+1)+'\n')

# Calculate the distance between node pairs and write MPCs
def MPC(NodeGroupList,Rot,Ori_Dist_Vect,J_RVE,NodeGroup_Label,MPC_Type,Eqns,RVE_Dim,n_macro_eles,n_macroele_GPs,N_GloDeriv,MacroEle_N):
    # Remove any volume scaling for 2D cases
    if RVE_Dim == 2:
        J_RVE = 1.0
        
    for n_nodegroups in range(len(NodeGroupList)):
        # Set up PBCs for each node pair in the group
        for n_nodepair_terms in range(len(NodeGroupList[n_nodegroups])-1):
            Dist_Vect = NodePairDist(Rot,Ori_Dist_Vect[n_nodepair_terms],J_RVE)
            NodePair_Label = [NodeGroup_Label[0],NodeGroup_Label[1],NodeGroup_Label[n_nodepair_terms+2]]
            
            # Solid-solid, displacement only
            if MPC_Type == 1: 
                for n_dofs in range(RVE_Dim):
                    Eqns.write('** Constraint: Ele'+str(n_macro_eles+1)+'-RVE'+str(n_macroele_GPs+1)+'-'+str(NodePair_Label[0])+str(NodePair_Label[1])+str(NodePair_Label[2])+str(n_nodegroups+1)+'-DOF'+str(n_dofs+1)+'\n') # Create an Equation type Constraint for the DOF
                    Eqns.write('*Equation'+'\n')
                    Eqns.write(str(2+len(N_GloDeriv))+'\n') # Addresses both linear and quadratic macroscale nodes in 2D and 3D
                    Eqns.write('Ele'+str(n_macro_eles+1)+'-RVE'+str(n_macroele_GPs+1)+'-'+str(NodePair_Label[0])+str(NodePair_Label[2])+'Node'+str(n_nodegroups+1)+', '+str(n_dofs+1)+', -1.0'+'\n')
                    Eqns.write('Ele'+str(n_macro_eles+1)+'-RVE'+str(n_macroele_GPs+1)+'-'+str(NodePair_Label[0])+str(NodePair_Label[1])+'Node'+str(n_nodegroups+1)+', '+str(n_dofs+1)+', 1.0'+'\n')
                    for n_macroele_nodes in range(len(N_GloDeriv)):
                        Coeff = 0.0
                        for n_coeff_terms in range(RVE_Dim):
                            Coeff += Dist_Vect[n_coeff_terms]*N_GloDeriv[n_macroele_nodes][n_coeff_terms]
                        Eqns.write('Ele'+str(n_macro_eles+1)+'-N'+str(n_macroele_nodes+1)+', '+str(n_dofs+1)+', '+str(Coeff)+'\n')
                        
        # Set up MPC for rigid body motion, applied on corner V1 
        if NodeGroup_Label[0] == 'V':
            Dist_Vect = NodePairDist(Rot,Ori_Dist_Vect[-1],J_RVE)
            NodePair_Label = [NodeGroup_Label[0],NodeGroup_Label[1]]
            
            # Solid-solid, displacement only
            if MPC_Type == 1: 
                for n_dofs in range(RVE_Dim):
                    Eqns.write('** Constraint: Ele'+str(n_macro_eles+1)+'-RVE'+str(n_macroele_GPs+1)+'-'+str(NodePair_Label[0])+str(NodePair_Label[1])+str(n_nodegroups+1)+'-RigidBody'+'-DOF'+str(n_dofs+1)+'\n') # Create an Equation type Constraint for the DOF
                    Eqns.write('*Equation'+'\n')
                    Eqns.write(str(1+len(N_GloDeriv))+'\n') # Addresses both linear and quadratic macroscale nodes in 2D and 3D
                    Eqns.write('Ele'+str(n_macro_eles+1)+'-RVE'+str(n_macroele_GPs+1)+'-'+str(NodePair_Label[0])+str(NodePair_Label[1])+'Node'+str(n_nodegroups+1)+', '+str(n_dofs+1)+', -1.0'+'\n')
                    for n_macroele_nodes in range(len(N_GloDeriv)):
                        Coeff = MacroEle_N[n_macroele_nodes]
                        for n_coeff_terms in range(RVE_Dim):
                            Coeff += Dist_Vect[n_coeff_terms]*N_GloDeriv[n_macroele_nodes][n_coeff_terms]
                        Eqns.write('Ele'+str(n_macro_eles+1)+'-N'+str(n_macroele_nodes+1)+', '+str(n_dofs+1)+', '+str(Coeff)+'\n')


### New input file modules
# Read lines from temporary data files and writes them to the new Direct FE2 input file
def DattoInp(Dat,Inp):
    f_dat = open(Dat,'r')
    
    # Read all lines in the .dat file and write it into the new Direct FE2 input file
    while 1:
        line = f_dat.readline()
        if not line:
            break
        Inp.write(line)
        
    f_dat.close()

# Write the new Direct FE2 input file    
def WriteDFE2(NewInpName,inp1,StartConst):
    f_inp = open(NewInpName,'w')
    
    # Write from the header in the macroscale input file until the end of the macroscale Part
    Mark1 = inp1.index('*End Part')
    for n_inp1_lines in range(0,Mark1+1):
        f_inp.write(inp1[n_inp1_lines]+'\n')
    
    # Write the information on the RVE Parts    
    DattoInp('RVEParts.dat',f_inp)
    
    # Write from the end of the macroscale Part in the macroscale input file until the end of the macroscale Instance
    Mark2 = inp1.index('*End Instance')
    for n_inp1_lines in range(Mark1+1,Mark2+2):
        f_inp.write(inp1[n_inp1_lines]+'\n')
    
    # Write the information on RVE Instances    
    DattoInp('Insts.dat',f_inp)
    
    # Write from the end of the macroscale Instance in the macroscale input file until the end of Assembly, including any macroscale Constraints 
    if StartConst == 'a': # If there are no macroscale Constraints
        StartConst = inp1.index('*End Assembly') # Index the end of Assembly in the macroscale input file      
    for n_inp1_lines in range(Mark2+2,StartConst):
        f_inp.write(inp1[n_inp1_lines]+'\n')
    
    # Write the information on Sets  
    DattoInp('Sets.dat',f_inp)
    
    # Write the information on MPCS
    DattoInp('Eqns.dat',f_inp)
    
    # Print from the end of Assembly in the macroscale input file until the end of Materials
    Mark3 = inp1.index('** ----------------------------------------------------------------')
    for n_inp1_lines in range(StartConst,Mark3):
        f_inp.write(inp1[n_inp1_lines]+'\n')
        
    # Write the information on RVE Materials
    DattoInp('RVEMats.dat',f_inp)
    
    # Print from the end of Materials in the macroscale input file until the end of the file
    for n_inp1_lines in range(Mark3,len(inp1)):
        f_inp.write(inp1[n_inp1_lines]+'\n')
        
    f_inp.close()
    
# Delete the temporary files
def DelTempFiles():
    os.remove('RVEParts.dat')
    os.remove('Insts.dat')
    os.remove('Sets.dat')
    os.remove('Eqns.dat')
    os.remove('RVEMats.dat')
            

'''
Revision log

250630 Original release
'''

