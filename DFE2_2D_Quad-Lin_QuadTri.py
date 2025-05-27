# -*- coding: utf-8 -*-
"""
Created on Fri Jul 14 15:07:06 2023

@author: Kirk Ming Yeoh (e0546208@u.nus.edu)

This script sets up a Direct FE2 input file for problems involving linear 2D quadrilateral elements (CPS4/CPS4R/CPE4/CPE4R) at the macroscale, and any 2D continuum elements at the microscale.
"""

### Importing required libraries
import os
import numpy as np
import math

### Obtain the user-defined inputs
# Read the user-defined inputs
execfile('DFE2_0_UserInput.py')

# Defining the macroscale integration points for full integration if not specified
if GP == '':
    GP = [[-3**-0.5,-3**-0.5],[3**-0.5,-3**-0.5],[3**-0.5,3**-0.5],[-3**-0.5,3**-0.5]]

# Defining the Gaussian weights based on the integration points
if len(GP) == 1:
    Weight = 4.0
elif len(GP) == 2:
    Weight = 2.0
elif len(GP) == 4:
    Weight = 1.0

# Defining the tolerance if not specified
if Tol == '':
    Tol = 1e-6

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

# isclose equivalent comparison for floating point numbers
# Some versions of Abaqus Python IDE have older numpy modules which do not have the standard isclose()
def FPisclose(FP1,FP2,tolerance):
    if (abs(FP1-FP2) <= tolerance):
        return 1
    else:
        return 0

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
                inp1.insert(n_inp1_lines+2+n_lines,Temp2) # Insert the current filled line into the input file lines
                n_lines = n_lines+1
                Temp2 = str(n_fulllist_terms) # Start a new line with the next term
                n_terms = 1
            else: # All other terms
                Temp2 = Temp2+', '+str(n_fulllist_terms)
                n_terms = n_terms+1
        inp1.insert(n_inp1_lines+2+n_lines,Temp2) # Insert the final line into the input file lines
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
    if (inp1[n_inp1_lines].count('** Section'))!=0: # Find the line containing the macroscale Section thickness, using the keyword '** Section'
        if inp1[n_inp1_lines+2] == ',': 
            Thickness = 1.0 # Set the thickness to 1.0 if it is not defined
        else:
            Thickness = float(inp1[n_inp1_lines+2].strip(',')) # Extract the thickness value if it is defined
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

RVEMats = open('RVEMats.dat','w') # Create a temporary file to store information on RVE Materials
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
# Sorting the nodal connectivity to match with Direct FE2 conventions
N_macro_eles = len(MacroNodalConnect) # Total number of macroscale elements
NodalConnect = [] # List to store sorted macroscale element nodal connectivities
NodalCoordX = [] # List to store nodal x coordinates based on sorted macroscale element nodal connectivities
NodalCoordY = [] # List to store nodal y coordinates based on sorted macroscale element nodal connectivities
for n_macro_eles in range(N_macro_eles): # Loop through all macroscale elements
    Nodes = [] # List to store node labels for all nodes of the current macroscale element
    X = [] # List to store x coordinates for all nodes of the current macroscale element
    Y = [] # List to store y coordinates for all nodes of the current macroscale element
    Ang = [] # List to store angles relative to the centroid for all nodes of the current macroscale element
    TempConnect = [] # List to store sorted nodal connectivity of the current macroscale element
    TempCoordX = [] # List to store nodal x coordinates based on sorted nodal connectivity of the current macroscale element
    TempCoordY = [] # List to store nodal y coordinates based on sorted nodal connectivity of the current macroscale element
    
    # Obtaining nodal information from this particular element
    for n_macroele_nodes in range(len(MacroNodalConnect[n_macro_eles])): # Loop through all nodes of the current macroscale element 
        Nodes.append(MacroNodalConnect[n_macro_eles][n_macroele_nodes]) # Store the node number 
        X.append(MacroNodalCoord[MacroNodalConnect[n_macro_eles][n_macroele_nodes]][0]) # Store the x coordinate
        Y.append(MacroNodalCoord[MacroNodalConnect[n_macro_eles][n_macroele_nodes]][1]) # Store the y coordinate
    
    # Obtaining centroid of the element
    X0 = sum(X)/float(len(MacroNodalConnect[n_macro_eles])) # x coordinate of centroid
    Y0 = sum(Y)/float(len(MacroNodalConnect[n_macro_eles])) # y coordinate of centroid
    
    # Obtaining orientation of each node relative to centroid
    for n_macroele_nodes in range(len(MacroNodalConnect[n_macro_eles])): # Loop through all nodes of the current macroscale element
        X1 = X[n_macroele_nodes]-X0 # Difference in x coordinate
        Y1 = Y[n_macroele_nodes]-Y0 # Difference in y coordinate
        if FPisclose(X1,0.0,Tol): # If the node of interest sits directly above/below the centroid, e.g., the macroscale element is rotated 45 deg
            if Y1>0: # If the node is above the centroid
                theta = 0.5*(math.pi)
            else: # If the node is below the centroid
                theta = 1.5*(math.pi)
        else:       
            theta = math.atan(Y1/X1) # Calculate the angle of the node wrt the centroid using arctan of the differences in coordinates
            if X1<0:
                theta = theta+math.pi
            if theta<0:
                theta = theta+2*(math.pi) # Convert any negative angles to their positive equivalent
        Ang.append(theta*360/(2*(math.pi))) # Convert the angle to degrees 
    SAng = sorted(Ang) # Sort the angles 

    # Sort the nodes based on the sorted angles
    # The sorted order is based on Direct FE2 conventions
    # The first node has natural coordinates [-1,-1], which is expected to have the 3rd sorted angle
    # The other nodes are found by going about the macroscale element in a counter-clockwise fashion
    Order = [2,3,0,1]
    for n_macroele_nodes in range(len(MacroNodalConnect[n_macro_eles])): # Loop through all nodes of the current macroscale element
        Ind = Ang.index(SAng[Order[n_macroele_nodes]]) # Index the sorted angle
        TempConnect.append(Nodes[Ind]) # Extract the corresponding node label
        TempCoordX.append(X[Ind]) # Extract the corresponding node x coordinate
        TempCoordY.append(Y[Ind]) # Extract the corresponding node y coordinate
        
    NodalConnect.append(TempConnect) # Add the sorted node for the current macroscale element to the main list
    NodalCoordX.append(TempCoordX) # Add the sorted x coordinates for the nodes of the current macroscale element to the main list
    NodalCoordY.append(TempCoordY) # Add the sorted y coordinates for the nodes of the current macroscale element to the main list


### Processing the RVE part information
# Finding the smallest and largest nodal coordinate in both directions
# Assumes a rectangular RVE aligned along the global x and y directions
RVE_ListX = [] # Temporary list to store x coordinates of all RVE nodes
RVE_ListY = [] # Temporary list to store y coordinates of all RVE nodes
for n_RVE_nodes in range(len(RVENodalCoord)): # Loop through all RVE nodes
    RVE_ListX.append(RVENodalCoord[n_RVE_nodes][0]) # Store the x coordinates
    RVE_ListY.append(RVENodalCoord[n_RVE_nodes][1]) # Store the y coordinates
xMin = min(RVE_ListX) # Smallest x coordinate
xMax = max(RVE_ListX) # Largest x coordinate
yMin = min(RVE_ListY) # Smallest y coordinate
yMax = max(RVE_ListY) # Largest x coordinate
del RVE_ListX # Remove the temporary list
del RVE_ListY # Remove the temporary list

# Sorting the RVE boundary nodes
FaceLNodes,FaceRNodes,FaceBNodes,FaceTNodes = [],[],[],[] # List to store the nodes on the RVE boundaries
for n_RVE_nodes in range(len(RVENodalCoord)): # Loop through all RVE nodes
    if FPisclose(RVENodalCoord[n_RVE_nodes][0],xMin,Tol): # If the nodal x coordinate matches the RVE smallest x coordinate
        FaceLNodes.append(n_RVE_nodes) # Store the node the left face list
    if FPisclose(RVENodalCoord[n_RVE_nodes][0],xMax,Tol): # If the nodal x coordinate matches the RVE largest x coordinate
        FaceRNodes.append(n_RVE_nodes) # Store the node the right face list
    if FPisclose(RVENodalCoord[n_RVE_nodes][1],yMin,Tol): # If the nodal y coordinate matches the RVE smallest y coordinate
        FaceBNodes.append(n_RVE_nodes) # Store the node the bottom face list
    if FPisclose(RVENodalCoord[n_RVE_nodes][1],yMax,Tol): # If the nodal y coordinate matches the RVE largest y coordinate
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
    
# Sorting RVE dimensions and offsets
B_RVE = xMax - xMin # RVE dimension along x direction
H_RVE = yMax - yMin # RVE dimension along y direction
OffsetX = (xMax + xMin)/2.0 # RVE centroid x coordinate
OffsetY = (yMax + yMin)/2.0 # RVE centroid y coordinate
Offset = [OffsetX,OffsetY] # RVE centroid coordinates

# Adjusting RVE nodal coordinates to correspond to a part centered at the origin
for n_RVE_node in RVENodalCoord: # Loop through all RVE nodes
    for n_nodal_coord in range(2): # Loop through all coordinates of each node
        n_RVE_node[n_nodal_coord] = n_RVE_node[n_nodal_coord] - Offset[n_nodal_coord] # Subtracting the offset from all RVE nodal coordinates
        

### Generating the RVE placement in the macroscale mesh
RVEParts = open('RVEParts.dat','w') # Create a temporary file to store information on RVE Parts
Insts = open('Insts.dat','w') # Create a temporary file to store information on RVE Instances
SF = [] # List to store required RVE volume scaling factors
for n_macro_eles in range(N_macro_eles): # Loop through all macroscale elements
    # Calculate the mapping function between natural and global coordinates
    # x = a0 + a1*tsi + a2*eta + a3*tsi*eta
    # y = b0 + b1*tsi + b2*eta + b3*tsi*eta
    # C is matrix where [C]*{a0;a1;a2;a3} = {x1;x2;x3;x4} and [C]*{b0;b1;b2;b3} = {y1;y2;y3;y4}, based on the equations above
    # 1,2,3,4 for x and y refer to the 4 macroscale nodes
    C = np.array([[1,-1,-1,1],[1,1,-1,-1],[1,1,1,1],[1,-1,1,-1]]) # Matrix of natural coordinates
    C_inv = np.linalg.inv(C) # Inverse of matrix C
    [a0,a1,a2,a3] = np.dot(C_inv,NodalCoordX[n_macro_eles]) # Coefficient a
    [b0,b1,b2,b3] = np.dot(C_inv,NodalCoordY[n_macro_eles]) # Coefficient b
    
    for n_macroele_GPs in range(len(GP)): # Loop through all integration points of the macroscale element
        [tsi,eta] = GP[n_macroele_GPs] # Natural coordinates of the current integration point
        J = np.array([[a1+a3*eta,b1+b3*eta],[a2+a3*tsi,b2+b3*tsi]]) # Jacobian matrix of the current integration point
        J_RVE = (Weight*abs(Weight*np.linalg.det(J))/(B_RVE*H_RVE)) # Scaling factor for RVE volume (to be applied as Section thickness) at the current integration point
        
        [N1,N2,N3,N4] = Bilin_Interpolation(tsi,eta) # Shape function values corresponding to the current integration point
        RVE_X = N1*NodalCoordX[n_macro_eles][0] + N2*NodalCoordX[n_macro_eles][1] + N3*NodalCoordX[n_macro_eles][2] + N4*NodalCoordX[n_macro_eles][3] # x coordinate of the current integration point
        RVE_Y = N1*NodalCoordY[n_macro_eles][0] + N2*NodalCoordY[n_macro_eles][1] + N3*NodalCoordY[n_macro_eles][2] + N4*NodalCoordY[n_macro_eles][3] # y coordinate of the current integration point

        if round(J_RVE,5) in SF: # Check if this volume scaling factor has been used before, reuse the same RVE Part if so
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
                print>>RVEParts,str(n_RVE_nodes+1)+', '+str(RVENodalCoord[n_RVE_nodes][0])+', '+str(RVENodalCoord[n_RVE_nodes][1])

            # Print the elements of the new Part
            for n_inp2_lines in range(StartEle,Sections[0]): # Loop through the RVE input file lines from the start of elements till the first Section
                print>>RVEParts,inp2[n_inp2_lines]

            # Print the Sections of the new Part
            for n_RVE_sections in range(len(Sections)): # Loop through the RVE Sections
                print>>RVEParts,inp2[Sections[n_RVE_sections]]+'-'+str(Ind_SF+1) # Print the Section name, with the SF index added
                print>>RVEParts,inp2[Sections[n_RVE_sections]+1] # Print the next line of the Section
                if inp2[Sections[n_RVE_sections]+1].count('*Solid Section')!=0: # If the Section is a Solid Section
                    print>>RVEParts,str(Thickness*J_RVE)+',' # Apply the SF to the thickness
                elif inp2[Sections[n_RVE_sections]+1].count('*Cohesive Section')!=0: # If the Section is a Cohesive Section
                    Line = inp2[Sections[n_RVE_sections]+2].split(',')
                    print>>RVEParts,str(Line[0])+','+str(Thickness*J_RVE) # Apply the SF to the thickness

            # Print the end of the new Part
            print>>RVEParts,'*End Part'

        # Print the Instance of the new Part
        print>>Insts,'*Instance, name=Ele'+str(n_macro_eles+1)+'-RVE'+str(n_macroele_GPs+1)+', part=RVE-'+str(Ind_SF+1) # Name of the Instance and the Part it refers to
        print>>Insts,str(RVE_X)+', '+str(RVE_Y)+', 0.' # Translation required for the Instance, from the origin to the macroscale integration point coordinates
        print>>Insts,'*End Instance'
        print>>Insts,'**'
            
RVEParts.close()
Insts.close()

            
### Setting up the MPCs
Sets = open('Sets.dat','w') # Create a temporary file to store information on Direct FE2 Sets
Eqns = open('Eqns.dat','w') # Create a temporary file to store information on Direct FE2 MPCs

# Pairing the nodes
# Assume the RVE mesh is perfectly periodic
# Left and right faces
FaceLNodes = TakeVertexOut(SortListofNodes1D(FaceLNodes,1)) # Sort the left face nodes based on their y coordinates and remove the two corner nodes
FaceRNodes = TakeVertexOut(SortListofNodes1D(FaceRNodes,1)) # Sort the right face nodes based on their y coordinates and remove the two corner nodes
PairingFacesLR = [] # List to store the paired left and right face nodes
for n_FaceL_nodes in range(len(FaceLNodes)): # Loop through the left face nodes
    Temp = [] # Temporary list to store each pair of left and right face nodes
    Temp.append(FaceLNodes[n_FaceL_nodes]) # Left face node number of the pair
    Temp.append(FaceRNodes[n_FaceL_nodes]) # Right face node number of the pair
    PairingFacesLR.append(Temp) # Store the pair of node numbers into the left-right list

# Bottom and top faces
FaceBNodes = TakeVertexOut(SortListofNodes1D(FaceBNodes,0)) # Sort the bottom face nodes based on their x coordinates and remove the two corner nodes
FaceTNodes = TakeVertexOut(SortListofNodes1D(FaceTNodes,0)) # Sort the top face nodes based on their x coordinates and remove the two corner nodes
PairingFacesBT = [] # List to store the paired bottom and top face nodes
for n_FaceB_nodes in range(len(FaceBNodes)): # Loop through the bottom face nodes
    Temp = [] # Temporary list to store each pair of bottom and top face nodes
    Temp.append(FaceBNodes[n_FaceB_nodes]) # Bottom face node number of the pair
    Temp.append(FaceTNodes[n_FaceB_nodes]) # Top face node number of the pair
    PairingFacesBT.append(Temp) # Store the pair of node numbers into the bottom-top list
    
# Calculating the coefficients and setting up the MPCs
# For each macroscale element
for n_macro_eles in range(N_macro_eles): # Loop through all macroscale elements
    # Same mapping function between natural and global coordinates
    C = np.array([[1,-1,-1,1],[1,1,-1,-1],[1,1,1,1],[1,-1,1,-1]]) # Same matrix of natural coordinates based on the mapping function between natural and global coordinates
    C_inv = np.linalg.inv(C) # Inverse of matrix C
    [a0,a1,a2,a3] = np.dot(C_inv,NodalCoordX[n_macro_eles]) # Coefficient a
    [b0,b1,b2,b3] = np.dot(C_inv,NodalCoordY[n_macro_eles]) # Coefficient b
    
    # Call macroscale nodes into Sets
    for n_macroele_nodes in range(4): # Loop through all nodes of the macroscale element
        print>>Sets,'*Nset, nset=Ele'+str(n_macro_eles+1)+'-N'+str(n_macroele_nodes+1)+', instance='+str(MacroInstName) # Create a Set for the macroscale node
        print>>Sets,str(NodalConnect[n_macro_eles][n_macroele_nodes]+1) # Node number of the macroscale node

    # For each macroscale integration point
    for n_macroele_GPs in range(len(GP)): # Loop through all integration points of the macroscale element
        [tsi,eta] = GP[n_macroele_GPs] # Natural coordinates of the current integration point
        J = np.array([[a1+a3*eta,b1+b3*eta],[a2+a3*tsi,b2+b3*tsi]]) # Jacobian matrix of the current integration point
        J_inv = np.linalg.inv(J) # Inverse of the Jacobian matrix of the current integration point

        # Shape function values at the current macroscale integration point
        Shape_fn = Bilin_Interpolation(tsi,eta)

        # Expressions for the shape function gradients
        # Derived by differentiating the shape functions wrt the natural coordinates
        dN1 = [-0.25*(1-eta),-0.25*(1-tsi)] # Gradients of the shape function wrt tsi and eta for macroscale node 1
        dN2 = [0.25*(1-eta),-0.25*(1+tsi)] # Gradients of the shape function wrt tsi and eta for macroscale node 2
        dN3 = [0.25*(1+eta),0.25*(1+tsi)] # Gradients of the shape function wrt tsi and eta for macroscale node 3
        dN4 = [-0.25*(1+eta),0.25*(1-tsi)] # Gradients of the shape function wrt tsi and eta for macroscale node 4
        N_NatDeriv = [dN1,dN2,dN3,dN4] # List of shape function gradients wrt tsi and eta
        N_GloDeriv = [[],[],[],[]] # List to store shape function gradients wrt to x and y
                
        # Calculate shape function gradients along x and y directions
        # Obtained by multiplying the inverse of the Jacobian matrix with the shape function gradients wrt to tsi and eta
        for n_macroele_nodes in range(4): # Loop through all nodes of the macroscale element
            N_GloDeriv[n_macroele_nodes] = np.dot(J_inv,np.transpose(np.array(N_NatDeriv[n_macroele_nodes]))) # Matrix multiplication between the inverse of the Jacobian matrix and shape function gradients wrt tsi and eta
        
        # Call Sets and set up the MPCs for left and right face nodes
        for n_FaceLR_nodepairs in range(len(PairingFacesLR)): # Loop through all left-right face node pairs
            print>>Sets,'*Nset, nset=Ele'+str(n_macro_eles+1)+'-RVE'+str(n_macroele_GPs+1)+'-FaceNodeL'+str(n_FaceLR_nodepairs+1)+', instance=Ele'+str(n_macro_eles+1)+'-RVE'+str(n_macroele_GPs+1) # Create a Set for the left face node of the pair
            print>>Sets,str(PairingFacesLR[n_FaceLR_nodepairs][0]+1) # Left face node number of the pair
            print>>Sets,'*Nset, nset=Ele'+str(n_macro_eles+1)+'-RVE'+str(n_macroele_GPs+1)+'-FaceNodeR'+str(n_FaceLR_nodepairs+1)+', instance=Ele'+str(n_macro_eles+1)+'-RVE'+str(n_macroele_GPs+1) # Create a Set for the right face node of the pair
            print>>Sets,str(PairingFacesLR[n_FaceLR_nodepairs][1]+1) # Right face node number of the pair
            
            for n_RVEnode_dofs in range(2): # Loop through all DOFs of the nodes
                print>>Eqns,'** Constraint: Ele'+str(n_macro_eles+1)+'-RVE'+str(n_macroele_GPs+1)+'-LR'+str(n_FaceLR_nodepairs+1)+'-DOF'+str(n_RVEnode_dofs+1) # Create an Equation type Constraint for the DOF
                print>>Eqns,'*Equation'
                print>>Eqns,'6' # Number of terms in the Constraint
                print>>Eqns,'Ele'+str(n_macro_eles+1)+'-RVE'+str(n_macroele_GPs+1)+'-FaceNodeR'+str(n_FaceLR_nodepairs+1)+', '+str(n_RVEnode_dofs+1)+', -1.0' # Right face node DOF, also the first DOF which will be removed
                print>>Eqns,'Ele'+str(n_macro_eles+1)+'-RVE'+str(n_macroele_GPs+1)+'-FaceNodeL'+str(n_FaceLR_nodepairs+1)+', '+str(n_RVEnode_dofs+1)+', 1.0' # Left face node DOF
                for n_macroele_nodes in range(4): # Loop through all macroscale nodes
                    # Coefficient of the macroscale node term obtained by multiplying RVE dimension and macroscale shape function gradient along x direction
                    print>>Eqns,'Ele'+str(n_macro_eles+1)+'-N'+str(n_macroele_nodes+1)+', '+str(n_RVEnode_dofs+1)+', '+str(B_RVE*N_GloDeriv[n_macroele_nodes][0]) # Macroscale node DOF
                    
        # Call Sets and set up the MPCs for bottom and top face nodes
        for n_FaceBT_nodepairs in range(len(PairingFacesBT)): # Loop through all bottom-top face node pairs
            print>>Sets,'*Nset, nset=Ele'+str(n_macro_eles+1)+'-RVE'+str(n_macroele_GPs+1)+'-FaceNodeB'+str(n_FaceBT_nodepairs+1)+', instance=Ele'+str(n_macro_eles+1)+'-RVE'+str(n_macroele_GPs+1) # Create a Set for the bottom face node of the pair
            print>>Sets,str(PairingFacesBT[n_FaceBT_nodepairs][0]+1) # Bottom face node number of the pair
            print>>Sets,'*Nset, nset=Ele'+str(n_macro_eles+1)+'-RVE'+str(n_macroele_GPs+1)+'-FaceNodeT'+str(n_FaceBT_nodepairs+1)+', instance=Ele'+str(n_macro_eles+1)+'-RVE'+str(n_macroele_GPs+1) # Create a Set for the top face node of the pair
            print>>Sets,str(PairingFacesBT[n_FaceBT_nodepairs][1]+1) # Top face node number of the pair
            
            for n_RVEnode_dofs in range(2): # Loop through all DOFs of the nodes
                print>>Eqns,'** Constraint: Ele'+str(n_macro_eles+1)+'-RVE'+str(n_macroele_GPs+1)+'-BT'+str(n_FaceBT_nodepairs+1)+'-DOF'+str(n_RVEnode_dofs+1) # Create an Equation type Constraint for the DOF
                print>>Eqns,'*Equation'
                print>>Eqns,'6' # Number of terms in the Constraint
                print>>Eqns,'Ele'+str(n_macro_eles+1)+'-RVE'+str(n_macroele_GPs+1)+'-FaceNodeT'+str(n_FaceBT_nodepairs+1)+', '+str(n_RVEnode_dofs+1)+', -1.0' # Top face node DOF, also the first DOF which will be removed
                print>>Eqns,'Ele'+str(n_macro_eles+1)+'-RVE'+str(n_macroele_GPs+1)+'-FaceNodeB'+str(n_FaceBT_nodepairs+1)+', '+str(n_RVEnode_dofs+1)+', 1.0' # Bottom face node DOF
                for n_macroele_nodes in range(4): # Loop through all macroscale nodes
                    # Coefficient of the macroscale node term obtained by multiplying RVE dimension and macroscale shape function gradient along y direction
                    print>>Eqns,'Ele'+str(n_macro_eles+1)+'-N'+str(n_macroele_nodes+1)+', '+str(n_RVEnode_dofs+1)+', '+str(H_RVE*N_GloDeriv[n_macroele_nodes][1]) # Macroscale node DOF
                    
        # Call Sets and set up the MPCs for vertices
        print>>Sets,'*Nset, nset=Ele'+str(n_macro_eles+1)+'-RVE'+str(n_macroele_GPs+1)+'-V1, instance=Ele'+str(n_macro_eles+1)+'-RVE'+str(n_macroele_GPs+1) # Create a Set for the node V1
        print>>Sets,str(V1+1) # Node number for V1
        print>>Sets,'*Nset, nset=Ele'+str(n_macro_eles+1)+'-RVE'+str(n_macroele_GPs+1)+'-V2, instance=Ele'+str(n_macro_eles+1)+'-RVE'+str(n_macroele_GPs+1) # Create a Set for the node V2
        print>>Sets,str(V2+1) # Node number for V2
        print>>Sets,'*Nset, nset=Ele'+str(n_macro_eles+1)+'-RVE'+str(n_macroele_GPs+1)+'-V3, instance=Ele'+str(n_macro_eles+1)+'-RVE'+str(n_macroele_GPs+1) # Create a Set for the node V3
        print>>Sets,str(V3+1) # Node number for V3
        print>>Sets,'*Nset, nset=Ele'+str(n_macro_eles+1)+'-RVE'+str(n_macroele_GPs+1)+'-V4, instance=Ele'+str(n_macro_eles+1)+'-RVE'+str(n_macroele_GPs+1) # Create a Set for the node V4
        print>>Sets,str(V4+1) # Node number for V4

        # Nodes V2 and V3
        for n_RVEnode_dofs in range(2): # Loop through all DOFs of the nodes
            print>>Eqns,'** Constraint: Ele'+str(n_macro_eles+1)+'-RVE'+str(n_macroele_GPs+1)+'-V23-DOF'+str(n_RVEnode_dofs+1) # Create an Equation type Constraint for the DOF
            print>>Eqns,'*Equation'
            print>>Eqns,'6' # Number of terms in the Constraint
            print>>Eqns,'Ele'+str(n_macro_eles+1)+'-RVE'+str(n_macroele_GPs+1)+'-V3, '+str(n_RVEnode_dofs+1)+', -1.0' # V3 node DOF, also the first DOF which will be removed
            print>>Eqns,'Ele'+str(n_macro_eles+1)+'-RVE'+str(n_macroele_GPs+1)+'-V2, '+str(n_RVEnode_dofs+1)+', 1.0' # V2 node DOF
            for n_macroele_nodes in range(4): # Loop through all macroscale nodes
                # Coefficient of the macroscale node term obtained by multiplying RVE dimension and macroscale shape function gradient along y direction
                print>>Eqns,'Ele'+str(n_macro_eles+1)+'-N'+str(n_macroele_nodes+1)+', '+str(n_RVEnode_dofs+1)+', '+str(H_RVE*N_GloDeriv[n_macroele_nodes][1]) # Macroscale node DOF

        # Nodes V1 and V4
        for n_RVEnode_dofs in range(2): # Loop through all DOFs of the nodes
            print>>Eqns,'** Constraint: Ele'+str(n_macro_eles+1)+'-RVE'+str(n_macroele_GPs+1)+'-V14-DOF'+str(n_RVEnode_dofs+1) # Create an Equation type Constraint for the DOF
            print>>Eqns,'*Equation'
            print>>Eqns,'6' # Number of terms in the Constraint
            print>>Eqns,'Ele'+str(n_macro_eles+1)+'-RVE'+str(n_macroele_GPs+1)+'-V4, '+str(n_RVEnode_dofs+1)+', -1.0' # V4 node DOF, also the first DOF which will be removed
            print>>Eqns,'Ele'+str(n_macro_eles+1)+'-RVE'+str(n_macroele_GPs+1)+'-V1, '+str(n_RVEnode_dofs+1)+', 1.0' # V1 node DOF
            for n_macroele_nodes in range(4): # Loop through all macroscale nodes
                # Coefficient of the macroscale node term obtained by multiplying RVE dimension and macroscale shape function gradient along y direction
                print>>Eqns,'Ele'+str(n_macro_eles+1)+'-N'+str(n_macroele_nodes+1)+', '+str(n_RVEnode_dofs+1)+', '+str(H_RVE*N_GloDeriv[n_macroele_nodes][1]) # Macroscale node DOF

        # Nodes V1 and V2
        for n_RVEnode_dofs in range(2): # Loop through all DOFs of the nodes
            print>>Eqns,'** Constraint: Ele'+str(n_macro_eles+1)+'-RVE'+str(n_macroele_GPs+1)+'-V12-DOF'+str(n_RVEnode_dofs+1) # Create an Equation type Constraint for the DOF
            print>>Eqns,'*Equation'
            print>>Eqns,'6' # Number of terms in the Constraint
            print>>Eqns,'Ele'+str(n_macro_eles+1)+'-RVE'+str(n_macroele_GPs+1)+'-V2, '+str(n_RVEnode_dofs+1)+', -1.0' # V2 node DOF, also the first DOF which will be removed
            print>>Eqns,'Ele'+str(n_macro_eles+1)+'-RVE'+str(n_macroele_GPs+1)+'-V1, '+str(n_RVEnode_dofs+1)+', 1.0' # V1 node DOF
            for n_macroele_nodes in range(4): # Loop through all macroscale nodes
                # Coefficient of the macroscale node term obtained by multiplying RVE dimension and macroscale shape function gradient along x direction
                print>>Eqns,'Ele'+str(n_macro_eles+1)+'-N'+str(n_macroele_nodes+1)+', '+str(n_RVEnode_dofs+1)+', '+str(B_RVE*N_GloDeriv[n_macroele_nodes][0]) # Macroscale node DOF

        # Set the rigid body constraint for V1
        for n_RVEnode_dofs in range(2): # Loop through all DOFs of node V1
            print>>Eqns,'** Constraint: Ele'+str(n_macro_eles+1)+'-RVE'+str(n_macroele_GPs+1)+'-V1-DOF'+str(n_RVEnode_dofs+1) # Create an Equation type Constraint for the DOF
            print>>Eqns,'*Equation'
            print>>Eqns,'5' # Number of terms in the Constraint
            print>>Eqns,'Ele'+str(n_macro_eles+1)+'-RVE'+str(n_macroele_GPs+1)+'-V1, '+str(n_RVEnode_dofs+1)+', -1.0' # V1 node DOF, also the first DOF which will be removed
            for n_macroele_nodes in range(4): # Loop through all macroscale nodes
                # Coefficient of the macroscale node term obtained by as shape function values at RVE centroid (or macroscale integration point) subtracted with half of RVE dimensions multiplied with macroscale shape function along x and y directions
                print>>Eqns,'Ele'+str(n_macro_eles+1)+'-N'+str(n_macroele_nodes+1)+', '+str(n_RVEnode_dofs+1)+', '+str(Shape_fn[n_macroele_nodes]-0.5*B_RVE*N_GloDeriv[n_macroele_nodes][0]-0.5*H_RVE*N_GloDeriv[n_macroele_nodes][1])

Sets.close()
Eqns.close()


### Writing the new Direct FE2 input file
# Read the information written previously
# Opening the files containing the information
RVEParts = open('RVEParts.dat','r') # RVE Parts
Insts = open('Insts.dat','r') # RVE Instances
Sets = open('Sets.dat','r') # Sets
Eqns = open('Eqns.dat','r') # MPCs
RVEMats = open('RVEMats.dat','r') # Materials

# Empty lists to store the information
RVEParts_lines = [] # RVE Parts
Insts_lines = [] # RVE Instances
Sets_lines = [] # Sets
Eqns_lines = [] # MPCs
RVEMats_lines = [] # Materials

# Read the lines containing the RVE Parts
while 1:
    line = RVEParts.readline()
    if not line:
        break
    line = line.strip()
    RVEParts_lines.append(line)

# Read the lines containing the RVE Instances
while 1:
    line = Insts.readline()
    if not line:
        break
    line = line.strip()
    Insts_lines.append(line)

# Read the lines containing the Sets
while 1:
    line = Sets.readline()
    if not line:
        break
    line = line.strip()
    Sets_lines.append(line)

# Read the lines containing the MPCs
while 1:
    line = Eqns.readline()
    if not line:
        break
    line = line.strip()
    Eqns_lines.append(line)

# Read the lines containing the Materials
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
        
f3 = open(NewInpName,'w') # Create the final Direct FE2 input file

# Print from the header in the macroscale input file until the end of the macroscale Part
Mark1 = inp1.index('*End Part')
for n_inp1_lines in range(0,Mark1+1):
    print>>f3,inp1[n_inp1_lines]

# Print the information on RVE Parts
for n_RVEParts_lines in range(len(RVEParts_lines)):
    print>>f3,RVEParts_lines[n_RVEParts_lines] 

# Print from the end of the macroscale Part in the macroscale input file until the end of the macroscale Instance
Mark2 = inp1.index('*End Instance')
for n_inp1_lines in range(Mark1+1,Mark2+2):
    print>>f3,inp1[n_inp1_lines]

# Print the information on RVE Instances
for n_RVEInsts_lines in range(len(Insts_lines)):
    print>>f3,Insts_lines[n_RVEInsts_lines]

# Print from the end of the macroscale Instance in the macroscale input file until the end of Assembly, including any macroscale Constraints 
if StartConst == 'a': # If there are no macroscale Constraints
    StartConst = inp1.index('*End Assembly') # Index the end of Assembly in the macroscale input file   
for n_inp1_lines in range(Mark2+2,StartConst):
    print>>f3,inp1[n_inp1_lines]

# Print the information on Sets
for n_Sets_lines in range(len(Sets_lines)):
    print>>f3,Sets_lines[n_Sets_lines]

# Print the information on MPCs
for n_MPCs_lines in range(len(Eqns_lines)):
    print>>f3,Eqns_lines[n_MPCs_lines]                    

# Print from the end of Assembly in the macroscale input file until the end of Materials
Mark3 = inp1.index('** ----------------------------------------------------------------')
for n_inp1_lines in range(StartConst,Mark3):
    print>>f3,inp1[n_inp1_lines]

# Print the information on RVE Materials
for n_RVEMats_lines in range(len(RVEMats_lines)):
    print>>f3,RVEMats_lines[n_RVEMats_lines]

# Print from the end of Materials in the macroscale input file until the end of the file
for n_inp1_lines in range(Mark3,len(inp1)):
    print>>f3,inp1[n_inp1_lines]
    
f3.close()

# Delete the temporary lists and files 
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

Added '+1' to the end of RVE Material definition if the RVE input file has no further information
When reading RVE Material definitions, previous code will miss the last line of the entire RVE input file

End of 240916 Revision

250519 Revision
Revisions for improved clarity:
Replaced one letter, non-descriptive variables with more explanatory variable names
Added additional comments to most lines to explain their functoions
Renamed the file to clarify that it is meant for linear quadrilateral macroscale elements

Replaced all == floating point comparisons wth isclose() functions
More portable and flexible comparison for floating point numbers

Added the cases where isclose(X1,0.0) when sorting the macroscale ndoes
Accounts for the cases where the node of interest is located directly above or below the centroid, where atan cannot be uesd

Revised the user-defined input section to read another Python script
Prevents users from modifying this current script unless necessary

Revised the GP and Gaussian weight section
Allows for more flexibility to use reduced integration if desired

End of 250519 Revision

Proposed Revisions (yet to be implemented)
Generalise RVE Sections portion to account for possible material orientations
Account for RVE sets and surfaces from the microscale input file, found from Part?
Account for RVE level interactions and constraints, found from Instance, including the rearrangement? adopt from code for Yuhao and ZB
Account for BCs applied at Initial step when writing RVE materials, which appears before the ------- line? adopt from code for Yuhao and ZB
Allow for (partial) reduced integration by changing the list of GP? adopt from code for ZB
Account for multiple macroscale parts, multiple types of RVE and different RVE orientations? adopt from code for Yuhao
Relabel J as J^T for better clarity
Merge with quadratic macroscale? Probably not a good idea, better to split up
Merge the RVE and MPC loops?
'''     
        
            
        
            
            


            





























