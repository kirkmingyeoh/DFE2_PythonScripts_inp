# -*- coding: utf-8 -*-
"""
Created on Thu Jul  3 11:36:42 2025

@author: Kirk Ming Yeoh (e0546208@u.nus.edu)

This script checks the selected outputs from the ODB against reference values for the Demo 2D case
"""

from abaqus import *
from abaqusConstants import *
from odbAccess import *

# Defining an isclose equivalent comparison for floating point numbers
# Some versions of Abaqus Python IDE have older numpy modules which do not have the standard isclose()
def FPisclose(FP1,FP2,tolerance):
    if (abs(FP1-FP2) <= tolerance): # If the absolute difference between the two floating point numbers are smaller than the tolerance
        return 1 # Return true
    else:
        return 0 # Return false

Tol = 1e-5 # Tolerance for floating point comparison

### Extracting results from the ODB
ODBFile = 'DFE2_2D.odb' # Name of the ODB file
odb = openOdb(path=ODBFile) # Open the ODB file

Inst = odb.rootAssembly.instances['MACRO-1'] # Call the macroscale Instance as Inst

# Extract y reaction force (RF2) of the loaded edge
RF2 = 0.0 # Variable to save RF2
NodeList_RF2 = [1,22,43] # Nodes on the loaded edge
for n_nodelist_nodes in range(len(NodeList_RF2)): # Loop through the nodes on the loaded edge
    RF2 += odb.steps['Step-1'].frames[-1].fieldOutputs['RF'].getSubset(region=Inst.nodes[NodeList_RF2[n_nodelist_nodes]-1]).values[0].data[1] # Sum up RF2 of the nodes on the loaded edge

# Extract the x displacement (U1) of selected macroscale nodes
U1 = [] # List to store U1 of selected macroscale nodes
NodeList_U1 = [6,29,56] # Selected macroscale nodes
for n_nodelist_nodes in range(len(NodeList_U1)): # Loop through selected macroscale nodes
    U1.append(odb.steps['Step-1'].frames[-1].fieldOutputs['U'].getSubset(region=Inst.nodes[NodeList_U1[n_nodelist_nodes]-1]).values[0].data[0]) # Store U1 of the selected macroscale nodes

# Extract the y displacement (U2) of selected macroscale nodes
U2 = [] # List to store U2 of selected macroscale nodes
NodeList_U2 = [14,33,51] # Selected macroscale nodes
for n_nodelist_nodes in range(len(NodeList_U2)): # Loop through selected macroscale nodes
    U2.append(odb.steps['Step-1'].frames[-1].fieldOutputs['U'].getSubset(region=Inst.nodes[NodeList_U2[n_nodelist_nodes]-1]).values[0].data[1]) # Store U2 of the selected macroscale nodes

odb.close() # Close the ODB file


### Checking the extracted rsults
# Checking RF2 of the loaded edge
if FPisclose(RF2,13.4743955135,Tol): # If the extracted RF2 is within tolerance from the reference value
    print('Check on RF2 of the loaded edge passed.')
else:
    print('Check on RF2 of the loaded edge failed.')

# Checking U1 of the selected macroscale nodes
U1_Check = 0 # Marker to indicate if U1 of all nodes are close to the reference value
if FPisclose(U1[0],-0.69768679,Tol): # If the extracted U1 of the first node is within tolerance from the reference value
    if FPisclose(U1[1],-2.7427836e-06,Tol): # If the extracted U1 of the second node is within tolerance from the reference value
        if FPisclose(U1[2],0.43040556,Tol): # If the extracted U1 of the third node is within tolerance from the reference value
            U1_Check = 1 # Update the marker
            print('Check on U1 of selected macroscale nodes passed.')
if U1_Check == 0: # If U1 of at least one node is far from the reference value
    print('Check on U1 of selected macroscale nodes failed.')
    
# Checking U2 of selected macroscale nodes
U2_Check = 0 # Marker to indicate if U2 of all nodes are close to the reference value
if FPisclose(U2[0],1.6428696,Tol): # If the extracted U2 of the first node is within tolerance from the reference value
    if FPisclose(U2[1],2.5954218,Tol): # If the extracted U2 of the second node is within tolerance from the reference value
        if FPisclose(U2[2],4.3374419,Tol): # If the extracted U2 of the third node is within tolerance from the reference value
            U2_Check = 1 # Update the marker
            print('Check on U2 of selected macroscale nodes passed.')
if U2_Check == 0: # If U2 of at least one node is far from the reference value
    print('Check on U2 of selected macroscale nodes failed.')            
