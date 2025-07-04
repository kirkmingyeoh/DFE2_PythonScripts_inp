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
ODBFile = 'DFE2_3D.odb' # Name of the ODB file
odb = openOdb(path=ODBFile) # Open the ODB file

Inst = odb.rootAssembly.instances['MACRO-1'] # Call the macroscale Instance as Inst

# Extract y reaction force (RF3) of the loaded edge
RF3 = 0.0 # Variable to save RF3
NodeList_RF3 = [91,92,93,94,95,96,97,98,99] # Nodes on the loaded edge
for n_nodelist_nodes in range(len(NodeList_RF3)): # Loop through the nodes on the loaded edge
    RF3 += odb.steps['Step-1'].frames[-1].fieldOutputs['RF'].getSubset(region=Inst.nodes[NodeList_RF3[n_nodelist_nodes]-1]).values[0].data[2] # Sum up RF2 of the nodes on the loaded edge

# Extract the x displacement (U1) of selected macroscale nodes
U1 = [] # List to store U1 of selected macroscale nodes
NodeList_U1 = [28,36,51,56,61] # Selected macroscale nodes
for n_nodelist_nodes in range(len(NodeList_U1)): # Loop through selected macroscale nodes
    U1.append(odb.steps['Step-1'].frames[-1].fieldOutputs['U'].getSubset(region=Inst.nodes[NodeList_U1[n_nodelist_nodes]-1]).values[0].data[0]) # Store U1 of the selected macroscale nodes

# Extract the y displacement (U2) of selected macroscale nodes
U2 = [] # List to store U2 of selected macroscale nodes
NodeList_U2 = [15,20,52,56,80] # Selected macroscale nodes
for n_nodelist_nodes in range(len(NodeList_U2)): # Loop through selected macroscale nodes
    U2.append(odb.steps['Step-1'].frames[-1].fieldOutputs['U'].getSubset(region=Inst.nodes[NodeList_U2[n_nodelist_nodes]-1]).values[0].data[1]) # Store U2 of the selected macroscale nodes

# Extract the y displacement (U2) of selected macroscale nodes
U3 = [] # List to store U2 of selected macroscale nodes
NodeList_U3 = [1,30,47,64,85] # Selected macroscale nodes 12
for n_nodelist_nodes in range(len(NodeList_U3)): # Loop through selected macroscale nodes
    U3.append(odb.steps['Step-1'].frames[-1].fieldOutputs['U'].getSubset(region=Inst.nodes[NodeList_U3[n_nodelist_nodes]-1]).values[0].data[2]) # Store U3 of the selected macroscale nodes

odb.close() # Close the ODB file


### Checking the extracted rsults
# Checking RF2 of the loaded edge
if FPisclose(RF3,2098.14198303,Tol): # If the extracted RF2 is within tolerance from the reference value
    print('Check on RF2 of the loaded edge passed.')
else:
    print('Check on RF2 of the loaded edge failed.')

# Checking U1 of the selected macroscale nodes
U1_Check = 0 # Marker to indicate if U1 of all nodes are close to the reference value
if FPisclose(U1[0],-0.36230746,Tol): # If the extracted U1 of the first node is within tolerance from the reference value
    if FPisclose(U1[1],0.36230746,Tol): # If the extracted U1 of the second node is within tolerance from the reference value
        if FPisclose(U1[2],0.00060640799,Tol): # If the extracted U1 of the third node is within tolerance from the reference value
            if FPisclose(U1[3],-0.59760743,Tol): # If the extracted U1 of the fourth node is within tolerance from the reference value
                if FPisclose(U1[4],0.59657943,Tol): # If the extracted U1 of the fifth node is within tolerance from the reference value
                    U1_Check = 1 # Update the marker
                    print('Check on U1 of selected macroscale nodes passed.')
if U1_Check == 0: # If U1 of at least one node is far from the reference value
    print('Check on U1 of selected macroscale nodes failed.')
    
# Checking U2 of selected macroscale nodes
U2_Check = 0 # Marker to indicate if U2 of all nodes are close to the reference value
if FPisclose(U2[0],-6.3585787e-05,Tol): # If the extracted U2 of the first node is within tolerance from the reference value
    if FPisclose(U2[1],-0.00010871163,Tol): # If the extracted U2 of the second node is within tolerance from the reference value
        if FPisclose(U2[2],0.010938519,Tol): # If the extracted U2 of the third node is within tolerance from the reference value
            if FPisclose(U2[3],-0.00094869977,Tol): # If the extracted U2 of the fourth node is within tolerance from the reference value
                if FPisclose(U2[4],-0.0015525675,Tol): # If the extracted U2 of the fifth node is within tolerance from the reference value
                    U2_Check = 1 # Update the marker
                    print('Check on U2 of selected macroscale nodes passed.')
if U2_Check == 0: # If U2 of at least one node is far from the reference value
    print('Check on U2 of selected macroscale nodes failed.')    

# Checking U3 of selected macroscale nodes
U3_Check = 0 # Marker to indicate if U3 of all nodes are close to the reference value
if FPisclose(U3[0],0.18848495,Tol): # If the extracted U3 of the first node is within tolerance from the reference value
    if FPisclose(U3[1],1.3033993,Tol): # If the extracted U3 of the second node is within tolerance from the reference value
        if FPisclose(U3[2],3.2233429,Tol): # If the extracted U3 of the third node is within tolerance from the reference value
            if FPisclose(U3[3],5.7022977,Tol): # If the extracted U3 of the fourth node is within tolerance from the reference value
                if FPisclose(U3[4],8.5283012,Tol): # If the extracted U3 of the fifth node is within tolerance from the reference value
                    U3_Check = 1 # Update the marker
                    print('Check on U3 of selected macroscale nodes passed.')
if U3_Check == 0: # If U2 of at least one node is far from the reference value
    print('Check on U3 of selected macroscale nodes failed.')         

