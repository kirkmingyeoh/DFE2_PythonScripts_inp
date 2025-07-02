# -*- coding: utf-8 -*-
"""
Created on Fri Jul 14 15:07:06 2023

@author: Kirk Ming Yeoh (e0546208@u.nus.edu)

This script sets up a Direct FE2 input file for problems involving linear 2D quadrilateral elements (CPS4/CPS4R/CPE4/CPE4R) at the macroscale, and any 2D continuum elements at the microscale.
"""

### Import the user-defined inputs
#from Test import DFE2_0_UserInput as Input
from DFE2_0_UserInput import MacroInpName,RVEInpName,NewInpName # Input file names
from DFE2_0_UserInput import GP,Tol,Sort # Optional/additional input information

# Define the macroscale integration points for full integration if not specified
if GP == '':
    GP = [[-3**-0.5,-3**-0.5],[3**-0.5,-3**-0.5],[3**-0.5,3**-0.5],[-3**-0.5,3**-0.5]]

# Define the Gaussian weights based on the integration points
if len(GP) == 1:
    Weight = 4.0
elif len(GP) == 2:
    Weight = 2.0
elif len(GP) == 4:
    Weight = 1.0

# Define the tolerance if not specified
if Tol == '':
    Tol = 1e-6

# Define script specific parameters
Macro_Dim = 2
MacroEle_Type = 'Quad-Lin'
RVE_Dim = 2


### Read and save the Macro and RVE input file lines
from DFE2_0_Modules import ReadInp
inp1,inp2 = [],[]
ReadInp(MacroInpName,inp1)
ReadInp(RVEInpName,inp2)    

# Extracting macroscale element info from old inp file
from DFE2_0_Modules import MacroInfo
MacroNodalConnect,MacroNodalCoord,MacroInstName,StartConst,Thickness = MacroInfo(inp1,RVE_Dim)

# Sort macroscale nodes
from DFE2_0_Modules import SortMacroNodes
NodalConnect,NodalCoordX,NodalCoordY = SortMacroNodes(MacroNodalConnect,MacroNodalCoord,Macro_Dim,Sort)

# Extract key RVE information
from DFE2_0_Modules import RVEInfo
RVENodalConnect,RVENodalCoord,Sections,Materials,StartEle,B_RVE,H_RVE,Offset,MaxCoord = RVEInfo(inp2,RVE_Dim)

# Sort RVE boundary nodes
from DFE2_0_Modules import RVENodes
FaceLNodes,FaceRNodes,FaceBNodes,FaceTNodes,V1,V2,V3,V4 = RVENodes(RVENodalCoord,RVE_Dim,MaxCoord,Tol,Offset)   


### Pairing the RVE boundary nodes
# Assumes the RVE mesh is perfectly periodic
from DFE2_0_Modules import SortListofNodes1D,ExcludeNodes

# Left and right faces
FaceLNodes = SortListofNodes1D(ExcludeNodes(FaceLNodes,[],[RVENodalCoord[V1][1],RVENodalCoord[V4][1]],[],RVENodalCoord,RVE_Dim),1,RVENodalCoord) # Sort the left face nodes based on their y coordinates and remove the two corner nodes
FaceRNodes = SortListofNodes1D(ExcludeNodes(FaceRNodes,[],[RVENodalCoord[V1][1],RVENodalCoord[V4][1]],[],RVENodalCoord,RVE_Dim),1,RVENodalCoord) # Sort the right face nodes based on their y coordinates and remove the two corner nodes
PairingFacesLR = [] # List to store the paired left and right face nodes
for n_FaceL_nodes in range(len(FaceLNodes)): # Loop through the left face nodes
    Temp = [] # Temporary list to store each pair of left and right face nodes
    Temp.append(FaceLNodes[n_FaceL_nodes]) # Left face node number of the pair
    Temp.append(FaceRNodes[n_FaceL_nodes]) # Right face node number of the pair
    PairingFacesLR.append(Temp) # Store the pair of node numbers into the left-right list

# Bottom and top faces
FaceBNodes = SortListofNodes1D(ExcludeNodes(FaceBNodes,[RVENodalCoord[V1][0],RVENodalCoord[V2][0]],[],[],RVENodalCoord,RVE_Dim),0,RVENodalCoord) # Sort the bottom face nodes based on their x coordinates and remove the two corner nodes
FaceTNodes = SortListofNodes1D(ExcludeNodes(FaceTNodes,[RVENodalCoord[V1][0],RVENodalCoord[V2][0]],[],[],RVENodalCoord,RVE_Dim),0,RVENodalCoord) # Sort the top face nodes based on their x coordinates and remove the two corner nodes
PairingFacesBT = [] # List to store the paired bottom and top face nodes
for n_FaceB_nodes in range(len(FaceBNodes)): # Loop through the bottom face nodes
    Temp = [] # Temporary list to store each pair of bottom and top face nodes
    Temp.append(FaceBNodes[n_FaceB_nodes]) # Bottom face node number of the pair
    Temp.append(FaceTNodes[n_FaceB_nodes]) # Top face node number of the pair
    PairingFacesBT.append(Temp) # Store the pair of node numbers into the bottom-top list


### Set up the Direct FE2 model
# Open temporary files to store information
RVEParts = open('RVEParts.dat','w') # Temporary file to store information on RVE Parts
Insts = open('Insts.dat','w') # Temporary file to store information on RVE Instances
Sets = open('Sets.dat','w') # Temporary file to store information on Direct FE2 Sets
Eqns = open('Eqns.dat','w') # Temporary file to store information on Direct FE2 MPCs

# Import required functions
from DFE2_0_Modules import NatGloCoord,MacroEleJ,ShapeFn,RVEPartInst,MacroNodeSet,ShapeFnDeriv,RVENodeSet,MPC

SF = [] # List to store required RVE volume scaling factors
for n_macro_eles in range(len(MacroNodalConnect)): # Loop through all macroscale elements
    # Calculate the mapping function between natural and global coordinates
    Coeff_a,Coeff_b = NatGloCoord(Macro_Dim,MacroEle_Type,NodalCoordX[n_macro_eles],NodalCoordY[n_macro_eles],'')
    
    # Call macroscale nodes into Sets
    MacroNodeSet(MacroNodalConnect,Sets,n_macro_eles,MacroInstName,NodalConnect)
    
    for n_macroele_GPs in range(len(GP)): # Loop through all integration points of the macroscale element
        # Natural coordinates and shape function values of the current integration point
        [tsi,eta] = GP[n_macroele_GPs]
        MacroEle_N = ShapeFn(MacroEle_Type,tsi,eta,'')
        
        # Calculate parameters related to the Jacobian of the current macroscale integration point 
        J_T_inv,J_det = MacroEleJ(MacroEle_Type,Coeff_a,Coeff_b,'',tsi,eta,'')
        J_RVE = (Weight*J_det/(B_RVE*H_RVE)) # Scaling factor for RVE volume (to be applied as Section thickness) at the current integration point
        
        # Generate RVE Part and Instance
        RVEPartInst(J_RVE,SF,RVEParts,RVENodalCoord,RVE_Dim,StartEle,Sections,inp2,Thickness,Insts,n_macro_eles,n_macroele_GPs,MacroEle_N,NodalCoordX[n_macro_eles],NodalCoordY[n_macro_eles],'')
        
        # Calculate shape function gradients along global coordinates
        N_GloDeriv = ShapeFnDeriv(MacroEle_Type,tsi,eta,'',J_T_inv)
        
        # Call Sets and set up MPCs for left and right faces
        RVENodeSet(PairingFacesLR,Sets,n_macro_eles,n_macroele_GPs,['Face','L','R'])
        MPC(PairingFacesLR,'',[[B_RVE,0.0]],J_RVE,['Face','L','R'],1,Eqns,RVE_Dim,n_macro_eles,n_macroele_GPs,N_GloDeriv,MacroEle_N)
        
        # Call Sets and set up MPCs for bottom and top faces
        RVENodeSet(PairingFacesBT,Sets,n_macro_eles,n_macroele_GPs,['Face','B','T'])
        MPC(PairingFacesLR,'',[[0.0,H_RVE]],J_RVE,['Face','B','T'],1,Eqns,RVE_Dim,n_macro_eles,n_macroele_GPs,N_GloDeriv,MacroEle_N)
        
        # Call Sets and set up MPCs for vertices
        NodeGroupList_V1V2V3V4 = [[V1,V2,V3,V4]]
        NodeGroup_Label_V1V2V3V4 = ['V','1','2','3','4']
        Ori_Dist_Vects_V1V2V3V4 = [[B_RVE,0.0],[B_RVE,H_RVE],[0.0,H_RVE]]
        RVENodeSet(NodeGroupList_V1V2V3V4,Sets,n_macro_eles,n_macroele_GPs,NodeGroup_Label_V1V2V3V4)
        MPC(NodeGroupList_V1V2V3V4,'',Ori_Dist_Vects_V1V2V3V4,J_RVE,NodeGroup_Label_V1V2V3V4,1,Eqns,RVE_Dim,n_macro_eles,n_macroele_GPs,N_GloDeriv,MacroEle_N)

RVEParts.close()
Insts.close()
Sets.close()
Eqns.close()


### Write the new Direct FE2 input file
from DFE2_0_Modules import WriteDFE2,DelTempFiles
WriteDFE2(NewInpName,inp1,StartConst)
DelTempFiles()


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

250630 Revision
** Major revision - modularisation
All major functions are now modularised and shared across scripts, accounting for variations
All modules are carried in another script 'DFE2_0_Modules.py'

Minor revisions
Relabel J as J^T for better clarity
Merged the RVE and MPC loops
Replaced print>> function with write() for cross-compatibility across different versions of Python

End of 250630 Revision

End of Revision


Proposed Revisions (yet to be implemented)
Generalise RVE Sections portion to account for possible material orientations
Account for RVE sets and surfaces from the microscale input file, found from Part?
Account for RVE level interactions and constraints, found from Instance, including the rearrangement? adopt from code for Yuhao and ZB
Account for BCs applied at Initial step when writing RVE materials, which appears before the ------- line? adopt from code for Yuhao and ZB
Account for multiple macroscale parts, multiple types of RVE and different RVE orientations? adopt from code for Yuhao
Merge with quadratic macroscale? Probably not a good idea, better to split up
Account for different elements within the same macroscale part
'''     
        
            
        
            
            


            





























