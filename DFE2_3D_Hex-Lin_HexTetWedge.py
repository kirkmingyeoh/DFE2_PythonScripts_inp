# -*- coding: utf-8 -*-
"""
Created on Sun Jul 16 14:40:46 2023

@author: Kirk Ming Yeoh (e0546208@u.nus.edu)

This script sets up a Direct FE2 input file for problems involving linear 3D hexahedral elements (C3D8/C3D8R) at the macroscale, and any 3D continuum elements at the microscale.
"""

### Import the user-defined inputs
from DFE2_0_UserInput import MacroInpName,RVEInpName,NewInpName # Input file names
from DFE2_0_UserInput import GP,Tol,Sort # Optional/additional input information

# Define the macroscale integration points for full integration if not specified
if GP == '':
    GP = [[-3**-0.5,-3**-0.5,-3**-0.5],[3**-0.5,-3**-0.5,-3**-0.5],[3**-0.5,3**-0.5,-3**-0.5],[-3**-0.5,3**-0.5,-3**-0.5],[-3**-0.5,-3**-0.5,3**-0.5],[3**-0.5,-3**-0.5,3**-0.5],[3**-0.5,3**-0.5,3**-0.5],[-3**-0.5,3**-0.5,3**-0.5]]

# Define the Gaussian weights based on the integration points
if len(GP) == 1:
    Weight = 8.0
elif len(GP) == 2:
    Weight = 4.0
elif len(GP) == 4:
    Weight = 2.0
elif len(GP) == 8:
    Weight = 1.0

# Define the tolerance if not specified
if Tol == '':
    Tol = 1e-6
    
# Define script specific parameters
Macro_Dim = 3
MacroEle_Type = 'Hex-Lin'
RVE_Dim = 3
 
   
### Read and process information from the Macro and RVE input files
# Read and save the Macro and RVE input file lines
from DFE2_0_Modules import ReadInp
inp1,inp2 = [],[]
ReadInp(MacroInpName,inp1)
ReadInp(RVEInpName,inp2)    
    
# Extract key macroscale information
from DFE2_0_Modules import MacroInfo
MacroNodalConnect,MacroNodalCoord,MacroInstName,StartConst = MacroInfo(inp1,RVE_Dim)

# Sort macroscale nodes
from DFE2_0_Modules import SortMacroNodes
NodalConnect,NodalCoordX,NodalCoordY,NodalCoordZ = SortMacroNodes(MacroNodalConnect,MacroNodalCoord,Macro_Dim,Sort)

# Extract key RVE information
from DFE2_0_Modules import RVEInfo
RVENodalConnect,RVENodalCoord,Sections,Materials,StartEle,B_RVE,H_RVE,T_RVE,Offset,MaxCoord = RVEInfo(inp2,RVE_Dim)

# Sort RVE boundary nodes
from DFE2_0_Modules import RVENodes
FaceLNodes,FaceRNodes,FaceBaNodes,FaceFNodes,FaceBNodes,FaceTNodes,EdgeLBaNodes,EdgeLFNodes,EdgeLBNodes,EdgeLTNodes,EdgeRBaNodes,EdgeRFNodes,EdgeRBNodes,EdgeRTNodes,EdgeBaBNodes,EdgeBaTNodes,EdgeFBNodes,EdgeFTNodes,V1,V2,V3,V4,V5,V6,V7,V8 = RVENodes(RVENodalCoord,RVE_Dim,MaxCoord,Tol,Offset)


### Pairing the RVE boundary nodes
# Assumes the RVE mesh is perfectly periodic
from DFE2_0_Modules import SortListofNodes1D,SortListofNodes2D,ExcludeNodes,PairFaceNodes

# Edges parallel to the z axis
EdgeLBaNodes = SortListofNodes1D(ExcludeNodes(EdgeLBaNodes,[],[],[RVENodalCoord[V1][2],RVENodalCoord[V5][2]],RVENodalCoord,RVE_Dim),2,RVENodalCoord)
EdgeLFNodes = SortListofNodes1D(ExcludeNodes(EdgeLFNodes,[],[],[RVENodalCoord[V1][2],RVENodalCoord[V5][2]],RVENodalCoord,RVE_Dim),2,RVENodalCoord)
EdgeRBaNodes = SortListofNodes1D(ExcludeNodes(EdgeRBaNodes,[],[],[RVENodalCoord[V1][2],RVENodalCoord[V5][2]],RVENodalCoord,RVE_Dim),2,RVENodalCoord)
EdgeRFNodes = SortListofNodes1D(ExcludeNodes(EdgeRFNodes,[],[],[RVENodalCoord[V1][2],RVENodalCoord[V5][2]],RVENodalCoord,RVE_Dim),2,RVENodalCoord)
PairingEdges_LBaLFRBaRF = []
for n_edgeZ_nodes in range(len(EdgeLBaNodes)):
    PairingEdges_LBaLFRBaRF.append([EdgeLBaNodes[n_edgeZ_nodes],EdgeLFNodes[n_edgeZ_nodes],EdgeRBaNodes[n_edgeZ_nodes],EdgeRFNodes[n_edgeZ_nodes]])

# Edges parallel to the y axis
EdgeLBNodes = SortListofNodes1D(ExcludeNodes(EdgeLBNodes,[],[RVENodalCoord[V1][1],RVENodalCoord[V4][1]],[],RVENodalCoord,RVE_Dim),1,RVENodalCoord)
EdgeLTNodes = SortListofNodes1D(ExcludeNodes(EdgeLTNodes,[],[RVENodalCoord[V1][1],RVENodalCoord[V4][1]],[],RVENodalCoord,RVE_Dim),1,RVENodalCoord)
EdgeRBNodes = SortListofNodes1D(ExcludeNodes(EdgeRBNodes,[],[RVENodalCoord[V1][1],RVENodalCoord[V4][1]],[],RVENodalCoord,RVE_Dim),1,RVENodalCoord)
EdgeRTNodes = SortListofNodes1D(ExcludeNodes(EdgeRTNodes,[],[RVENodalCoord[V1][1],RVENodalCoord[V4][1]],[],RVENodalCoord,RVE_Dim),1,RVENodalCoord)
PairingEdges_LBLTRBRT = []
for n_edgeY_nodes in range(len(EdgeLBNodes)):
    PairingEdges_LBLTRBRT.append([EdgeLBNodes[n_edgeY_nodes],EdgeLTNodes[n_edgeY_nodes],EdgeRBNodes[n_edgeY_nodes],EdgeRTNodes[n_edgeY_nodes]])

# Edges parallel to the x axis
EdgeBaBNodes = SortListofNodes1D(ExcludeNodes(EdgeBaBNodes,[RVENodalCoord[V1][0],RVENodalCoord[V2][0]],[],[],RVENodalCoord,RVE_Dim),0,RVENodalCoord)
EdgeBaTNodes = SortListofNodes1D(ExcludeNodes(EdgeBaTNodes,[RVENodalCoord[V1][0],RVENodalCoord[V2][0]],[],[],RVENodalCoord,RVE_Dim),0,RVENodalCoord)
EdgeFBNodes = SortListofNodes1D(ExcludeNodes(EdgeFBNodes,[RVENodalCoord[V1][0],RVENodalCoord[V2][0]],[],[],RVENodalCoord,RVE_Dim),0,RVENodalCoord)
EdgeFTNodes = SortListofNodes1D(ExcludeNodes(EdgeFTNodes,[RVENodalCoord[V1][0],RVENodalCoord[V2][0]],[],[],RVENodalCoord,RVE_Dim),0,RVENodalCoord)
PairingEdges_BaBBaTFBFT = []
for n_edgeX_nodes in range(len(EdgeBaBNodes)):
    PairingEdges_BaBBaTFBFT.append([EdgeBaBNodes[n_edgeX_nodes],EdgeBaTNodes[n_edgeX_nodes],EdgeFBNodes[n_edgeX_nodes],EdgeFTNodes[n_edgeX_nodes]])

# Left and right faces
FaceLNodes = SortListofNodes2D(ExcludeNodes(FaceLNodes,[],[RVENodalCoord[V3][1],RVENodalCoord[V6][1]],[RVENodalCoord[V3][2],RVENodalCoord[V6][2]],RVENodalCoord,RVE_Dim),1,2,RVENodalCoord) # Sort the left face nodes based on their y and z coordinates, excluding nodes on the edges using the y and z coordinates of V3 and V6
FaceRNodes = SortListofNodes2D(ExcludeNodes(FaceRNodes,[],[RVENodalCoord[V3][1],RVENodalCoord[V6][1]],[RVENodalCoord[V3][2],RVENodalCoord[V6][2]],RVENodalCoord,RVE_Dim),1,2,RVENodalCoord) # Sort the right face nodes based on their y and z coordinates, excluding nodes on the edges using the y and z coordinates of V3 and V6
PairingFacesLR = PairFaceNodes(FaceLNodes,FaceRNodes,1,2,RVENodalCoord)

# Bottom and top faces
FaceBNodes = SortListofNodes2D(ExcludeNodes(FaceBNodes,[RVENodalCoord[V3][0],RVENodalCoord[V8][0]],[RVENodalCoord[V3][1],RVENodalCoord[V8][1]],[],RVENodalCoord,RVE_Dim),0,1,RVENodalCoord) # Sort the bottom face nodes based on their x and y coordinates, excluding nodes on the edges using the x and y coordinates of V3 and V8
FaceTNodes = SortListofNodes2D(ExcludeNodes(FaceTNodes,[RVENodalCoord[V3][0],RVENodalCoord[V8][0]],[RVENodalCoord[V3][1],RVENodalCoord[V8][1]],[],RVENodalCoord,RVE_Dim),0,1,RVENodalCoord) # Sort the top face nodes based on their x and y coordinates, excluding nodes on the edges using the x and y coordinates of V3 and V8
PairingFacesBT = PairFaceNodes(FaceBNodes,FaceTNodes,0,1,RVENodalCoord)

# Front and back faces
FaceBaNodes = SortListofNodes2D(ExcludeNodes(FaceBaNodes,[RVENodalCoord[V1][0],RVENodalCoord[V6][0]],[],[RVENodalCoord[V1][2],RVENodalCoord[V6][2]],RVENodalCoord,RVE_Dim),0,2,RVENodalCoord) # Sort the back face nodes based on their x and z coordinates, excluding nodes on the edges using the x and z coordinates of V1 and V6
FaceFNodes = SortListofNodes2D(ExcludeNodes(FaceFNodes,[RVENodalCoord[V1][0],RVENodalCoord[V6][0]],[],[RVENodalCoord[V1][2],RVENodalCoord[V6][2]],RVENodalCoord,RVE_Dim),0,2,RVENodalCoord) # Sort the front face nodes based on their x and z coordinates, excluding nodes on the edges using the x and z coordinates of V1 and V6
PairingFacesBaF = PairFaceNodes(FaceBaNodes,FaceFNodes,0,2,RVENodalCoord)


### Set up the Direct FE2 model
# Open temporary files to store information
RVEParts = open('RVEParts.dat','w') # Temporary file to store information on RVE Parts
Insts = open('Insts.dat','w') # Temporary file to store information on RVE Instances
Sets = open('Sets.dat','w') # Temporary file to store information on Direct FE2 Sets
Eqns = open('Eqns.dat','w') # Temporary file to store information on Direct FE2 MPCs

from DFE2_0_Modules import NatGloCoord,MacroEleJ,ShapeFn,RVEPartInst,MacroNodeSet,ShapeFnDeriv,RVENodeSet,MPC

SF = [] # List to store required RVE volume scaling factors
for n_macro_eles in range(len(MacroNodalConnect)): # Loop through all macroscale elements
    # Calculate the mapping function between natural and global coordinates
    Coeff_a,Coeff_b,Coeff_d = NatGloCoord(Macro_Dim,MacroEle_Type,NodalCoordX[n_macro_eles],NodalCoordY[n_macro_eles],NodalCoordZ[n_macro_eles])
    
    # Call macroscale nodes into Sets
    MacroNodeSet(MacroNodalConnect,Sets,n_macro_eles,MacroInstName,NodalConnect)
    
    for n_macroele_GPs in range(len(GP)): # Loop through all integration points of the macroscale element
        # Natural coordinates and shape function values of the current integration point
        [tsi,eta,zeta] = GP[n_macroele_GPs] 
        MacroEle_N = ShapeFn(MacroEle_Type,tsi,eta,zeta) 
        
        # Calculate parameters related to the Jacobian of the current macroscale integration point 
        J_T_inv,J_det = MacroEleJ(MacroEle_Type,Coeff_a,Coeff_b,Coeff_d,tsi,eta,zeta)
        J_RVE = (Weight*J_det/(B_RVE*H_RVE*T_RVE))**(1.0/3.0) # Scaling factor for RVE volume at the current integration point
        
        # Generate RVE Part and Instance
        RVEPartInst(J_RVE,SF,RVEParts,RVENodalCoord,RVE_Dim,StartEle,Sections,inp2,'',Insts,n_macro_eles,n_macroele_GPs,MacroEle_N,NodalCoordX[n_macro_eles],NodalCoordY[n_macro_eles],NodalCoordZ[n_macro_eles])
        
        # Calculate shape function gradients along global coordinates
        N_GloDeriv = ShapeFnDeriv(MacroEle_Type,tsi,eta,zeta,J_T_inv)
        
        # Call Sets and set up MPCs for left and right faces
        RVENodeSet(PairingFacesLR,Sets,n_macro_eles,n_macroele_GPs,['Face','L','R'])
        MPC(PairingFacesLR,'',[[B_RVE,0.0,0.0]],J_RVE,['Face','L','R'],1,Eqns,RVE_Dim,n_macro_eles,n_macroele_GPs,N_GloDeriv,MacroEle_N)
        
        # Call Sets and set up MPCs for back and front faces
        RVENodeSet(PairingFacesBaF,Sets,n_macro_eles,n_macroele_GPs,['Face','Ba','F'])
        MPC(PairingFacesBaF,'',[[0.0,H_RVE,0.0]],J_RVE,['Face','Ba','F'],1,Eqns,RVE_Dim,n_macro_eles,n_macroele_GPs,N_GloDeriv,MacroEle_N)
        
        # Call Sets and set up MPCs for left and right faces
        RVENodeSet(PairingFacesBT,Sets,n_macro_eles,n_macroele_GPs,['Face','B','T'])
        MPC(PairingFacesBT,'',[[0.0,0.0,T_RVE]],J_RVE,['Face','B','T'],1,Eqns,RVE_Dim,n_macro_eles,n_macroele_GPs,N_GloDeriv,MacroEle_N)
        
        # Call Sets and set up MPCs for edges parallel to the z axis
        NodeGroup_Label_LBaLFRBaRF = ['Edge','LBa','LF','RBa','RF']
        Ori_Dist_Vects_LBaLFRBaRF = [[0.0,H_RVE,0.0],[B_RVE,0.0,0.0],[B_RVE,H_RVE,0.0]]
        RVENodeSet(PairingEdges_LBaLFRBaRF,Sets,n_macro_eles,n_macroele_GPs,NodeGroup_Label_LBaLFRBaRF)
        MPC(PairingEdges_LBaLFRBaRF,'',Ori_Dist_Vects_LBaLFRBaRF,J_RVE,NodeGroup_Label_LBaLFRBaRF,1,Eqns,RVE_Dim,n_macro_eles,n_macroele_GPs,N_GloDeriv,MacroEle_N)
        
        # Call Sets and set up MPCs for edges parallel to the y axis
        NodeGroup_Label_LBLTRBRT = ['Edge','LB','LT','RB','RT']
        Ori_Dist_Vects_LBLTRBRT = [[0.0,0.0,T_RVE],[B_RVE,0.0,0.0],[B_RVE,0.0,T_RVE]]
        RVENodeSet(PairingEdges_LBLTRBRT,Sets,n_macro_eles,n_macroele_GPs,NodeGroup_Label_LBLTRBRT)
        MPC(PairingEdges_LBLTRBRT,'',Ori_Dist_Vects_LBLTRBRT,J_RVE,NodeGroup_Label_LBLTRBRT,1,Eqns,RVE_Dim,n_macro_eles,n_macroele_GPs,N_GloDeriv,MacroEle_N)
        
        # Call Sets and set up MPCs for edges parallel to the x axis
        NodeGroup_Label_BaBBaTFBFT = ['Edge','BaB','BaT','FB','FT']
        Ori_Dist_Vects_BaBBaTFBFT = [[0.0,0.0,T_RVE],[0.0,H_RVE,0.0],[0.0,H_RVE,T_RVE]]
        RVENodeSet(PairingEdges_BaBBaTFBFT,Sets,n_macro_eles,n_macroele_GPs,NodeGroup_Label_BaBBaTFBFT)
        MPC(PairingEdges_BaBBaTFBFT,'',Ori_Dist_Vects_BaBBaTFBFT,J_RVE,NodeGroup_Label_BaBBaTFBFT,1,Eqns,RVE_Dim,n_macro_eles,n_macroele_GPs,N_GloDeriv,MacroEle_N)
        
        # Call Sets and set up MPCs for the vertices
        NodeGroupList_V1V2V3V4V5V6V7V8 = [[V1,V2,V3,V4,V5,V6,V7,V8]]
        NodeGroup_Label_V1V2V3V4V5V6V7V8 = ['V','1','2','3','4','5','6','7','8']
        Ori_Dist_Vects_V1V2V3V4V5V6V7V8 = [[B_RVE,0.0,0.0,0.0],[B_RVE,H_RVE,0.0],[0.0,H_RVE,0.0],[0.0,0.0,T_RVE],[B_RVE,0.0,T_RVE],[B_RVE,H_RVE,T_RVE],[0.0,H_RVE,T_RVE],[-0.5*B_RVE,-0.5*H_RVE,-0.5*T_RVE]]
        RVENodeSet(NodeGroupList_V1V2V3V4V5V6V7V8,Sets,n_macro_eles,n_macroele_GPs,NodeGroup_Label_V1V2V3V4V5V6V7V8)
        MPC(NodeGroupList_V1V2V3V4V5V6V7V8,'',Ori_Dist_Vects_V1V2V3V4V5V6V7V8,J_RVE,NodeGroup_Label_V1V2V3V4V5V6V7V8,1,Eqns,RVE_Dim,n_macro_eles,n_macroele_GPs,N_GloDeriv,MacroEle_N)
        
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

230714 Original release

240916 Revision
Replaced 'remove' function with 'del' function
'remove' function searches and deletes the first match, while 'del' function deletes the specific line as intended 

250519 Revision
Revisions for improved clarity:
Replaced one letter, non-descriptive variables with more explanatory variable names
Added additional comments to most lines to explain their functoions
Renamed the file to clarify that it is meant for linear quadrilateral macroscale elements

Replaced all == floating point comparisons wth isclose() functions
More portable and flexible comparison for floating point numbers

Revised the user-defined input section to read another Python script
Prevents users from modifying this current script unless necessary

Revised the GP and Gaussian weight section
Allows for more flexibility to use reduced integration if desired

End of 250519 Revision

End of Revision

250630 Revision
** Major revision - modularisation
All major functions are now modularised and shared across scripts, accounting for variations
All modules are carried in another script 'DFE2_0_Modules.py'

Minor revisions
Relabel J as J^T for better clarity
Merged the RVE and MPC loops
Replaced print>> function with write() for cross-compatibility across different versions of Python

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

















