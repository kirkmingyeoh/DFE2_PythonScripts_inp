# Instructions 

This document details the general instructions to run the Python scripts in this repository to generate Direct FE2 input files for Abaqus. These are general instructions which are applicable across all scripts. Specific instructions unique to certain scripts should be detailed within the respective scripts. 

-----
User-provided Abaqus input files
-----
Two Abaqus input files need to be provided by the user - one each to describe the macroscale problem as well as the microscale RVE. This can be generated through Abaqus' CAE interface. 

1) Macroscale input file - This input file should contain all the information detailing the macroscale problem to be analysed. The macroscale Part and Instance should be meshed and assigned an Elastic material with arbitrary property values. The script will replace the arbitrary material with null material properties. The boundary conditions should also be imposed on the macroscale Instance, and all other analysis settings such as time increment should be included. 

2) Microscale input file - This input file should contain all the information detailing the microscale RVE of the macroscale structure. The RVE should be rectangular (2D) or cuboidal (3D) and aligned along the global x, y and z axes. The RVE should also be periodic with the appropriate materials assigned to the corresponding regions. An Instance should be included in the Assembly such that the Part information will be captured when the input file is generated. 

----
User-provided information to the Python script
-----
1) 

2)

-----
Abaqus input file returned by the Python script
-----
The Python script will return 
