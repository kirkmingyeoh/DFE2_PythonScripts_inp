# UserInput_Notes

This document provides explanatory notes for the user-specified inputs to be provided in DFE2_0_UserInput.py

-----
Input file names
-----
MacroInpName
- Name of macroscale input file
- Input file should contain all the information detailing the macroscale problem to be analysed

RVEInpName
- Name of RVE input file
- Input file should contain all the information detailing the microscale RVE of the macroscale structure

NewInpName
- Name of new Direct FE2 input file
- This input file will be returned by the script

-----
Additional information
-----
GP
- Specify the GPs as a list in terms of natural coordinates of the macroscale element (optional)
- E.g., GP = [[-3**-0.5,-3**-0.5],[3**-0.5,-3**-0.5],[3**-0.5,3**-0.5],[-3**-0.5,3**-0.5]] for 2D quadrilateral elements with full integration
- E.g., GP = [[0.0,0.0]] for 2D quadrilateral elements with fuly reduced integration
- Leave as '' if not intending to specify, where full integration will be used

Tol 
- Specify the tolerance when performing floating point comparisons (optional)
- Leave as default of 1e-6 or '' if not intending to specify

Sort
- Specify how the macroscale nodes are sorted
- '' will use Abaqus' default nodal connectivity order; 'global' will align natural coordinates tsi,eta,zeta to global coordinates x,y,z respectively as much as possible;
  'user' will follow a user-defined method which the user would need to specify in SortMacroNodes module in DFE2_0_Modules.py. The latter two options will be added in a 
  future update
