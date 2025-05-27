# -*- coding: utf-8 -*-
"""
Created on Mon May 19 14:00:01 2025

@author: Kirk Ming Yeoh (e0546208@u.nus.edu)

This user input file is to specify the user-defined inputs, such as the names of the macroscale, microscale RVE and output Direct FE2 input files.
Users can also specify additional information such as the number of integration points per macroscale element, as well as tolerance for floating point comparisons.
This input file is used with and will be called by the other Python scripts to read the user-defined inputs.
"""

### Input file names
MacroInpName = 'Demo_2D_Macro.inp' # Name of macroscale input file
RVEInpName = 'Demo_2D_RVE.inp' # Name of RVE input file Demo_2D_RVE
NewInpName = 'DFE2_2D.inp' # Name of new Direct FE2 input file


### Additional information
GP = ''
'''
- Specify the GPs as a list in terms of natural coordinates of the macroscale element (optional)
- E.g., GP = [[-3**-0.5,-3**-0.5],[3**-0.5,-3**-0.5],[3**-0.5,3**-0.5],[-3**-0.5,3**-0.5]] for 2D quadrilateral elements with full integration
- E.g., GP = [[0.0,0.0]] for 2D quadrilateral elements with fuly reduced integration
- Leave as '' if not intending to specify, where full integration will be used 
'''

Tol = 1e-6
'''
- Specify the tolerance when performing floating point comparisons (optional)
- Leave as default of 1e-6 or '' if not intending to specify
'''










