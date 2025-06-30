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
Tol = 1e-6
Sort = ''


'''
Revision log

250526 Original release

250630 Revision
Moved all explanatory notes to another document for conciseness 
Added Sort as new input for users to define the order to sort the macroscale nodes

End of 250630 Revision

End of Revision
'''








