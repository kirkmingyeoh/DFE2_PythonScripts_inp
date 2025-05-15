# DFE2_PythonScripts_inp

Direct FE2 is a novel multiscale modelling approach that was first proposed in 2020 [1]. In contrast to the staggered solution schemes of conventional FE2 approaches, Direct FE2 combines both the macroscale and microscale analyses into a single FE problem, resulting in a monolithic solution scheme. Such an approach has been shown to be more computationally efficient [2]. This is done using multi-point constraints and volume scaling, which are readily available functions in most commercial FE codes. Thus, it does not require any user-defined subroutines. 

Moreover, Direct FE2 is also much simpler to implement and use, as it does not require any user intervention to handle the data exchange between scales, nor does it require the user to specify how the tangent modulus across scales should be calculated for different problems. As of this writing, Direct FE2 has been extended for various analyses such as structural beam, plate and shell elements, strain rate dependent problems, thermal problems, composite damage, micromorphic continua, crack propagation, sandwich structures, and topology optimization problems. 

In this repository, Python scripts to help users set up a Direct FE2 analysis input file for use in Abaqus are shared. The script takes as input two Abaqus input files, one each for the macroscale and the microscale representative volume element (RVE), and then returns a separate Direct FE2 input file that incorporates the RVE into the macroscale analysis. Detailed descriptions of the input files as well as instructions to run the scripts are presented in (2) Instructions.md below. The user can then immediately submit the Direct FE2 input file as a job in Abaqus. Two examples, one each for 2D quadrilateral and 3D hexahedral elements, are also presented.

This repository contains the following files:

(1) Python scripts - These are the Python scripts used to set up Direct FE2 input files, which are further detailed below.

(2) Instructions.md - Instructions on how to prepare the user-provided input files as well as the Python script for setting up a Direct FE2 input file. 

(3) Demo.zip - Zip folders containing simple examples for the Python scripts. Each zip folder contains the user-provided macroscale and microscale input files, as well as the resulting Direct FE2 input files returned by the Python script. 'Demo_2D_Quad_Quad.zip' runs with 'DFE2_2D_Quad_QuadTri.py', while 'Demo_3D_Hex_Hex.zip' runs with 'DFE2_3D_Hex_HexTetWedge.py'. 

Kindly note that this repository is preserved to be the same as described in the manuscript 'A Python script repository for multiscale modelling with Direct FE2 in Abaqus' (DOI: xxx) for the reader's reference. As such, it will not be updated further. The live repository which will be continuously updated as improvements are made or when new user-friendly scripts for various Direct FE2 models are developed can be found at: https://github.com/kirkmingyeoh/DFE2_PythonScripts_inp_live.

-----
Python scripts
-----

1) **DFE2_2D_Quad_QuadTri.py** - This script sets up a Direct FE2 input file for problems involving 2D quadrilateral elements (CPS4 or CPS4R) with a 2x2 integration scheme at the macroscale, and any 2D continuum elements at the microscale. 

2) **DFE2_3D_Hex_HexTetWedge.py** - This script sets up a Direct FE2 input file for problems involving 3D hexahedral elements (C3D8 or C3D8R) with a 2x2x2 integration scheme at the macroscale, and any 3D continuum elements at the microscale.

-----
Instructions for examples
-----
xxx
