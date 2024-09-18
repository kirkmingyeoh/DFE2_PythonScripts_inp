# DFE2_PythonScripts_inp

Direct FE2 is a novel multiscale modelling approach that was first proposed in 2020 [1]. In contrast to the staggered solution schemes of traditional FE2 approaches, Direct FE2 combines both the macroscale and microscale analyses into one, resulting in a monolithic solution scheme. Such an approach has been shown to be more computationally efficient [x]]. This is done using multi-point constraints, which are readily available functions in most commercial FE codes. Thus, it does not require any special subroutines to be developed. 

Moreover, Direct FE2 is also much simpler to implement and use, as it does not require any user intervention to handle the data exchange between scales, nor does it require the user to specify how the tangent modulus across scales should be calculated for different problems. As of this writing, Direct FE2 has been extended for various analyses such as structural beam, plate and shell elements, strain rate dependent problems, thermal problems, composite damage, micromorphic continua, crack propagation, sandwich structures, and topology optimization problems. A full list of Direct FE2 references are presented at the end of this document. 

In this repository, Python scripts to help users set up a Direct FE2 analysis input file for use in Abaqus are shared. The script takes as input two Abaqus input files, one each for the macroscale and the microscale representative volume element (RVE), and then returns a separate Direct FE2 input file that incorporates the RVE into the macroscale analysis. Detailed descriptions of the input files as well as instructions to run the scripts are presented in (2) Instructions.docx below. The user can then immediately submit the Direct FE2 input file as a job in Abaqus. Two examples, one each for 2D and 3D quadrilateral elements, are also presented.

This repository contains the following files:

(1) Python scripts - These are the Python scripts used to set up Direct FE2 input files, which are further detailed below.

(2) Instructions.docx - Instructions on how to prepare the user-provided input files as well as the Python script for setting up a Direct FE2 input file. 

-----
Python scripts
-----

1) **** - explain 

-----
References
-----

1) 
