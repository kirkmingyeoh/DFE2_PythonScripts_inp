# DFE2_PythonScripts_inp

Direct FE2 is a multiscale modelling approach that was first proposed in 2020 [1]. In contrast to the staggered solution schemes of conventional FE2 approaches, Direct FE2 combines both the macroscale and microscale analyses into a single FE problem, resulting in a monolithic solution scheme. Such an approach has been shown to be more computationally efficient [2]. This is done using multi-point constraints and volume scaling, which are readily available functions in most commercial FE codes. Thus, it does not require any user-defined subroutines. 

Moreover, Direct FE2 is also much simpler to implement and use, as it does not require any user intervention to handle the data exchange between scales, nor does it require the user to specify how the tangent modulus across scales should be calculated for different problems. As of this writing, Direct FE2 has been extended for various analyses such as structural beam, plate and shell elements, strain rate dependent problems, thermal problems, composite damage, micromorphic continua, crack propagation, sandwich structures, and topology optimization problems. 

In this repository, Python scripts to help users set up a Direct FE2 analysis input file for use in Abaqus are shared. The script takes as input two Abaqus input files, one each for the macroscale and the microscale representative volume element (RVE), and then returns a separate Direct FE2 input file that incorporates the RVE into the macroscale analysis. Detailed descriptions of the input files as well as instructions to run the scripts are presented in (2) Instructions.md below. The user can then immediately submit the Direct FE2 input file as a job in Abaqus. Three examples, one each for 2D quadrilateral, one for 3D hexahedral elements, as well as one for a 2D verification case, are also presented.

This repository contains the following files:

(1) Python scripts - These are the Python scripts used to set up Direct FE2 input files, which are further detailed below.

(2) Instructions.md - Instructions on how to prepare the user-provided input files and run the Python script for setting up a Direct FE2 input file. Instructions on validating the generated Direct FE2 input file are also provided. 

(3) Demo_xx.zip - Zip folders containing simple examples to demonstrate the functionalities of the Python scripts, which are further detailed below. 

(4) Refences.md - A list of Direct FE2 references.

Kindly note that this repository is preserved to be the same as described in the manuscript 'A Python script repository for multiscale modelling with Direct FE2 in Abaqus' (DOI: xxx) for the reader's reference. As such, it will not be updated further. The live repository which will be continuously updated as improvements are made or when new user-friendly scripts for various Direct FE2 models are developed can be found at: https://github.com/kirkmingyeoh/DFE2_PythonScripts_inp_live.

-----
Python scripts
-----
The following are the Python scripts presented in this repository. 

0) **DFE2_0_UserInput.py** - This user input file is to specify the user-defined inputs, such as the names of the macroscale, microscale RVE and output Direct FE2 input files. Users can also specify additional information such as the number of integration points per macroscale element, as well as tolerance for floating point comparisons. This input file is used with and will be called by the other Python scripts to read the user-defined inputs. 

1) **DFE2_2D_Quad-Lin_QuadTri.py** - This script sets up a Direct FE2 input file for problems involving linear 2D quadrilateral elements (CPS4, CPS4R, CPE4, CPE4R) at the macroscale, and any 2D continuum elements at the microscale. 

2) **DFE2_3D_Hex-Lin_HexTetWedge.py** - This script sets up a Direct FE2 input file for problems involving linear 3D hexahedral elements (C3D8 or C3D8R) at the macroscale, and any 3D continuum elements at the microscale.

It is noted that the type of macroscale element specified in the macroscale input file, whether with reduced or full integration, does not affect the final results. This is because Direct FE2 only utilises the degrees of freedom of the macroscale element's nodes, and not the macroscale element's integration points. 

-----
Examples
-----
Three examples are presented in this repository.

1) **Demo_2D_Quad-Lin_Quad-Lin.zip** - This is a 2D example which uses linear quadrilateral elements (CPS4) at both scales. 
   
2) **Demo_3D_Hex-Lin_Hex_Lin.zip** - This is a 3D example which uses linear hexahedral elements (C3D8) at both scales. 

3) **Demon_2D_Homogeneous_Validation** - This is a 2D homogeneous example using linear quadrilateral elements (CPS4) which demonstrates a validation case for Direct FE2. 

Examples 1) and 2) are demonstrative examples to show how the Python scripts presented in this repository work, one each for 2D and 3D linear continuum macroscale elements. They each contain the user-provided macroscale and microscale input files, the scripts required, the resulting Direct FE2 input files returned by the Python script, as well as a list briefly describing each file within. 

To run these examples, extract the files from the zip folder and place them in a folder/repository. In the user input file 'DFE2_0_UserInput.py', specify a new name for the output Direct FE2 input file in line 11. Change the working directory of Abaqus CAE or the Python IDE to the folder/repository where these files are and execute the main Python script. Instructions on how to change the working directory and execute the script are provided in (2) Instructions.md. A new Direct FE2 input file with the name specified above will be generated. If done correctly, this generated input file should be identical to the provided reference. This generated input file can be submitted directly to Abaqus for solving.

Example 3) is demonstration case to show how a Direct FE2 model generated by the Python scripts can be validated. It contains the user-provided macroscale and microscale inputs files, the scripts, the resulting Direct FE2 input file returned by the Python scripts as well as the job files for the Direct FE2 model. It also contains the input file and job files for the DNS reference case. An accompanying documentation containing the model descriptions, results, comparison and a brief discussion is also provided. Detailed instructions on how to perform such a validation can be found in (2) Instructions.md in this repository. 

