# DFE2_PythonScripts_inp

-----
Introduction
-----

Direct FE2 is a multiscale modelling approach that was first proposed in 2020 [1]. In contrast to the staggered solution schemes of conventional FE2 approaches, Direct FE2 combines both the macroscale and microscale analyses into a single FE problem, resulting in a monolithic solution scheme. Such an approach has been shown to be more computationally efficient [2]. This is done using multi-point constraints and volume scaling, which are readily available functions in most commercial FE codes. Thus, it does not require any user-defined subroutines. 

Moreover, Direct FE2 is also much simpler to implement and use, as it does not require any user intervention to handle the data exchange between scales, nor does it require the user to specify how the tangent modulus across scales should be calculated for different problems. As of this writing, Direct FE2 has been extended for various analyses such as structural beam, plate and shell elements, strain rate dependent problems, thermal problems, composite damage, micromorphic continua, crack propagation, sandwich structures, and topology optimization problems. 

In this repository, Python scripts to help users set up a Direct FE2 analysis input file for use in Abaqus are shared. The script takes as input two Abaqus input files, one each for the macroscale and the microscale representative volume element (RVE), and then returns a separate Direct FE2 input file that incorporates the RVE into the macroscale analysis. Detailed descriptions of the input files as well as instructions to run the scripts are presented in (2) Instructions.md below. The user can then immediately submit the Direct FE2 input file as a job in Abaqus. Three examples, one each for 2D quadrilateral, one for 3D hexahedral elements, as well as one for a 2D verification case, are also presented.

This repository contains the following files and folder:

(1) Python scripts - These are the Python scripts used to set up Direct FE2 input files, which are detailed below.

(2) Instructions.md - Instructions on how to prepare the user-provided input files and run the Python script for setting up a Direct FE2 input file. Instructions on validating the generated Direct FE2 input file are also provided. 

(3) Demo - A folder containing simple examples to demonstrate the functionalities of the Python scripts, which are further detailed below. 

(4) Refences.md - A list of Direct FE2 references.

(5) UserInput_Notes.md - Notes on specifying user-defined inputs in DFE2_0_UserInput.py. 

-----
Examples
-----
Three examples are presented in the Demo folder of this repository.

1) **Demo_2D_Quad-Lin_Quad-Lin** - This is a 2D example which uses linear quadrilateral elements (CPS4) at both scales. 
   
2) **Demo_3D_Hex-Lin_Hex_Lin** - This is a 3D example which uses linear hexahedral elements (C3D8) at both scales. 

3) **Demo_2D_Homogeneous_Validation** - This is a 2D homogeneous example using linear quadrilateral elements (CPS4) which demonstrates a validation case for Direct FE2. 

Examples 1) and 2) are demonstrative examples to show how the Python scripts presented in this repository work, one each for 2D and 3D linear continuum macroscale elements. They each contain the user-provided macroscale and microscale input files, the scripts required, the resulting Direct FE2 input files returned by the Python script, as well as a list briefly describing each file within. 

To run these examples, download the input files from the respective subfolder in Demo and place them in a local folder/repository. Download 'DFE2_0_UserInput.py', 'DFE2_0_Modules.py' and the appropriate main Python script ('DFE2_2D_Quad-Lin_QuadTri.py' for Example 1) and 'DFE2_3D_Hex-Lin_HexTetWedge.py' for example 2)) from the main folder of the repository and place them in the same local folder/repository. In the user input file 'DFE2_0_UserInput.py', specify the names of the macroscale input file, RVE input file and final Direct FE2 input file in lines 13, 14 and 15 respectively. Change the working directory of Abaqus CAE or the Python IDE to the folder/repository where these files are and execute the main Python script. Instructions on how to change the working directory and execute the script are provided in (2) Instructions.md. A new Direct FE2 input file with the name specified above will be generated. If done correctly, this generated input file should be identical to the provided reference. It is noted that differences in the level of precision for floating point numbers may be observed depending on whether a Python IDE or Abaqus CAE used to execute the main script. This generated Direct FE2 input file can be submitted directly to Abaqus for solving. 

In each of these examples, a checking script is also provided in the respective folders. The checking script extracts the outputs at specific nodes from the final Abaqus ODB file and compares them against reference values. This provides a check to ensure that the main Python scripts generated the Direct FE2 input file correctly, and that it is solved correctly by Abaqus. Note that this checking script needs to be placed in the same local folder/repository as the ODB file of the respective demonstrative example and run using either Abaqus CAE or Abaqus CAE nogui. Abaqus CAE or Abaqus CAE nogui is required as the script uses Abaqus in-built modules to read the ODB files. If the name of the generated Direct FE2 input file or the final ODB file was changed, line 25 of the checking script needs to be changed from the default accordingly. 

Example 3) is demonstration case to show how a Direct FE2 model generated by the Python scripts can be validated. It contains the user-provided macroscale and microscale inputs files, the scripts, the resulting Direct FE2 input file returned by the Python scripts as well as the job files for the Direct FE2 model. It also contains the input file and job files for the DNS reference case. An accompanying documentation containing the model descriptions, results, comparison and a brief discussion is also provided. Detailed instructions on how to perform such a validation can be found in (2) Instructions.md in this repository. 

-----
Python scripts
-----
Unless otherwise stated, these scripts assume the following:  
1) Macroscale - all macroscale elements are homogenised with the same RVE and same number of integration points. It is noted that the type of macroscale element specified in the user-defined macroscale input file, whether with reduced or full integration, does not affect the final results. This is because Direct FE2 only utilises the degrees of freedom of the macroscale element's nodes, and not the macroscale element's integration points.
2) RVE - rectangular (2D) or cuboidal (3D) geometry; perfectly periodic mesh  
  
  
The Python scripts presented in this repository are detailed below:  

1) **DFE2_0_UserInput.py** - User input file to specify user-defined inputs. To be used with and wil be called by all other Python scripts.  
2) **DFE2_0_Modules.py** - Contains all the modules which will be called by the main scripts as listed below to set up the Direct FE2 input files.  
3) **Main Python scripts** - Reads the two user-provided Abaqus input files, calls the necessary modules and performs the necessary calculations to set up the Direct FE2 model, and writes the final Direct FE2 input file. The scripts are generally named in the format 'DFE2_Dimension_MacroEle-FurtherInfo_RVEEle-FurtherInfo_FurtherModelInfo.py'.

| Name | Dimension | Macro element | RVE element | Additional details |
| :-----: | :-----: | :-----: | :-----: | :-----: |
| DFE2_2D_Quad-Lin_QuadTri.py | 2D | Quadrilateral - linear | Quadrilateral or triangle | - |
| DFE2_3D_Hex-Lin_HexTetWedge.py | 3D | Hexahedral - linear | Hexahedral, wedge or tetrahedral | - |




