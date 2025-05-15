# DFE2_PythonScripts_inp

Direct FE2 is a novel multiscale modelling approach that was first proposed in 2020 [1]. In contrast to the staggered solution schemes of conventional FE2 approaches, Direct FE2 combines both the macroscale and microscale analyses into a single FE problem, resulting in a monolithic solution scheme. Such an approach has been shown to be more computationally efficient [2]. This is done using multi-point constraints and volume scaling, which are readily available functions in most commercial FE codes. Thus, it does not require any user-defined subroutines. 

Moreover, Direct FE2 is also much simpler to implement and use, as it does not require any user intervention to handle the data exchange between scales, nor does it require the user to specify how the tangent modulus across scales should be calculated for different problems. As of this writing, Direct FE2 has been extended for various analyses such as structural beam, plate and shell elements, strain rate dependent problems, thermal problems, composite damage, micromorphic continua, crack propagation, sandwich structures, and topology optimization problems. 

In this repository, Python scripts to help users set up a Direct FE2 analysis input file for use in Abaqus are shared. The script takes as input two Abaqus input files, one each for the macroscale and the microscale representative volume element (RVE), and then returns a separate Direct FE2 input file that incorporates the RVE into the macroscale analysis. Detailed descriptions of the input files as well as instructions to run the scripts are presented in (2) Instructions.md below. The user can then immediately submit the Direct FE2 input file as a job in Abaqus. Two examples, one each for 2D quadrilateral and 3D hexahedral elements, are also presented.

This repository is preserved to be the same as described in the manuscript 'A Python script repository for multiscale modelling with Direct FE2 in Abaqus' (DOI: xxx) for the reader's reference. As such, it will not be updated further. The live repository which will be continuously updated as improvements are made or when new user-friendly scripts for various Direct FE2 models are developed can be found at: https://github.com/kirkmingyeoh/DFE2_PythonScripts_inp_live.

This repository contains the following files:

(1) Python scripts - These are the Python scripts used to set up Direct FE2 input files, which are further detailed below.

(2) Instructions.md - Instructions on how to prepare the user-provided input files as well as the Python script for setting up a Direct FE2 input file. 

(3) Demo.zip - Zip folders containing simple examples for the Python scripts. Each zip folder contains the user-provided macroscale and microscale input files, as well as the resulting Direct FE2 input files returned by the Python script. 'Demo_2D_Quad_Quad.zip' runs with 'DFE2_2D_Quad_QuadTri.py', while 'Demo_3D_Hex_Hex.zip' runs with 'DFE2_3D_Hex_HexTetWedge.py'. 

-----
Python scripts
-----

1) **DFE2_2D_Quad_QuadTri.py** - This script sets up a Direct FE2 input file for problems involving 2D quadrilateral elements (CPS4 or CPS4R) with a 2x2 integration scheme at the macroscale, and any 2D continuum elements at the microscale. 

2) **DFE2_3D_Hex_HexTetWedge.py** - This script sets up a Direct FE2 input file for problems involving 3D hexahedral elements (C3D8 or C3D8R) with a 2x2x2 integration scheme at the macroscale, and any 3D continuum elements at the microscale.

-----
References
-----
1) Tan, V.B.C., Raju, K. and Lee, H.P., 2020. Direct FE2 for concurrent multilevel modelling of heterogeneous structures. Computer Methods in Applied Mechanics and Engineering, 360, p.112694.
2) Raju, K., Zhi, J., Su, Z.C., Tay, T.E. and Tan, V.B.C., 2021. Analysis of nonlinear shear and damage behaviour of angle-ply laminates with Direct FE2. Composites Science and Technology, 216, p.109050.
3) Koyanagi, J., Kawamoto, K., Higuchi, R., Tan, V.B.C. and Tay, T.E., 2021. Direct FE2 for simulating strain-rate dependent compressive failure of cylindrical CFRP. Composites Part C: Open Access, 5, p.100165.
4) Zhi, J., Raju, K., Tay, T.E. and Tan, V.B.C., 2021. Multiscale analysis of thermal problems in heterogeneous materials with Direct FE2 method. International Journal for Numerical Methods in Engineering, 122(24), pp.7482-7503.
5) Xu, J., Li, P., Poh, L.H., Zhang, Y. and Tan, V.B.C., 2022. Direct FE2 for concurrent multilevel modeling of heterogeneous thin plate structures. Computer Methods in Applied Mechanics and Engineering, 392, p.114658.
6) Zhi, J., Poh, L.H., Tay, T.E. and Tan, V.B.C., 2022. Direct FE2 modeling of heterogeneous materials with a micromorphic computational homogenization framework. Computer Methods in Applied Mechanics and Engineering, 393, p.114837.
7) Yeoh, K.M., Poh, L.H., Tay, T.E. and Tan, V.B.C., 2022. Multiscale computational homogenisation of shear-flexible beam elements: a Direct FE2 approach. Computational Mechanics, 70(5), pp.891-910.
8) Zhi, J., Yang, B., Li, Y., Tay, T.E. and Tan, V.B.C., 2023. Multiscale thermo-mechanical analysis of cure-induced deformation in composite laminates using Direct FE2. Composites Part A: Applied Science and Manufacturing, 173, p.107704.
9) Yeoh, K.M., Poh, L.H., Tay, T.E. and Tan, V.B.C., 2023. Multiscale modelling of sandwich structured composites using direct FE2. Composites Science and Technology, 239, p.110066.
10) Zhi, J., Leong, K.H., Yeoh, K.M., Tay, T.E. and Tan, V.B.C., 2023. Multiscale modeling of laminated thin-shell structures with Direct FE2. Computer Methods in Applied Mechanics and Engineering, 407, p.115942.
11) Liu, K., Meng, L., Zhao, A., Wang, Z., Chen, L. and Li, P., 2023. A hybrid direct FE2 method for modeling of multiscale materials and structures with strain localization. Computer Methods in Applied Mechanics and Engineering, 412, p.116080.
12) Zhao, A., Tan, V.B.C., Li, P., Liu, K. and Hu, Z., 2023. A Reconstruction Approach for Concurrent Multiscale Topology Optimization Based on Direct FE2 Method. Mathematics, 11(12), p.2779.
13) Zhao, A., Li, P., Cui, Y., Hu, Z. and Tan, V.B.C., 2024. Multiscale topology optimization with Direct FE2. Computer Methods in Applied Mechanics and Engineering, 419, p.116662.
14) Li, H., Chen, L., Zhi, G., Meng, L., Lian, H., Liu, Z., Yu, T. and Li, P., 2024. A direct FE2 method for concurrent multilevel modeling of piezoelectric materials and structures. Computer Methods in Applied Mechanics and Engineering, 420, p.116696.
15) Meng, L., Zhang, H., Liu, Z., Shu, X. and Li, P., 2024. A direct FE2 method for concurrent multilevel modeling of a coupled thermoelectric problem–Joule heating effect–in multiscale materials and structures. Thin-Walled Structures, 203, p.112166.
16) Christoff, B.G., Almeida Jr, J.H.S., Ribeiro, M.L., Maciel, M.M., Guedes, R.M. and Tita, V., 2024. Multiscale modelling of composite laminates with voids through the direct FE2 method. Theoretical and Applied Fracture Mechanics, 131, p.104424.
17) Zhao, A., Liu, K., Li, P. and Cui, Y., 2024. Tunable deformation design of porous Al2O3 based on the Direct FE2 method. Modelling and Simulation in Materials Science and Engineering, 32(5), p.055015.
18) Zhao, H., Yeoh, K.M., Zhi, J. and Tan, V.B.C., 2024. Direct FE2 multiscale modeling of hydrogen-induced cracking in reactor pressure vessels. International Journal of Mechanical Sciences, 274, p.109285.
19) Dhari, R.S., Hall, W., Feih, S. and Javanbakht, Z., 2024. Direct FE2 failure prediction of fused-filament fabricated parts. Engineering Failure Analysis, 163, p.108408.
20) Zhao, H., Zheng, X., Yang, S., Yang, X. and Li, W., 2024. Direct FE2 multiscale simulation of hydrogen diffusion in Zircaloy cladding. Acta Mechanica Sinica, 40(12), pp.1-14.
21) Lan, Y., Ma, L., Du, X. and Zhou, W., 2024. Multiscale Simulation of the Coupled Chemo-Mechanical Behavior of Porous Electrode Materials by Direct FE2 Method. International Journal of Applied Mechanics, 16(3), p.2450039.
22) Cui, Y. and Zhang, Z., 2025. A novel concurrent multiscale method based on the coupling of Direct FE2 and CPFEM. Thin-Walled Structures, 206, p.112610.
