# Entropy-stable-OFDG
This repository contains software and instructions to reproduce two numerical experiments in the paper
> An entropy stable oscillation-free discontinuous Galerkin method for hyperbolic conservation laws 
>
> * authors: Yong Liu, Jianfang Lu, Chi-Wang Shu
> * Academy of Mathematics and Systems Science, Chinese Academy of Sciences; South China University of Technology; Brown University.

# Get Started
The codes are written in Fortran language. The intel fortran compile and the intel mkl package are required.

# Burgers 1D
To run the code, you need to go to the burgers1d file directory and do the following steps:
(1) first step, you need to compile it and generate an executable file "all". The command is: csh creat.txt
(2) Second step, you need to run the executable file "all". The command is: ./all 

# Double Mach Reflection 2D
This is an MPI code. You need the MPI environment. To run the MPI code (with 8 nodes and every node has 64 cores), you need to go to the DMR file directory and do the following steps:
(1) first step, you need to compile it and get an executable file "all". The command is: csh creat.txt 
(2) Second step, to run the executable file "all" you need to submit the "job.sh". 

