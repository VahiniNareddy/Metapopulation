# Metapopulation 
Ising representation of metapopulation models
---------------------------------------------

We studied metapopulation models using coupled lattice maps in two-dimensions. The metapopulation lattice configurations are represented using dynamical Ising model with memory as both the models exist in Ising universality class. The parameters for the metapopulation model are growth parameter, dispersal and noise whereas the parameters of the Ising model with memory are memory and coupling. We showed that working with Ising parameters rather than metapopulation parameters is useful in many ways in our paper.  

Getting started
---------------
These instructions will help you simulate metapopulation Model A (see our paper for details) and obtain the lattice configurations which are later used to infer the Ising parameters.

Prerequisites
-------------
You need to be able to run C++ code for metapopulation simulations and Mathematica code for Ising inference.

Metapopulation simulation
-------------------------
The dynamical details of the code (ModelA.cpp) is given in the paper. The code requires GNU scientific library for the random numbers.

To compile the code,

        g++ ModelA.cpp -lgsl -lgslcblas -o ModelA.out

then use the generated executable (as an example),

        ./ModelA.out 937595 60 0.1162 0.1

Here 937595 is a random number, 60 is the lattice size , 0.1162 is the noise and 0.1 is the dispersal

Once the execution is done, there should be five output files in the folder.
observables_..txt
state_size_….mcs_0.txt
state_size_….mcs_2999997.txt
state_size_….mcs_2999998.txt
state_size_….mcs_2999999.txt

The last three files are the consecutive lattice configurations used to obtain the Ising parameters using maximum likelihood inference methods using Mathematica.

Inference methods
-----------------
The details of the methods used in the Wolfram Notebook (inference_methods) are given in appendix D. 

To obtain Ising parameters, the inputs are the three consecutive lattice configurations. Import the files by entering the appropriate filepath.
