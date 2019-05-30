# KineticMonteCarlo
HPC 2019 Final Project on Kinetic Monte Carlo

Anya Katsevich and Terrence Alsup

## Running the code.
The main files are

1. kmc_serial.c
2. kmc_parallel.c

and rely on the files

3. parameters_serial.txt
4. parameters.txt

respectively.  The parameters_serial.txt file contains

0 L K 1 1 T

0 is the is_dH parameter (which has to do with using symmetric rates and doesn't need to be changed). L is the length of the crystal ~400.  K is the inverse temperature ~0.5.  1 corresponds to the number of samples to take.  The second 1 correspond to the number of times to run until (just take 1 to run until some final time).  T is the final time we run until.

The parameters.txt file contains

L K T 0.1

L, K, and T are as before and c = 0.1 corresponds to the average number of events we want to happen in each section. This is the recommended value for this parameter.

To run the code first load the following modules:

1) gcc-8.1
2) mpi/openmpi-x86_64

The command to run kmc_serial is just

./kmc_serial 1

and the 1 corresponds to the random seed.  To run kmc_parallel type the command

mpirun -n 4 ./kmc_parallel 1

where 4 is the number of processors and 1 again corresponds to the random seed. (Note that each processor has a different random seed still).

## Output of code

kmc_serial will output a file called "h.txt" which contains all of the heights for each lattice site in the crystal.

kmc_parallel will give similar ouput but instead for each processor.  These must be combined when we plot the results.

## Plotting the results.
The Matlab file KMC_plotter.m can be used to plot the results.  Note that the data must be in the same directory.  For convenience, data from earlier runs is included in the data directory of this repository.
