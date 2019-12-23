# README

This directory contains the SBP staggered interpolation and differentiation
operators as well as example code used in *Energy conservative SBP
discretizations of the acoustic wave equation in covariant form on staggered
curvilinear grids*.


## Operators

The operators have been implemented in the package
[sbp.jl](github.com/ooreilly/sbp.jl). A matlab/octave implementation can be
found in this directory.

The script `test_operators.m` tests that the SBP properties are satisfied to rounding
error and reports accuracy information. 

The directory `resources/` contains the upper block and interior stencil of each
operator stored in .json format. This directory also contains binary matlab
files that are used by the provided matlab/octave scripts. 

To load and test the operators used in all of the numerical experiments, 
```bash
$ octave test_operators.m build3
```
As
can be seen from the output, the interpolation operators have a norm that is
close to unity,
```
||PpPm||_2 =  1.03625
```
Here, `Pm` denotes `\hat{P}` in the paper and `Pp` denotes `P` in the paper. 
We also include the operators with large norms,
```
$ octave test_operators.m build7
||PpPm||_2 =  8.52851
```


## Simulation of acoustic waves over a Gaussian hill

The final part of the paper includes a numerical examples that features acoustic
waves interacting with a curved boundary. By following these instructions, it
should be possible to reproduce this simulation. 

1. You need to install the package [sbp.jl](github.com/ooreilly/sbp.jl). See
that repository for how to install the package.

2. Run the script `simulation.jl` to perform the simulation. This script takes
   several command-line arguments, as can be seen from its header:

        usage: 
        simulation.jl <outputdir> <builddir> <refine> <scheme> <metric_tensor> <run> <vtk>
        
             builddir path      Directory to write discretization.
             outputdir path     Directory to write output data to.
             refine int         Level of grid refinement to use. refine = 0 (no grid
                                                                          refinement)
             scheme string      Choose between 'staggered' or 'collocated'
             metric_tensor string  Set to 'modified' to use the modified metric tensor
                                   discretization.      
             run     int        1 (Run simulation in Julia) 0 (Otherwise) 
             vtk     int        How often to write vtk files
        
        This example code solves the acoustic wave equation in covariant form on a
        curvilinear staggered grid
        
        The governing equations are discretized in space and the zero pressure boundary
        condition is applied on all boundaries. The resulting discretization is stored
        in a sparse matrix `A`. The solution is advanced in time by solving the linear
        system
        
        dq/dt = A*q
        
        using a low-storage Runge-Kutta time integration scheme.

 Perhaps the easiest way to get started is to use the provided `Makefile` that
 wraps this script and other scripts. Feel free to open up the `Makefile` in
 your favorite editor to learn more about what it can do. The following commands
 run a collocated and staggered grid simulation at the coarsest resolution in
 Julia.

        make init build run_w_julia=1 scheme=collocated
        make init build run_w_julia=1 scheme=staggered

3. After the simulation has completed, the directory `default` should have been
   created for you and this directory should contain the subdirectories
   ```
   collocated/output_modified_0
   staggered/output_modified_0
   ```
   These directories contains a sequence of .vtk files that can be opened up in
   [Paraview](https://paraview.org) for visualization purposes.


### Matrix-based GPU code
While that should be pretty much all you need to do run the example code, you
will notice that simulations run extremely slow once the grid refinement number
`refine` increases. For this reason, we wrote a CUDA C++ code that runs on the
GPU. 

1. To use this code, you need to have the CUDA toolkit installed, and make
sure that you install it with both CUBLAS and CUSparse. 

2. Download the header-only library [lmv](https://github.com/ooreilly/lmv) and create an
   environment variable called `LMV_INCLUDE_DIR` set it to the lmv include dir,
   for example, 

   ```
   export LMV_INCLUDE_DIR=~/codes/lmv/include      
   ```

3. Run cmake to build the GPU solvers
        
   ```
   mkdir build
   cd build
   cmake ..
   ```

4. If you ran the Julia script `simulation.jl`, it will have written several
   files that are needed to run the GPU solver. Otherwise, you will need to
   rerun the Julia script. Using the provided `Makefile`, we can accomplish this
   and disable running with Julia. For example,

   ```
   make init build run run_w_julia=0 refine=4 scheme=collocated
   ```

   The command `run` is responsible for executing the solver.  After completion,
   the output directory will again contain vtk files, but also a .csv file that
   holds the pressure field at the receiver location for each time step. The
   first column in this file is the time, and the second is the pressure. Some
   runtime diagnostics are written to the file `log.txt`.


### Stencil-based GPU code
The stencil-based GPU code does only include an implementation of the interior
computation. Hence, the SBP boundary stencil computations have not implemented.
The original intent of this code is to investigate the computational efficiency
of the collocated and staggered schemes. In most cases, the boundary makes up a
small fraction of the computational domain therefore its computational cost is
insignificant compared to the interior computation. The stencil-based solver can
be run by using the provided makefile. As for the matrix-based solver, this
solver also requires that a preprocessing step is run by calling a Julia script.
The script `stencil.jl` takes care of this part. 


Since this code is meant for timing the kernels, you will probably want to
switch off the vtk output option to avoid writing snapshots of the wavefield to
disk. The following command handles both preprocessing and running the solver.
```
 make init build-stencil run-stencil vtk_stride=0 scheme=collocated
```
or 
```
 make init build-stencil run-stencil vtk_stride=0 scheme=staggered
```
for running the staggered scheme instead of the collocated scheme.

The output file `log.txt` will be written to the output directory after
completion. This file contains some statistics about the simulation.
