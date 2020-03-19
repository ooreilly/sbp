# Energy conservation vs Energy dissipation

This is the source code used to produce the numerical experiment titled "Energy conservation vs
energy dissipation. Unlike the previous numerical experiments, which are 1D dimensional, this
experiment is 2D dimensional. The source code is written in C++ and CUDA and requires a NVIDIA CUDA
compatible graphics card to run.

## Requirements
* CUDA (11.3 tested)
* CMake 3.10 or higher
* ParaView (for visualization)
* Python 3 (for generating figures)

## Installation
Use CMake to compile the code
```bash
$ mkdir -p build
$ cd build
$ cd cmake ..
$ make
```
The above commands should create the executable `acoustic`.

## Usage 
To display usage information, type
```
$ ./acoustic
usage: ./acoustic <name> <refinement> <time_stepping_scheme> <nx> <ny> <h> <nu> <vtk_stride> <T> <fp> <t0> <sx> <sy> 
Options
 name:                  Output filename 
 refinement:            Grid refinement number 
 time_stepping_scheme:  0 = leap-frog, 1 = staggered RK4, 2 = RK4 
 coupling method:       2 = energy conserving, 3 = energy dissipating, 4 = energy conserving (naive)
 nx:                    Grid points x 
 ny:                    Grid points y 
 Lx:                    Domain size x 
 Ly:                    Domain size y 
 nu:                    Time step scaling factor, dt = nu * h
 vtk_stride:            VTK output stride factor. 
                        Output every `vtk_stride` time steps (also influenced by refinement)
                        Use negative number to disable
 stride:            Time series output stride factor. 
  T:            Simulation end time. 
  fp:            Source central frequency. 
  t0:            Source start time. 
  sx:            Source x-coordinate. 
  sy:            Source y-coordinate. 
```
## Pre-configured simulations
The provided `Makefile` can be used to run some of the preconfigured simulations presented in the
paper. Start by calling the initialization command
```bash
$ make init
```
This command will create directories to place output, logs, and figures in.

Open the `Makefile` to explore more options. For example,
```
$ make rk4_ed refine=1
```
runs the RK4 Energy dissipation (ed) experiment at grid refinement level `1`.

## Visualization
Use any of the `vtk_*` commands to output data in the VTK file format. For example, to output
simulation data for the RK4 Energy Dissipation experiment, use
```bash
$ make vtk_ed
```
The resulting output gets dumped to the directory `output`. Open one of the ".vtk" files in
ParaView. The naming convention should be self-explanatory. Nonetheless, `step_0` denotes the
initial condition and `step_1` is the resulting output after completing one time step. Block
`block_0` is to the left (top in the paper) and `block_1` is to the left (bottom in the paper).

## Error in time figures
Use script `error_in_time.py` to regenerate the error in time plots for a given grid refinement
level, e.g.,
```
$ python3 error_in_time.py <grid>
```
You need to generate the output before you can use this script. Otherwise, it will complain.

To reproduce the figures found in the paper, use 
```
$ python3 error_in_time.py 0
$ python3 error_in_time.py 1
```
These figures get exported as pdf files and placed in the directory `figures`.

## Computational efficiency study
The computational efficiency should be reproducible by running the bash script
```bash
$ ./efficiency.sh
```
Expect it to take a few hours to complete. Please keep in mind that timings are platform and compiler version-dependent. 

Once the output data has been generated, you can reproduce the computational efficiency and modeled speedup figures using
```bash
$ python3 efficiency.py
```




