# Problem 2 - Superconvergence

This directory contains the source code used to generate the results presented in the numerical
experiment titled "Superconvergence".

## Requirements

 * [sbp.jl](https://github.com/ooreilly/sbp.jl) (home-made Julia package needed to run the Julia
   codes)
 * [Julia](https://julialang.org/) (tested on v1.7)

## Initialization
Once installed, you can create the necessary data directory by calling
```bash
$ make init
```
This command creates the directory `data` used for storing output and `figures` for storing the pdf
figures found in the paper.

## Energy dissipative penalty parameter choice
Reproduce the data, tables, and figures for the energy dissipative penalty parameter choice, use
```bash
$ make dissipative dissipative-table
```
The command `make dissipative` calls the Julia script `convergence.jl` which is quite slow despite
running a 1D convergence test. To speed it up, you can change the hard-coded parameter
`nrefinements` that specifies the number of grid refinements to use.

The latter command will output Latex tables to stdout that can be copied into a Latex document.
To 
 
## Energy conservative penalty parameter choice
Much of the same instructions for reproducing the convergence study for the energy dissipative
treatment also applied to the energy conservative treatment. Use:
```bash
$ make conservative conservative-table 
```                        

