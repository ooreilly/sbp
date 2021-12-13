# Problem 1 - Matrix Norm behavior

This directory contains the source code used to generate the results presented in the numerical
experiment titled "Matrix norm behavior".

## Requirements
 * [sbp.jl](https://github.com/ooreilly/sbp.jl) (home-made Julia package needed to run the Julia
   codes)
 * [Julia](https://julialang.org/) (tested on 1.4, 1.7)
 * [Python3](https://www.python.org/), [Numpy](https://numpy.org/) [Matplotlib](https://matplotlib.org/) (for plotting)
 * [hdf5storage](https://pypi.org/project/hdf5storage/) (for loading MATLAB data files)
 * [pdfcrop](https://ctan.org/pkg/pdfcrop?lang=en) (to remove white space around figures)


## Initialization
Once installed, you can create the necessary data and figure directories by calling
```bash
$ make init
```
This command creates the directory `data` used for storing output and `figures` for storing the pdf
figures found in the paper.

## Naive penalty parameter choice
Reproduce the naive penalty parameter choice figures via
```bash
$ make naive naive_fig
```
 
## Optimal penalty parameter choice
Reproduce the optimal penalty parameter choice figures via
```bash
$ make optimal optimal_fig
```                        

