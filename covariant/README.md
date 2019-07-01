# README

This directory contains the SBP staggered interpolation and differentiation
operators used in *Energy conservative SBP discretizations of the acoustic wave equation in covariant
form on staggered curvilinear grids*.


The script `test_operators.m` tests that the SBP properties are satisfied to rounding
error and reports accuracy information. 

The directory `resources/` contains the upper block and interior stencil of each
operator stored in .json format. This directory also contains binary matlab
files that are used by the provided octave scripts. 

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



