# Matrix bounds and the Projection matrix formula
This directory contains example code to demonstrate how to manually bound penalty parameters using
matrix norms and how to use the Projection Matrix Formula (PMF) to derive penalty terms.

## Requirements
 * [Python3](https://www.python.org/)
 * [NumPy](https://numpy.org/)
 * [SymPy](https://www.sympy.org/en/index.html)

## Matrix bounds
The script `matrix_bounds.py` can be used to reproduce the bounds found in the Section titled
"Suboptimal scaling". To obtain the Frobenius norm bound for the Naive penalty parameter choice, run
the script as is
```bash
$ python3 matrix_bounds.py
\| . \|_2 = max(
⎡   ___________________                     ⎤
⎢  ╱   2   2     2   2                      ⎥
⎢╲╱  c₁ ⋅r₁  + c₂ ⋅r₂    √2⋅│c₁│            ⎥
⎢──────────────────────, ───────, 0, 0, 0, 0⎥
⎣        2⋅│r₁│             2               ⎦
)
\| . \|_F^2 = 
    2     2   2
3⋅c₁    c₂ ⋅r₂ 
───── + ───────
  4          2 
         4⋅r₁  

```
To obtain the bound for the non-stiff penalty parameter choice, comment out
the Naive coupling parameter
```python
# Naive
# k1 = sp.Rational(1, 2)
# l1 = sp.Rational(1, 2)

```
and run the script again

```bash
$ python3 matrix_bounds.py
\| . \|_2 = max(
⎡   ___________________                                    ⎤
⎢  ╱   2   2     2   2                                     ⎥
⎢╲╱  c₁ ⋅r₁  + c₂ ⋅r₂  ⋅│c₁│  √2⋅│c₁│⋅│c₂│⋅│r₂│            ⎥
⎢───────────────────────────, ─────────────────, 0, 0, 0, 0⎥
⎣      │c₁⋅r₁ + c₂⋅r₂│         │c₁⋅r₁ + c₂⋅r₂│             ⎦
)
\| . \|_F^2 = 
  2 ⎛  2   2       2   2⎞
c₁ ⋅⎝c₁ ⋅r₁  + 3⋅c₂ ⋅r₂ ⎠
─────────────────────────
                    2    
     (c₁⋅r₁ + c₂⋅r₂)
```
This is the bound presented in the end of Section "Revisiting the motivating example" 

## Projection Matrix Formula

Run the script to obtain the penalty terms for the energy conservative coupling
```bash
$ python3 penalty_terms.py
...
Penalty term
⎡Z₁⋅Z₂⋅(vx₁ - vx₂)⎤
⎢─────────────────⎥
⎢   ρ₁⋅(Z₁ + Z₂)  ⎥
⎢                 ⎥
⎢  Z₁⋅(p₁ - p₂)   ⎥
⎢  ────────────   ⎥
⎢  ρ₁⋅(Z₁ + Z₂)   ⎥
⎢                 ⎥
⎢        0        ⎥
⎢                 ⎥
⎢Z₁⋅Z₂⋅(vx₁ - vx₂)⎥
⎢─────────────────⎥
⎢   ρ₂⋅(Z₁ + Z₂)  ⎥
⎢                 ⎥
⎢  Z₂⋅(p₁ - p₂)   ⎥
⎢  ────────────   ⎥
⎢  ρ₂⋅(Z₁ + Z₂)   ⎥
⎢                 ⎥
⎣        0        ⎦
```
To obtain the penalty terms for the energy dissipative coupling, modify the script by commenting out the energy conserving choice:
```bash
# Energy dissipating penalty parameter choice
penalty_parameter = lambda R : 0 * R

# Energy conserving penalty parameter choice
# penalty_parameter = lambda R : -R.T
```
Feel free to explore other options.


