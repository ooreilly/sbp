import sympy as sp
import pmf

# Energy dissipating penalty parameter choice
penalty_parameter = lambda R : 0 * R

# Energy conserving penalty parameter choice
 penalty_parameter = lambda R : -R.T

rho1, c1, rho2, c2 = sp.symbols("rho1 c1 rho2 c2", real=True, positive=True)
Z1, Z2 = sp.symbols("Z1 Z2", real=True, positive=True)
nx, ny = sp.symbols("nx ny", real=True)
p1, vx1, vy1 = sp.symbols("p1 vx1 vy1")
p2, vx2, vy2 = sp.symbols("p2 vx2 vy2")


nx = 1
ny = 0
# A-hat presented in the energy-method section (before Proposition 1)
A = sp.zeros(6)
A[0,1] = Z1 * nx 
A[0,2] = Z1 * ny
A[1,0] = Z1 * nx 
A[2,0] = Z1 * ny

A[3,4] = -Z2 * nx
A[3,5] = -Z2 * ny
A[4,3] = -Z2 * nx
A[5,3] = -Z2 * ny

# Solution vector
u = sp.Matrix([p1 / Z1, vx1, vy1, p2 / Z2, vx2, vy2])

# Reflection matrix
Rf = lambda x1, x2:  sp.Matrix([[x2 - x1, 2 * sp.sqrt(x1*x2)], 
               [2 * sp.sqrt(x1*x2), x1 - x2]]) / (x1 + x2)

R = Rf(Z1, Z2)
Im = sp.eye(R.shape[0])

dR = penalty_parameter(R)

# Stability condition
print("The coupling is stable if the eigenvalues of the following matrix are bounded by 1")
sp.pprint(sp.simplify(dR*dR.T))

# Eigendecomposition
lm = sp.Matrix([[sp.sqrt(Z1), 0],[0, sp.sqrt(Z2)]])
lp = sp.Matrix([[sp.sqrt(Z1), 0],[0, sp.sqrt(Z2)]])
Xp = sp.Matrix([[1, nx, ny,0,0,0], [0,0,0,1,-nx,-ny]]).T / sp.sqrt(
        nx**2 + ny**2 + 1)
Xm = sp.Matrix([[1, -nx, -ny,0,0,0], [0,0,0,1,nx,ny]]).T / sp.sqrt(
        nx**2 + ny**2 + 1)
wm = lm * Xm.T * u
wp = lp * Xp.T * u
assert -Xm * lm**2 * Xm.T + Xp * lp**2 * Xp.T == A

assert  u.T * A * u == sp.simplify(wp.T * wp  - wm.T * wm)
print("Transformed coupling conditions")
sp.pprint(sp.simplify(wm - R*wp))

print("Energy rate for the weak primal condition")
Es = wp.T * wp - wp.T * R.T * R * wp
sp.pprint(sp.simplify(Es))

print("Energy rate for the weak dual condition")
Ed =  wm.T * wm - wm.T * dR.T * dR * wm
sp.pprint(sp.simplify(Ed))

# Projection matrix obtained using the Projection matrix formula
P = pmf.projection_matrix_formula(lm, lp, Xp, Xm, R, dR)

print("Complete energy rate")
sp.pprint(sp.simplify(u.T * (A - 2 * A * P) * u))


print("Projection matrix")
sp.pprint(sp.simplify(P))


# Penalty term
print("Penalty term")

Rhoi = 0 * A
for i in range(3):
    Rhoi[i,i] = 1/rho1
for i in range(3):
    Rhoi[3+i,3+i] = 1/rho2

sp.pprint(sp.simplify(Rhoi * A * P * u))

