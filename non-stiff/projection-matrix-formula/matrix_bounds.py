import sympy as sp

nx, ny, c1, c2, r1, r2, k1, l1 = sp.symbols("nx ny c1 c2 r1 r2 k1 l1", real=True)
# Pick a particular side to focus on to simplify the Symbolic computations
nx = 1
ny = 0

a, b = sp.symbols("a b", real=True, positive=True)

LT = sp.Matrix([[r1 * c1, 0, 0, -r2 * c2, 0, 0],
                [0, nx, ny, 0, -nx, -ny]])


# Optimal with polynomial power
p = sp.symbols("p", integer=True, positive=True)
p = 1
f1 = (c1 * r1) ** p
f2 = (c2 * r2) ** p
k1 = f2 / (f1 + f2)
l1 = f1 / (f1 + f2)

# Naive
 k1 = sp.Rational(1, 2)
 l1 = sp.Rational(1, 2)

s11 = 0
s22 = 0
s32 = 0
s12 = c1 * k1
s21 = nx * l1 / r1
s31 = ny * l1 / r1

S = sp.Matrix([[s11, s12],
              [s21, s22],
              [s31, s32]])
V = S * LT
print("\| . \|_2 = max(")
sp.pprint(sp.factor(V.singular_values()))
print(")")
print("\| . \|_F^2 = ")
v = sum([vi**2 for vi in V[:]])
sp.pprint(sp.simplify(v))
