import numpy as np
import matplotlib.pyplot as plt
import helper
import sys

grid = int(sys.argv[1])

print("Grid:", grid)

prec = np.float64

solutions = ["rk4_naive", "srk4_ec", "srk4_opt", "rk4_ed"]
labels=["RK4-NA", "SRK4-EC1", "SRK4-EC2", "RK4-ED"]
p0 = {}
p1 = {}
max_err = {}
t = {}
dlabels = {}


for lbl, soln in zip(labels, solutions):
    p0[soln], p1[soln], t[soln] = helper.load(soln, grid, prec)
    max_err[soln] = np.max(np.abs(p0[soln]-p1[soln]))
    dlabels[soln] = lbl

for i, soln in enumerate(solutions):
    plt.plot(t[soln], p0[soln] - p1[soln], "C%d" % i, label=dlabels[soln])
plt.xlabel("t")
plt.ylabel("Error")
plt.legend()
filename = "figures/problem1_grid_%d.pdf" % grid
plt.savefig(filename)
print("Wrote:", filename)

print("Maximum error")
print(max_err)
print("Maximum error (SRK4_EC / RK4)")
print(max_err["srk4_ec"] / max_err["rk4_ed"])
