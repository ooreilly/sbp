import numpy as np
import matplotlib.pyplot as plt
import helper

show_plot = False
poly_exp=3
rk4="rk4_ed"
srk4="srk4_ec"
srk4_ec="srk4_opt"
naive="rk4_naive"
num_grids = 6
grids = range(0, num_grids)
h0 = 0.01
h = [0] * num_grids
lam0 = 0.05
lam_h = [0] * num_grids
labels=["RK4-NA", "SRK4-EC1", "SRK4-EC2", "RK4-ED"]

def show():
    if show_plot:
        plt.show()

prec = np.float64

def plot_solution(name, r0=0, r1=1, color="C0"):
    for ref in range(r0, r1 + 1):
        plt.plot(t, p, color)


soln = {
        naive : helper.Struct(),
        srk4 : helper.Struct(), 
        srk4_ec : helper.Struct(),
        rk4 : helper.Struct()
        }

for si in soln:
    soln[si].t = []
    soln[si].p = []
    soln[si].err = []
    soln[si].time = []
    soln[si].coeff = []

for si in soln:
    for grid in grids:
        p0, p1, t = helper.load(si, grid, prec)
        h[grid] = h0 * 0.5 ** grid
        lam_h[grid] = lam0 / h[grid]
        soln[si].p.append(p0-p1)
        soln[si].t.append(t)
        soln[si].time.append(helper.parse_log("logs/%s" % si, grid))
        soln[si].err.append(np.max(np.abs(p0-p1)))

print(soln[srk4_ec].err)
print(soln[naive].err)
print(soln[rk4].err)
table = [h[:grid]]
print(labels)
for si in soln:
    table.append(soln[si].err)
    rate = helper.convergence_rates(soln[si].err)
    table.append(rate)

print(helper.latex_table(table))
 
fig, ax = plt.subplots()
ax.set_xscale('log', basex=2)
ax.set_yscale('log', basey=10)


for si in soln:
    ax.plot(lam_h, soln[si].err,"o-")
plt.xlabel("Minimum grid points per wavelength ($\lambda / h$)")
plt.ylabel("Error")
ax.legend(labels)
plt.savefig("figures/convergence.pdf")
show()

fig, ax = plt.subplots()
ax.set_xscale('log', basex=2)
ax.set_yscale('log', basey=10)

for si in soln:
    ax.plot(soln[si].time, soln[si].err,"o-")
plt.xlabel("Computational time (s)")
plt.ylabel("Error")
ax.legend(labels)
plt.savefig("figures/efficiency.pdf")
show()


def poly(c, x):
    return sum([ci * x ** p for p, ci in enumerate(c[::-1])])

for si in soln:
    d = np.polyfit(np.log10(soln[si].err),
            np.log10(soln[si].time), poly_exp, full=True)
    soln[si].coeff = d[0]


time1= soln[rk4].time
time2= soln[srk4_ec].time
err1 = soln[rk4].err
err2 = soln[srk4_ec].err
plt.clf()
x = np.linspace(-1.0, -6.0, 100)
y1 = poly(soln[rk4].coeff, x)
y2 = poly(soln[srk4_ec].coeff, x)
plt.plot(np.log10(err1), np.log10(time1), "C0-x")
plt.plot(np.log10(err2), np.log10(time2), "C1-x")
plt.plot(x, y1)
plt.plot(x, y2)
show()
plt.figure()
plt.semilogx(10**x, 10**y1/10**y2)
plt.xlabel("Relative error")
plt.ylabel("Speedup")
plt.savefig("figures/speedup.pdf")
show()

