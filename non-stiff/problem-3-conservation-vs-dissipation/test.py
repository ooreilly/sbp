import numpy as np
import matplotlib.pyplot as plt
import helper

prec = np.float64

solutions = ["rk4_ed", "rk4_ec3", "srk4_opt"]
labels = ["PMF-ED", "EC3", "EC2"]
p0 = {}
p1 = {}
t = {}
dlabels = {}

for lbl, soln in zip(labels, solutions):
    p0[soln], p1[soln], t[soln] = helper.load(soln, 2, prec)
    dlabels[soln] = lbl

for i, soln in enumerate(solutions):
    plt.plot(t[soln], p0[soln] - p1[soln], "C%d" % i, label=dlabels[soln])
    #plt.plot(t[soln], , "C%d" % i, label=dlabels[soln])
#plt.plot(t, q0-q1, label="EC-OPT")
#plt.plot(t, w0-w1, label="ED")
plt.legend()
plt.show()

