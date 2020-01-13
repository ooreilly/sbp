import helper
import matplotlib.pyplot as plt
import numpy as np
import sys

data = []
grid0=0
gridn=6
ngrid = gridn - grid0 + 1
ctask="timing/collocated/output_default"
smtask="timing/staggered/output_modified"

ctime="timing_stencil/collocated/output_modified"
smtime="timing_stencil/staggered/output_modified"

# collocated, staggered
c = []
sm = []

# Configs for efficiency plot
ecc = []
ecm = []

for num in range(grid0, gridn):
    c.append(helper.load_csv("%s_%d" % (ctask, num), "u.csv"))
    sm.append(helper.load_csv("%s_%d" % (smtask, num), "u.csv"))

    ecc.append(helper.load_config("%s_%d" % (ctime, num), "log.txt"))
    ecm.append(helper.load_config("%s_%d" % (smtime, num), "log.txt"))

for num in range(gridn, gridn+1):
    ref = helper.load_csv("%s_%d" % (smtask, num), "u.csv")
    c.append(ref)
    sm.append(ref)

pc = helper.compute_errors(c, 1)
pm = helper.compute_errors(sm, 1)

tc = []
tm = []

for num in range(grid0, gridn):
    tc.append(ecc[num].avg)
    tm.append(ecm[num].avg)

fig, ax = plt.subplots()
ax.loglog(tc, pc.norm_rel_err, 'o-', label="collocated")
ax.loglog(tm, pm.norm_rel_err, 'o-', label="staggered")
plt.xlabel('Avg. time / time step (ms)')
plt.ylabel('Relative error')
ax.set_xscale('log', basex=2)
ax.set_yscale('log', basey=10)
plt.legend()
plt.savefig("efficiency.pdf")
plt.show()
