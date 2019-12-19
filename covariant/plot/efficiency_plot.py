import helper
import matplotlib.pyplot as plt
import numpy as np
import sys

data = []
grid0=0
gridn=6
ngrid = gridn - grid0 + 1
ctask="timing/collocated/output_default"
sdtask="timing/staggered/output_default"
smtask="timing/staggered/output_modified"

ctime="timing_stencil/collocated/output_modified"
smtime="timing_stencil/staggered/output_modified"

# collocated, staggered modified, staggered default data
c = []
sm = []
sd = []

# Configs for efficiency plot
ecc = []
ecm = []

# Configs for convergence plot
cc = []
cm = []
cd = []

for num in range(grid0, gridn):
    c.append(helper.load_csv("%s_%d" % (ctask, num), "u.csv"))
    sm.append(helper.load_csv("%s_%d" % (smtask, num), "u.csv"))

    sd.append(helper.load_csv("%s_%d" % (sdtask, num), "u.csv"))
    ecc.append(helper.load_config("%s_%d" % (ctime, num), "log.txt"))
    ecm.append(helper.load_config("%s_%d" % (smtime, num), "log.txt"))

    cc.append(helper.load_config("%s_%d" % (ctask, num), "log.txt"))
    cm.append(helper.load_config("%s_%d" % (smtask, num), "log.txt"))
    cd.append(helper.load_config("%s_%d" % (sdtask, num), "log.txt"))

for num in range(gridn, gridn+1):
    ref = helper.load_csv("%s_%d" % (smtask, num), "u.csv")
    c.append(ref)
    sm.append(ref)
    sd.append(ref)

pc = helper.compute_errors(c, 1)
pm = helper.compute_errors(sm, 1)
pd = helper.compute_errors(sd, 1)

xc = []
xm = []
xd = []
tc = []
tm = []
td = []
err_c = []
err_m = []
err_d = []

for num in range(grid0, gridn):
    err_c.append(np.linalg.norm(pc.err[num]))
    err_m.append(np.linalg.norm(pm.err[num]))
    err_d.append(np.linalg.norm(pd.err[num]))
    tc.append(ecc[num].avg)
    tm.append(ecm[num].avg)

    xc.append(cc[num].m)
    xm.append(cm[num].m)
    xd.append(cd[num].m)

print("collocated", tc)
print("staggered", tm)
print("ratio", [b/a for a, b in zip(tc, tm)])

print("collocated", err_c)
print("staggered", err_m)

fig, ax = plt.subplots()
ax.loglog(tc, err_c, 'o-', label="collocated")
ax.loglog(tm, err_m, 'o-', label="staggered-modified")
plt.xlabel('Avg. time / time step (ms)')
plt.ylabel('Error')
ax.set_xscale('log', basex=2)
ax.set_yscale('log', basey=2)
plt.legend()
plt.savefig("efficiency.pdf")
plt.show()

fig, ax = plt.subplots()
ax.loglog(xc, err_c, 'o-', label="collocated")
ax.loglog(xm, err_m, 'o-', label="staggered-modified")
ax.loglog(xd, err_d, 'o-', label="staggered-energy-pp")
plt.xlabel('Number of grid points')
plt.ylabel('Error')
ax.set_xscale('log', basex=2)
ax.set_yscale('log', basey=2)
plt.legend()
plt.savefig("convergence.pdf")
