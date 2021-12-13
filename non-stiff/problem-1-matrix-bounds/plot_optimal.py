import hdf5storage
import matplotlib.pyplot as plt
import sys

filename = sys.argv[1]
show= int(sys.argv[2])
out = sys.argv[3]

spec = hdf5storage.loadmat(filename)

dr = spec["dR"]
Mnorm = spec["Mnorm"]
rd = spec["rd"]
rl = spec["rl"]
sd = spec["sd"]
sl = spec["sl"]
lower = rl * abs(spec["sigma_norm_lower"])
upper = sd + rl * spec["sigma_norm_upper"]
plt.plot(dr, Mnorm, 'k-', label="$ h || M_h ||_2$")
plt.plot(dr, lower, 'b--', label="LB")
plt.plot(dr, upper, 'r--', label="UB")
plt.xlabel("$\\delta R$")
plt.legend()
plt.savefig(out)
if show:
    plt.show()
