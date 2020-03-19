import hdf5storage
import matplotlib.pyplot as plt
import sys

filename = sys.argv[1]
show= int(sys.argv[2])
out = sys.argv[3]

spec = hdf5storage.loadmat(filename)

alpha = spec["alpha"]
Mnorm = spec["Mnorm"]
rd = spec["rd"]
rl = spec["rl"]
sd = spec["sd"]
sl = spec["sl"]
lower = rl * abs(spec["sigma_norm_lower"])
upper = sd + rl * spec["sigma_norm_upper"]
plt.plot(alpha, Mnorm, 'k-', label="$ h || M_h ||_2$")
plt.plot(alpha, lower, 'b--', label="LB")
plt.plot(alpha, upper, 'r--', label="UB")
#plt.ylabel("$h \|M\|_2 / c$")
plt.xlabel("$\\alpha$")
plt.legend()
plt.savefig(out)
if show:
    plt.show()
