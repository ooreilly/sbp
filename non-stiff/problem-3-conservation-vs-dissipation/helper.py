import numpy as np

class Struct(dict):
  def __init__(self,**kw):
        dict.__init__(self,kw)
        self.__dict__ = self

def load(name, ref, prec=np.float32):
        p0 = np.fromfile("output/%s_%d_p0.bin" % (name, ref), dtype=prec)
        p1 = np.fromfile("output/%s_%d_p1.bin" % (name, ref), dtype=prec)
        t = np.fromfile("output/%s_%d_t.bin" % (name, ref), dtype=prec)
        return p0, p1, t

def parse_log(name, ref):
    import re
    lines = open("%s_%d.txt" % (name ,ref), "r").readlines()
    pattern = "Simulation took:\s([\d\.s]+)"
    for line in lines:
        match = re.findall(pattern, line)
        if match:
            return float(match[0])


def convergence_rates(errors):
    """

    Compute convergence rates assuming factor of two refinement and that the
    error on the finest grid should be discarded.

    """
    return ["-"] +  list(np.log2(np.array(errors[:-2]) /
        np.array(errors[1:-1])))

def error(u, v, dt):
    """
    l2 error

    """
    return np.linalg.norm(u - v) * dt

def normalized_error(u, v, dt):
    return error(u, v, dt) / np.linalg.norm(u) / dt

def latex_table(cols):
    strout = ""
    row_vals = []
    for row in range(len(cols[0])):
        vals = []
        for col in cols:
            vals.append(texformat(col[row]))
        row_vals.append(" & ".join(vals))
    return " \\\ \n".join(row_vals) + " \\\ \n"

def texformat(val):
    if isinstance(val, float):
        exp = int(np.round(np.log10(val)))
        man = val / 10 ** exp
        if man < 1:
            man = man * 10
            exp = exp - 1
        if exp == 0:
            return "$%2.3f$" % val
        else:
            return "$%2.3f \\times 10^{%d}$" % (man, exp)
    else:
        return val





    return strout

