import numpy as np

class Struct(dict):
    """
    Make a dict behave as a struct.

    Example:
    
        test = Struct(a=1, b=2, c=3)

    """
    def __init__(self,**kw):
        dict.__init__(self,kw)
        self.__dict__ = self

def load_csv(path, name):
    d = np.genfromtxt('%s/%s'%(path, name), delimiter=' ')
    data = Struct()
    data.t = d[:,0]
    data.d = d
    return data

def load_config(path, name):
    d = np.genfromtxt('%s/%s'%(path, name), delimiter='=', dtype=[('string',
    'S8'), ('float', 'f8')])
    f = open("%s/%s"%(path, name))
    d = Struct()
    for line in f.readlines():
        key, value = line.strip('\n').split("=")
        d[key] = value
    if 'm' in d:
        d.m = int(d.m)
    if 'n' in d:
        d.n = int(d.n)
    if 'nnz' in d:
        d.nnz = int(d.nnz)
    if 'nt' in d:
        d.nt = int(d.nt)
    if 'elapsed' in d:
        d.elapsed = float(d.elapsed)
        d.avg = 1e3 * d.elapsed / d.nt
    if 'nx' in d:
        d.nx = int(d.nx)
    if 'ny' in d:
        d.ny = int(d.ny)
    return d

def norm_error(err, h):
    return [np.linalg.norm(h*erri) for erri in err]

def convergence(err):
    import numpy as np
    q = []
    for i, e in enumerate(err): 
        if i == len(err) - 1:
            break
        q.append(np.log2(err[i]/err[i+1]))
    return q

def rate(e0, e1, n0, n1, nr, p):
    """
    This equation is solved for p to determine the convergence rate when the
    reference solution is used to measure the error.

    y(p) = 0 determines the convergence rate.

    e0, n0 : Error and grid number for grid "0"
    e1, n1 : Error and grid number for grid "1"
    nr : Reference grid number 
    p : Convergence rate to solve for

    """
    h0 = 0.5 ** (n0 * p)
    h1 = 0.5 ** (n1 * p)
    hr = 0.5 ** (nr * p)
    y = e0 * (h1 - hr) - e1 * (h0 - hr)
    return y

def scaling_constant(n0, n1, nr, p):
    """
    Determine the scaling constant.

    """

    h0 = 0.5 ** (n0 * p)
    h1 = 0.5 ** (n1 * p)
    hr = 0.5 ** (nr * p)
    y = (h0 - hr) / (h1 - hr)
    return y


def adjusted_rate(err, nr, pmax=6):
    from scipy import optimize
 

    import numpy as np
    q = []
    for i, e in enumerate(err): 
        if i == len(err) - 1:
            break

        e0 = err[i]
        e1 = err[i+1]
        n0 = i
        n1 = i + 1
        f = lambda p : rate(e0, e1, n0, n1, nr, p)
        sol = optimize.brentq(f, 1, pmax)
        q.append(sol)
    return q


def error(data, field):
    n = len(data)
    ref = data[-1].d[:,field]
    err = []
    for di in data[:-1]:
        err.append(di.d[:,field] - ref)
    return err

def compute_errors(data, field):
    out = Struct()
    out.err = error(data, field)
    out.field = []
    ngrid = len(data)
    for num in range(ngrid):
        out.field.append(data[num].d[:,field])
    out.t = data[0].t
    h = data[0].t[1]
    out.norm_err = norm_error(out.err, h)
    ref = data[-1].d[:,field]
    out.norm_rel_err = np.array(out.norm_err)/norm_error([ref], h)[0]
    #out.conv = adjusted_rate(out.norm_err, ngrid)
    return out

def error_table(p, u1, u2, pts_wl=2.5):

    fmt = '%.2e'
    pe = [fmt%ei for ei in p.norm_rel_err]
    pq = ['-'] + ['%3.2f'%q for q in p.conv]

    u1e = [fmt%ei for ei in u1.norm_rel_err]
    u1q = ['-'] + ['%3.2f'%q for q in u1.conv]
    u2e = [fmt%ei for ei in u2.norm_rel_err]
    u2q = ['-'] + ['%3.2f'%q for q in u2.conv]
    out = ''
    N = [pts_wl * 2 ** i for i in range(len(pe))]
    for Ni, pi, qi, u1i, u1qi, u2i, u2qi  in zip(N, pe, pq, u1e, u1q, u2e, u2q):
        out += '%d & %s & %s & %s & %s & %s & %s \\\ \n' % (Ni, pi, qi, u1i, u1qi, u2i, u2qi)
    return out
