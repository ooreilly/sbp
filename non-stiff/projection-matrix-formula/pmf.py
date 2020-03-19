import sympy as sp

def projection_matrix_formula(lm, lp, Xp, Xm, R, dR):
    L = Xm * lm - Xp * lp * R.T
    dL = Xm * (lm ** - 1) + Xp * lp ** - 1 * dR
    #V = (sp.simplify(L.T * dL)) ** -1
    Im = sp.eye(Xm.shape[1])
    V = (Im - R*dR) ** -1
    P = dL * V * L.T
    return P

