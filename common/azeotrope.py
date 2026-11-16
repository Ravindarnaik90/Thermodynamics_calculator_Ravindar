"""Binary azeotrope search: solve y1 - x1 = 0 (with VLE) after scanning for sign changes."""
import numpy as np
from scipy.optimize import brentq
from .bubble_dew import bubble_P, bubble_T

def _scan(g, n=600):
    xs = np.linspace(0.001, 0.999, n); v = np.array([g(x) for x in xs])
    roots = [brentq(g, xs[i], xs[i + 1], xtol=1e-14) for i in range(n - 1) if v[i] * v[i + 1] < 0]
    k = int(np.argmin(np.abs(v)))
    return roots, dict(min_abs_y1_minus_x1=float(abs(v[k])), at_x1=float(xs[k]))

def azeotrope_T(m, T):
    def g(x1):
        P, y = bubble_P(m, T, [x1, 1 - x1]); return y[0] - x1
    roots, diag = _scan(g)
    return [dict(x1=r, P=bubble_P(m, T, [r, 1 - r])[0]) for r in roots], diag

def azeotrope_P(m, P):
    def g(x1):
        T, y = bubble_T(m, P, [x1, 1 - x1]); return y[0] - x1
    roots, diag = _scan(g, 300)
    return [dict(x1=r, T=bubble_T(m, P, [r, 1 - r])[0]) for r in roots], diag
