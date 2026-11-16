"""(P,T) flash: Rachford-Rice + gamma successive substitution.
Single-phase cases detected from P_bubble / P_dew at the flash T."""
import numpy as np
from scipy.optimize import brentq
from .bubble_dew import bubble_P, dew_P

def flash_PT(m, P, T, z, tol=1e-13, maxit=2000):
    z = np.asarray(z, float)
    Pb, _ = bubble_P(m, T, z); Pd, _, _ = dew_P(m, T, z)
    out = dict(P_bubble_at_T=Pb, P_dew_at_T=Pd)
    if P >= Pb:
        return dict(out, phase="subcooled liquid", V=0.0, x=z, y=None, converged=True, iters=0)
    if P <= Pd:
        return dict(out, phase="superheated vapour", V=1.0, x=None, y=z, converged=True, iters=0)
    x = z.copy(); ps = m.psat(T)
    for it in range(maxit):
        K = m.gamma(x, T) * ps / P
        f = lambda V: np.sum(z * (K - 1) / (1 + V * (K - 1)))
        if not (K.max() > 1 > K.min()):
            return dict(out, phase="two-phase (FAILED: all K on one side of 1)", V=np.nan, x=x, y=K * x, converged=False, iters=it)
        # negative-flash bracket: RR is monotone between its poles (V may leave [0,1] while gamma iterates)
        lo = 1.0 / (1 - K.max()) + 1e-10; hi = 1.0 / (1 - K.min()) - 1e-10
        V = brentq(f, lo, hi, xtol=1e-15)
        xn = z / (1 + V * (K - 1)); xn /= xn.sum()
        if np.max(np.abs(xn - x)) < tol:
            x = xn; y = K * x; y /= y.sum()
            ok = 0 <= V <= 1
            return dict(out, phase="two-phase" if ok else "two-phase (FAILED: V outside [0,1])", V=V, x=x, y=y, converged=ok, iters=it + 1)
        x = xn
    return dict(out, phase="two-phase (not converged)", V=V, x=x, y=K * x, converged=False, iters=maxit)
