"""Independent cross-checks with scipy.optimize.fsolve (a different algorithm from the main solvers).
Several starting guesses (ideal-Raoult based, then a grid) are tried; only CONVERGED, PHYSICAL solutions
(x_i > 0, sum x = 1) are accepted. If none converges the function returns ok=False (reported as a failed check)."""
import numpy as np
from scipy.optimize import fsolve, brentq, least_squares
from .vapour_pressure import tsat

def _fs(f, u0):
    with np.errstate(all="ignore"):
        s, info, ier, _ = fsolve(f, u0, xtol=1e-13, full_output=True)
        res = float(np.max(np.abs(f(s))))
    return s, (ier == 1 and np.isfinite(res) and res < 1e-8)   # never trust the solver flag alone: check the residual

def _starts_x(n):
    g = np.linspace(0.05, 0.9, 5); out = []
    for a in g:
        for b in g:
            x = np.array([a, b, 1 - a - b]) if n == 3 else np.array([a, 1 - a])
            if np.all(x > 0.01): out.append(x)
            if n == 2: break
    return out

def _raoult_dew(m, T, y):
    ps = m.psat(T); P = 1 / np.sum(y / ps); x = y * P / ps; return P, x / x.sum()

def xc_dew_P(m, T, y, **_):
    n = len(y); y = np.asarray(y); ps = m.psat(T); P0, xr = _raoult_dew(m, T, y)
    f = lambda u: [u[n] * y[i] - u[i] * m.gamma(u[:n], T)[i] * ps[i] for i in range(n)] + [u[:n].sum() - 1]
    for x0 in [xr] + _starts_x(n):
        s, ok = _fs(f, list(x0) + [P0])
        if ok and np.all(s[:n] > 0) and abs(s[:n].sum() - 1) < 1e-9 and s[n] > 0: return s[n], s[:n], True
    return float("nan"), np.full(n, np.nan), False

def xc_bubble_T(m, P, x, **_):
    x = np.asarray(x); Tt = np.array([tsat(c, P) for c in m.comps])
    f = lambda T: np.sum(x * m.gamma(x, T[0]) * m.psat(T[0])) - P
    for T0 in (float(np.sum(x * Tt)), float(Tt.min()), float(Tt.max())):
        s, ok = _fs(f, [T0])
        if ok: return s[0], True
    return float("nan"), False

def xc_dew_T(m, P, y, **_):
    n = len(y); y = np.asarray(y); Tt = np.array([tsat(c, P) for c in m.comps])
    def f(u):
        x, T = u[:n], u[n]; ps = m.psat(T); g = m.gamma(x, T)
        return [P * y[i] - x[i] * g[i] * ps[i] for i in range(n)] + [x.sum() - 1]
    for T0 in (float(np.sum(y * Tt)), float(Tt.max()), float(Tt.min()) + 10):
        _, xr = _raoult_dew(m, T0, y)
        for x0 in [xr] + _starts_x(n):
            s, ok = _fs(f, list(x0) + [T0])
            if ok and np.all(s[:n] > 0) and abs(s[:n].sum() - 1) < 1e-9 and 250 < s[n] < 450: return s[n], s[:n], True
    return float("nan"), np.full(n, np.nan), False

def xc_flash(m, P, T, z, **_):
    """unknowns V, x_i. Equations: x_i (1 + V (K_i - 1)) = z_i (n eqs) and sum (K_i - 1) x_i = 0 (i.e. sum y = sum x = 1;
    this closure removes the trivial root V = 0, x = z). K_i = gamma_i Psat_i / P. Accepts only 0<=V<=1, x_i>0."""
    n = len(z); z = np.asarray(z); ps = m.psat(T)
    def f(u):
        V, x = u[0], u[1:]; K = m.gamma(x, T) * ps / P
        return [x[i] * (1 + V * (K[i] - 1)) - z[i] for i in range(n)] + [np.sum((K - 1) * x)]
    K0 = ps / P; starts = []
    if K0.max() > 1 > K0.min():                       # ideal-Raoult Rachford-Rice start
        rr = lambda V: np.sum(z * (K0 - 1) / (1 + V * (K0 - 1)))
        V0 = brentq(rr, 1 / (1 - K0.max()) + 1e-9, 1 / (1 - K0.min()) - 1e-9); starts.append((V0, z / (1 + V0 * (K0 - 1))))
    starts += [(V0, z) for V0 in (0.5, 0.2, 0.8)]
    for V0, x0 in starts:
        s, ok = _fs(f, [V0] + list(x0))
        if ok and 0 <= s[0] <= 1 and np.all(s[1:] > 0) and abs(s[1:].sum() - 1) < 1e-9: return s[0], s[1:], True
    rng = np.random.default_rng(12345)             # seeded random multi-start (reproducible)
    for _ in range(300):
        x0 = rng.dirichlet(np.ones(n)); V0 = rng.uniform(0.05, 0.95); s, ok = _fs(f, [V0] + list(x0))
        if ok and 0 <= s[0] <= 1 and np.all(s[1:] > 0) and abs(s[1:].sum() - 1) < 1e-9: return s[0], s[1:], True
    # second independent algorithm: bounded trust-region least squares, grid of starting vapour fractions
    for V0 in (0.1, 0.3, 0.5, 0.7, 0.9):
        with np.errstate(all="ignore"):
            r = least_squares(f, [V0] + list(z), bounds=([0] + [1e-9] * n, [1] + [1] * n), xtol=1e-15, ftol=1e-15, gtol=1e-15)
        if np.max(np.abs(r.fun)) < 1e-8 and abs(r.x[1:].sum() - 1) < 1e-9: return r.x[0], r.x[1:], True
    return float("nan"), np.full(n, np.nan), False
