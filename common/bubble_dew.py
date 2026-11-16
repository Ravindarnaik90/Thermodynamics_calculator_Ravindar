"""BUBL P, DEW P, BUBL T, DEW T from Eq. (10.5)."""
import numpy as np
from scipy.optimize import brentq

def bubble_P(m, T, x):
    x = np.asarray(x, float); t = x * m.gamma(x, T) * m.psat(T)
    P = t.sum(); return P, t / P

def dew_P(m, T, y, tol=1e-13, maxit=5000):
    """x_i = y_i P/(gamma_i Psat_i); P = 1/sum(y_i/(gamma_i Psat_i)); damped successive substitution."""
    y = np.asarray(y, float); ps = m.psat(T); x = y.copy()
    for _ in range(maxit):
        g = m.gamma(x, T); P = 1.0 / np.sum(y / (g * ps))
        xn = y * P / (g * ps); xn /= xn.sum()
        if np.max(np.abs(xn - x)) < tol:
            g = m.gamma(xn, T); P = 1.0 / np.sum(y / (g * ps))
            return P, xn, True
        x = 0.3 * x + 0.7 * xn
    return P, x, False

def bubble_T(m, P, x, lo=250.0, hi=430.0):
    x = np.asarray(x, float)
    f = lambda T: np.sum(x * m.gamma(x, T) * m.psat(T)) - P
    T = brentq(f, lo, hi, xtol=1e-12)
    t = x * m.gamma(x, T) * m.psat(T)
    return T, t / t.sum()

def dew_T(m, P, y, lo=250.0, hi=430.0):
    y = np.asarray(y, float)
    f = lambda T: dew_P(m, T, y)[0] - P
    T = brentq(f, lo, hi, xtol=1e-10)
    _, x, conv = dew_P(m, T, y)
    return T, x
