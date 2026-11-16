"""VLE model wrapper: modified Raoult's law  y_i P = x_i gamma_i Psat_i  (Eq. 10.5)."""
import numpy as np
from .vapour_pressure import psat
from .wilson import gamma_wilson
from .nrtl import gamma_nrtl

class VLEModel:
    def __init__(self, comps, model):
        self.comps, self.model = list(comps), model
        self._g = {"Wilson": gamma_wilson, "NRTL": gamma_nrtl}[model]
    def psat(self, T): return np.array([psat(c, T) for c in self.comps])
    def gamma(self, x, T): return self._g(x, T, self.comps)

def check_vle(m, T, P, x, y):
    x, y = np.asarray(x), np.asarray(y)
    resid = np.abs(y * P - x * m.gamma(x, T) * m.psat(T))
    return dict(sum_x=float(x.sum()), sum_y=float(y.sum()), max_resid_kPa=float(resid.max()))
