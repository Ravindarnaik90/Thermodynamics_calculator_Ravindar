"""Antoine equation: ln(Psat/kPa) = A - B/(t/degC + C). T given in K."""
import numpy as np
from .params import COMPONENTS

def psat(name, T):
    c = COMPONENTS[name]
    return float(np.exp(c["A"] - c["B"] / (T - 273.15 + c["C"])))

def tsat(name, P):
    """Inverse Antoine, returns K."""
    c = COMPONENTS[name]
    return c["B"] / (c["A"] - np.log(P)) - c["C"] + 273.15
