"""Chemical-reaction equilibrium helpers (SVNA Chapter 13).
Ideal-gas heat capacity:  Cp/R = A + B T + C T^2 + D T^-2   (T in K)           [Table C.1]
Standard enthalpy / Gibbs energy of formation at 298.15 K, J/mol               [Table C.4]
Equilibrium constant from Eq. (13.18):
    dG/RT = (dG0 - dH0)/(R T0) + dH0/(R T) + (1/T) Int(dCp/R dT) - Int((dCp/R)/T dT)
"""
import numpy as np
from scipy.integrate import quad

R, T0, P0 = 8.314, 298.15, 1.0          # J/(mol K), K, bar (standard-state pressure)

# name: (A, B, C, D)
CP = {
    "H2":           (3.249, 0.422e-3,  0.0,        0.083e5),   # Table C.1
    "N2":           (3.280, 0.593e-3,  0.0,        0.040e5),   # Table C.1
    "NH3":          (3.578, 3.020e-3,  0.0,       -0.186e5),   # Table C.1
    "H2O":          (3.470, 1.450e-3,  0.0,        0.121e5),   # Table C.1 (gas)
    "H2S":          (3.931, 1.490e-3,  0.0,       -0.232e5),   # Table C.1
    "SO2":          (5.699, 0.801e-3,  0.0,       -1.015e5),   # Table C.1
    "styrene":      (2.050, 50.192e-3, -16.662e-6, 0.0),       # Table C.1
    "ethylbenzene": (1.124, 55.380e-3, -18.476e-6, 0.0),       # Table C.1
    "1-butene":     (1.967, 31.630e-3, -9.873e-6,  0.0),       # Table C.1
    "butadiene":    (2.734, 26.786e-3, -8.882e-6,  0.0),       # Table C.1 (1,3-butadiene)
    "S(s)":         (4.114, -1.728e-3, 0.0,       -0.783e5),   # Table C.2, rhombic sulfur (check Cp(298) = 22.6 J/mol/K)
}
HF = {"H2": 0.0, "N2": 0.0, "NH3": -46110.0, "H2O": -241818.0, "H2S": -20630.0, "SO2": -296830.0,
      "styrene": 147360.0, "ethylbenzene": 29920.0, "1-butene": -540.0, "butadiene": 109240.0, "S(s)": 0.0}
GF = {"H2": 0.0, "N2": 0.0, "NH3": -16450.0, "H2O": -228572.0, "H2S": -33560.0, "SO2": -300194.0,
      "styrene": 213900.0, "ethylbenzene": 130890.0, "1-butene": 70340.0, "butadiene": 149795.0, "S(s)": 0.0}

def delta_coeffs(nu):
    """(dA, dB, dC, dD) = sum nu_i (A_i, B_i, C_i, D_i)"""
    return tuple(sum(v * CP[s][k] for s, v in nu.items()) for k in range(4))

def dH298(nu): return sum(v * HF[s] for s, v in nu.items())
def dG298(nu): return sum(v * GF[s] for s, v in nu.items())

def J_int(nu, T):
    """Int_{T0}^{T} dCp/R dT  (analytic)"""
    A, B, C, D = delta_coeffs(nu)
    return A * (T - T0) + B / 2 * (T**2 - T0**2) + C / 3 * (T**3 - T0**3) - D * (1 / T - 1 / T0)

def I_int(nu, T):
    """Int_{T0}^{T} (dCp/R)/T dT  (analytic)"""
    A, B, C, D = delta_coeffs(nu)
    return A * np.log(T / T0) + B * (T - T0) + C / 2 * (T**2 - T0**2) - D / 2 * (1 / T**2 - 1 / T0**2)

def dH(nu, T):
    """Standard reaction enthalpy at T, J/mol"""
    return dH298(nu) + R * J_int(nu, T)

def lnK(nu, T):
    """ln K(T) from Eq. (13.18)"""
    dG_RT = (dG298(nu) - dH298(nu)) / (R * T0) + dH298(nu) / (R * T) + J_int(nu, T) / T - I_int(nu, T)
    return -dG_RT

def lnK_numeric(nu, T):
    """Same quantity with the two integrals done by numerical quadrature (independent check)."""
    f = lambda t: sum(v * (CP[s][0] + CP[s][1] * t + CP[s][2] * t**2 + CP[s][3] / t**2) for s, v in nu.items())
    J = quad(f, T0, T, epsabs=1e-13, epsrel=1e-13)[0]; I = quad(lambda t: f(t) / t, T0, T, epsabs=1e-13, epsrel=1e-13)[0]
    return -((dG298(nu) - dH298(nu)) / (R * T0) + dH298(nu) / (R * T) + J / T - I)

def describe_reaction(rep, nu, T, label="T"):
    """Print the thermodynamic data and intermediate quantities for K(T)."""
    rep.table(["species", "nu_i", "A", "10^3 B", "10^6 C", "10^-5 D", "dHf298 [J/mol]", "dGf298 [J/mol]"],
              [[s, v, CP[s][0], CP[s][1] * 1e3, CP[s][2] * 1e6, CP[s][3] / 1e5, HF[s], GF[s]] for s, v in nu.items()])
    A, B, C, D = delta_coeffs(nu)
    rep.kv("dA, dB, dC, dD", f"{A:.5g}, {B:.5g}, {C:.5g}, {D:.5g}")
    rep.kv("dH298 = sum nu_i dHf_i", dH298(nu) / 1e3, "kJ/mol"); rep.kv("dG298 = sum nu_i dGf_i", dG298(nu) / 1e3, "kJ/mol")
    rep.kv("K(298.15 K) = exp(-dG298/RT0)", float(np.exp(-dG298(nu) / (R * T0))))
    rep.kv(f"dH({T} K)", dH(nu, T) / 1e3, "kJ/mol"); rep.kv(f"ln K({T} K)", lnK(nu, T)); rep.kv(f"K({T} K)", float(np.exp(lnK(nu, T))))

def thermo_checks(rep, nu, Ts):
    """Independent checks of the K(T) machinery."""
    rep.check("ln K(298.15) = -dG298/(R T0)", abs(lnK(nu, T0) - (-dG298(nu) / (R * T0))) < 1e-12)
    worst = max(abs(lnK(nu, T) - lnK_numeric(nu, T)) for T in Ts)
    rep.check("analytic integrals agree with numerical quadrature", worst < 1e-9, f"(max |diff in ln K| = {worst:.1e})")
    wv = 0
    for T in Ts:
        h = 1e-3; wv = max(wv, abs((lnK(nu, T + h) - lnK(nu, T - h)) / (2 * h) - dH(nu, T) / (R * T**2)))
    rep.check("van't Hoff relation d ln K/dT = dH/(R T^2)", wv < 1e-8, f"(max abs deviation {wv:.1e})")
