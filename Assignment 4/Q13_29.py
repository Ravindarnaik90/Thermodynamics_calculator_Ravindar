"""Q13_29 -- Problem 13.29 (simultaneous removal of H2S and SO2).   Run:  python Q13_29.py"""
import os, sys
import numpy as np
from scipy.optimize import brentq, fsolve
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from common.report import Report
from common.reaction import R, T0, P0, CP, dH298, dG298, dH, lnK, describe_reaction, thermo_checks

QID = "Q13_29"
SOURCE = "Smith, Van Ness & Abbott, 6th ed.: Table C.1 (gas heat capacities), Table C.2 (solid sulfur), Table C.4 (formation properties)"
STATEMENT = ("Refineries often have both H2S and SO2 to dispose of. The reaction 2 H2S(g) + SO2(g) -> 3 S(s) + 2 H2O(g) suggests a "
             "means of getting rid of both. For reactants in the stoichiometric proportion, estimate the percent conversion of "
             "each reactant if the reaction comes to equilibrium at 723.15 K (450 degC) and 8 bar.")

def main():
    rep = Report(QID, "Problem 13.29 - Simultaneous removal of H2S and SO2, 723.15 K, 8 bar", HERE, source=SOURCE)
    rep.h("Problem statement"); rep.p(STATEMENT)
    rep.h("Governing equations")
    rep.eq("Reaction: 2 H2S + SO2 -> 3 S(s) + 2 H2O(g),  nu = (-2, -1, +3, +2)")
    rep.eq("K = exp(-dG0(T)/RT) from Eq. (13.18)")
    rep.eq("Solid sulfur is a pure phase (activity = 1), so it does not appear in K:")
    rep.eq("K = y_H2O^2 / ( y_H2S^2 y_SO2 (P/P0) )          (ideal gas, gas-phase species only)")
    rep.h("Given data")
    T, P = 723.15, 8.0
    rep.kv("T", T, "K"); rep.kv("P", P, "bar"); rep.kv("feed basis", "2 mol H2S + 1 mol SO2 (stoichiometric)")
    nu = {"H2S": -2, "SO2": -1, "S(s)": 3, "H2O": 2}
    rep.h("Step 1 - equilibrium constant at T")
    rep.kv("Cp(S, solid) at 298.15 K = R (A + B T + D/T^2)", R * (CP["S(s)"][0] + CP["S(s)"][1] * T0 + CP["S(s)"][3] / T0**2), "J/(mol K)  (literature: ~22.7)")
    describe_reaction(rep, nu, T); K = float(np.exp(lnK(nu, T)))
    rep.h("Step 2 - mole balance (basis 2 mol H2S + 1 mol SO2, extent e)")
    rep.eq("gas phase: n_H2S = 2 - 2e,  n_SO2 = 1 - e,  n_H2O = 2e,   N_gas = 3 - e;   S(s) = 3e")
    rep.eq("K = [2e/N]^2 / ( [(2-2e)/N]^2 [(1-e)/N] P ),  conversion of H2S = 2e/2 = e,  of SO2 = e/1 = e")
    def resid(ee, TT=T):
        nH2S, nSO2, nW = 2 - 2 * ee, 1 - ee, 2 * ee; Ng = nH2S + nSO2 + nW
        return lnK(nu, TT) - np.log((nW / Ng)**2 / ((nH2S / Ng)**2 * (nSO2 / Ng) * P / P0))
    rep.h("Step 3 - solve for e (scipy brentq, bracket 1e-12 < e < 1-1e-12)")
    e = brentq(resid, 1e-12, 1 - 1e-12, xtol=1e-15)
    nH2S, nSO2, nW = 2 - 2 * e, 1 - e, 2 * e; Ng = nH2S + nSO2 + nW
    y = np.array([nH2S, nSO2, nW]) / Ng
    rep.kv("extent e", e, "mol"); rep.kv("gas moles N", Ng, "mol"); rep.kv("solid S formed", 3 * e, "mol")
    rep.kv("y_H2S, y_SO2, y_H2O", ", ".join(f"{v:.5f}" for v in y))
    rep.h("Final answers")
    st = "converged (brentq)"
    rep.result("-", "K at 723.15 K", K, "-", "explicit")
    rep.result("-", "conversion of H2S", 100 * e, "%", st); rep.result("-", "conversion of SO2", 100 * e, "%", st)
    rep.h("Verification")
    thermo_checks(rep, nu, [450.0, 723.15, 900.0])
    rep.check("sum of gas mole fractions = 1", abs(y.sum() - 1) < 1e-12)
    rep.check("equilibrium relation y_W^2/(y_H2S^2 y_SO2 P) = K", abs(y[2]**2 / (y[0]**2 * y[1] * P) - K) / K < 1e-10)
    rep.check("S balance: S in = S out", abs((2 + 1) - (nH2S + nSO2 + 3 * e)) < 1e-12)
    rep.check("O balance: O in = O out", abs((2 * 1) - (2 * nSO2 + nW)) < 1e-12)
    rep.check("H balance: H in = H out", abs((2 * 2) - (2 * nH2S + 2 * nW)) < 1e-12)
    sol = fsolve(lambda u: resid(u[0]), [0.5], xtol=1e-12, full_output=True)
    rep.check("independent fsolve gives the same extent", abs(resid(sol[0][0])) < 1e-9 and abs(sol[0][0] - e) < 1e-9, f"(e = {sol[0][0]:.10f}, residual {abs(resid(sol[0][0])):.1e})")
    rep.check("0 <= e <= 1", 0 <= e <= 1)
    rep.note("Elemental sulfur is treated as a pure solid, as stated in the problem. At 723 K real sulfur is a liquid/vapour; the solid heat capacity (Table C.2) is therefore extrapolated beyond its stated range.")
    # ---------------- graph
    Ts = np.linspace(450, 900, 120)
    cv = [100 * brentq(lambda ee: resid(ee, t), 1e-14, 1 - 1e-14) for t in Ts]
    fig, ax = plt.subplots(figsize=(7.2, 4.5), dpi=300)
    ax.plot(Ts - 273.15, cv, lw=2.2, color="#1f4e79", label="equilibrium conversion, 8 bar")
    ax.plot(T - 273.15, 100 * e, "o", color="#c0392b", ms=8, label=f"450 degC: {100*e:.2f} %")
    ax.set(xlabel="Temperature, t [°C]", ylabel="Equilibrium conversion of H$_2$S and SO$_2$ [%]", title="Problem 13.29: 2 H$_2$S + SO$_2$ $\\rightarrow$ 3 S(s) + 2 H$_2$O at 8 bar")
    ax.title.set_fontsize(10); ax.grid(alpha=.35); ax.legend(); fig.tight_layout()
    fig.savefig(os.path.join(rep.dirs["plots"], QID + ".png"), dpi=300); plt.close(fig)
    rep.p(f"Graph saved as plots/{QID}/{QID}.png")
    ok = rep.save(); sys.exit(0 if ok else 1)

if __name__ == "__main__": main()
