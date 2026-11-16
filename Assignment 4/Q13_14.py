"""Q13_14 -- Problem 13.14 (styrene hydrogenation).   Run:  python Q13_14.py"""
import os, sys
import numpy as np
from scipy.optimize import brentq
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from common.report import Report
from common.reaction import R, T0, P0, dH298, dG298, dH, lnK, describe_reaction, thermo_checks

QID = "Q13_14"
SOURCE = "Smith, Van Ness & Abbott, 6th ed.: Table C.1 (ideal-gas heat capacities), Table C.4 (formation properties)"
STATEMENT = ("The reaction C6H5CH:CH2(g) + H2(g) -> C6H5.C2H5(g) reaches equilibrium at 923.15 K (650 degC) and "
             "atmospheric pressure. If the system initially contains 1.5 mol H2 for each mole of styrene, what is the "
             "composition of the system at equilibrium? Assume ideal gases.")

def main():
    rep = Report(QID, "Problem 13.14 - Styrene hydrogenation at 923.15 K and 1 atm", HERE, source=SOURCE)
    rep.h("Problem statement"); rep.p(STATEMENT)
    rep.h("Governing equations")
    rep.eq("Reaction: styrene + H2 -> ethylbenzene,  nu = (-1, -1, +1)")
    rep.eq("K = exp(-dG0(T)/RT),  dG0/RT = (dG298-dH298)/(R T0) + dH298/(R T) + (1/T) Int(dCp/R dT) - Int((dCp/R)/T dT)   [Eq. 13.18]")
    rep.eq("Ideal gas:  K = prod (y_i P/P0)^nu_i = y_EB / (y_S y_H2 P/P0),   P0 = 1 bar")
    rep.h("Given data")
    T = 923.15; P = 1.01325; n_S, n_H2 = 1.0, 1.5          # K, bar (1 atm), mol
    rep.kv("T", T, "K"); rep.kv("P", P, "bar (1 atm)"); rep.kv("feed basis", "1 mol styrene + 1.5 mol H2")
    nu = {"styrene": -1, "H2": -1, "ethylbenzene": 1}
    rep.h("Step 1 - equilibrium constant at T")
    describe_reaction(rep, nu, T); K = float(np.exp(lnK(nu, T)))
    rep.h("Step 2 - mole balance with reaction extent e (basis 1 mol styrene)")
    rep.eq("n_S = 1 - e,  n_H2 = 1.5 - e,  n_EB = e,  N = 2.5 - e")
    rep.eq("K = [e/N] / ( [(1-e)/N] [(1.5-e)/N] P )  =>  K P (1-e)(1.5-e) = e (2.5-e)")
    comp = lambda e: (np.array([1 - e, 1.5 - e, e]), 2.5 - e)
    resid = lambda e: lnK(nu, T) - np.log((e / (2.5 - e)) / (((1 - e) / (2.5 - e)) * ((1.5 - e) / (2.5 - e)) * P / P0))
    rep.h("Step 3 - solve for e (scipy brentq on ln K - ln Q = 0, bracket 1e-9 < e < 1-1e-9)")
    e = brentq(resid, 1e-9, 1 - 1e-9, xtol=1e-14)
    n, N = comp(e); y = n / N
    rep.kv("extent e", e, "mol"); rep.kv("total moles N", N, "mol")
    for s_, ni, yi in zip(["styrene", "hydrogen", "ethylbenzene"], n, y): rep.kv(f"n_{s_} / y_{s_}", f"{ni:.5f} mol / {yi:.5f}")
    rep.h("Final answers")
    st = "converged (brentq)"
    rep.result("-", "K at 923.15 K", K, "-", "explicit")
    rep.result("-", "y styrene", y[0], "-", st); rep.result("-", "y hydrogen", y[1], "-", st); rep.result("-", "y ethylbenzene", y[2], "-", st)
    rep.result("-", "styrene conversion", 100 * e, "%", st)
    rep.h("Verification")
    thermo_checks(rep, nu, [500.0, 700.0, 923.15, 1100.0])
    rep.check("sum of mole fractions = 1", abs(y.sum() - 1) < 1e-12)
    rep.check("equilibrium relation y_EB/(y_S y_H2 P) = K", abs(y[2] / (y[0] * y[1] * P) - K) / K < 1e-10)
    q = K * P; a = 1 + q; e_exact = (2.5 * a - np.sqrt(6.25 * a**2 - 6 * q * a)) / (2 * a)
    rep.check("independent closed-form quadratic root agrees", abs(e - e_exact) < 1e-10, f"(e = {e_exact:.10f})")
    rep.check("0 <= e <= 1 (limiting reactant styrene)", 0 <= e <= 1)
    # ---------------- graph: equilibrium conversion vs temperature
    Ts = np.linspace(600, 1100, 101); conv = []
    for t in Ts:
        r_ = lambda ee: lnK(nu, t) - np.log((ee / (2.5 - ee)) / (((1 - ee) / (2.5 - ee)) * ((1.5 - ee) / (2.5 - ee)) * P / P0))
        conv.append(100 * brentq(r_, 1e-12, 1 - 1e-12))
    fig, ax = plt.subplots(figsize=(7.2, 4.5), dpi=300)
    ax.plot(Ts - 273.15, conv, color="#1f4e79", lw=2.2, label="equilibrium conversion, 1 atm")
    ax.plot(T - 273.15, 100 * e, "o", color="#c0392b", ms=8, label=f"650 degC: {100*e:.1f} %")
    ax.set(xlabel="Temperature, t [°C]", ylabel="Equilibrium styrene conversion [%]", title="Problem 13.14: styrene + H$_2$ $\\rightarrow$ ethylbenzene, 1 atm, H$_2$/styrene = 1.5")
    ax.title.set_fontsize(10); ax.grid(alpha=.35); ax.legend(); fig.tight_layout()
    png = os.path.join(rep.dirs["plots"], QID + ".png"); fig.savefig(png, dpi=300); plt.close(fig)
    rep.p(f"Graph saved as plots/{QID}/{QID}.png")
    ok = rep.save(); sys.exit(0 if ok else 1)

if __name__ == "__main__": main()
