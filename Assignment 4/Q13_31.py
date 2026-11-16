"""Q13_31 -- Problem 13.31 (liquid-phase isomerisation A -> B).   Run:  python Q13_31.py"""
import os, sys
import numpy as np
from scipy.optimize import brentq, minimize_scalar
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from common.report import Report

QID = "Q13_31"
SOURCE = "Smith, Van Ness & Abbott, 6th ed., Chapter 13 (liquid-phase reaction equilibrium); data as given in the problem"
STATEMENT = ("The isomerisation reaction A -> B occurs in the liquid phase. A and B are miscible in all proportions, and the "
             "excess Gibbs energy of the mixture is G^E/RT = 0.1 xA xB. If dG0(298 K) = -1000 J/mol, what is the equilibrium "
             "composition at 298.15 K? How much error is introduced if A and B are assumed to form an ideal solution?")

def main():
    rep = Report(QID, "Problem 13.31 - Liquid-phase isomerisation A -> B at 298.15 K", HERE, source=SOURCE)
    rep.h("Problem statement"); rep.p(STATEMENT)
    rep.h("Governing equations")
    rep.eq("K = exp(-dG0/RT)  with standard states of pure liquid A and B at T (and P);  nu = (-1, +1)")
    rep.eq("K = (x_B gamma_B) / (x_A gamma_A)")
    rep.eq("G^E/RT = 0.1 xA xB  =>  ln gamma_A = 0.1 xB^2,   ln gamma_B = 0.1 xA^2        (symmetric Margules, A12 = A21 = 0.1)")
    rep.eq("Ideal solution (gamma = 1):  K = xB/(1 - xB)  =>  xB = K/(1+K)")
    rep.h("Given data")
    R, T, dG0 = 8.314, 298.15, -1000.0
    rep.kv("T", T, "K"); rep.kv("dG0(298)", dG0, "J/mol"); rep.kv("R", R, "J/(mol K)")
    rep.h("Step 1 - equilibrium constant")
    K = float(np.exp(-dG0 / (R * T))); rep.kv("K = exp(1000/(R T))", K); rep.kv("ln K", np.log(K))
    rep.h("Step 2 - non-ideal solution: solve ln K = ln(xB/xA) + ln gamma_B - ln gamma_A for xB (brentq)")
    lgA = lambda xB: 0.1 * xB**2; lgB = lambda xB: 0.1 * (1 - xB)**2
    f = lambda xB: np.log(K) - (np.log(xB / (1 - xB)) + lgB(xB) - lgA(xB))
    xB = brentq(f, 1e-9, 1 - 1e-9, xtol=1e-15); xA = 1 - xB
    rep.kv("xB", xB); rep.kv("xA", xA); rep.kv("gamma_A, gamma_B", f"{np.exp(lgA(xB)):.5f}, {np.exp(lgB(xB)):.5f}")
    rep.h("Step 3 - ideal solution")
    xB_id = K / (1 + K); rep.kv("xB (ideal)", xB_id); rep.kv("xA (ideal)", 1 - xB_id)
    err = xB_id - xB
    rep.h("Final answers")
    st = "converged (brentq)"
    rep.result("-", "K", K, "-", "explicit")
    rep.result("-", "xB, non-ideal solution (G^E/RT = 0.1 xA xB)", xB, "-", st); rep.result("-", "xA, non-ideal solution", xA, "-", st)
    rep.result("-", "xB, ideal solution", xB_id, "-", "explicit"); rep.result("-", "xA, ideal solution", 1 - xB_id, "-", "explicit")
    rep.result("-", "error in xB from assuming ideal solution", err, "mole fraction", "explicit")
    rep.result("-", "relative error in xB", 100 * err / xB, "%", "explicit")
    rep.h("Verification")
    rep.check("equilibrium relation xB gB / (xA gA) = K", abs(xB * np.exp(lgB(xB)) / (xA * np.exp(lgA(xB))) - K) / K < 1e-12)
    rep.check("xA + xB = 1 and 0 < xB < 1", 0 < xB < 1)
    # activity coefficients from G^E/RT by differentiation: ln gamma_i = d(n G^E/RT)/dn_i
    nGE = lambda nA, nB: 0.1 * nA * nB / (nA + nB); h = 1e-6
    lgA_num = (nGE(xA + h, xB) - nGE(xA - h, xB)) / (2 * h); lgB_num = (nGE(xA, xB + h) - nGE(xA, xB - h)) / (2 * h)
    rep.check("ln gamma_i = d(nG^E/RT)/dn_i agrees with the closed forms", abs(lgA_num - lgA(xB)) < 1e-8 and abs(lgB_num - lgB(xB)) < 1e-8)
    Gmix = lambda x: x * dG0 / (R * T) + (1 - x) * np.log(1 - x) + x * np.log(x) + 0.1 * x * (1 - x)     # G/RT, reference pure A
    mn = minimize_scalar(Gmix, bounds=(1e-9, 1 - 1e-9), method="bounded", options={"xatol": 1e-12})
    rep.check("Gibbs-energy minimum of the mixture coincides with the K-equation root", abs(mn.x - xB) < 1e-6, f"(minimum at xB = {mn.x:.8f})")
    rep.check("Gibbs-energy minimum (ideal) coincides with K/(1+K)", abs(minimize_scalar(lambda x: x * dG0 / (R * T) + (1 - x) * np.log(1 - x) + x * np.log(x), bounds=(1e-9, 1 - 1e-9), method="bounded", options={"xatol": 1e-12}).x - xB_id) < 1e-6)
    # ---------------- graph
    x = np.linspace(1e-4, 1 - 1e-4, 800)
    G_id = x * dG0 / (R * T) + x * np.log(x) + (1 - x) * np.log(1 - x); G_re = G_id + 0.1 * x * (1 - x)
    fig, ax = plt.subplots(figsize=(7.2, 4.5), dpi=300)
    ax.plot(x, G_id, lw=2.2, color="#2e8b57", label="ideal solution"); ax.plot(x, G_re, lw=2.2, color="#1f4e79", label="$G^E/RT = 0.1\\,x_Ax_B$")
    ax.plot(xB_id, Gmix(xB_id) - 0.1 * xB_id * (1 - xB_id), "o", color="#2e8b57", ms=7, label=f"ideal minimum, $x_B$ = {xB_id:.4f}")
    ax.plot(xB, Gmix(xB), "o", color="#c0392b", ms=7, label=f"non-ideal minimum, $x_B$ = {xB:.4f}")
    ax.set(xlabel="$x_B$ [-]", ylabel="G / (RT)   (reference: pure A)", title="Problem 13.31: Gibbs energy of the A/B liquid mixture at 298.15 K")
    ax.title.set_fontsize(10); ax.grid(alpha=.35); ax.legend(); fig.tight_layout()
    fig.savefig(os.path.join(rep.dirs["plots"], QID + ".png"), dpi=300); plt.close(fig)
    rep.p(f"Graph saved as plots/{QID}/{QID}.png")
    ok = rep.save(); sys.exit(0 if ok else 1)

if __name__ == "__main__": main()
