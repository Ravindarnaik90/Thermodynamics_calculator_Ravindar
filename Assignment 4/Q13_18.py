"""Q13_18 -- Problem 13.18 (1,3-butadiene from 1-butene with steam).   Run:  python Q13_18.py"""
import os, sys
import numpy as np
from scipy.optimize import brentq
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from common.report import Report
from common.reaction import R, T0, P0, dH298, dG298, dH, lnK, describe_reaction, thermo_checks

QID = "Q13_18"
SOURCE = "Smith, Van Ness & Abbott, 6th ed.: Table C.1 (ideal-gas heat capacities), Table C.4 (formation properties)"
STATEMENT = ("1,3-Butadiene is produced by the dehydrogenation of 1-butene: C2H5CH:CH2(g) -> CH2:CHCH:CH2(g) + H2(g). "
             "Side reactions are suppressed by steam. If equilibrium is attained at 950 K and 1 bar and the reactor product "
             "contains 10 mol-% 1,3-butadiene, find (a) the mole fractions of the other species in the product gas, "
             "(b) the mole fraction of steam required in the feed.")

def main():
    rep = Report(QID, "Problem 13.18 - 1,3-Butadiene from 1-butene with steam, 950 K, 1 bar", HERE, source=SOURCE)
    rep.h("Problem statement"); rep.p(STATEMENT)
    rep.h("Governing equations")
    rep.eq("Reaction: 1-butene -> 1,3-butadiene + H2,  nu = (-1, +1, +1);  steam is inert")
    rep.eq("K = exp(-dG0(T)/RT) from Eq. (13.18);   ideal gas:  K = y_B y_H2 (P/P0) / y_C4")
    rep.h("Given data")
    T, P, yB = 950.0, 1.0, 0.10
    rep.kv("T", T, "K"); rep.kv("P", P, "bar"); rep.kv("y_butadiene in product", yB)
    nu = {"1-butene": -1, "butadiene": 1, "H2": 1}
    rep.h("Step 1 - equilibrium constant at 950 K")
    describe_reaction(rep, nu, T); K = float(np.exp(lnK(nu, T)))
    rep.h("Step 2 - part (a): product composition")
    rep.eq("Steam is inert, so butadiene and H2 are formed in equal amounts:  y_H2 = y_B = 0.10")
    rep.eq("K = y_B y_H2 P / y_C4   =>   y_C4 = y_B^2 P / K,      y_steam = 1 - y_B - y_H2 - y_C4")
    yH2 = yB; yC4 = yB * yH2 * (P / P0) / K; ys = 1 - yB - yH2 - yC4
    rep.kv("y_1-butene", yC4); rep.kv("y_steam", ys)
    rep.h("Step 3 - part (b): steam in the feed (basis 1 mol 1-butene fed, e mol reacted, s mol steam)")
    rep.eq("Product: butene (1-e), butadiene e, H2 e, steam s;  N = 1 + s + e")
    rep.eq("y_B / y_C4 = e / (1-e)  =>  e = y_B/(y_B + y_C4);   N = e / y_B;   s = N - 1 - e;   feed fraction = s/(1+s)")
    e = yB / (yB + yC4); N = e / yB; s = N - 1 - e; y_feed = s / (1 + s)
    rep.kv("extent e", e, "mol"); rep.kv("total product moles N", N, "mol"); rep.kv("steam fed s", s, "mol per mol butene")
    rep.h("Final answers")
    rep.result("a", "y 1,3-butadiene (given)", yB, "-", "given"); rep.result("a", "y hydrogen", yH2, "-", "explicit")
    rep.result("a", "y 1-butene", yC4, "-", "explicit"); rep.result("a", "y steam", ys, "-", "explicit")
    rep.result("b", "steam fed per mol 1-butene", s, "mol/mol", "explicit"); rep.result("b", "mole fraction of steam in feed", y_feed, "-", "explicit")
    rep.h("Verification")
    thermo_checks(rep, nu, [600.0, 950.0, 1200.0])
    rep.check("sum of product mole fractions = 1", abs(yB + yH2 + yC4 + ys - 1) < 1e-12)
    rep.check("equilibrium relation y_B y_H2 P/y_C4 = K", abs(yB * yH2 * P / yC4 - K) / K < 1e-12)
    rep.check("steam feed is positive (s > 0)", s > 0)
    # forward check: take the computed feed (1 butene + s steam), solve the equilibrium extent from K, recover y_B = 0.10
    res = lambda ee: lnK(nu, T) - np.log(((ee / (1 + s + ee))**2 * P / P0) / ((1 - ee) / (1 + s + ee)))
    e_fw = brentq(res, 1e-12, 1 - 1e-12, xtol=1e-14); yB_fw = e_fw / (1 + s + e_fw)
    rep.check("forward check: that feed gives y_butadiene = 0.10 at equilibrium (independent brentq)", abs(yB_fw - yB) < 1e-9, f"(y_B = {yB_fw:.10f})")
    # ---------------- graph: steam required vs butadiene content
    yBs = np.linspace(0.02, 0.24, 120); feed = []
    for yb in yBs:
        yc = yb * yb * P / K; ee = yb / (yb + yc); NN = ee / yb; ss = NN - 1 - ee
        feed.append(ss / (1 + ss) if ss > 0 else np.nan)
    fig, ax = plt.subplots(figsize=(7.2, 4.5), dpi=300)
    ax.plot(100 * yBs, feed, color="#1f4e79", lw=2.2, label="950 K, 1 bar")
    ax.plot(100 * yB, y_feed, "o", color="#c0392b", ms=8, label=f"10 mol-% butadiene: steam {y_feed:.3f}")
    ax.set(xlabel="1,3-Butadiene in product gas [mol-%]", ylabel="Steam mole fraction in feed [-]", title="Problem 13.18: steam required for a given butadiene content (950 K, 1 bar)")
    ax.title.set_fontsize(10); ax.grid(alpha=.35); ax.legend(); fig.tight_layout()
    fig.savefig(os.path.join(rep.dirs["plots"], QID + ".png"), dpi=300); plt.close(fig)
    rep.p(f"Graph saved as plots/{QID}/{QID}.png")
    ok = rep.save(); sys.exit(0 if ok else 1)

if __name__ == "__main__": main()
