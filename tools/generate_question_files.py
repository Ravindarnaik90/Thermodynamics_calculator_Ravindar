"""Writes Q12_12.py ... Q12_23.py from six templates. The generated files are plain, standalone scripts."""
import os
HERE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")

HEAD = '''"""@@QID@@  --  Question @@QNUM@@  (@@MODEL@@).  Run:  python @@QID@@.py
Self-contained: uses only the helper package common/ (thermodynamic functions); does NOT run other question files."""
import os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from common.models import VLEModel, check_vle
from common.bubble_dew import bubble_P, dew_P, bubble_T, dew_T
from common.flash import flash_PT
from common.azeotrope import azeotrope_T, azeotrope_P
from common.report import Report, describe_model
from common.diagrams import pxy_data, txy_data, plot_pxy, plot_txy, curve_checks
from common.crosscheck import xc_dew_P, xc_bubble_T, xc_dew_T, xc_flash
from common.vapour_pressure import psat, tsat

QID, MODEL = "@@QID@@", "@@MODEL@@"
STATEMENT = "@@STATEMENT@@"
EQ_RAOULT = "Eq. (10.5), modified Raoult's law:  y_i P = x_i gamma_i P_i^sat      (ideal vapour phase)"
EQ_ANTOINE = "Antoine:  ln(P_i^sat/kPa) = A_i - B_i / (t/degC + C_i),  t = T - 273.15"
EQ_MODEL = @@EQMODEL@@
EQ_W = ["Wilson:  Lambda_ij = (V_j/V_i) exp(-a_ij/RT),  Lambda_ii = 1",
        "         ln gamma_i = 1 - ln(sum_j x_j Lambda_ij) - sum_k [ x_k Lambda_ki / sum_j x_j Lambda_kj ]"]
EQ_N = ["NRTL:  tau_ij = b_ij/RT,  G_ij = exp(-alpha_ij tau_ij),  tau_ii = 0, G_ii = 1",
        "       ln gamma_i = sum_j x_j tau_ji G_ji / sum_k x_k G_ki",
        "                    + sum_j [ x_j G_ij / sum_k x_k G_kj ] ( tau_ij - sum_m x_m tau_mj G_mj / sum_k x_k G_kj )"]

def header(rep):
    rep.h("Problem statement"); rep.p(STATEMENT)
    rep.h("Model and governing equations")
    for e in (EQ_RAOULT, EQ_ANTOINE): rep.eq(e)
    for e in (EQ_W if MODEL == "Wilson" else EQ_N): rep.eq(e)
    rep.eq("Units: T in K inside all functions, P in kPa, R = 8.314 J/(mol K), a_ij and b_ij in J/mol")
'''

CURVE_P = HEAD + '''
def main():
    rep = Report(QID, "Question @@QNUM@@ - @@MODEL@@ P-x-y diagram, acetone(1)/water(2), T = 333.15 K", HERE)
    header(rep)
    rep.h("Given data")
    T = 333.15; COMPS = ["acetone", "water"]
    rep.kv("System", "acetone(1) / water(2)  (binary chosen by the student from the table)"); rep.kv("T", T, "K (60 degC)")
    rep.kv("Composition grid", "x1 = 0, 0.02, ..., 1.0 (51 points)")
    m = VLEModel(COMPS, MODEL); rep.range_warnings(COMPS, [T])
    rep.h("Step 1 - pure-component vapour pressures and model parameters at T")
    describe_model(rep, m, T, [0.5, 0.5])
    rep.h("Step 2 - for every x1: BUBL P  (P = sum x_i gamma_i Psat_i,  y_i = x_i gamma_i Psat_i / P)")
    rep.h("Step 3 - DEW P computed independently for the same compositions used as y1 (dew curve)")
    d = pxy_data(m, T)
    rep.h("Computed diagram data (every 5th point; full data in the CSV)")
    rep.table(["x1", "P_bubble [kPa]", "y1", "P_dew(y1=x1 col) [kPa]", "x1 on dew curve"],
              [[d["x1_or_y1"][i], d["P_bubble_kPa"][i], d["y1_vapour_in_equilibrium_with_x1"][i], d["P_dew_kPa"][i], d["x1_liquid_in_equilibrium_with_y1"][i]] for i in range(0, 51, 5)])
    import pandas as pd
    pd.DataFrame(d).to_csv(os.path.join(rep.dirs["results"], QID + "_diagram_data.csv"), index=False)
    png = os.path.join(rep.dirs["plots"], QID + "_" + MODEL + "_Pxy.png"); plot_pxy(d, MODEL, T, QID, png)
    rep.h("Final answers")
    rep.result("diagram", "P_sat water", d["P_bubble_kPa"][0], "kPa", "explicit (Antoine)")
    rep.result("diagram", "P_sat acetone", d["P_bubble_kPa"][-1], "kPa", "explicit (Antoine)")
    rep.result("diagram", "diagram points", 51, "-", "all points converged")
    rep.p(f"Graph saved as plots/{QID}/{os.path.basename(png)}")
    rep.h("Verification"); curve_checks(rep, m, d, "P", T)
    ok = rep.save(); sys.exit(0 if ok else 1)

if __name__ == "__main__": main()
'''

CURVE_T = HEAD + '''
def main():
    rep = Report(QID, "Question @@QNUM@@ - @@MODEL@@ t-x-y diagram, acetone(1)/water(2), P = 101.33 kPa", HERE)
    header(rep)
    rep.h("Given data")
    P = 101.33; COMPS = ["acetone", "water"]
    rep.kv("System", "acetone(1) / water(2)  (binary chosen by the student from the table)"); rep.kv("P", P, "kPa")
    rep.kv("Composition grid", "x1 = 0, 0.02, ..., 1.0 (51 points)")
    m = VLEModel(COMPS, MODEL)
    rep.range_warnings(COMPS, [tsat("acetone", P), 363.0], "(bubble/dew range at this P)")
    rep.h("Step 1 - pure-component saturation temperatures from inverse Antoine: T_sat = B/(A - ln P) - C + 273.15")
    rep.kv("T_sat acetone", tsat("acetone", P), "K"); rep.kv("T_sat water", tsat("water", P), "K")
    rep.h("Step 2 - for every x1: BUBL T. Solve  sum x_i gamma_i(T,x) Psat_i(T) - P = 0  (scipy brentq, bracket 250-430 K); then y_i = x_i gamma_i Psat_i / P")
    rep.h("Step 3 - DEW T for the same compositions used as y1: solve  P_dew(T,y) - P = 0  (brentq around damped successive substitution for x)")
    d = txy_data(m, P)
    rep.h("Computed diagram data (every 5th point; full data in the CSV)")
    rep.table(["x1", "T_bubble [K]", "T_bubble [degC]", "y1", "T_dew(y1=x1 col) [K]", "x1 on dew curve"],
              [[d["x1_or_y1"][i], d["T_bubble_K"][i], d["T_bubble_K"][i] - 273.15, d["y1_vapour_in_equilibrium_with_x1"][i], d["T_dew_K"][i], d["x1_liquid_in_equilibrium_with_y1"][i]] for i in range(0, 51, 5)])
    import pandas as pd
    pd.DataFrame(d).to_csv(os.path.join(rep.dirs["results"], QID + "_diagram_data.csv"), index=False)
    png = os.path.join(rep.dirs["plots"], QID + "_" + MODEL + "_txy.png"); plot_txy(d, MODEL, P, QID, png)
    rep.h("Final answers")
    rep.result("diagram", "T_sat water", d["T_bubble_K"][0], "K", "explicit (inverse Antoine)")
    rep.result("diagram", "T_sat acetone", d["T_bubble_K"][-1], "K", "explicit (inverse Antoine)")
    rep.result("diagram", "diagram points", 51, "-", "all brentq solves converged")
    rep.p(f"Graph saved as plots/{QID}/{os.path.basename(png)}")
    rep.h("Verification"); curve_checks(rep, m, d, "T", P)
    ok = rep.save(); sys.exit(0 if ok else 1)

if __name__ == "__main__": main()
'''

# ---------- point-calculation bodies ----------
COMMON_FLASH = '''
    # ------------------------------------------------------------------ (c) flash
    rep.h("(c) P,T-flash")
    rep.p("Material balance z_i = V y_i + (1-V) x_i with y_i = K_i x_i,  K_i = gamma_i Psat_i / P  =>  x_i = z_i/(1+V(K_i-1)).  Rachford-Rice: sum z_i(K_i-1)/(1+V(K_i-1)) = 0, solved with brentq; gamma updated with x until x converges (tol 1e-13).")
    fl = flash_PT(m, P_f, T_f, z)
    rep.kv("flash P", P_f, "kPa"); rep.kv("flash T", T_f, "K"); rep.kv("z", list(z))
    rep.kv("P_bubble(z) and P_dew(z) at flash T", f"{fl['P_bubble_at_T']:.6g} kPa, {fl['P_dew_at_T']:.6g} kPa")
    rep.kv("phase state", fl["phase"]); rep.kv("outer iterations", fl["iters"])
    st = "converged" if fl["converged"] else "NOT CONVERGED"
    rep.result("c", "flash pressure", P_f, "kPa", st); rep.result("c", "flash temperature", T_f, "K", st)
    rep.result("c", "vapour fraction V", fl["V"], "mol/mol", st)
    if fl["x"] is not None:
        for c, v in zip(COMPS, fl["x"]): rep.result("c", f"x_{c}", v, "-", st)
    if fl["y"] is not None:
        for c, v in zip(COMPS, fl["y"]): rep.result("c", f"y_{c}", v, "-", st)
    rep.check("flash converged", fl["converged"]); rep.check("0 <= V <= 1", 0 <= fl["V"] <= 1)
    if fl["phase"] == "two-phase":
        x_f, y_f, V = fl["x"], fl["y"], fl["V"]
        K = m.gamma(x_f, T_f) * m.psat(T_f) / P_f
        rep.vec("K_i = gamma_i Psat_i / P", COMPS, K)
        rep.check("sum x = 1 and sum y = 1", abs(x_f.sum() - 1) < 1e-12 and abs(y_f.sum() - 1) < 1e-12)
        rep.check("0 <= x_i, y_i <= 1", bool(np.all((x_f >= 0) & (x_f <= 1) & (y_f >= 0) & (y_f <= 1))))
        rep.check("component material balance z_i = V y_i + (1-V) x_i", np.max(np.abs(z - (V * y_f + (1 - V) * x_f))) < 1e-10, f"(max err {np.max(np.abs(z - (V*y_f + (1-V)*x_f))):.1e})")
        rep.check("equilibrium y_i = K_i x_i  (Eq. 10.5 residual)", check_vle(m, T_f, P_f, x_f, y_f)["max_resid_kPa"] < 1e-8)
        Vx, xx, okx = xc_flash(m, P_f, T_f, z, V0=0.5)
        rep.check("independent fsolve flash (different algorithm/start) agrees", okx and abs(Vx - V) < 1e-7 and np.allclose(xx, x_f, atol=1e-7), f"(V {Vx:.8f} vs {V:.8f})")
    else:
        rep.note("Specified condition is outside the two-phase region; reported as a single phase (no forced two-phase solution).")
'''

BIN_P = HEAD + '''
def main():
    rep = Report(QID, "Question @@QNUM@@ - @@MODEL@@ BUBL P / DEW P / flash / azeotrope, acetone(1)/water(2), T = 333.15 K", HERE)
    header(rep)
    rep.h("Given data")
    COMPS = ["acetone", "water"]; T = 333.15
    x = np.array([0.3, 0.7]); y_d = np.array([0.3, 0.7]); z = np.array([0.3, 0.7]); T_f = T
    rep.kv("System", "acetone(1)/water(2)"); rep.kv("T", T, "K"); rep.kv("(a) x1", 0.3); rep.kv("(b) y1", 0.3); rep.kv("(c) z1", 0.3, "and P = (P_bubble + P_dew)/2")
    m = VLEModel(COMPS, MODEL); rep.range_warnings(COMPS, [T])
    # ------------------------------------------------------------------ (a)
    rep.h("(a) BUBL P:  P = sum x_i gamma_i Psat_i ;  y_i = x_i gamma_i Psat_i / P")
    describe_model(rep, m, T, x)
    Pb, yb = bubble_P(m, T, x)
    rep.vec("partial pressures x_i gamma_i Psat_i", COMPS, x * m.gamma(x, T) * m.psat(T), "kPa")
    rep.result("a", "P_bubble", Pb, "kPa", "explicit"); [rep.result("a", f"y_{c}", v, "-", "explicit") for c, v in zip(COMPS, yb)]
    c = check_vle(m, T, Pb, x, yb)
    rep.check("sum x = 1, sum y = 1", abs(c["sum_x"] - 1) < 1e-12 and abs(c["sum_y"] - 1) < 1e-12)
    rep.check("Eq.(10.5) residual y_i P - x_i gamma_i Psat_i = 0", c["max_resid_kPa"] < 1e-9, f"({c['max_resid_kPa']:.1e} kPa)")
    Pr, xr, okr = dew_P(m, T, yb); rep.check("reverse: DEW P at y_bubble returns same P and x", okr and abs(Pr - Pb) < 1e-8 and np.allclose(xr, x, atol=1e-8))
    # ------------------------------------------------------------------ (b)
    rep.h("(b) DEW P:  P = 1 / sum(y_i/(gamma_i Psat_i)) ;  x_i = y_i P / (gamma_i Psat_i)  (iterate: gamma depends on x)")
    Pd, xd, ok = dew_P(m, T, y_d)
    rep.vec("liquid composition at dew point", COMPS, xd); rep.vec("gamma at dew liquid", COMPS, m.gamma(xd, T))
    st = "converged (damped successive substitution)" if ok else "NOT CONVERGED"
    rep.result("b", "P_dew", Pd, "kPa", st); [rep.result("b", f"x_{c}", v, "-", st) for c, v in zip(COMPS, xd)]
    c = check_vle(m, T, Pd, xd, y_d)
    rep.check("DEW P converged", ok); rep.check("sum x = 1, sum y = 1", abs(c["sum_x"] - 1) < 1e-10 and abs(c["sum_y"] - 1) < 1e-12)
    rep.check("Eq.(10.5) residual", c["max_resid_kPa"] < 1e-8, f"({c['max_resid_kPa']:.1e} kPa)")
    Px, xx, okx = xc_dew_P(m, T, y_d); rep.check("independent fsolve DEW P agrees", okx and abs(Px - Pd) < 1e-6 and np.allclose(xx, xd, atol=1e-7), f"({Px:.6f} vs {Pd:.6f} kPa)")
    rep.check("P_dew < P_bubble (two-phase region exists for z = 0.3)", Pd < Pb)
    P_f = 0.5 * (Pb + Pd); rep.kv("P for flash = (P_bubble + P_dew)/2", P_f, "kPa")
''' + COMMON_FLASH + '''
    # ------------------------------------------------------------------ (d)
    rep.h("(d) Azeotrope at T:  condition x1 = y1 together with Eq.(10.5)  <=>  gamma1 Psat1 = gamma2 Psat2 at a single x1")
    roots, diag = azeotrope_T(m, T)
    xs = np.linspace(1e-6, 1 - 1e-6, 4000); alpha = np.array([m.gamma([a, 1 - a], T)[0] * m.psat(T)[0] / (m.gamma([a, 1 - a], T)[1] * m.psat(T)[1]) for a in xs])
    rep.kv("relative volatility alpha12 = K1/K2: min / max over 0<x1<1", f"{alpha.min():.6g} / {alpha.max():.6g}")
    if roots:
        for k, r in enumerate(roots):
            rep.result("d", f"x1_az (root {k+1})", r["x1"], "-", "brentq on y1-x1"); rep.result("d", f"P_az (root {k+1})", r["P"], "kPa", "brentq on y1-x1")
            xa = np.array([r["x1"], 1 - r["x1"]]); Pa, ya = bubble_P(m, T, xa)
            rep.check(f"azeotrope root {k+1}: y1 = x1 within 1e-10", abs(ya[0] - xa[0]) < 1e-10)
    else:
        rep.result("d", "azeotrope exists", 0, "(0 = no)", "scan of 600 points: no sign change of y1-x1")
        rep.kv("closest approach |y1-x1| (excluding pure-component ends)", diag["min_abs_y1_minus_x1"]); rep.kv("at x1", diag["at_x1"])
        rep.check("no azeotrope: alpha12 > 1 for all x1 (independent criterion)", alpha.min() > 1.0, f"(min alpha12 = {alpha.min():.4f})")
        rep.p("Interpretation: an azeotrope requires alpha12 = 1 at some composition. alpha12 stays above 1 over the whole range, so acetone is always enriched in the vapour and no azeotrope exists at this temperature for the model.")
    ok = rep.save(); sys.exit(0 if ok else 1)

if __name__ == "__main__": main()
'''

BIN_T = HEAD + '''
def main():
    rep = Report(QID, "Question @@QNUM@@ - @@MODEL@@ BUBL T / DEW T / flash / azeotrope, acetone(1)/water(2), P = 101.33 kPa", HERE)
    header(rep)
    rep.h("Given data")
    COMPS = ["acetone", "water"]; P = 101.33
    x = np.array([0.3, 0.7]); y_d = np.array([0.3, 0.7]); z = np.array([0.3, 0.7]); P_f = P
    rep.kv("System", "acetone(1)/water(2)"); rep.kv("P", P, "kPa"); rep.kv("(a) x1", 0.3); rep.kv("(b) y1", 0.3); rep.kv("(c) z1", 0.3, "and T = (T_bubble + T_dew)/2")
    m = VLEModel(COMPS, MODEL)
    # ------------------------------------------------------------------ (a)
    rep.h("(a) BUBL T:  solve  f(T) = sum x_i gamma_i(T) Psat_i(T) - P = 0  (brentq, 250-430 K) ;  y_i = x_i gamma_i Psat_i / P")
    Tb, yb = bubble_T(m, P, x)
    describe_model(rep, m, Tb, x)
    rep.vec("partial pressures at T_bubble", COMPS, x * m.gamma(x, Tb) * m.psat(Tb), "kPa")
    rep.result("a", "T_bubble", Tb, "K", "converged (brentq)"); rep.result("a", "T_bubble", Tb - 273.15, "degC", "converged (brentq)")
    [rep.result("a", f"y_{c}", v, "-", "converged (brentq)") for c, v in zip(COMPS, yb)]
    c = check_vle(m, Tb, P, x, yb); rep.range_warnings(COMPS, [Tb])
    rep.check("sum y = 1", abs(c["sum_y"] - 1) < 1e-12); rep.check("Eq.(10.5) residual", c["max_resid_kPa"] < 1e-8, f"({c['max_resid_kPa']:.1e} kPa)")
    rep.check("bubble equation residual sum x gamma Psat - P = 0", abs(np.sum(x * m.gamma(x, Tb) * m.psat(Tb)) - P) < 1e-8)
    Tx, okx = xc_bubble_T(m, P, x); rep.check("independent fsolve BUBL T agrees", okx and abs(Tx - Tb) < 1e-6, f"({Tx:.6f} vs {Tb:.6f} K)")
    # ------------------------------------------------------------------ (b)
    rep.h("(b) DEW T:  solve  P_dew(T, y) - P = 0  (brentq; inner DEW P by damped successive substitution)")
    Td, xd = dew_T(m, P, y_d); Pchk, xchk, okd = dew_P(m, Td, y_d)
    rep.vec("gamma at dew liquid", COMPS, m.gamma(xd, Td))
    rep.result("b", "T_dew", Td, "K", "converged (brentq)"); rep.result("b", "T_dew", Td - 273.15, "degC", "converged (brentq)")
    [rep.result("b", f"x_{c}", v, "-", "converged (brentq)") for c, v in zip(COMPS, xd)]
    c = check_vle(m, Td, P, xd, y_d); rep.range_warnings(COMPS, [Td])
    rep.check("sum x = 1, sum y = 1", abs(c["sum_x"] - 1) < 1e-10 and abs(c["sum_y"] - 1) < 1e-12)
    rep.check("Eq.(10.5) residual", c["max_resid_kPa"] < 1e-7, f"({c['max_resid_kPa']:.1e} kPa)")
    rep.check("dew equation residual sum y/(gamma Psat) - 1/P = 0", abs(np.sum(y_d / (m.gamma(xd, Td) * m.psat(Td))) - 1 / P) < 1e-9)
    Tx, xx, okx = xc_dew_T(m, P, y_d); rep.check("independent fsolve DEW T agrees", okx and abs(Tx - Td) < 1e-5 and np.allclose(xx, xd, atol=1e-6), f"({Tx:.6f} vs {Td:.6f} K)")
    rep.check("T_bubble < T_dew", Tb < Td)
    T_f = 0.5 * (Tb + Td); rep.kv("T for flash = (T_bubble + T_dew)/2", T_f, "K")
''' + COMMON_FLASH + '''
    # ------------------------------------------------------------------ (d)
    rep.h("(d) Azeotrope at P:  condition x1 = y1 together with Eq.(10.5)  <=>  alpha12 = gamma1 Psat1/(gamma2 Psat2) = 1 at one x1 (along the bubble curve)")
    roots, diag = azeotrope_P(m, P)
    xs = np.linspace(1e-4, 1 - 1e-4, 400); al = []
    for a in xs:
        Ta, _ = bubble_T(m, P, [a, 1 - a]); g = m.gamma([a, 1 - a], Ta); ps = m.psat(Ta); al.append(g[0] * ps[0] / (g[1] * ps[1]))
    al = np.array(al); rep.kv("relative volatility alpha12 along bubble curve: min / max", f"{al.min():.6g} / {al.max():.6g}")
    if roots:
        for k, r in enumerate(roots):
            rep.result("d", f"x1_az (root {k+1})", r["x1"], "-", "brentq on y1-x1"); rep.result("d", f"T_az (root {k+1})", r["T"], "K", "brentq on y1-x1")
            xa = np.array([r["x1"], 1 - r["x1"]]); Ta, ya = bubble_T(m, P, xa)
            rep.check(f"azeotrope root {k+1}: y1 = x1 within 1e-10", abs(ya[0] - xa[0]) < 1e-10)
    else:
        rep.result("d", "azeotrope exists", 0, "(0 = no)", "scan of 300 points: no sign change of y1-x1")
        rep.kv("closest approach |y1-x1| (excluding pure-component ends)", diag["min_abs_y1_minus_x1"]); rep.kv("at x1", diag["at_x1"])
        rep.check("no azeotrope: alpha12 > 1 for all x1 (independent criterion)", al.min() > 1.0, f"(min alpha12 = {al.min():.4f})")
        rep.p("Interpretation: an azeotrope requires alpha12 = 1 at some composition. alpha12 stays above 1 along the whole bubble curve, so no azeotrope exists at this pressure for the model.")
    ok = rep.save(); sys.exit(0 if ok else 1)

if __name__ == "__main__": main()
'''

TER_P = HEAD + '''
def main():
    rep = Report(QID, "Question @@QNUM@@ - @@MODEL@@ BUBL P / DEW P / flash, acetone(1)/methanol(2)/water(3), T = 338.15 K", HERE)
    header(rep)
    rep.h("Given data")
    COMPS = ["acetone", "methanol", "water"]; T = 338.15
    x = np.array([0.3, 0.4, 0.3]); y_d = np.array([0.3, 0.4, 0.3]); z = np.array([0.3, 0.4, 0.3]); T_f = T
    rep.kv("System", "acetone(1)/methanol(2)/water(3)"); rep.kv("T", T, "K (65 degC)")
    rep.kv("(a) x1, x2, x3 = 1 - x1 - x2", list(x)); rep.kv("(b) y1, y2, y3 = 1 - y1 - y2", list(y_d)); rep.kv("(c) z1, z2, z3 = 1 - z1 - z2", list(z), "and P = (P_bubble + P_dew)/2")
    m = VLEModel(COMPS, MODEL); rep.range_warnings(COMPS, [T])
    rep.h("(a) BUBL P:  P = sum x_i gamma_i Psat_i ;  y_i = x_i gamma_i Psat_i / P   (complete ternary model, all three binary pairs)")
    describe_model(rep, m, T, x)
    Pb, yb = bubble_P(m, T, x)
    rep.vec("partial pressures x_i gamma_i Psat_i", COMPS, x * m.gamma(x, T) * m.psat(T), "kPa")
    rep.result("a", "P_bubble", Pb, "kPa", "explicit"); [rep.result("a", f"y_{c}", v, "-", "explicit") for c, v in zip(COMPS, yb)]
    c = check_vle(m, T, Pb, x, yb)
    rep.check("sum x = 1, sum y = 1", abs(c["sum_x"] - 1) < 1e-12 and abs(c["sum_y"] - 1) < 1e-12)
    rep.check("Eq.(10.5) residual", c["max_resid_kPa"] < 1e-9, f"({c['max_resid_kPa']:.1e} kPa)")
    Pr, xr, okr = dew_P(m, T, yb); rep.check("reverse: DEW P at y_bubble returns same P and x", okr and abs(Pr - Pb) < 1e-7 and np.allclose(xr, x, atol=1e-7))
    rep.h("(b) DEW P:  P = 1/sum(y_i/(gamma_i Psat_i)) ;  x_i = y_i P/(gamma_i Psat_i)  (iterated)")
    Pd, xd, ok = dew_P(m, T, y_d); st = "converged (damped successive substitution)" if ok else "NOT CONVERGED"
    rep.vec("gamma at dew liquid", COMPS, m.gamma(xd, T))
    rep.result("b", "P_dew", Pd, "kPa", st); [rep.result("b", f"x_{c}", v, "-", st) for c, v in zip(COMPS, xd)]
    c = check_vle(m, T, Pd, xd, y_d)
    rep.check("DEW P converged", ok); rep.check("sum x = 1, sum y = 1", abs(c["sum_x"] - 1) < 1e-10 and abs(c["sum_y"] - 1) < 1e-12)
    rep.check("Eq.(10.5) residual", c["max_resid_kPa"] < 1e-8, f"({c['max_resid_kPa']:.1e} kPa)")
    Px, xx, okx = xc_dew_P(m, T, y_d); rep.check("independent fsolve DEW P agrees", okx and abs(Px - Pd) < 1e-5 and np.allclose(xx, xd, atol=1e-6), f"({Px:.6f} vs {Pd:.6f} kPa)")
    rep.check("P_dew < P_bubble", Pd < Pb)
    P_f = 0.5 * (Pb + Pd); rep.kv("P for flash = (P_bubble + P_dew)/2", P_f, "kPa")
''' + COMMON_FLASH + '''
    rep.p("(No azeotrope part is asked in this question.)")
    ok = rep.save(); sys.exit(0 if ok else 1)

if __name__ == "__main__": main()
'''

TER_T = HEAD + '''
def main():
    rep = Report(QID, "Question @@QNUM@@ - @@MODEL@@ BUBL T / DEW T / flash, acetone(1)/methanol(2)/water(3), P = 101.33 kPa", HERE)
    header(rep)
    rep.h("Given data")
    COMPS = ["acetone", "methanol", "water"]; P = 101.33
    x = np.array([0.3, 0.4, 0.3]); y_d = np.array([0.3, 0.4, 0.3]); z = np.array([0.3, 0.2, 0.5]); P_f = P
    rep.kv("System", "acetone(1)/methanol(2)/water(3)"); rep.kv("P", P, "kPa")
    rep.kv("(a) x1, x2, x3 = 1 - x1 - x2", list(x)); rep.kv("(b) y1, y2, y3 = 1 - y1 - y2", list(y_d))
    rep.kv("(c) z1, z2, z3 = 1 - z1 - z2", list(z), "and T = (T_bubble + T_dew)/2")
    rep.note("In (c) the assignment image prints z2 = 0.2 (whereas Q12.20(c) uses z2 = 0.4). Used as printed: z = (0.3, 0.2, 0.5). Confirm with your handout.")
    m = VLEModel(COMPS, MODEL)
    rep.h("(a) BUBL T:  solve  sum x_i gamma_i(T) Psat_i(T) - P = 0  (brentq) ;  y_i = x_i gamma_i Psat_i/P")
    Tb, yb = bubble_T(m, P, x); describe_model(rep, m, Tb, x)
    rep.vec("partial pressures at T_bubble", COMPS, x * m.gamma(x, Tb) * m.psat(Tb), "kPa")
    rep.result("a", "T_bubble", Tb, "K", "converged (brentq)"); rep.result("a", "T_bubble", Tb - 273.15, "degC", "converged (brentq)")
    [rep.result("a", f"y_{c}", v, "-", "converged (brentq)") for c, v in zip(COMPS, yb)]
    c = check_vle(m, Tb, P, x, yb); rep.range_warnings(COMPS, [Tb])
    rep.check("sum y = 1", abs(c["sum_y"] - 1) < 1e-12); rep.check("Eq.(10.5) residual", c["max_resid_kPa"] < 1e-8, f"({c['max_resid_kPa']:.1e} kPa)")
    rep.check("bubble equation residual", abs(np.sum(x * m.gamma(x, Tb) * m.psat(Tb)) - P) < 1e-8)
    Tx, okx = xc_bubble_T(m, P, x); rep.check("independent fsolve BUBL T agrees", okx and abs(Tx - Tb) < 1e-6, f"({Tx:.6f} vs {Tb:.6f} K)")
    rep.h("(b) DEW T:  solve  P_dew(T, y) - P = 0  (brentq; inner DEW P by successive substitution)")
    Td, xd = dew_T(m, P, y_d); rep.vec("gamma at dew liquid", COMPS, m.gamma(xd, Td))
    rep.result("b", "T_dew", Td, "K", "converged (brentq)"); rep.result("b", "T_dew", Td - 273.15, "degC", "converged (brentq)")
    [rep.result("b", f"x_{c}", v, "-", "converged (brentq)") for c, v in zip(COMPS, xd)]
    c = check_vle(m, Td, P, xd, y_d); rep.range_warnings(COMPS, [Td])
    rep.check("sum x = 1, sum y = 1", abs(c["sum_x"] - 1) < 1e-10 and abs(c["sum_y"] - 1) < 1e-12)
    rep.check("Eq.(10.5) residual", c["max_resid_kPa"] < 1e-7, f"({c['max_resid_kPa']:.1e} kPa)")
    rep.check("dew equation residual", abs(np.sum(y_d / (m.gamma(xd, Td) * m.psat(Td))) - 1 / P) < 1e-9)
    Tx, xx, okx = xc_dew_T(m, P, y_d); rep.check("independent fsolve DEW T agrees", okx and abs(Tx - Td) < 1e-5 and np.allclose(xx, xd, atol=1e-6), f"({Tx:.6f} vs {Td:.6f} K)")
    rep.check("T_bubble < T_dew", Tb < Td)
    T_f = 0.5 * (Tb + Td); rep.kv("T for flash = (T_bubble + T_dew)/2", T_f, "K")
''' + COMMON_FLASH + '''
    rep.p("(No azeotrope part is asked in this question.)")
    ok = rep.save(); sys.exit(0 if ok else 1)

if __name__ == "__main__": main()
'''

BINSTAT = "(system: acetone(1)/water(2), chosen by the student from the table) "
TERSTAT = "Acetone(1)/methanol(2)/water(3). "
ST = {
 "12.12": (CURVE_P, "Wilson", "Based on Eq. (10.5) and the Wilson equation prepare a P-x-y diagram for T = 333.15 K (60 degC). " + BINSTAT),
 "12.13": (CURVE_T, "Wilson", "Based on Eq. (10.5) and the Wilson equation prepare a t-x-y diagram for P = 101.33 kPa. " + BINSTAT),
 "12.14": (CURVE_P, "NRTL", "Based on Eq. (10.5) and the NRTL equation prepare a P-x-y diagram for T = 333.15 K (60 degC). " + BINSTAT),
 "12.15": (CURVE_T, "NRTL", "Based on Eq. (10.5) and the NRTL equation prepare a t-x-y diagram for P = 101.33 kPa. " + BINSTAT),
 "12.16": (BIN_P, "Wilson", "Wilson equation: (a) BUBL P: T = 333.15 K, x1 = 0.3; (b) DEW P: T = 333.15 K, y1 = 0.3; (c) P,T-flash: T = 333.15 K, P = (P_bubble+P_dew)/2, z1 = 0.3; (d) if an azeotrope exists at T = 333.15 K find P_az and x1_az = y1_az. " + BINSTAT),
 "12.17": (BIN_P, "NRTL", "Work 12.16 for the NRTL equation. " + BINSTAT),
 "12.18": (BIN_T, "Wilson", "Wilson equation: (a) BUBL T: P = 101.33 kPa, x1 = 0.3; (b) DEW T: P = 101.33 kPa, y1 = 0.3; (c) P,T-flash: P = 101.33 kPa, T = (T_bubble+T_dew)/2, z1 = 0.3; (d) if an azeotrope exists at P = 101.33 kPa find T_az and x1_az = y1_az. " + BINSTAT),
 "12.19": (BIN_T, "NRTL", "Work 12.18 for the NRTL equation. " + BINSTAT),
 "12.20": (TER_P, "Wilson", TERSTAT + "Wilson equation: (a) BUBL P: T = 338.15 K, x1 = 0.3, x2 = 0.4; (b) DEW P: T = 338.15 K, y1 = 0.3, y2 = 0.4; (c) P,T-flash: T = 338.15 K, P = (P_bubble+P_dew)/2, z1 = 0.3, z2 = 0.4."),
 "12.21": (TER_P, "NRTL", "Work 12.20 for the NRTL equation. " + TERSTAT),
 "12.22": (TER_T, "Wilson", TERSTAT + "Wilson equation: (a) BUBL T: P = 101.33 kPa, x1 = 0.3, x2 = 0.4; (b) DEW T: P = 101.33 kPa, y1 = 0.3, y2 = 0.4; (c) P,T-flash: P = 101.33 kPa, T = (T_bubble+T_dew)/2, z1 = 0.3, z2 = 0.2 (as printed)."),
 "12.23": (TER_T, "NRTL", "Work 12.22 for the NRTL equation. " + TERSTAT),
}
for q, (tpl, model, stmt) in ST.items():
    code = tpl.replace("@@QID@@", "Q" + q.replace(".", "_")).replace("@@QNUM@@", q).replace("@@MODEL@@", model).replace("@@STATEMENT@@", stmt).replace("@@EQMODEL@@", "None")
    open(os.path.join(HERE, "Q" + q.replace(".", "_") + ".py"), "w", encoding="utf-8").write(code)
print("written", len(ST), "files")
