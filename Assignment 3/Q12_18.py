"""Q12_18  --  Question 12.18  (Wilson).  Run:  python Q12_18.py
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

QID, MODEL = "Q12_18", "Wilson"
STATEMENT = "Wilson equation: (a) BUBL T: P = 101.33 kPa, x1 = 0.3; (b) DEW T: P = 101.33 kPa, y1 = 0.3; (c) P,T-flash: P = 101.33 kPa, T = (T_bubble+T_dew)/2, z1 = 0.3; (d) if an azeotrope exists at P = 101.33 kPa find T_az and x1_az = y1_az. (system: acetone(1)/water(2), chosen by the student from the table) "
EQ_RAOULT = "Eq. (10.5), modified Raoult's law:  y_i P = x_i gamma_i P_i^sat      (ideal vapour phase)"
EQ_ANTOINE = "Antoine:  ln(P_i^sat/kPa) = A_i - B_i / (t/degC + C_i),  t = T - 273.15"
EQ_MODEL = None
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

def main():
    rep = Report(QID, "Question 12.18 - Wilson BUBL T / DEW T / flash / azeotrope, acetone(1)/water(2), P = 101.33 kPa", HERE)
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
