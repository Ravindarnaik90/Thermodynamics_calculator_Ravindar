"""Independent audit of the existing ../src code. Run: python audit_existing_code.py
Compares the generic matrix implementations with INDEPENDENT closed-form binary equations,
checks component-order invariance, Gibbs-Duhem, Antoine normal boiling points, gamma-infinity."""
import os, sys, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))
from src.params import COMPONENTS, PAIRS, R, pair
from src.models import VLEModel
from src.vapour_pressure import psat, tsat
from src import bubble_dew as bd
from src.flash import flash_PT
from scipy.optimize import fsolve

L = ["# Audit of existing project (src/, main.py results)", ""]
bad = []
def rec(name, ok, detail=""):
    L.append(f"- [{'PASS' if ok else 'FAIL'}] {name} {detail}"); 
    if not ok: bad.append(name)

# 1 Antoine / normal boiling points (literature NBP: acetone 56.05 C, methanol 64.7 C, water 100.0 C)
L.append("## 1. Antoine equation form, units")
for n, nbp in (("acetone", 56.05), ("methanol", 64.70), ("water", 100.0)):
    t = tsat(n, 101.325) - 273.15
    rec(f"{n}: T_sat(101.325 kPa) = {t:.2f} C vs literature {nbp} C", abs(t - nbp) < 0.6)
    rec(f"{n}: psat round trip", abs(psat(n, tsat(n, 101.325)) - 101.325) < 1e-9)

# 2 Wilson closed-form binary (SVNA eq. 12.?)
def wilson_bin(x1, T, a12, a21, V1, V2):
    x2 = 1 - x1; L12 = V2 / V1 * np.exp(-a12 / (R * T)); L21 = V1 / V2 * np.exp(-a21 / (R * T))
    d = L12 / (x1 + x2 * L12) - L21 / (x2 + x1 * L21)
    return np.exp(-np.log(x1 + x2 * L12) + x2 * d), np.exp(-np.log(x2 + x1 * L21) - x1 * d)
def nrtl_bin(x1, T, b12, b21, al):
    x2 = 1 - x1; t12, t21 = b12 / (R * T), b21 / (R * T); G12, G21 = np.exp(-al * t12), np.exp(-al * t21)
    l1 = x2**2 * (t21 * (G21 / (x1 + x2 * G21))**2 + t12 * G12 / (x2 + x1 * G12)**2)
    l2 = x1**2 * (t12 * (G12 / (x2 + x1 * G12))**2 + t21 * G21 / (x1 + x2 * G21)**2)
    return np.exp(l1), np.exp(l2)
L += ["", "## 2/3. Wilson & NRTL vs independent closed-form binary equations (incl. component-order swap)"]
for (a, b), p in PAIRS.items():
    for T in (333.15, 338.15, 350.0):
        for x1 in (0.05, 0.3, 0.7, 0.95):
            gw = VLEModel([a, b], "Wilson").gamma([x1, 1 - x1], T); gn = VLEModel([a, b], "NRTL").gamma([x1, 1 - x1], T)
            cw = wilson_bin(x1, T, p["a12"], p["a21"], COMPONENTS[a]["V"], COMPONENTS[b]["V"]); cn = nrtl_bin(x1, T, p["b12"], p["b21"], p["alpha"])
            gws = VLEModel([b, a], "Wilson").gamma([1 - x1, x1], T)[::-1]; gns = VLEModel([b, a], "NRTL").gamma([1 - x1, x1], T)[::-1]
            ok = np.allclose(gw, cw, rtol=1e-10) and np.allclose(gn, cn, rtol=1e-10) and np.allclose(gw, gws, rtol=1e-10) and np.allclose(gn, gns, rtol=1e-10)
            if not ok: rec(f"{a}/{b} T={T} x1={x1}", False)
    rec(f"{a}/{b}: generic = closed form and order-invariant (12 states x 2 models)", not any(a in s and b in s for s in bad))

# 4 Gibbs-Duhem ternary
L += ["", "## 4. Gibbs-Duhem (ternary, finite difference): sum x_i dln(gamma_i) = 0"]
for mod in ("Wilson", "NRTL"):
    m = VLEModel(["acetone", "methanol", "water"], mod); worst = 0
    for x in ([0.3, 0.4, 0.3], [0.1, 0.2, 0.7], [0.6, 0.1, 0.3]):
        x = np.array(x); h = 1e-6
        for d in ([1, -1, 0], [1, 0, -1], [0, 1, -1]):
            d = np.array(d, float); dl = (np.log(m.gamma(x + h * d, 338.15)) - np.log(m.gamma(x - h * d, 338.15))) / (2 * h)
            worst = max(worst, abs(np.dot(x, dl)))
    rec(f"{mod}: max |sum x dlnγ| = {worst:.2e}", worst < 1e-6)

# 5 infinite dilution
L += ["", "## 5. Parameter plausibility: gamma-infinity (Wilson vs NRTL fitted to same data should agree)"]
for (a, b) in PAIRS:
    e = 1e-9; T = 338.15
    w = (VLEModel([a, b], "Wilson").gamma([e, 1 - e], T)[0], VLEModel([a, b], "Wilson").gamma([1 - e, e], T)[1])
    n = (VLEModel([a, b], "NRTL").gamma([e, 1 - e], T)[0], VLEModel([a, b], "NRTL").gamma([1 - e, e], T)[1])
    ratio = max(max(w[i] / n[i], n[i] / w[i]) for i in (0, 1))
    rec(f"{a}/{b}: γ∞ Wilson {w[0]:.3g}/{w[1]:.3g} vs NRTL {n[0]:.3g}/{n[1]:.3g}", ratio < 2.0,
        "" if ratio < 2 else " <-- NRTL and Wilson parameters are NOT mutually consistent; one set is probably wrong/mis-transcribed")

# 6 independent solver cross-check
L += ["", "## 6. Root-finding cross-check with scipy.fsolve (independent algorithm, different starting guess)"]
for mod in ("Wilson", "NRTL"):
    m = VLEModel(["acetone", "water"], mod); y = np.array([0.3, 0.7]); T = 333.15
    P, x, ok = bd.dew_P(m, T, y)
    f = lambda u: [u[2] * y[0] - u[0] * m.gamma(u[:2], T)[0] * m.psat(T)[0], u[2] * y[1] - u[1] * m.gamma(u[:2], T)[1] * m.psat(T)[1], u[0] + u[1] - 1]
    sol = fsolve(f, [0.1, 0.9, 50], xtol=1e-13)
    rec(f"{mod} DEW P (binary): successive substitution {P:.6f} vs fsolve {sol[2]:.6f} kPa", abs(P - sol[2]) < 1e-6)
    mt = VLEModel(["acetone", "methanol", "water"], mod); yt = np.array([0.3, 0.4, 0.3]); Pt, xt, ok = bd.dew_P(mt, 338.15, yt)
    ft = lambda u: [u[3] * yt[i] - u[i] * mt.gamma(u[:3], 338.15)[i] * mt.psat(338.15)[i] for i in range(3)] + [u[0] + u[1] + u[2] - 1]
    s = fsolve(ft, [0.1, 0.2, 0.7, 60], xtol=1e-13)
    rec(f"{mod} DEW P (ternary): {Pt:.6f} vs fsolve {s[3]:.6f} kPa", abs(Pt - s[3]) < 1e-5)

L += ["", "## 7. Code-review notes (manual)",
 "- Component ordering: pair(i,j) returns a_ij for ordered names; verified by swap-invariance test above.",
 "- Wilson indexing: Λ_ij = (V_j/V_i) exp(-a_ij/RT) uses a_ij of row species i - consistent with SVNA footnote.",
 "- Temperatures are K inside all functions; Antoine converts with T-273.15; pressure in kPa throughout.",
 "- Antoine validity ranges: the pasted ranges (acetone -13..55 C, methanol 15..84 C) are NARROWER than needed for several state points (acetone/water dew T ~ 90 C, ternary dew T ~ 76 C). Antoine is extrapolated there; flagged in each question report.",
 "- Wilson cannot predict liquid-liquid splitting and no stability test is done in the flash; acceptable for these state points but unverified.",
 "- Issue found in earlier main.py: nothing numerically wrong was found in the equations; the main weakness is INPUT DATA provenance (typed text, not book pages).",
 "- 12.22(c) uses z2 = 0.2 as printed on the image (differs from 12.20(c) z2 = 0.4); keep as printed unless the instructor says otherwise.", ""]
L += ["## Summary", f"- Failed audit items: {len(bad)}" + ("" if not bad else " -> " + "; ".join(bad))]
open(os.path.join(HERE, "reports", "AUDIT_REPORT.md"), "w").write("\n".join(L)); print("\n".join(L))
