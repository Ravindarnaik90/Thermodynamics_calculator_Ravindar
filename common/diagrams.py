"""P-x-y and t-x-y diagram data (from actual VLE calculations), plots (300 dpi) and curve checks."""
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from . import bubble_dew as bd
from .vapour_pressure import psat, tsat
BIN = ["acetone", "water"]   # component 1 = acetone, component 2 = water

def pxy_data(m, T, n=51):
    xs = np.linspace(0, 1, n); Pb, Yb, Pd, Xd = [], [], [], []
    for x1 in xs:
        if x1 in (0.0, 1.0):
            p = psat(BIN[0] if x1 == 1 else BIN[1], T); Pb.append(p); Yb.append(x1); Pd.append(p); Xd.append(x1); continue
        P, y = bd.bubble_P(m, T, [x1, 1 - x1]); Pb.append(P); Yb.append(y[0])
        Pdw, xd, ok = bd.dew_P(m, T, [x1, 1 - x1]); assert ok, "DEW P did not converge"; Pd.append(Pdw); Xd.append(xd[0])
    return dict(x1_or_y1=xs, P_bubble_kPa=np.array(Pb), y1_vapour_in_equilibrium_with_x1=np.array(Yb),
                P_dew_kPa=np.array(Pd), x1_liquid_in_equilibrium_with_y1=np.array(Xd))

def txy_data(m, P, n=51):
    xs = np.linspace(0, 1, n); Tb, Yb, Td, Xd = [], [], [], []
    for x1 in xs:
        if x1 in (0.0, 1.0):
            t = tsat(BIN[0] if x1 == 1 else BIN[1], P); Tb.append(t); Yb.append(x1); Td.append(t); Xd.append(x1); continue
        T, y = bd.bubble_T(m, P, [x1, 1 - x1]); Tb.append(T); Yb.append(y[0])
        Tdw, xd = bd.dew_T(m, P, [x1, 1 - x1]); Td.append(Tdw); Xd.append(xd[0])
    return dict(x1_or_y1=xs, T_bubble_K=np.array(Tb), y1_vapour_in_equilibrium_with_x1=np.array(Yb),
                T_dew_K=np.array(Td), x1_liquid_in_equilibrium_with_y1=np.array(Xd))

def plot_pxy(d, model, T, qid, path):
    fig, ax = plt.subplots(figsize=(7.5, 5.5), dpi=300)
    ax.plot(d["x1_or_y1"], d["P_bubble_kPa"], "b-", lw=2, label="Bubble-point curve: P vs $x_1$ (liquid)")
    ax.plot(d["x1_or_y1"], d["P_dew_kPa"], "r--", lw=2, label="Dew-point curve: P vs $y_1$ (vapour)")
    ax.plot([0, 1], [d["P_bubble_kPa"][0], d["P_bubble_kPa"][-1]], "ko", ms=7, label="Pure-component $P^{sat}$ endpoints")
    ax.annotate(f"water  {d['P_bubble_kPa'][0]:.2f} kPa", (0, d["P_bubble_kPa"][0]), xytext=(0.04, d["P_bubble_kPa"][0] + 6), fontsize=8)
    ax.annotate(f"acetone  {d['P_bubble_kPa'][-1]:.2f} kPa", (1, d["P_bubble_kPa"][-1]), xytext=(0.62, d["P_bubble_kPa"][-1] - 14), fontsize=8)
    ax.set(xlim=(0, 1), ylim=(0, 1.1 * d["P_bubble_kPa"].max()), xlabel="Mole fraction of acetone, $x_1$ or $y_1$ [-]", ylabel="Pressure, P [kPa]",
           title=f"Question {qid.replace('Q','').replace('_','.')} – {model}: P–x–y, acetone(1)/water(2), T = {T:.2f} K ({T-273.15:.2f} °C)")
    ax.title.set_fontsize(10); ax.grid(alpha=.4); ax.legend(loc="upper left", fontsize=8); fig.tight_layout(); fig.savefig(path, dpi=300); plt.close(fig)

def plot_txy(d, model, P, qid, path):
    Tc = lambda a: a - 273.15
    fig, ax = plt.subplots(figsize=(7.5, 5.5), dpi=300)
    ax.plot(d["x1_or_y1"], Tc(d["T_bubble_K"]), "b-", lw=2, label="Bubble-point curve: t vs $x_1$ (liquid)")
    ax.plot(d["x1_or_y1"], Tc(d["T_dew_K"]), "r--", lw=2, label="Dew-point curve: t vs $y_1$ (vapour)")
    ax.plot([0, 1], [Tc(d["T_bubble_K"][0]), Tc(d["T_bubble_K"][-1])], "ko", ms=7, label="Pure-component $T^{sat}$ endpoints")
    ax.annotate(f"water  {Tc(d['T_bubble_K'][0]):.2f} °C", (0, Tc(d["T_bubble_K"][0])), xytext=(0.03, Tc(d["T_bubble_K"][0]) + 3.5), fontsize=8)
    ax.annotate(f"acetone  {Tc(d['T_bubble_K'][-1]):.2f} °C", (1, Tc(d["T_bubble_K"][-1])), xytext=(0.70, Tc(d["T_bubble_K"][-1]) - 5), fontsize=8)
    lo, hi = Tc(d["T_bubble_K"]).min(), Tc(d["T_bubble_K"]).max()
    ax.set(xlim=(0, 1), ylim=(lo - 8, hi + 8), xlabel="Mole fraction of acetone, $x_1$ or $y_1$ [-]", ylabel="Temperature, t [°C]",
           title=f"Question {qid.replace('Q','').replace('_','.')} – {model}: t–x–y, acetone(1)/water(2), P = {P:.2f} kPa")
    ax.title.set_fontsize(10); ax.grid(alpha=.4); ax.legend(loc="upper right", fontsize=8); fig.tight_layout(); fig.savefig(path, dpi=300); plt.close(fig)

def curve_checks(rep, m, d, kind, cond):
    """Independent checks on a computed diagram. kind 'P' (cond=T) or 'T' (cond=P)."""
    xs = d["x1_or_y1"]; y = d["y1_vapour_in_equilibrium_with_x1"]
    rep.check("0 <= y1 <= 1 everywhere", np.all((y >= 0) & (y <= 1)))
    rep.check("acetone is the more volatile species: y1 > x1 on 0<x1<1", np.all(y[1:-1] > xs[1:-1]))
    if kind == "P":
        e = 1e-12; Pw, _ = bd.bubble_P(m, cond, [e, 1 - e]); Pa, _ = bd.bubble_P(m, cond, [1 - e, e])
        rep.check("pure-component endpoints: bubble P(x1->0) = Psat(water), P(x1->1) = Psat(acetone)",
                  abs(Pw - psat("water", cond)) < 1e-6 and abs(Pa - psat("acetone", cond)) < 1e-6, f"(P_sat water {psat('water', cond):.4f}, acetone {psat('acetone', cond):.4f} kPa)")
        # dew pressure computed independently for y = y_bubble must reproduce the bubble pressure
        dd = [abs(bd.dew_P(m, cond, [yy, 1 - yy])[0] - P) for yy, P in zip(y[1:-1], d["P_bubble_kPa"][1:-1])]
        rep.check("dew P at y(bubble) reproduces bubble P (independent algorithm)", max(dd) < 1e-6, f"(max diff {max(dd):.2e} kPa)")
    else:
        rep.check("pure-component endpoints: T_bubble(0) = Tsat(water), T_bubble(1) = Tsat(acetone)",
                  abs(d["T_bubble_K"][0] - tsat("water", cond)) < 1e-9 and abs(d["T_bubble_K"][-1] - tsat("acetone", cond)) < 1e-9)
        dd = [abs(bd.dew_T(m, cond, [yy, 1 - yy])[0] - Tb) for yy, Tb in zip(y[1:-1:5], d["T_bubble_K"][1:-1:5])]
        rep.check("dew T at y(bubble) reproduces bubble T (independent algorithm)", max(dd) < 1e-5, f"(max diff {max(dd):.2e} K)")
    # Gibbs-Duhem: x1 dln g1 + x2 dln g2 = 0
    worst = 0
    for x1 in (0.1, 0.3, 0.5, 0.7, 0.9):
        Tg = cond if kind == "P" else float(np.interp(x1, xs, d["T_bubble_K"])); h = 1e-6
        dl = (np.log(m.gamma([x1 + h, 1 - x1 - h], Tg)) - np.log(m.gamma([x1 - h, 1 - x1 + h], Tg))) / (2 * h)
        worst = max(worst, abs(x1 * dl[0] + (1 - x1) * dl[1]))
    rep.check("Gibbs-Duhem consistency (model identity)", worst < 1e-6, f"(max |sum x dlnγ| = {worst:.1e})")
