"""Q13_20 -- Problem 13.20 (ammonia synthesis equilibrium).   Run:  python Q13_20.py"""
import os, sys
import numpy as np
from scipy.optimize import brentq
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from common.report import Report
from common.reaction import R, T0, P0, dH298, dG298, dH, lnK, describe_reaction, thermo_checks

QID = "Q13_20"
SOURCE = "Smith, Van Ness & Abbott, 6th ed.: Table C.1 (heat capacities), Table C.4 (formation properties), Table B.1 (critical constants)"
STATEMENT = ("For the ammonia synthesis reaction 1/2 N2(g) + 3/2 H2(g) -> NH3(g), with a stoichiometric feed of H2 and N2: "
             "(a) What is the equilibrium mole fraction of ammonia at 1 bar and 300 K? (b) At what temperature does the "
             "equilibrium mole fraction of ammonia equal 0.50 for a pressure of 1 bar? (c) At what temperature is y_NH3 = 0.50 "
             "for a pressure of 100 bar, assuming the equilibrium mixture is an ideal gas? (d) At what temperature is "
             "y_NH3 = 0.50 for 100 bar, assuming the equilibrium mixture is an ideal solution of gases?")
# critical constants and acentric factors (Table B.1): Tc [K], Pc [bar], omega
CRIT = {"NH3": (405.7, 112.8, 0.253), "N2": (126.2, 34.00, 0.038), "H2": (33.19, 13.13, -0.216)}
C_STOICH = 3**1.5 / 16          # (1/4)^(1/2) (3/4)^(3/2) = 0.32476

def ln_phi(sp, T, P):
    """Pure-species fugacity coefficient, Pitzer virial correlation: ln phi = (Pr/Tr)(B0 + omega B1)"""
    Tc, Pc, w = CRIT[sp]; Tr, Pr = T / Tc, P / Pc
    B0 = 0.083 - 0.422 / Tr**1.6; B1 = 0.139 - 0.172 / Tr**4.2
    return Pr / Tr * (B0 + w * B1)

def y_from_KP(KP):
    """solve y / (C (1-y)^2) = K P for y (monotonic in y -> unique root in (0,1))"""
    return brentq(lambda y: y / (C_STOICH * (1 - y)**2) - KP, 1e-15, 1 - 1e-15, xtol=1e-15)

def main():
    rep = Report(QID, "Problem 13.20 - Ammonia synthesis equilibrium", HERE, source=SOURCE)
    rep.h("Problem statement"); rep.p(STATEMENT)
    rep.h("Governing equations")
    rep.eq("Reaction: 1/2 N2 + 3/2 H2 -> NH3,  nu = (-1/2, -3/2, +1);  K = exp(-dG0(T)/RT)  [Eq. 13.18]")
    rep.eq("Stoichiometric feed (3 H2 : 1 N2):  y_N2 = (1-y)/4,  y_H2 = 3(1-y)/4,  y_NH3 = y")
    rep.eq("Ideal gas:  K = y / [ (y_N2)^(1/2) (y_H2)^(3/2) (P/P0) ] = y / [ C (1-y)^2 (P/P0) ],   C = 3^1.5/16 = 0.32476")
    rep.eq("Ideal solution of gases (Lewis-Randall): K = y / [C (1-y)^2 P] * Phi,   Phi = phi_NH3 / (phi_N2^(1/2) phi_H2^(3/2)),  phi_i = pure-species values at (T, P)")
    rep.eq("Pure-gas phi_i from the Pitzer virial correlation: ln phi = (Pr/Tr)(B0 + omega B1),  B0 = 0.083 - 0.422/Tr^1.6,  B1 = 0.139 - 0.172/Tr^4.2")
    nu = {"N2": -0.5, "H2": -1.5, "NH3": 1.0}
    rep.h("Step 1 - reaction thermodynamics")
    describe_reaction(rep, nu, 300.0)
    # ---------------- (a)
    rep.h("(a) y_NH3 at 1 bar and 300 K")
    K300 = float(np.exp(lnK(nu, 300.0))); ya = y_from_KP(K300 * 1.0)
    rep.kv("K(300 K)", K300); rep.kv("K P = y/(C (1-y)^2)", K300); rep.kv("y_NH3", ya)
    # ---------------- (b),(c): required K for y = 0.5
    ytarget = 0.5; target = lambda P: np.log(ytarget / (C_STOICH * (1 - ytarget)**2 * P))
    rep.h("(b) temperature for y_NH3 = 0.50 at 1 bar (ideal gas)")
    rep.eq("ln K(T) = ln[ 0.5 / (C x 0.25 x P) ]; solved for T with brentq on 250-1500 K")
    Tb = brentq(lambda T: lnK(nu, T) - target(1.0), 250.0, 1500.0, xtol=1e-12)
    rep.kv("required K", float(np.exp(target(1.0)))); rep.kv("T", Tb, "K"); rep.kv("t", Tb - 273.15, "degC")
    rep.h("(c) temperature for y_NH3 = 0.50 at 100 bar (ideal gas)")
    Tc_ = brentq(lambda T: lnK(nu, T) - target(100.0), 250.0, 1500.0, xtol=1e-12)
    rep.kv("required K", float(np.exp(target(100.0)))); rep.kv("T", Tc_, "K"); rep.kv("t", Tc_ - 273.15, "degC")
    rep.h("(d) temperature for y_NH3 = 0.50 at 100 bar (ideal solution of gases)")
    lnPhi = lambda T, P: ln_phi("NH3", T, P) - 0.5 * ln_phi("N2", T, P) - 1.5 * ln_phi("H2", T, P)
    Td = brentq(lambda T: lnK(nu, T) - (target(100.0) + lnPhi(T, 100.0)), 250.0, 1500.0, xtol=1e-12)
    rep.kv("ln Phi at solution", lnPhi(Td, 100.0)); rep.kv("Phi", float(np.exp(lnPhi(Td, 100.0))))
    rep.kv("phi_NH3, phi_N2, phi_H2", ", ".join(f"{np.exp(ln_phi(s, Td, 100.0)):.4f}" for s in ("NH3", "N2", "H2")))
    rep.kv("T", Td, "K"); rep.kv("t", Td - 273.15, "degC")
    rep.h("Final answers")
    rep.result("a", "y NH3 at 1 bar, 300 K", ya, "-", "converged (brentq)")
    rep.result("b", "T for y = 0.5 at 1 bar", Tb, "K", "converged (brentq)"); rep.result("b", "T for y = 0.5 at 1 bar", Tb - 273.15, "degC", "converged (brentq)")
    rep.result("c", "T for y = 0.5 at 100 bar, ideal gas", Tc_, "K", "converged (brentq)"); rep.result("c", "T for y = 0.5 at 100 bar, ideal gas", Tc_ - 273.15, "degC", "converged (brentq)")
    rep.result("d", "T for y = 0.5 at 100 bar, ideal solution", Td, "K", "converged (brentq)"); rep.result("d", "T for y = 0.5 at 100 bar, ideal solution", Td - 273.15, "degC", "converged (brentq)")
    rep.h("Verification")
    thermo_checks(rep, nu, [300.0, 400.0, 600.0, 800.0])
    def y_of(T, P, ideal=True):
        l = lnK(nu, T) + np.log(1.0) - (0 if ideal else lnPhi(T, P)); return y_from_KP(np.exp(l) * P)
    rep.check("(a) y satisfies K P = y/(C (1-y)^2)", abs(ya / (C_STOICH * (1 - ya)**2) - K300) / K300 < 1e-10)
    rep.check("(b) forward: y(T_b, 1 bar) = 0.50", abs(y_of(Tb, 1.0) - 0.5) < 1e-9)
    rep.check("(c) forward: y(T_c, 100 bar, ideal gas) = 0.50", abs(y_of(Tc_, 100.0) - 0.5) < 1e-9)
    rep.check("(d) forward: y(T_d, 100 bar, ideal solution) = 0.50", abs(y_of(Td, 100.0, ideal=False) - 0.5) < 1e-9)
    rep.check("composition sums to 1 at (a): y_N2 + y_H2 + y_NH3 = 1", abs((1 - ya) / 4 + 3 * (1 - ya) / 4 + ya - 1) < 1e-14)
    rep.check("higher pressure raises the temperature at which y = 0.5 (Le Chatelier): T_c > T_b", Tc_ > Tb)
    rep.check("virial correlation applicable: all Pr/Tr values moderate at solution", all(100.0 / CRIT[s][1] / (Td / CRIT[s][0]) < 1.5 for s in CRIT))
    # ---------------- graph
    Ts = np.linspace(350, 1000, 260)
    y1 = [y_of(t, 1.0) for t in Ts]; y100 = [y_of(t, 100.0) for t in Ts]; y100s = [y_of(t, 100.0, ideal=False) for t in Ts]
    fig, ax = plt.subplots(figsize=(7.2, 4.6), dpi=300)
    ax.plot(Ts, y1, lw=2.2, color="#1f4e79", label="1 bar (ideal gas)"); ax.plot(Ts, y100, lw=2.2, color="#2e8b57", label="100 bar, ideal gas")
    ax.plot(Ts, y100s, lw=2.2, ls="--", color="#e67e22", label="100 bar, ideal solution of gases"); ax.axhline(0.5, color="grey", lw=1, ls=":")
    for t, c in ((Tb, "#1f4e79"), (Tc_, "#2e8b57"), (Td, "#e67e22")): ax.plot(t, 0.5, "o", color=c, ms=7)
    ax.set(xlabel="Temperature, T [K]", ylabel="Equilibrium mole fraction NH$_3$ [-]", ylim=(0, 1.02), title="Problem 13.20: ammonia equilibrium, stoichiometric N$_2$/H$_2$ feed")
    ax.title.set_fontsize(10); ax.grid(alpha=.35); ax.legend(); fig.tight_layout()
    fig.savefig(os.path.join(rep.dirs["plots"], QID + ".png"), dpi=300); plt.close(fig)
    rep.p(f"Graph saved as plots/{QID}/{QID}.png")
    ok = rep.save(); sys.exit(0 if ok else 1)

if __name__ == "__main__": main()
