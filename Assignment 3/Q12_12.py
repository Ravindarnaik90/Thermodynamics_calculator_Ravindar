"""Q12_12  --  Question 12.12  (Wilson).  Run:  python Q12_12.py
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

QID, MODEL = "Q12_12", "Wilson"
STATEMENT = "Based on Eq. (10.5) and the Wilson equation prepare a P-x-y diagram for T = 333.15 K (60 degC). (system: acetone(1)/water(2), chosen by the student from the table) "
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
    rep = Report(QID, "Question 12.12 - Wilson P-x-y diagram, acetone(1)/water(2), T = 333.15 K", HERE)
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
