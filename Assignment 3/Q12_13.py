"""Q12_13  --  Question 12.13  (Wilson).  Run:  python Q12_13.py
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

QID, MODEL = "Q12_13", "Wilson"
STATEMENT = "Based on Eq. (10.5) and the Wilson equation prepare a t-x-y diagram for P = 101.33 kPa. (system: acetone(1)/water(2), chosen by the student from the table) "
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
    rep = Report(QID, "Question 12.13 - Wilson t-x-y diagram, acetone(1)/water(2), P = 101.33 kPa", HERE)
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
