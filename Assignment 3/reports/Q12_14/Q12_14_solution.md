# Question 12.14 - NRTL P-x-y diagram, acetone(1)/water(2), T = 333.15 K

*Generated 2026-11-11 by `Q12_14.py`*

Parameter source: Smith, Van Ness & Abbott, 6th ed. (Antoine constants, molar volumes, Wilson and NRTL parameters as supplied by the student)
> **Parameter provenance:** Smith, Van Ness & Abbott, 6th ed. (Antoine constants, molar volumes, Wilson and NRTL parameters as supplied by the student). Verified against book: **False**.


## Problem statement

Based on Eq. (10.5) and the NRTL equation prepare a P-x-y diagram for T = 333.15 K (60 degC). (system: acetone(1)/water(2), chosen by the student from the table) 


## Model and governing equations

    Eq. (10.5), modified Raoult's law:  y_i P = x_i gamma_i P_i^sat      (ideal vapour phase)
    Antoine:  ln(P_i^sat/kPa) = A_i - B_i / (t/degC + C_i),  t = T - 273.15
    NRTL:  tau_ij = b_ij/RT,  G_ij = exp(-alpha_ij tau_ij),  tau_ii = 0, G_ii = 1
           ln gamma_i = sum_j x_j tau_ji G_ji / sum_k x_k G_ki
                        + sum_j [ x_j G_ij / sum_k x_k G_kj ] ( tau_ij - sum_m x_m tau_mj G_mj / sum_k x_k G_kj )
    Units: T in K inside all functions, P in kPa, R = 8.314 J/(mol K), a_ij and b_ij in J/mol

## Given data

- System = acetone(1) / water(2)  (binary chosen by the student from the table) 
- T = 333.15 K (60 degC)
- Composition grid = x1 = 0, 0.02, ..., 1.0 (51 points) 
- **NOTE:** Antoine for acetone used at 60.00 C , outside stated range -13.0..55.0 C (extrapolation).

## Step 1 - pure-component vapour pressures and model parameters at T

- T = 333.15 K
- t = 60 deg C
- P_sat (Antoine) = [acetone: 115.155691, water: 20.007045] kPa
NRTL tau_ij = b_ij/RT:

|  | acetone | water |
|---|---|---|
| acetone | 0 | 0.953892 |
| water | 1.80998 | 0 |

NRTL G_ij = exp(-alpha tau_ij):

|  | acetone | water |
|---|---|---|
| acetone | 1 | 0.600697 |
| water | 0.380194 | 1 |

- gamma at x=[0.5, 0.5] = [acetone: 1.434729, water: 1.641448] 

## Step 2 - for every x1: BUBL P  (P = sum x_i gamma_i Psat_i,  y_i = x_i gamma_i Psat_i / P)


## Step 3 - DEW P computed independently for the same compositions used as y1 (dew curve)


## Computed diagram data (every 5th point; full data in the CSV)

| x1 | P_bubble [kPa] | y1 | P_dew(y1=x1 col) [kPa] | x1 on dew curve |
|---|---|---|---|---|
| 0 | 20.007 | 0 | 20.007 | 0 |
| 0.1 | 75.0452 | 0.750535 | 22.1902 | 0.00181113 |
| 0.2 | 86.7608 | 0.791049 | 24.9069 | 0.00416292 |
| 0.3 | 91.7294 | 0.807222 | 28.3793 | 0.00734234 |
| 0.4 | 95.4815 | 0.820324 | 32.9717 | 0.0118859 |
| 0.5 | 99.0288 | 0.834187 | 39.3245 | 0.0189327 |
| 0.6 | 102.548 | 0.850172 | 48.6689 | 0.0314352 |
| 0.7 | 106.065 | 0.869739 | 63.6475 | 0.0605849 |
| 0.8 | 109.556 | 0.895681 | 89.5343 | 0.249785 |
| 0.9 | 112.846 | 0.934105 | 110.026 | 0.813687 |
| 1 | 115.156 | 1 | 115.156 | 1 |


## Final answers

- [RESULT diagram] P_sat water = 20.007 kPa
- [RESULT diagram] P_sat acetone = 115.156 kPa
- [RESULT diagram] diagram points = 51 -
Graph saved as plots/Q12_14/Q12_14_NRTL_Pxy.png


## Verification

- [PASS] 0 <= y1 <= 1 everywhere 
- [PASS] acetone is the more volatile species: y1 > x1 on 0<x1<1 
- [PASS] pure-component endpoints: bubble P(x1->0) = Psat(water), P(x1->1) = Psat(acetone) (P_sat water 20.0070, acetone 115.1557 kPa)
- [PASS] dew P at y(bubble) reproduces bubble P (independent algorithm) (max diff 4.26e-14 kPa)
- [PASS] Gibbs-Duhem consistency (model identity) (max |sum x dlnγ| = 4.2e-11)

## Verification summary

Checks passed: 5/5  -> overall: ALL PASSED