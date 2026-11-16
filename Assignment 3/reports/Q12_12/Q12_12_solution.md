# Question 12.12 - Wilson P-x-y diagram, acetone(1)/water(2), T = 333.15 K

*Generated 2026-11-11 by `Q12_12.py`*

Parameter source: Smith, Van Ness & Abbott, 6th ed. (Antoine constants, molar volumes, Wilson and NRTL parameters as supplied by the student)
> **Parameter provenance:** Smith, Van Ness & Abbott, 6th ed. (Antoine constants, molar volumes, Wilson and NRTL parameters as supplied by the student). Verified against book: **False**.


## Problem statement

Based on Eq. (10.5) and the Wilson equation prepare a P-x-y diagram for T = 333.15 K (60 degC). (system: acetone(1)/water(2), chosen by the student from the table) 


## Model and governing equations

    Eq. (10.5), modified Raoult's law:  y_i P = x_i gamma_i P_i^sat      (ideal vapour phase)
    Antoine:  ln(P_i^sat/kPa) = A_i - B_i / (t/degC + C_i),  t = T - 273.15
    Wilson:  Lambda_ij = (V_j/V_i) exp(-a_ij/RT),  Lambda_ii = 1
             ln gamma_i = 1 - ln(sum_j x_j Lambda_ij) - sum_k [ x_k Lambda_ki / sum_j x_j Lambda_kj ]
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
Wilson Lambda_ij = (Vj/Vi) exp(-a_ij/RT):

|  | acetone | water |
|---|---|---|
| acetone | 1 | 0.157116 |
| water | 0.45919 | 1 |

- gamma at x=[0.5, 0.5] = [acetone: 1.445291, water: 1.639140] 

## Step 2 - for every x1: BUBL P  (P = sum x_i gamma_i Psat_i,  y_i = x_i gamma_i Psat_i / P)


## Step 3 - DEW P computed independently for the same compositions used as y1 (dew curve)


## Computed diagram data (every 5th point; full data in the CSV)

| x1 | P_bubble [kPa] | y1 | P_dew(y1=x1 col) [kPa] | x1 on dew curve |
|---|---|---|---|---|
| 0 | 20.007 | 0 | 20.007 | 0 |
| 0.1 | 74.0826 | 0.747261 | 22.1905 | 0.00179955 |
| 0.2 | 86.9645 | 0.792286 | 24.9075 | 0.00414779 |
| 0.3 | 92.6433 | 0.810403 | 28.3801 | 0.00734225 |
| 0.4 | 96.3984 | 0.823113 | 32.9724 | 0.0119447 |
| 0.5 | 99.6139 | 0.835393 | 39.3233 | 0.0191618 |
| 0.6 | 102.772 | 0.849651 | 48.6588 | 0.0321603 |
| 0.7 | 106.052 | 0.868144 | 63.5949 | 0.0629659 |
| 0.8 | 109.459 | 0.89416 | 89.391 | 0.235773 |
| 0.9 | 112.757 | 0.933731 | 110.07 | 0.817899 |
| 1 | 115.156 | 1 | 115.156 | 1 |


## Final answers

- [RESULT diagram] P_sat water = 20.007 kPa
- [RESULT diagram] P_sat acetone = 115.156 kPa
- [RESULT diagram] diagram points = 51 -
Graph saved as plots/Q12_12/Q12_12_Wilson_Pxy.png


## Verification

- [PASS] 0 <= y1 <= 1 everywhere 
- [PASS] acetone is the more volatile species: y1 > x1 on 0<x1<1 
- [PASS] pure-component endpoints: bubble P(x1->0) = Psat(water), P(x1->1) = Psat(acetone) (P_sat water 20.0070, acetone 115.1557 kPa)
- [PASS] dew P at y(bubble) reproduces bubble P (independent algorithm) (max diff 4.26e-14 kPa)
- [PASS] Gibbs-Duhem consistency (model identity) (max |sum x dlnγ| = 7.6e-11)

## Verification summary

Checks passed: 5/5  -> overall: ALL PASSED