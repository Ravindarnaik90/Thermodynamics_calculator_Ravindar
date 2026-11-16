# Question 12.15 - NRTL t-x-y diagram, acetone(1)/water(2), P = 101.33 kPa

*Generated 2026-11-11 by `Q12_15.py`*

Parameter source: Smith, Van Ness & Abbott, 6th ed. (Antoine constants, molar volumes, Wilson and NRTL parameters as supplied by the student)
> **Parameter provenance:** Smith, Van Ness & Abbott, 6th ed. (Antoine constants, molar volumes, Wilson and NRTL parameters as supplied by the student). Verified against book: **False**.


## Problem statement

Based on Eq. (10.5) and the NRTL equation prepare a t-x-y diagram for P = 101.33 kPa. (system: acetone(1)/water(2), chosen by the student from the table) 


## Model and governing equations

    Eq. (10.5), modified Raoult's law:  y_i P = x_i gamma_i P_i^sat      (ideal vapour phase)
    Antoine:  ln(P_i^sat/kPa) = A_i - B_i / (t/degC + C_i),  t = T - 273.15
    NRTL:  tau_ij = b_ij/RT,  G_ij = exp(-alpha_ij tau_ij),  tau_ii = 0, G_ii = 1
           ln gamma_i = sum_j x_j tau_ji G_ji / sum_k x_k G_ki
                        + sum_j [ x_j G_ij / sum_k x_k G_kj ] ( tau_ij - sum_m x_m tau_mj G_mj / sum_k x_k G_kj )
    Units: T in K inside all functions, P in kPa, R = 8.314 J/(mol K), a_ij and b_ij in J/mol

## Given data

- System = acetone(1) / water(2)  (binary chosen by the student from the table) 
- P = 101.33 kPa
- Composition grid = x1 = 0, 0.02, ..., 1.0 (51 points) 
- **NOTE:** Antoine for acetone used at 56.20 C (bubble/dew range at this P), outside stated range -13.0..55.0 C (extrapolation).

## Step 1 - pure-component saturation temperatures from inverse Antoine: T_sat = B/(A - ln P) - C + 273.15

- T_sat acetone = 329.35 K
- T_sat water = 373.149 K

## Step 2 - for every x1: BUBL T. Solve  sum x_i gamma_i(T,x) Psat_i(T) - P = 0  (scipy brentq, bracket 250-430 K); then y_i = x_i gamma_i Psat_i / P


## Step 3 - DEW T for the same compositions used as y1: solve  P_dew(T,y) - P = 0  (brentq around damped successive substitution for x)


## Computed diagram data (every 5th point; full data in the CSV)

| x1 | T_bubble [K] | T_bubble [degC] | y1 | T_dew(y1=x1 col) [K] | x1 on dew curve |
|---|---|---|---|---|---|
| 0 | 373.149 | 99.9991 | 0 | 373.149 | 0 |
| 0.1 | 342.112 | 68.9616 | 0.724441 | 370.314 | 0.00346276 |
| 0.2 | 337.67 | 64.52 | 0.780687 | 367.207 | 0.0076617 |
| 0.3 | 336.019 | 62.8691 | 0.801362 | 363.763 | 0.012946 |
| 0.4 | 334.858 | 61.7082 | 0.817104 | 359.891 | 0.019954 |
| 0.5 | 333.81 | 60.6599 | 0.833037 | 355.454 | 0.0300143 |
| 0.6 | 332.806 | 59.6562 | 0.850716 | 350.237 | 0.0465268 |
| 0.7 | 331.834 | 58.6838 | 0.87155 | 343.89 | 0.0824308 |
| 0.8 | 330.891 | 57.7408 | 0.898153 | 336.125 | 0.29191 |
| 0.9 | 330.008 | 56.8579 | 0.936237 | 330.836 | 0.805901 |
| 1 | 329.35 | 56.2002 | 1 | 329.35 | 1 |


## Final answers

- [RESULT diagram] T_sat water = 373.149 K
- [RESULT diagram] T_sat acetone = 329.35 K
- [RESULT diagram] diagram points = 51 -
Graph saved as plots/Q12_15/Q12_15_NRTL_txy.png


## Verification

- [PASS] 0 <= y1 <= 1 everywhere 
- [PASS] acetone is the more volatile species: y1 > x1 on 0<x1<1 
- [PASS] pure-component endpoints: T_bubble(0) = Tsat(water), T_bubble(1) = Tsat(acetone) 
- [PASS] dew T at y(bubble) reproduces bubble T (independent algorithm) (max diff 6.20e-12 K)
- [PASS] Gibbs-Duhem consistency (model identity) (max |sum x dlnγ| = 5.6e-11)

## Verification summary

Checks passed: 5/5  -> overall: ALL PASSED