# Question 12.13 - Wilson t-x-y diagram, acetone(1)/water(2), P = 101.33 kPa

*Generated 2026-11-11 by `Q12_13.py`*

Parameter source: Smith, Van Ness & Abbott, 6th ed. (Antoine constants, molar volumes, Wilson and NRTL parameters as supplied by the student)
> **Parameter provenance:** Smith, Van Ness & Abbott, 6th ed. (Antoine constants, molar volumes, Wilson and NRTL parameters as supplied by the student). Verified against book: **False**.


## Problem statement

Based on Eq. (10.5) and the Wilson equation prepare a t-x-y diagram for P = 101.33 kPa. (system: acetone(1)/water(2), chosen by the student from the table) 


## Model and governing equations

    Eq. (10.5), modified Raoult's law:  y_i P = x_i gamma_i P_i^sat      (ideal vapour phase)
    Antoine:  ln(P_i^sat/kPa) = A_i - B_i / (t/degC + C_i),  t = T - 273.15
    Wilson:  Lambda_ij = (V_j/V_i) exp(-a_ij/RT),  Lambda_ii = 1
             ln gamma_i = 1 - ln(sum_j x_j Lambda_ij) - sum_k [ x_k Lambda_ki / sum_j x_j Lambda_kj ]
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
| 0.1 | 342.592 | 69.4424 | 0.718339 | 370.308 | 0.00326884 |
| 0.2 | 337.738 | 64.5877 | 0.780572 | 367.196 | 0.0073084 |
| 0.3 | 335.821 | 62.6713 | 0.804285 | 363.748 | 0.0125011 |
| 0.4 | 334.628 | 61.4777 | 0.82004 | 359.874 | 0.0195543 |
| 0.5 | 333.653 | 60.5031 | 0.834454 | 355.438 | 0.0299571 |
| 0.6 | 332.736 | 59.5858 | 0.850325 | 350.23 | 0.0475415 |
| 0.7 | 331.818 | 58.6682 | 0.869927 | 343.909 | 0.08655 |
| 0.8 | 330.893 | 57.7432 | 0.896373 | 336.162 | 0.277229 |
| 0.9 | 330.01 | 56.8604 | 0.935402 | 330.79 | 0.811293 |
| 1 | 329.35 | 56.2002 | 1 | 329.35 | 1 |


## Final answers

- [RESULT diagram] T_sat water = 373.149 K
- [RESULT diagram] T_sat acetone = 329.35 K
- [RESULT diagram] diagram points = 51 -
Graph saved as plots/Q12_13/Q12_13_Wilson_txy.png


## Verification

- [PASS] 0 <= y1 <= 1 everywhere 
- [PASS] acetone is the more volatile species: y1 > x1 on 0<x1<1 
- [PASS] pure-component endpoints: T_bubble(0) = Tsat(water), T_bubble(1) = Tsat(acetone) 
- [PASS] dew T at y(bubble) reproduces bubble T (independent algorithm) (max diff 2.11e-11 K)
- [PASS] Gibbs-Duhem consistency (model identity) (max |sum x dlnγ| = 3.7e-11)

## Verification summary

Checks passed: 5/5  -> overall: ALL PASSED