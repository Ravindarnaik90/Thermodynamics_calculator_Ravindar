# Question 12.20 - Wilson BUBL P / DEW P / flash, acetone(1)/methanol(2)/water(3), T = 338.15 K

*Generated 2026-11-11 by `Q12_20.py`*

Parameter source: Smith, Van Ness & Abbott, 6th ed. (Antoine constants, molar volumes, Wilson and NRTL parameters as supplied by the student)
> **Parameter provenance:** Smith, Van Ness & Abbott, 6th ed. (Antoine constants, molar volumes, Wilson and NRTL parameters as supplied by the student). Verified against book: **False**.


## Problem statement

Acetone(1)/methanol(2)/water(3). Wilson equation: (a) BUBL P: T = 338.15 K, x1 = 0.3, x2 = 0.4; (b) DEW P: T = 338.15 K, y1 = 0.3, y2 = 0.4; (c) P,T-flash: T = 338.15 K, P = (P_bubble+P_dew)/2, z1 = 0.3, z2 = 0.4.


## Model and governing equations

    Eq. (10.5), modified Raoult's law:  y_i P = x_i gamma_i P_i^sat      (ideal vapour phase)
    Antoine:  ln(P_i^sat/kPa) = A_i - B_i / (t/degC + C_i),  t = T - 273.15
    Wilson:  Lambda_ij = (V_j/V_i) exp(-a_ij/RT),  Lambda_ii = 1
             ln gamma_i = 1 - ln(sum_j x_j Lambda_ij) - sum_k [ x_k Lambda_ki / sum_j x_j Lambda_kj ]
    Units: T in K inside all functions, P in kPa, R = 8.314 J/(mol K), a_ij and b_ij in J/mol

## Given data

- System = acetone(1)/methanol(2)/water(3) 
- T = 338.15 K (65 degC)
- (a) x1, x2, x3 = 1 - x1 - x2 = [np.float64(0.3), np.float64(0.4), np.float64(0.3)] 
- (b) y1, y2, y3 = 1 - y1 - y2 = [np.float64(0.3), np.float64(0.4), np.float64(0.3)] 
- (c) z1, z2, z3 = 1 - z1 - z2 = [np.float64(0.3), np.float64(0.4), np.float64(0.3)] and P = (P_bubble + P_dew)/2
- **NOTE:** Antoine for acetone used at 65.00 C , outside stated range -13.0..55.0 C (extrapolation).

## (a) BUBL P:  P = sum x_i gamma_i Psat_i ;  y_i = x_i gamma_i Psat_i / P   (complete ternary model, all three binary pairs)

- T = 338.15 K
- t = 65 deg C
- P_sat (Antoine) = [acetone: 135.575914, methanol: 102.530309, water: 25.101354] kPa
Wilson Lambda_ij = (Vj/Vi) exp(-a_ij/RT):

|  | acetone | methanol | water |
|---|---|---|---|
| acetone | 1 | 0.699994 | 0.158142 |
| methanol | 0.762903 | 1 | 0.220474 |
| water | 0.474294 | 1.92089 | 1 |

- gamma at x=[0.3, 0.4, 0.3] = [acetone: 1.539365, methanol: 0.977854, water: 1.431242] 
- partial pressures x_i gamma_i Psat_i = [acetone: 62.610247, methanol: 40.103866, water: 10.777834] kPa
- [RESULT a] P_bubble = 113.492 kPa
- [RESULT a] y_acetone = 0.551671 -
- [RESULT a] y_methanol = 0.353363 -
- [RESULT a] y_water = 0.0949656 -
- [PASS] sum x = 1, sum y = 1 
- [PASS] Eq.(10.5) residual (7.1e-15 kPa)
- [PASS] reverse: DEW P at y_bubble returns same P and x 

## (b) DEW P:  P = 1/sum(y_i/(gamma_i Psat_i)) ;  x_i = y_i P/(gamma_i Psat_i)  (iterated)

- gamma at dew liquid = [acetone: 3.668866, methanol: 1.048805, water: 1.077313] 
- [RESULT b] P_dew = 64.8647 kPa
- [RESULT b] x_acetone = 0.0391215 -
- [RESULT b] x_methanol = 0.24128 -
- [RESULT b] x_water = 0.719599 -
- [PASS] DEW P converged 
- [PASS] sum x = 1, sum y = 1 
- [PASS] Eq.(10.5) residual (4.3e-12 kPa)
- [PASS] independent fsolve DEW P agrees (64.864676 vs 64.864676 kPa)
- [PASS] P_dew < P_bubble 
- P for flash = (P_bubble + P_dew)/2 = 89.1783 kPa

## (c) P,T-flash

Material balance z_i = V y_i + (1-V) x_i with y_i = K_i x_i,  K_i = gamma_i Psat_i / P  =>  x_i = z_i/(1+V(K_i-1)).  Rachford-Rice: sum z_i(K_i-1)/(1+V(K_i-1)) = 0, solved with brentq; gamma updated with x until x converges (tol 1e-13).

- flash P = 89.1783 kPa
- flash T = 338.15 K
- z = [np.float64(0.3), np.float64(0.4), np.float64(0.3)] 
- P_bubble(z) and P_dew(z) at flash T = 113.492 kPa, 64.8647 kPa 
- phase state = two-phase 
- outer iterations = 29 
- [RESULT c] flash pressure = 89.1783 kPa
- [RESULT c] flash temperature = 338.15 K
- [RESULT c] vapour fraction V = 0.611082 mol/mol
- [RESULT c] x_acetone = 0.118795 -
- [RESULT c] x_methanol = 0.375558 -
- [RESULT c] x_water = 0.505648 -
- [RESULT c] y_acetone = 0.415327 -
- [RESULT c] y_methanol = 0.415556 -
- [RESULT c] y_water = 0.169117 -
- [PASS] flash converged 
- [PASS] 0 <= V <= 1 
- K_i = gamma_i Psat_i / P = [acetone: 3.496174, methanol: 1.106504, water: 0.334456] 
- [PASS] sum x = 1 and sum y = 1 
- [PASS] 0 <= x_i, y_i <= 1 
- [PASS] component material balance z_i = V y_i + (1-V) x_i (max err 5.6e-17)
- [PASS] equilibrium y_i = K_i x_i  (Eq. 10.5 residual) 
- [PASS] independent fsolve flash (different algorithm/start) agrees (V 0.61108169 vs 0.61108169)
(No azeotrope part is asked in this question.)


## Verification summary

Checks passed: 15/15  -> overall: ALL PASSED