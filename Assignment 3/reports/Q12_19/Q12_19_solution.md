# Question 12.19 - NRTL BUBL T / DEW T / flash / azeotrope, acetone(1)/water(2), P = 101.33 kPa

*Generated 2026-11-11 by `Q12_19.py`*

Parameter source: Smith, Van Ness & Abbott, 6th ed. (Antoine constants, molar volumes, Wilson and NRTL parameters as supplied by the student)
> **Parameter provenance:** Smith, Van Ness & Abbott, 6th ed. (Antoine constants, molar volumes, Wilson and NRTL parameters as supplied by the student). Verified against book: **False**.


## Problem statement

Work 12.18 for the NRTL equation. (system: acetone(1)/water(2), chosen by the student from the table) 


## Model and governing equations

    Eq. (10.5), modified Raoult's law:  y_i P = x_i gamma_i P_i^sat      (ideal vapour phase)
    Antoine:  ln(P_i^sat/kPa) = A_i - B_i / (t/degC + C_i),  t = T - 273.15
    NRTL:  tau_ij = b_ij/RT,  G_ij = exp(-alpha_ij tau_ij),  tau_ii = 0, G_ii = 1
           ln gamma_i = sum_j x_j tau_ji G_ji / sum_k x_k G_ki
                        + sum_j [ x_j G_ij / sum_k x_k G_kj ] ( tau_ij - sum_m x_m tau_mj G_mj / sum_k x_k G_kj )
    Units: T in K inside all functions, P in kPa, R = 8.314 J/(mol K), a_ij and b_ij in J/mol

## Given data

- System = acetone(1)/water(2) 
- P = 101.33 kPa
- (a) x1 = 0.3 
- (b) y1 = 0.3 
- (c) z1 = 0.3 and T = (T_bubble + T_dew)/2

## (a) BUBL T:  solve  f(T) = sum x_i gamma_i(T) Psat_i(T) - P = 0  (brentq, 250-430 K) ;  y_i = x_i gamma_i Psat_i / P

- T = 336.019 K
- t = 62.8691 deg C
- P_sat (Antoine) = [acetone: 126.551077, water: 22.809922] kPa
NRTL tau_ij = b_ij/RT:

|  | acetone | water |
|---|---|---|
| acetone | 0 | 0.945748 |
| water | 1.79453 | 0 |

NRTL G_ij = exp(-alpha tau_ij):

|  | acetone | water |
|---|---|---|
| acetone | 1 | 0.603317 |
| water | 0.383347 | 1 |

- gamma at x=[0.3, 0.7] = [acetone: 2.138846, water: 1.260605] 
- partial pressures at T_bubble = [acetone: 81.201996, water: 20.128004] kPa
- [RESULT a] T_bubble = 336.019 K
- [RESULT a] T_bubble = 62.8691 degC
- [RESULT a] y_acetone = 0.801362 -
- [RESULT a] y_water = 0.198638 -
- **NOTE:** Antoine for acetone used at 62.87 C , outside stated range -13.0..55.0 C (extrapolation).
- [PASS] sum y = 1 
- [PASS] Eq.(10.5) residual (5.7e-14 kPa)
- [PASS] bubble equation residual sum x gamma Psat - P = 0 
- [PASS] independent fsolve BUBL T agrees (336.019119 vs 336.019119 K)

## (b) DEW T:  solve  P_dew(T, y) - P = 0  (brentq; inner DEW P by damped successive substitution)

- gamma at dew liquid = [acetone: 8.132988, water: 1.000708] 
- [RESULT b] T_dew = 363.763 K
- [RESULT b] T_dew = 90.6131 degC
- [RESULT b] x_acetone = 0.012946 -
- [RESULT b] x_water = 0.987054 -
- **NOTE:** Antoine for acetone used at 90.61 C , outside stated range -13.0..55.0 C (extrapolation).
- [PASS] sum x = 1, sum y = 1 
- [PASS] Eq.(10.5) residual (1.2e-11 kPa)
- [PASS] dew equation residual sum y/(gamma Psat) - 1/P = 0 
- [PASS] independent fsolve DEW T agrees (363.763096 vs 363.763096 K)
- [PASS] T_bubble < T_dew 
- T for flash = (T_bubble + T_dew)/2 = 349.891 K

## (c) P,T-flash

Material balance z_i = V y_i + (1-V) x_i with y_i = K_i x_i,  K_i = gamma_i Psat_i / P  =>  x_i = z_i/(1+V(K_i-1)).  Rachford-Rice: sum z_i(K_i-1)/(1+V(K_i-1)) = 0, solved with brentq; gamma updated with x until x converges (tol 1e-13).

- flash P = 101.33 kPa
- flash T = 349.891 K
- z = [np.float64(0.3), np.float64(0.7)] 
- P_bubble(z) and P_dew(z) at flash T = 160.046 kPa, 58.7618 kPa 
- phase state = two-phase 
- outer iterations = 29 
- [RESULT c] flash pressure = 101.33 kPa
- [RESULT c] flash temperature = 349.891 K
- [RESULT c] vapour fraction V = 0.451689 mol/mol
- [RESULT c] x_acetone = 0.0478979 -
- [RESULT c] x_water = 0.952102 -
- [RESULT c] y_acetone = 0.60603 -
- [RESULT c] y_water = 0.39397 -
- [PASS] flash converged 
- [PASS] 0 <= V <= 1 
- K_i = gamma_i Psat_i / P = [acetone: 12.652534, water: 0.413789] 
- [PASS] sum x = 1 and sum y = 1 
- [PASS] 0 <= x_i, y_i <= 1 
- [PASS] component material balance z_i = V y_i + (1-V) x_i (max err 1.1e-16)
- [PASS] equilibrium y_i = K_i x_i  (Eq. 10.5 residual) 
- [PASS] independent fsolve flash (different algorithm/start) agrees (V 0.45168861 vs 0.45168861)

## (d) Azeotrope at P:  condition x1 = y1 together with Eq.(10.5)  <=>  alpha12 = gamma1 Psat1/(gamma2 Psat2) = 1 at one x1 (along the bubble curve)

- relative volatility alpha12 along bubble curve: min / max = 1.15847 / 32.7587 
- [RESULT d] azeotrope exists = 0 (0 = no)
- closest approach |y1-x1| (excluding pure-component ends) = 0.000139542 
- at x1 = 0.999 
- [PASS] no azeotrope: alpha12 > 1 for all x1 (independent criterion) (min alpha12 = 1.1585)
Interpretation: an azeotrope requires alpha12 = 1 at some composition. alpha12 stays above 1 along the whole bubble curve, so no azeotrope exists at this pressure for the model.


## Verification summary

Checks passed: 17/17  -> overall: ALL PASSED