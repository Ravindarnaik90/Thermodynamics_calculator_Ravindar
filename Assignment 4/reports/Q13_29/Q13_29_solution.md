# Problem 13.29 - Simultaneous removal of H2S and SO2, 723.15 K, 8 bar

*Generated 2026-11-11 by `Q13_29.py`*

Parameter source: Smith, Van Ness & Abbott, 6th ed. (Antoine constants, molar volumes, Wilson and NRTL parameters as supplied by the student)
> **Parameter provenance:** Smith, Van Ness & Abbott, 6th ed. (Antoine constants, molar volumes, Wilson and NRTL parameters as supplied by the student). Verified against book: **False**.


## Problem statement

Refineries often have both H2S and SO2 to dispose of. The reaction 2 H2S(g) + SO2(g) -> 3 S(s) + 2 H2O(g) suggests a means of getting rid of both. For reactants in the stoichiometric proportion, estimate the percent conversion of each reactant if the reaction comes to equilibrium at 723.15 K (450 degC) and 8 bar.


## Governing equations

    Reaction: 2 H2S + SO2 -> 3 S(s) + 2 H2O(g),  nu = (-2, -1, +3, +2)
    K = exp(-dG0(T)/RT) from Eq. (13.18)
    Solid sulfur is a pure phase (activity = 1), so it does not appear in K:
    K = y_H2O^2 / ( y_H2S^2 y_SO2 (P/P0) )          (ideal gas, gas-phase species only)

## Given data

- T = 723.15 K
- P = 8 bar
- feed basis = 2 mol H2S + 1 mol SO2 (stoichiometric) 

## Step 1 - equilibrium constant at T

- Cp(S, solid) at 298.15 K = R (A + B T + D/T^2) = 22.5972 J/(mol K)  (literature: ~22.7)
| species | nu_i | A | 10^3 B | 10^6 C | 10^-5 D | dHf298 [J/mol] | dGf298 [J/mol] |
|---|---|---|---|---|---|---|---|
| H2S | -2 | 3.931 | 1.49 | 0 | -0.232 | -20630 | -33560 |
| SO2 | -1 | 5.699 | 0.801 | 0 | -1.015 | -296830 | -300194 |
| S(s) | 3 | 4.114 | -1.728 | 0 | -0.783 | 0 | 0 |
| H2O | 2 | 3.47 | 1.45 | 0 | 0.121 | -241818 | -228572 |

- dA, dB, dC, dD = 5.721, -0.006065, 0, -62800 
- dH298 = sum nu_i dHf_i = -145.546 kJ/mol
- dG298 = sum nu_i dGf_i = -89.83 kJ/mol
- K(298.15 K) = exp(-dG298/RT0) = 5.47534e+15 
- dH(723.15 K) = -137.304 kJ/mol
- ln K(723.15 K) = 2.55853 
- K(723.15 K) = 12.9169 

## Step 2 - mole balance (basis 2 mol H2S + 1 mol SO2, extent e)

    gas phase: n_H2S = 2 - 2e,  n_SO2 = 1 - e,  n_H2O = 2e,   N_gas = 3 - e;   S(s) = 3e
    K = [2e/N]^2 / ( [(2-2e)/N]^2 [(1-e)/N] P ),  conversion of H2S = 2e/2 = e,  of SO2 = e/1 = e

## Step 3 - solve for e (scipy brentq, bracket 1e-12 < e < 1-1e-12)

- extent e = 0.766668 mol
- gas moles N = 2.23333 mol
- solid S formed = 2.30001 mol
- y_H2S, y_SO2, y_H2O = 0.20895, 0.10448, 0.68657 

## Final answers

- [RESULT -] K at 723.15 K = 12.9169 -
- [RESULT -] conversion of H2S = 76.6668 %
- [RESULT -] conversion of SO2 = 76.6668 %

## Verification

- [PASS] ln K(298.15) = -dG298/(R T0) 
- [PASS] analytic integrals agree with numerical quadrature (max |diff in ln K| = 3.6e-15)
- [PASS] van't Hoff relation d ln K/dT = dH/(R T^2) (max abs deviation 3.1e-12)
- [PASS] sum of gas mole fractions = 1 
- [PASS] equilibrium relation y_W^2/(y_H2S^2 y_SO2 P) = K 
- [PASS] S balance: S in = S out 
- [PASS] O balance: O in = O out 
- [PASS] H balance: H in = H out 
- [PASS] independent fsolve gives the same extent (e = 0.7666683385, residual 4.4e-16)
- [PASS] 0 <= e <= 1 
- **NOTE:** Elemental sulfur is treated as a pure solid, as stated in the problem. At 723 K real sulfur is a liquid/vapour; the solid heat capacity (Table C.2) is therefore extrapolated beyond its stated range.
Graph saved as plots/Q13_29/Q13_29.png


## Verification summary

Checks passed: 10/10  -> overall: ALL PASSED