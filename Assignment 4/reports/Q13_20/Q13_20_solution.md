# Problem 13.20 - Ammonia synthesis equilibrium

*Generated 2026-11-11 by `Q13_20.py`*

Parameter source: Smith, Van Ness & Abbott, 6th ed. (Antoine constants, molar volumes, Wilson and NRTL parameters as supplied by the student)
> **Parameter provenance:** Smith, Van Ness & Abbott, 6th ed. (Antoine constants, molar volumes, Wilson and NRTL parameters as supplied by the student). Verified against book: **False**.


## Problem statement

For the ammonia synthesis reaction 1/2 N2(g) + 3/2 H2(g) -> NH3(g), with a stoichiometric feed of H2 and N2: (a) What is the equilibrium mole fraction of ammonia at 1 bar and 300 K? (b) At what temperature does the equilibrium mole fraction of ammonia equal 0.50 for a pressure of 1 bar? (c) At what temperature is y_NH3 = 0.50 for a pressure of 100 bar, assuming the equilibrium mixture is an ideal gas? (d) At what temperature is y_NH3 = 0.50 for 100 bar, assuming the equilibrium mixture is an ideal solution of gases?


## Governing equations

    Reaction: 1/2 N2 + 3/2 H2 -> NH3,  nu = (-1/2, -3/2, +1);  K = exp(-dG0(T)/RT)  [Eq. 13.18]
    Stoichiometric feed (3 H2 : 1 N2):  y_N2 = (1-y)/4,  y_H2 = 3(1-y)/4,  y_NH3 = y
    Ideal gas:  K = y / [ (y_N2)^(1/2) (y_H2)^(3/2) (P/P0) ] = y / [ C (1-y)^2 (P/P0) ],   C = 3^1.5/16 = 0.32476
    Ideal solution of gases (Lewis-Randall): K = y / [C (1-y)^2 P] * Phi,   Phi = phi_NH3 / (phi_N2^(1/2) phi_H2^(3/2)),  phi_i = pure-species values at (T, P)
    Pure-gas phi_i from the Pitzer virial correlation: ln phi = (Pr/Tr)(B0 + omega B1),  B0 = 0.083 - 0.422/Tr^1.6,  B1 = 0.139 - 0.172/Tr^4.2

## Step 1 - reaction thermodynamics

| species | nu_i | A | 10^3 B | 10^6 C | 10^-5 D | dHf298 [J/mol] | dGf298 [J/mol] |
|---|---|---|---|---|---|---|---|
| N2 | -0.5 | 3.28 | 0.593 | 0 | 0.04 | 0 | 0 |
| H2 | -1.5 | 3.249 | 0.422 | 0 | 0.083 | 0 | 0 |
| NH3 | 1 | 3.578 | 3.02 | 0 | -0.186 | -46110 | -16450 |

- dA, dB, dC, dD = -2.9355, 0.0020905, 0, -33050 
- dH298 = sum nu_i dHf_i = -46.11 kJ/mol
- dG298 = sum nu_i dGf_i = -16.45 kJ/mol
- K(298.15 K) = exp(-dG298/RT0) = 762.212 
- dH(300.0 K) = -46.1512 kJ/mol
- ln K(300.0 K) = 6.52146 
- K(300.0 K) = 679.572 

## (a) y_NH3 at 1 bar and 300 K

- K(300 K) = 679.572 
- K P = y/(C (1-y)^2) = 679.572 
- y_NH3 = 0.934914 

## (b) temperature for y_NH3 = 0.50 at 1 bar (ideal gas)

    ln K(T) = ln[ 0.5 / (C x 0.25 x P) ]; solved for T with brentq on 250-1500 K
- required K = 6.1584 
- T = 399.499 K
- t = 126.349 degC

## (c) temperature for y_NH3 = 0.50 at 100 bar (ideal gas)

- required K = 0.061584 
- T = 577.639 K
- t = 304.489 degC

## (d) temperature for y_NH3 = 0.50 at 100 bar (ideal solution of gases)

- ln Phi at solution = -0.125904 
- Phi = 0.881699 
- phi_NH3, phi_N2, phi_H2 = 0.9251, 1.0335, 1.0213 
- T = 584.529 K
- t = 311.379 degC

## Final answers

- [RESULT a] y NH3 at 1 bar, 300 K = 0.934914 -
- [RESULT b] T for y = 0.5 at 1 bar = 399.499 K
- [RESULT b] T for y = 0.5 at 1 bar = 126.349 degC
- [RESULT c] T for y = 0.5 at 100 bar, ideal gas = 577.639 K
- [RESULT c] T for y = 0.5 at 100 bar, ideal gas = 304.489 degC
- [RESULT d] T for y = 0.5 at 100 bar, ideal solution = 584.529 K
- [RESULT d] T for y = 0.5 at 100 bar, ideal solution = 311.379 degC

## Verification

- [PASS] ln K(298.15) = -dG298/(R T0) 
- [PASS] analytic integrals agree with numerical quadrature (max |diff in ln K| = 8.9e-16)
- [PASS] van't Hoff relation d ln K/dT = dH/(R T^2) (max abs deviation 1.8e-12)
- [PASS] (a) y satisfies K P = y/(C (1-y)^2) 
- [PASS] (b) forward: y(T_b, 1 bar) = 0.50 
- [PASS] (c) forward: y(T_c, 100 bar, ideal gas) = 0.50 
- [PASS] (d) forward: y(T_d, 100 bar, ideal solution) = 0.50 
- [PASS] composition sums to 1 at (a): y_N2 + y_H2 + y_NH3 = 1 
- [PASS] higher pressure raises the temperature at which y = 0.5 (Le Chatelier): T_c > T_b 
- [PASS] virial correlation applicable: all Pr/Tr values moderate at solution 
Graph saved as plots/Q13_20/Q13_20.png


## Verification summary

Checks passed: 10/10  -> overall: ALL PASSED