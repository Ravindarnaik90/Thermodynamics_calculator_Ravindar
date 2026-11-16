# Problem 13.31 - Liquid-phase isomerisation A -> B at 298.15 K

*Generated 2026-11-11 by `Q13_31.py`*

Parameter source: Smith, Van Ness & Abbott, 6th ed. (Antoine constants, molar volumes, Wilson and NRTL parameters as supplied by the student)
> **Parameter provenance:** Smith, Van Ness & Abbott, 6th ed. (Antoine constants, molar volumes, Wilson and NRTL parameters as supplied by the student). Verified against book: **False**.


## Problem statement

The isomerisation reaction A -> B occurs in the liquid phase. A and B are miscible in all proportions, and the excess Gibbs energy of the mixture is G^E/RT = 0.1 xA xB. If dG0(298 K) = -1000 J/mol, what is the equilibrium composition at 298.15 K? How much error is introduced if A and B are assumed to form an ideal solution?


## Governing equations

    K = exp(-dG0/RT)  with standard states of pure liquid A and B at T (and P);  nu = (-1, +1)
    K = (x_B gamma_B) / (x_A gamma_A)
    G^E/RT = 0.1 xA xB  =>  ln gamma_A = 0.1 xB^2,   ln gamma_B = 0.1 xA^2        (symmetric Margules, A12 = A21 = 0.1)
    Ideal solution (gamma = 1):  K = xB/(1 - xB)  =>  xB = K/(1+K)

## Given data

- T = 298.15 K
- dG0(298) = -1000 J/mol
- R = 8.314 J/(mol K)

## Step 1 - equilibrium constant

- K = exp(1000/(R T)) = 1.49693 
- ln K = 0.403418 

## Step 2 - non-ideal solution: solve ln K = ln(xB/xA) + ln gamma_B - ln gamma_A for xB (brentq)

- xB = 0.604517 
- xA = 0.395483 
- gamma_A, gamma_B = 1.03722, 1.01576 

## Step 3 - ideal solution

- xB (ideal) = 0.599509 
- xA (ideal) = 0.400491 

## Final answers

- [RESULT -] K = 1.49693 -
- [RESULT -] xB, non-ideal solution (G^E/RT = 0.1 xA xB) = 0.604517 -
- [RESULT -] xA, non-ideal solution = 0.395483 -
- [RESULT -] xB, ideal solution = 0.599509 -
- [RESULT -] xA, ideal solution = 0.400491 -
- [RESULT -] error in xB from assuming ideal solution = -0.00500826 mole fraction
- [RESULT -] relative error in xB = -0.828473 %

## Verification

- [PASS] equilibrium relation xB gB / (xA gA) = K 
- [PASS] xA + xB = 1 and 0 < xB < 1 
- [PASS] ln gamma_i = d(nG^E/RT)/dn_i agrees with the closed forms 
- [PASS] Gibbs-energy minimum of the mixture coincides with the K-equation root (minimum at xB = 0.60451683)
- [PASS] Gibbs-energy minimum (ideal) coincides with K/(1+K) 
Graph saved as plots/Q13_31/Q13_31.png


## Verification summary

Checks passed: 5/5  -> overall: ALL PASSED