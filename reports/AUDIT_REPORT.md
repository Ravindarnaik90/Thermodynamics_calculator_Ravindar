# Audit of existing project (src/, main.py results)

## 1. Antoine equation form, units
- [PASS] acetone: T_sat(101.325 kPa) = 56.20 C vs literature 56.05 C 
- [PASS] acetone: psat round trip 
- [PASS] methanol: T_sat(101.325 kPa) = 64.70 C vs literature 64.7 C 
- [PASS] methanol: psat round trip 
- [PASS] water: T_sat(101.325 kPa) = 100.00 C vs literature 100.0 C 
- [PASS] water: psat round trip 

## 2/3. Wilson & NRTL vs independent closed-form binary equations (incl. component-order swap)
- [PASS] acetone/water: generic = closed form and order-invariant (12 states x 2 models) 
- [PASS] methanol/water: generic = closed form and order-invariant (12 states x 2 models) 
- [PASS] acetone/methanol: generic = closed form and order-invariant (12 states x 2 models) 

## 4. Gibbs-Duhem (ternary, finite difference): sum x_i dln(gamma_i) = 0
- [PASS] Wilson: max |sum x dlnγ| = 1.23e-10 
- [PASS] NRTL: max |sum x dlnγ| = 6.94e-11 

## 5. Parameter plausibility: gamma-infinity (Wilson vs NRTL fitted to same data should agree)
- [PASS] acetone/water: γ∞ Wilson 10.7/4.89 vs NRTL 10.5/5.09 
- [FAIL] methanol/water: γ∞ Wilson 1.81/1.14 vs NRTL 27.1/11.4  <-- NRTL and Wilson parameters are NOT mutually consistent; one set is probably wrong/mis-transcribed
- [PASS] acetone/methanol: γ∞ Wilson 1.81/1.77 vs NRTL 1.79/1.78 

## 6. Root-finding cross-check with scipy.fsolve (independent algorithm, different starting guess)
- [PASS] Wilson DEW P (binary): successive substitution 28.380141 vs fsolve 28.380141 kPa 
- [PASS] Wilson DEW P (ternary): 64.864676 vs fsolve 64.864676 kPa 
- [PASS] NRTL DEW P (binary): successive substitution 28.379278 vs fsolve 28.379278 kPa 
- [PASS] NRTL DEW P (ternary): 80.837806 vs fsolve 80.837806 kPa 

## 7. Code-review notes (manual)
- Component ordering: pair(i,j) returns a_ij for ordered names; verified by swap-invariance test above.
- Wilson indexing: Λ_ij = (V_j/V_i) exp(-a_ij/RT) uses a_ij of row species i - consistent with SVNA footnote.
- Temperatures are K inside all functions; Antoine converts with T-273.15; pressure in kPa throughout.
- Antoine validity ranges: the pasted ranges (acetone -13..55 C, methanol 15..84 C) are NARROWER than needed for several state points (acetone/water dew T ~ 90 C, ternary dew T ~ 76 C). Antoine is extrapolated there; flagged in each question report.
- Wilson cannot predict liquid-liquid splitting and no stability test is done in the flash; acceptable for these state points but unverified.
- Issue found in earlier main.py: nothing numerically wrong was found in the equations; the main weakness is INPUT DATA provenance (typed text, not book pages).
- 12.22(c) uses z2 = 0.2 as printed on the image (differs from 12.20(c) z2 = 0.4); keep as printed unless the instructor says otherwise.

## Summary
- Failed audit items: 1 -> methanol/water: γ∞ Wilson 1.81/1.14 vs NRTL 27.1/11.4