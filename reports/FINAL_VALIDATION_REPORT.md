# Final validation report

Each question file was executed in a separate Python process (`python Q12_xx.py`).

| Question | Exit code | Checks passed | Notes |
|---|---|---|---|
| Q12.12 | 0 | 5/5 | Antoine for acetone used at 60.00 C , outside stated range -13.0..55.0 C (extrapolation). |
| Q12.13 | 0 | 5/5 | Antoine for acetone used at 56.20 C (bubble/dew range at this P), outside stated range -13 |
| Q12.14 | 0 | 5/5 | Antoine for acetone used at 60.00 C , outside stated range -13.0..55.0 C (extrapolation). |
| Q12.15 | 0 | 5/5 | Antoine for acetone used at 56.20 C (bubble/dew range at this P), outside stated range -13 |
| Q12.16 | 0 | 16/16 | Antoine for acetone used at 60.00 C , outside stated range -13.0..55.0 C (extrapolation). |
| Q12.17 | 0 | 16/16 | Antoine for acetone used at 60.00 C , outside stated range -13.0..55.0 C (extrapolation). |
| Q12.18 | 0 | 17/17 | Antoine for acetone used at 62.67 C , outside stated range -13.0..55.0 C (extrapolation).; Antoine for acetone used at 90.60 C , outside stated range -13.0..55.0 C (extrapolation). |
| Q12.19 | 0 | 17/17 | Antoine for acetone used at 62.87 C , outside stated range -13.0..55.0 C (extrapolation).; Antoine for acetone used at 90.61 C , outside stated range -13.0..55.0 C (extrapolation). |
| Q12.20 | 0 | 15/15 | Antoine for acetone used at 65.00 C , outside stated range -13.0..55.0 C (extrapolation). |
| Q12.21 | 0 | 15/15 | Antoine for acetone used at 65.00 C , outside stated range -13.0..55.0 C (extrapolation). |
| Q12.22 | 0 | 16/16 | In (c) the assignment image prints z2 = 0.2 (whereas Q12.20(c) uses z2 = 0.4). Used as pri; Antoine for acetone used at 61.80 C , outside stated range -13.0..55.0 C (extrapolation).; Antoine for acetone used at 75.94 C , outside stated range -13.0..55.0 C (extrapolation). |
| Q12.23 | 0 | 16/16 | In (c) the assignment image prints z2 = 0.2 (whereas Q12.20(c) uses z2 = 0.4). Used as pri; Antoine for acetone used at 57.77 C , outside stated range -13.0..55.0 C (extrapolation).; Antoine for acetone used at 70.24 C , outside stated range -13.0..55.0 C (extrapolation). |
| Q13.14 | 0 | 7/7 |  |
| Q13.18 | 0 | 7/7 |  |
| Q13.20 | 0 | 10/10 |  |
| Q13.29 | 0 | 10/10 | Elemental sulfur is treated as a pure solid, as stated in the problem. At 723 K real sulfu |
| Q13.31 | 0 | 5/5 |  |

Total failed checks: **0**

## Unresolved issues / limitations
- **Input data NOT verified against the book.** Antoine constants, molar volumes, Wilson a_ij and NRTL b_ij/alpha were typed text supplied by the student. Passing numerical tests does not prove they match the book. Compare every entry in `common/params.py` with the book pages.
- **NRTL methanol/water parameters are inconsistent with the Wilson set** (gamma-infinity of methanol in water: Wilson 1.8, NRTL 27). Q12.21 and Q12.23 (ternary NRTL) are therefore NOT reliable until that row is verified. See `reports/AUDIT_REPORT.md`.
- Antoine equation is used outside the validity ranges given in the pasted table for acetone (above 55 degC) in several questions (flagged in each report).
- Q12.22(c)/Q12.23(c) use z2 = 0.2 as printed on the assignment image (Q12.20(c) uses 0.4).
- No liquid-liquid stability analysis is done (Wilson cannot predict LLE).
