# Question-wise mapping (verified against the assignment image before coding)

System for 12.12-12.19: **acetone(1)/water(2)** (chosen by the student). Ternary 12.20-12.23: acetone(1)/methanol(2)/water(3).
All use Eq. (10.5) (modified Raoult's law) + Antoine; kPa, K.

| Q | Given | Model | Calculations | Graph | Output files (in results/, plots/, reports/ of that question) |
|---|---|---|---|---|---|
| 12.12 | T = 333.15 K | Wilson | BUBL P and DEW P over x1 = 0..1 | P-x-y | `*_results.csv`, `*_diagram_data.csv`, `*_Wilson_Pxy.png`, `*_solution.md` |
| 12.13 | P = 101.33 kPa | Wilson | BUBL T and DEW T over x1 = 0..1 | t-x-y | same, `*_Wilson_txy.png` |
| 12.14 | T = 333.15 K | NRTL | as 12.12 | P-x-y | same, `*_NRTL_Pxy.png` |
| 12.15 | P = 101.33 kPa | NRTL | as 12.13 | t-x-y | same, `*_NRTL_txy.png` |
| 12.16 | T = 333.15 K; (a) x1 = 0.3; (b) y1 = 0.3; (c) P = (Pb+Pd)/2, z1 = 0.3 | Wilson | (a) BUBL P (b) DEW P (c) P,T-flash (d) azeotrope P_az, x1_az | none | csv, md |
| 12.17 | as 12.16 | NRTL | as 12.16 | none | csv, md |
| 12.18 | P = 101.33 kPa; (a) x1 = 0.3; (b) y1 = 0.3; (c) T = (Tb+Td)/2, z1 = 0.3 | Wilson | (a) BUBL T (b) DEW T (c) P,T-flash (d) azeotrope T_az, x1_az | none | csv, md |
| 12.19 | as 12.18 | NRTL | as 12.18 | none | csv, md |
| 12.20 | T = 338.15 K; (a) x1 = 0.3, x2 = 0.4; (b) y1 = 0.3, y2 = 0.4; (c) P = (Pb+Pd)/2, z1 = 0.3, z2 = 0.4 | Wilson (ternary) | BUBL P, DEW P, P,T-flash | none (tables) | csv, md |
| 12.21 | as 12.20 | NRTL (ternary) | as 12.20 | none | csv, md |
| 12.22 | P = 101.33 kPa; (a) x1 = 0.3, x2 = 0.4; (b) y1 = 0.3, y2 = 0.4; (c) T = (Tb+Td)/2, z1 = 0.3, **z2 = 0.2 (as printed)** | Wilson (ternary) | BUBL T, DEW T, P,T-flash | none (tables) | csv, md |
| 12.23 | as 12.22 | NRTL (ternary) | as 12.22 | none | csv, md |
