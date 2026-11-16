"""Runs every Q12_xx.py in its own process (independence test) and writes reports/FINAL_VALIDATION_REPORT.md
and reports/ALL_FINAL_ANSWERS.csv.  Run: python validate_all.py"""
import os, sys, json, csv, subprocess, glob
HERE = os.path.dirname(os.path.abspath(__file__))
Q = [f"12_{n}" for n in range(12, 24)] + ["13_14", "13_18", "13_20", "13_29", "13_31"]
L = ["# Final validation report", "", "Each question file was executed in a separate Python process (`python Q12_xx.py`).", "",
     "| Question | Exit code | Checks passed | Notes |", "|---|---|---|---|"]
allrows, unresolved, total_fail = [], [], 0
for q in Q:
    r = subprocess.run([sys.executable, os.path.join(HERE, f"Q{q}.py")], cwd=HERE, capture_output=True, text=True)
    open(os.path.join(HERE, "reports", f"Q{q}", f"Q{q}_console_output.txt"), "w", encoding="utf-8").write(r.stdout + r.stderr)
    js = json.load(open(os.path.join(HERE, "results", f"Q{q}", f"Q{q}_checks.json")))
    n = len(js["checks"]); ok = sum(c["passed"] for c in js["checks"]); fails = [c for c in js["checks"] if not c["passed"]]
    total_fail += len(fails)
    L.append(f"| Q{q.replace('_', '.')} | {r.returncode} | {ok}/{n} | {'; '.join(x[:90] for x in js['notes']) if js['notes'] else ''} |")
    for c in fails: unresolved.append(f"Q{q}: FAILED CHECK - {c['name']} {c['detail']}")
    if r.returncode != 0 and not fails: unresolved.append(f"Q{q}: non-zero exit with no failed check - see console output")
    allrows += list(csv.DictReader(open(os.path.join(HERE, "results", f"Q{q}", f"Q{q}_results.csv"))))
L += ["", f"Total failed checks: **{total_fail}**", "", "## Unresolved issues / limitations"]
L += [f"- {u}" for u in unresolved]
L += ["- **Input data NOT verified against the book.** Antoine constants, molar volumes, Wilson a_ij and NRTL b_ij/alpha were typed text supplied by the student. Passing numerical tests does not prove they match the book. Compare every entry in `common/params.py` with the book pages.",
      "- **NRTL methanol/water parameters are inconsistent with the Wilson set** (gamma-infinity of methanol in water: Wilson 1.8, NRTL 27). Q12.21 and Q12.23 (ternary NRTL) are therefore NOT reliable until that row is verified. See `reports/AUDIT_REPORT.md`.",
      "- Antoine equation is used outside the validity ranges given in the pasted table for acetone (above 55 degC) in several questions (flagged in each report).",
      "- Q12.22(c)/Q12.23(c) use z2 = 0.2 as printed on the assignment image (Q12.20(c) uses 0.4).",
      "- No liquid-liquid stability analysis is done (Wilson cannot predict LLE).", ""]
open(os.path.join(HERE, "reports", "FINAL_VALIDATION_REPORT.md"), "w", encoding="utf-8").write("\n".join(L))
with open(os.path.join(HERE, "reports", "ALL_FINAL_ANSWERS.csv"), "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(allrows[0].keys())); w.writeheader(); w.writerows(allrows)
print("\n".join(L))
