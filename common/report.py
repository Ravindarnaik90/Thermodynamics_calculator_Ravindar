"""Console + Markdown + CSV + JSON reporting used by every Q12_xx.py."""
import os, csv, json, datetime
import numpy as np
from . import params

class Report:
    def __init__(self, qid, title, here, source=None):
        self.qid, self.title, self.here = qid, title, here
        self.md, self.rows, self.checks, self.notes = [], [], [], []
        self.dirs = {k: os.path.join(here, k, qid) for k in ("results", "plots", "reports")}
        for d in self.dirs.values(): os.makedirs(d, exist_ok=True)
        self.out("=" * 78); self.out(title); self.out("=" * 78)
        self.md = ["# " + title, "", f"*Generated {datetime.date.today()} by `{qid}.py`*", ""]
        self.out(f"Parameter source: {params.PARAMETER_SOURCE}")
        self.md += [f"> **Parameter provenance:** {params.PARAMETER_SOURCE}. Verified against book: **{params.PARAMETERS_VERIFIED_AGAINST_BOOK}**.", ""]
    def out(self, s=""): print(s); self.md.append(s)
    def h(self, s): print("\n--- " + s + " ---"); self.md += ["", "## " + s, ""]
    def p(self, s): print(s); self.md += [s, ""]
    def eq(self, s): print("    " + s); self.md.append("    " + s)
    def kv(self, name, val, unit=""):
        v = f"{val:.6g}" if isinstance(val, (float, np.floating)) else str(val)
        print(f"  {name} = {v} {unit}"); self.md.append(f"- {name} = {v} {unit}")
    def vec(self, name, comps, v, unit=""):
        self.kv(name, "[" + ", ".join(f"{c}: {x:.6f}" for c, x in zip(comps, v)) + "]", unit)
    def result(self, part, quantity, value, unit, status):
        self.rows.append(dict(question=self.qid, part=part, quantity=quantity,
                              value="" if value is None else f"{float(value):.8g}", unit=unit, solver_status=status))
        self.kv(f"[RESULT {part}] {quantity}", value if value is not None else "n/a", unit)
    def check(self, name, ok, detail=""):
        self.checks.append(dict(name=name, passed=bool(ok), detail=detail))
        line = f"[{'PASS' if ok else 'FAIL'}] {name} {detail}"; print("  " + line); self.md.append("- " + line)
    def note(self, s): self.notes.append(s); print("  NOTE: " + s); self.md.append("- **NOTE:** " + s)
    def range_warnings(self, comps, T_list_K, label=""):
        for c in comps:
            lo, hi = params.RANGES_C[c]
            for T in T_list_K:
                if not (lo <= T - 273.15 <= hi):
                    self.note(f"Antoine for {c} used at {T-273.15:.2f} C {label}, outside stated range {lo}..{hi} C (extrapolation).")
                    break
    def table(self, header, rows):
        print("  " + " | ".join(header)); self.md += ["| " + " | ".join(header) + " |", "|" + "---|" * len(header)]
        for r in rows:
            cells = [f"{x:.6g}" if isinstance(x, (float, np.floating)) else str(x) for x in r]
            print("  " + " | ".join(cells)); self.md.append("| " + " | ".join(cells) + " |")
        self.md.append("")
    def save(self):
        with open(os.path.join(self.dirs["results"], f"{self.qid}_results.csv"), "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=["question", "part", "quantity", "value", "unit", "solver_status"]); w.writeheader(); w.writerows(self.rows)
        allok = all(c["passed"] for c in self.checks)
        json.dump(dict(question=self.qid, all_checks_passed=allok, checks=self.checks, notes=self.notes,
                       parameters_verified=params.PARAMETERS_VERIFIED_AGAINST_BOOK),
                  open(os.path.join(self.dirs["results"], f"{self.qid}_checks.json"), "w"), indent=1)
        self.h("Verification summary")
        self.out(f"Checks passed: {sum(c['passed'] for c in self.checks)}/{len(self.checks)}  -> overall: {'ALL PASSED' if allok else 'SOME FAILED'}")
        open(os.path.join(self.dirs["reports"], f"{self.qid}_solution.md"), "w", encoding="utf-8").write("\n".join(self.md))
        print(f"\nSaved: results/{self.qid}/, reports/{self.qid}/" + (f", plots/{self.qid}/" if os.listdir(self.dirs['plots']) else ""))
        return allok

def describe_model(rep, m, T, x):
    """Print the intermediate quantities for a hand check."""
    from .wilson import lambda_matrix
    from .nrtl import tau_G
    comps = m.comps
    rep.kv("T", T, "K"); rep.kv("t", T - 273.15, "deg C")
    rep.vec("P_sat (Antoine)", comps, m.psat(T), "kPa")
    if m.model == "Wilson":
        Lm = lambda_matrix(comps, T)
        rep.p("Wilson Lambda_ij = (Vj/Vi) exp(-a_ij/RT):")
        rep.table([""] + comps, [[comps[i]] + list(Lm[i]) for i in range(len(comps))])
    else:
        tau, G = tau_G(comps, T)
        rep.p("NRTL tau_ij = b_ij/RT:"); rep.table([""] + comps, [[comps[i]] + list(tau[i]) for i in range(len(comps))])
        rep.p("NRTL G_ij = exp(-alpha tau_ij):"); rep.table([""] + comps, [[comps[i]] + list(G[i]) for i in range(len(comps))])
    rep.vec("gamma at x=" + str([round(float(v), 4) for v in x]), comps, m.gamma(x, T))
