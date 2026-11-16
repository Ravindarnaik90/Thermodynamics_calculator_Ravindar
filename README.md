# individual_solutions - one independent Python program per question

## Read this first (accuracy)
* All 12 programs run and pass their numerical checks, **but the input data (Antoine constants, molar volumes, Wilson a_ij,
  NRTL b_ij and alpha) were typed text, not book pages.** Compare every number in `common/params.py` with your book.
* The NRTL **methanol/water** parameters are inconsistent with the Wilson ones (gamma-infinity 27 vs 1.8). Q12.21 and Q12.23
  depend on them and should not be trusted until that row is checked (see `reports/AUDIT_REPORT.md`).
* Q12.22(c)/12.23(c) use z2 = 0.2 exactly as printed on the assignment image.

## Setup (Windows PowerShell)
```powershell
cd path\to\assignment_3_solution\individual_solutions
python -m venv venv
.\venv\Scripts\Activate.ps1
pip install -r requirements.txt
```
If PowerShell blocks the activate script run once: `Set-ExecutionPolicy -Scope Process -ExecutionPolicy Bypass`

## Run one question
```powershell
python Q12_12.py
python Q12_16.py
```
(each file works alone; replace the number for other questions)

## Run all questions + validation report
```powershell
python validate_all.py
```

## Re-run the audit of the earlier project
```powershell
python audit_existing_code.py
```

## Where results go (per question, e.g. Q12_12)
* `results\Q12_12\` : `Q12_12_results.csv` (final answers), `Q12_12_checks.json`, diagram data CSV
* `plots\Q12_12\`   : PNG at 300 dpi (only for 12.12-12.15)
* `reports\Q12_12\` : `Q12_12_solution.md` (step-by-step) and the console output text
* `reports\FINAL_VALIDATION_REPORT.md`, `reports\ALL_FINAL_ANSWERS.csv`, `reports\AUDIT_REPORT.md`, `reports\QUESTION_MAPPING.md`

## Folder layout
`Q12_12.py ... Q12_23.py`, `common/` (Antoine, Wilson, NRTL, bubble/dew, flash, azeotrope, plots, cross-checks, reporting),
`tools/generate_question_files.py` (template that wrote the Q files; not needed to run them).
To use your own book values: edit `common/params.py`, set `PARAMETERS_VERIFIED_AGAINST_BOOK = True` only after checking, rerun.
