"""All thermodynamic parameters in ONE place.
SOURCE: supplied by the student as text (claimed to be Smith, Van Ness & Abbott 6th ed.,
Table B.2 / Chapter 12 table). NOT yet verified against photographed book pages.
Edit this file if your book's values differ; everything else recomputes.
"""
R = 8.314  # J/(mol K)

# Antoine: ln(Psat/kPa) = A - B/(t/degC + C);  V = molar volume, cm3/mol (25 degC)
COMPONENTS = {
    "acetone":  dict(A=14.3145, B=2756.22, C=228.06, V=74.05),
    "methanol": dict(A=16.5785, B=3638.27, C=239.50, V=40.73),
    "water":    dict(A=16.3872, B=3885.70, C=230.17, V=18.07),
}

# key=(species1, species2); a12,a21 Wilson (J/mol); b12,b21 NRTL (J/mol); alpha NRTL
PAIRS = {
    ("acetone", "water"):    dict(a12=1219.5, a21=6062.5, b12=2642.1, b21=5013.3, alpha=0.5343),
    ("methanol", "water"):   dict(a12=1965.9, a21=449.6,  b12=3538.7, b21=6852.0, alpha=0.2994),
    ("acetone", "methanol"): dict(a12=-677.8, a21=2441.4, b12=773.3,  b21=932.1,  alpha=0.3084),
}

def pair(i, j):
    """Return (a_ij, a_ji, b_ij, b_ji, alpha) for species names i, j (any order)."""
    if (i, j) in PAIRS:
        p = PAIRS[(i, j)]; return p["a12"], p["a21"], p["b12"], p["b21"], p["alpha"]
    p = PAIRS[(j, i)]; return p["a21"], p["a12"], p["b21"], p["b12"], p["alpha"]

# ---- provenance flags (read by every question report) ---------------------------------
PARAMETERS_VERIFIED_AGAINST_BOOK = False   # set True only after you compare each number with the book
PARAMETER_SOURCE = "Smith, Van Ness & Abbott, 6th ed. (Antoine constants, molar volumes, Wilson and NRTL parameters as supplied by the student)"
# Antoine validity ranges (deg C) as given in the pasted table (themselves unverified)
RANGES_C = {"acetone": (-13.0, 55.0), "methanol": (15.0, 84.0), "water": (0.0, 200.0)}
