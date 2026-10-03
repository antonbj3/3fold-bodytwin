# Mature work in the old project that nothing in the workspace calls — measured 2026-10-02

Written in response to Anton's question *"is there more important stuff like that you missed before"*, after the
adaptive resolution he asked for turned out to have been built and graded since 2026-07-04 without
any of our workspaces calling it.

## The measure
The names of everything carrying a result in `~/projects/bodytwin` — 821 reports in `reports/`, 374
JSON in `scripts/physics_exp/`, 417 distilled memories in `docs/inherited_memory/` — swept against the entire
workspace in one pass. Log files, `night_rounds/`, `results/MAP/` and `.tsv` do NOT count as a
hit: an agent that mentioned a name in a log has not connected anything.

| pool | okallade | av |
|---|---|---|
| `docs/inherited_memory/*.md` (distilled, graded lessons) | **410** | 417 |
| `reports/*.json` | 339 | 821 |
| `scripts/physics_exp/*.json` | 147 | 374 |
| **totalt** | **896** | 1 612 |

The densest loss is therefore not the reports but the distilled memories: **7 of 417 are
mentioned**. 117 of the uncalled ones are `type: project`. Many carry the text `RECALL-ONLY — index
saturated` themselves, so they were written by a session that knew they would be lost. The name list is in
`OLD_WORK_UNCALLED_NAMES.txt`.

The natural objection is that much of it is certification physics from the motion track and not
body simulation. That is true for most. But four hits bear directly on work run
this evening, and two of them defeat something we have already done.

## 1. How a scaling exponent is certified — DEFEATS OUR OWN WAY OF READING THEM
`exponent-cert-by-self-similar-collapse-and-causal-coupling-not-3pt-fit`

A three-point log-log slope is **almost tautological**, and a prefactor that itself depends on
the control parameter gives the right slope but is the wrong mechanism. Two independent certificates beat the fit:
**self-similar collapse** (y/Rⁿ against t/R^a falls on ONE master curve; a parameter-dependent prefactor gives
the slope still but BREAKS the collapse) and **a causal coupling knob** (scale the mechanism's
coupling and the exponent should move by a PREDICTED amount). Proven on pipe transient
growth: G_max/Re² constant 7,2e-5, peak deviation of the collapse 0,0068 over Re 1k–4k, and C₀∝√Re gave
Re exponent 3,000 = 2+2·0,5.

And sharper still: **a stable exponent can be a stably WRONG exponent.** An xmin plateau can
coexist with a poor absolute fit — one quantity was xmin-stable with drift 0,12 but KS=0,37.
So BOTH plateau AND KS goodness are required.

**What this defeats for us:** the fibril lane's segment exponents 0,2353 / 0,3558 / 0,1877 use few points,
and the network's bending slope 1,186 → 1,033 → 0,980 I recorded as converging toward the
published 1,0. That is precisely the reading the memory forbids. Convergence of a fitted exponent over
three runs is no certificate. The lane has been redirected.

## 2. Abstention has TWO independent axes, not one — refines the adaptive resolution
`graded-abstention-paydownable-cliff-fundamental-action`

The type of abstention determines the action. **GRADUAL** (λ<0): the band decreases monotonically with
refinement effort and crosses below the tolerance at a FINITE effort_crit → pay it down,
meaning collect more resolution. **CLIFF** (λ>0): the band GROWS under refinement and never gets below
the tolerance → no amount of data helps, change readout. Mixing them wastes effort: refining a cliff
never certifies.

And the own correction in the same file, which is the important part: **smoothness and paydownability are TWO
INDEPENDENT axes.** The first version mixed them because its only control was
mesh refinement. On a real LPBF forward model, identifiability decreases SMOOTHLY with laser power but
is still STRUCTURAL and not paydownable — a *gradual but fundamental* case. **The control axis can
be the physical OPERATING POINT, not just resolution.**

**What this gives us:** my hemostasis certificate says that bleeding rate is structurally
unidentifiable below f = 0,5. That is the paydownability axis. I have not measured the smoothness axis, and
the operating point — arterial pressure, vessel caliber — is a control axis I did not test. Class 1/2/3 is
therefore too coarse; the right map is 2×2.

## 3. Which spectral quantity certifies distance to instability depends on the operator's normality
`criticality-governor-hierarchy-normality-selects-sigmamin-abscissa-kreiss` — three levels:
self-adjoint (the smallest eigenvalue IS the margin, no transient), normal nonsymmetric, and
nonnormal (where the spectrum does NOT license absence of transient growth). The pipe case is
subcritical for all Re and still G_max reaches 288.
**Relevance:** the immune lane's quasi-stationary region and class changes are a stability question about a
nonsymmetric operator. If it reads stability from eigenvalues, it can miss transient growth.

## 4. Contact force margin with model-error-σ anchored to omitted force RMS
`p4a-contact-force-margin-overdet-independent-depth-gradients` — 2D peg-in-hole insertion force
with model-error-σ from friction uncertainty, anchored against omitted force RMS with ratio ~1,0.
**Relevance:** needle and tool contact with tissue has the same problem form, and we have chased
contact-age artifacts for nine rounds.

## The resulting rule
Before a lane derives anything: search ALL FOUR pools, not just the `data/` index (which is the
only one `notes/OLD_CELL_INDEX.json` covers). Four pools: `data/` (1 620 cell records), `reports/`
(821), `scripts/physics_exp/` (374) and `docs/inherited_memory/` (417). My cell index covered one of
four, and that is why I could write that something was "missing" three times in one evening.
