# MECHANISM LIGAMENT STIFFNESS — Blankevoort1991Ligament `linear_stiffness` vs cadaveric tensile-test literature (2026-07-22)

Validates the twin's knee-ligament **STIFFNESS** parameter — distinct from the reference-strain
fragility already swept in `docs/MECHANISM_LIGAMENT_EPSR_SENSITIVITY.md`. Extracts the per-bundle
`Blankevoort1991Ligament.linear_stiffness` for all 84 ACL/PCL/MCL/LCL bundles (`docs/MECHANISM_KNEE_LIGAMENTS.md`),
derives the whole-ligament **STRUCTURAL** stiffness (N/mm, force per unit elongation — what a
materials-testing machine reports) two independent ways, and compares against 8 decorrelated
cadaveric tensile-test studies, all PMIDs verified live this session. Script:
`scripts/msk/ligament_stiffness_extract.py`. Run:
`source_repository/.venv-msk/bin/python3 scripts/msk/ligament_stiffness_extract.py`.

**Isolation respected:** read-only against the model/config; nothing pushed/committed; this doc,
the script, and the two evidence JSONs are the only new files.

---

## 0. Headline

**ACL and PCL PASS; MCL and LCL FAIL, in opposite directions.** The model's whole-ligament
structural stiffness (derived from the per-bundle `linear_stiffness`, side `r`) is **ACL 283 N/mm**
(literature range 124–286, young/healthy range 218–286 → **inside**), **PCL 169 N/mm** (Trent1976
whole-PCL 179.5, Harner1995-bundle-sum 177 → **inside**, within 6% of 2 of 3 anchors), **MCL 172
N/mm** (literature 50–85 → **2.0–2.7x too stiff**), **LCL 39 N/mm** (literature 47–71 → **~35–45%
too soft**). The single most decisive test is the **MCL:LCL ratio**: two decorrelated cadaveric
studies 36 years apart (Trent 1976, Wilson 2012), on totally different donor cohorts, both find
MCL and LCL in near **parity** (ratio 1.07–1.18); this model's ratio is **4.43** — a **3.7–4.1x**
departure from a pattern confirmed twice independently. Tracing this to its root cause (not just
its magnitude): Lenhart2015's own methods (fetched live, PMC4886716) state every bundle's stiffness
is `cross-sectional-area × an ASSUMED UNIFORM 125 MPa modulus` — 125 MPa itself sits squarely
inside the task's own cited ~100–150 MPa collagen range (the **material** input is fine in
isolation), but applying **one modulus to both MCL and LCL** mechanically predicts exactly this
kind of imbalance, because Wilson 2012's own measurements show MCL is ~680% wider than LCL yet the
two carry near-identical real structural stiffness — meaning their **true** moduli must differ
substantially, which a uniform-modulus model cannot capture. **A genuine, geometrically-traced
falsifier hit for MCL/LCL, alongside a genuine PASS for ACL/PCL — reported symmetrically.**

**Also flagged up front, matching this repo's established pattern**: the task's own supplied PMID
for Woo 1991 (**1962712**) is **WRONG** — verified live to resolve to an unrelated 1991 clinical
paper (Shelbourne & Nitz, "The O'Donoghue triad revisited"). The correct PMID is **1867330**. This
is the **3rd** instance of this exact failure mode against this same body of knee-ligament work (2
priors already documented in `docs/MECHANISM_KNEE_LIGAMENTS.md`).

---

## 1. What's actually in the model (machine-verified, not eyeballed)

`Blankevoort1991Ligament.linear_stiffness` is **N per unit (dimensionless) strain** — despite the
"N" units this is a **material-like** parameter, not directly a structural N/mm number (the task's
own framing). Read from `scripts/msk/knee_ligament_config.json` (84 bundles: ACLam×6, ACLpl×6,
PCLal×5, PCLpm×5, MCLd×5, MCLs×6, MCLp×5, LCL×4, ×2 sides), machine-audited for uniformity rather
than assumed:

| bundle group | ligament | n fibers (side r) | `linear_stiffness` (N/strain) | uniform within group? |
|---|---|---:|---:|:-:|
| ACLam, ACLpl | ACL | 12 | **820** | PASS (12/12 identical) |
| PCLal | PCL | 5 | **900** | PASS |
| PCLpm | PCL | 5 | **300** | PASS |
| MCLd, MCLs, MCLp | MCL | 16 | **500** | PASS |
| LCL | LCL | 4 | **600** | PASS |

`transition_strain` = 0.06 and `damping_coefficient` = 0.003 are **identical across all 84
bundles** (machine-checked, not assumed). All values trace to the vendored Lenhart2015
`lenhart2015.osim` companion model, per `docs/MECHANISM_KNEE_LIGAMENTS.md` §2.

---

## 2. Material vs structural — the distinction the task asked to enforce

The Blankevoort force law's linear region is `F = k·(strain − εₜ/2)`, so `dF/d(strain) = k` — a
force-per-unit-strain ("material-like") quantity. A cadaveric tensile test instead reports
`dF/dx` (force per unit **elongation** in mm — "structural", depends on cross-section AND length).
Since `strain = x/L₀`, the chain rule gives the conversion: **`dF/dx = k / L₀`** (L₀ = slack
length). This is not a heuristic — it is the same relation a linear elastic bar obeys
(`k_structural = E·A/L₀`), and, as found in §6, it is *exactly* how the model's own `k` was built.

**Two independent derivations**, cross-checked against each other:

- **(A) Naive all-linear upper bound**: closed-form `Σ(k_bundle / slack_length_bundle)` per group
  — what the group's stiffness would be if every bundle were simultaneously past its
  `transition_strain`. Rote but exact.
- **(B) Finite-difference virtual tensile test**: a genuinely **geometric** derivation directly on
  the live model. At a fixed knee angle, the femur is held fixed and the tibia's real attachment
  points (via `findStationLocationInGround`) are rigidly translated along the ligament group's own
  anatomical pull axis (centroid-to-centroid, in ground), from 0 to 16 mm elongation in 0.1 mm
  steps. Each bundle's length/strain/force is recomputed via the closed-form Blankevoort law — the
  *same* form already machine-cross-checked to `0.00e+00 N` against `getTotalForce()` in
  `add_knee_ligaments.py` — and the group's linear-region slope is a least-squares fit over the
  last 30% of the sweep, mimicking exactly what a materials-testing analyst does with real
  force-elongation curve data (fit the post-toe-region line, not a near-zero-displacement tangent).

**Sign-convention self-catch**: the first version of test (B) had the pull direction backwards —
every `dF/dx` came out *negative* even though the engaged-bundle count clearly rose with the
nominally "positive" perturbation. This was an OODA-loop catch, not a silent trust: the
`sign_check_pass` diagnostic (engaged-count comparison between the two perturbation directions)
flagged it, the sign was re-derived and fixed, and the corrected version passes **20/20**
(group × angle) sign checks.

**Result**: (A) and (B) agree to **<1%** everywhere full bundle recruitment is achieved by 16 mm of
pull (ACL/PCL/LCL at all 5 tested angles 0°/15°/30°/60°/90°; MCL at 0°/15°/30°/60°). Only MCL at
90° is not fully recruited even at 16 mm (13/16 bundles engaged-linear, R²=0.9992 vs 1.0000
elsewhere) — honestly flagged, not chased further (diminishing returns; the asymptote is visibly
converging toward the same ~172 N/mm from below).

| ligament | (A) naive upper bound (N/mm) | (B) deep-pull FD, angle range (N/mm) | R² |
|---|---:|---:|---:|
| ACL | 283.14 | 282.4 – 282.9 (angle-invariant) | 1.0000 |
| PCL | 168.55 | 167.3 – 167.9 (angle-invariant) | 1.0000 |
| MCL | 171.74 | 134.8 (90°, under-recruited) – 172 (0-60°) | 0.9992–1.0000 |
| LCL | 38.79 | 38.48 – 38.78 | 0.9999–1.0000 |

For comparison, the **"as-configured" local stiffness** (small perturbation around the model's
actual current pose, i.e. what the ligament contributes to load-sharing *right now*, not its
nominal capacity) is far lower almost everywhere — e.g. ACL at 0° is only 18.1 N/mm (of its 283
N/mm capacity) because just 2 of 12 bundles are even in the toe region there; ACL/LCL both drop to
**0** by 60–90° flexion (fully slack, consistent with the already-published flexion-sweep pattern
in `docs/MECHANISM_KNEE_LIGAMENTS.md` §4). This is a distinct, and arguably more operationally
relevant, fact from the "what's the ligament's rated capacity" question this task centers on — both
are reported, not conflated.

**External anchor, checked before trusting anything else**: the freshly-loaded
`model_with_knee_ligaments.osim`'s baseline (elongation=0, knee_angle_r=0°) total forces per
group — ACL=6.833N, PCL=13.497N, MCL=91.441N, LCL=72.002N — match
`docs/MECHANISM_KNEE_LIGAMENTS.md`'s already-published flexion-sweep table at 0° **exactly** (ACL
6.83, PCL 13.50, MCL 91.42, LCL 72.00), confirming this task is using the same, already-verified
model artifact, not a stale or diverged copy.

---

## 3. Literature (all PMIDs verified live this session, per the task's instruction)

### 3.0 Honesty flag — the task's own PMID for Woo 1991 is wrong

| given PMID | what it actually is (verified live via NCBI eutils) |
|---|---|
| **1962712** | Shelbourne KD, Nitz PA. "The O'Donoghue triad revisited. Combined knee injuries involving anterior cruciate and medial collateral ligament tears." *Am J Sports Med.* 1991. An unrelated clinical/epidemiological paper — not a tensile-testing study. |

Correct PMID, found via NCBI ESearch (author + title terms) and confirmed via ESummary/EFetch:

**Woo SL, Hollis JM, Adams DJ, Lyon RM, Takai S.** "Tensile properties of the human
femur-anterior cruciate ligament-tibia complex. The effects of specimen age and orientation."
*Am J Sports Med.* 1991;19(3):217-25. **PMID [1867330](https://pubmed.ncbi.nlm.nih.gov/1867330/).**
This is the **3rd** instance of a wrong task-supplied PMID against this exact body of work (2 prior
corrections already in `docs/MECHANISM_KNEE_LIGAMENTS.md` for Lenhart2015 and Blankevoort&Huiskes1991).

### 3.1 Per-ligament literature (structural stiffness, N/mm)

All studies cross-tabulated against **Peters AE et al., "Tissue material properties and
computational modelling of the human tibiofemoral joint: a critical review," 2018** (PMID
[29379690](https://pubmed.ncbi.nlm.nih.gov/29379690/), open access, **PMC5787350**, Table 3) — a
peer-reviewed systematic review, itself an over-determination anchor, not a single cherry-picked
number. Where a first extraction pass from this table looked suspect, it was independently
re-verified (see §4).

**ACL:**

| study | PMID | stiffness (N/mm) | donor age |
|---|---|---|---|
| Trent, Walker & Wolf 1976 | — | 138.3 | 29–55 |
| Noyes & Grood 1976 | — | young 182.0 / old 129.0 | 16–86 |
| **Woo et al. 1991** | **1867330** | young(22-35) 218–242 [242±28 confirmed direct from abstract]; middle(40-50) 192–220; old(60-97) 124–180 | 22–97 |
| Butler et al. 1992 | — | AM 238.1 / AL 285.9 / posterior 154.9 | 26±4 |

Aggregate range **124–285.9 N/mm**; young/healthy-donor range **218–285.9 N/mm**.

**PCL** (task-named "Race & Amis"):

| study | PMID | stiffness (N/mm) | donor age |
|---|---|---|---|
| Trent, Walker & Wolf 1976 (whole PCL) | — | 179.5 | 29–55 |
| **Race & Amis 1994** | **8106532** | AL(aPC) 347.0 / PM(pPC) 77.0 | 53–98 |
| Harner et al. 1995 | **8600743** | AL 120.0 / PM 57.0 | 48–77 |

Race & Amis's table-derived failure loads (1620N AL / 258N PM) match their own abstract's directly-
stated "aPC mean strength of 1.6 kN" to <2% — corroborating that table row. Whole-PCL
bundle-summed estimates: Trent(direct)=179.5, Harner-summed=177.0, Race&Amis-summed=424.0 — 2 of 3
anchors converge near 177–180; Race&Amis is a ~2.4x outlier, itself **not** explained by an
age-decline story (its donors, 53–98, are *older* than Harner's 48–77, which should predict *lower*
not *higher* stiffness if age-decline applied — flagged as an unresolved inter-study disagreement,
not adjudicated here). Aggregate range **57–347 N/mm**.

**MCL:**

| study | PMID | stiffness (N/mm) | donor age |
|---|---|---|---|
| Trent, Walker & Wolf 1976 | — | 70.6 | 29–55 |
| Wilson et al. 2012 | **22030378** | 63±14 | 81±11 |
| Robinson et al. 2005 | **15797588** | *not reported* (only failure loads: sMCL 534N/10.2mm, dMCL 194N/7.1mm, PMC 425N/12.0mm) — derived secant ≈52 N/mm (sMCL), an **approximation**, likely an underestimate | — |

Aggregate range **50–85 N/mm**.

**LCL:**

| study | PMID | stiffness (N/mm) | donor age |
|---|---|---|---|
| Trent, Walker & Wolf 1976 | — | 59.8 | 29–55 |
| Wilson et al. 2012 | **22030378** | 59±12 | 81±11 |
| Sugita & Amis 2001 | **11476388** | *not reported* (tensile strength only: LCL 309N, popliteofibular 186N) | — |

Aggregate range **47–71 N/mm**. Trent1976 and Wilson2012 converge remarkably (59.8 vs 59±12)
despite a **36-year gap** and vastly different donor-age cohorts (29–55 vs 81±11) — a strong,
decorrelated anchor.

---

## 4. Symmetric QC — catching the tool's OWN errors before they became the report's errors

Per the watertight mandate ("MACHINE CROSS-CHECK... never eyeballing"), every literature number
that entered the falsifier verdict was treated as a claim needing verification, **including** the
secondary review table's own extraction, not just the primary PMIDs:

1. **Robinson 2005 "stiffness" mislabel, caught.** A first table-read pass reported Robinson2005's
   MCL "stiffness" as 534/194/425 N/mm. This conflicted directly with an independently-fetched copy
   of the paper's *own* abstract, which explicitly calls these numbers "maximum loads... mm mean
   extension" (i.e. failure load in N, not stiffness in N/mm). A targeted re-fetch of the table's
   column structure confirmed: the Stiffness (N/mm) column is **blank** for this study; 534/194/425
   belong under Failure load (N) only. Not propagated as stiffness in this report.
2. **Harner 1995 "no stiffness data," caught and reversed.** A second extraction pass claimed the
   table had no stiffness column for this row. A third, maximally explicit pass (quoting the row
   cell-by-cell with column names) found Stiffness=120.0/57.0 N/mm **and** Failure load=1120.0/419.0
   N as two separate, both-populated columns — reproducing the *first* pass's original numbers
   exactly. 2-of-3 agreement trusted.
3. **Trent 1976, independently re-verified.** Because its numbers support 4 separate ligament
   claims (ACL/PCL/MCL/LCL all from one 7-knee cohort), it was re-fetched a second time,
   independently, specifically in response to catching issues 1–2 above. 2-of-2 extractions agree
   exactly (138.3/179.5/70.6/59.8 N/mm stiffness; 620.8/658.0/515.8/376.6 N failure load).

This is the SAME discipline already logged in this repo's memory
(`subagent-claimed-live-check-was-fabricated`, `orient-by-measurement`) applied to a tool
(`WebFetch`'s own summarization step) rather than a sub-agent — a hit is a hit regardless of which
process produced it.

---

## 5. Falsifier verdict

| ligament | model structural stiffness (N/mm) | literature range (N/mm) | central estimate | model/central ratio | verdict |
|---|---:|---:|---:|---:|:-:|
| ACL | **283.1** | 124.0 – 285.9 | 251.9 | 1.12 | **PASS** — inside range, near young/healthy top |
| PCL | **168.6** | 57.0 – 347.0 | 178.2 | 0.95 | **PASS** — 2 of 3 anchors within 6% |
| MCL | **171.7** | 50.0 – 85.0 | 66.8 | **2.57** | **FAIL** — ~2.6x too stiff |
| LCL | **38.8** | 47.0 – 71.0 | 59.4 | **0.65** | **FAIL** — ~35% too soft |

**The single most decisive test — the MCL:LCL ratio, not either absolute number alone** (a ratio
cancels out any shared systematic bias from the parallel-bundle-summing approximation, so it is a
cleaner test than either magnitude on its own):

| | MCL:LCL ratio |
|---|---:|
| Model (this task) | **4.43** |
| Trent 1976 (same cohort, both ligaments) | 1.18 |
| Wilson 2012 (same cohort, both ligaments; own abstract: "no significant difference") | 1.07 |

Two decorrelated cadaveric studies, 36 years and completely different donor cohorts apart, both
find MCL and LCL structural stiffness in **near parity**. The model's ratio is **3.7–4.1x** larger
than either — a robust, cross-validated departure.

**Root cause, traced (not merely observed)**: Lenhart et al. 2015 (Lenhart RL, Kaiser J, Smith CR,
Thelen DG. "Prediction and Validation of Load-Dependent Behavior of the Tibiofemoral and
Patellofemoral Joints During Movement." Ann Biomed Eng. 2015;43(11):2675-2685. PMID
[25917122](https://pubmed.ncbi.nlm.nih.gov/25917122/) — already the source model for all 4
ligaments per `docs/MECHANISM_KNEE_LIGAMENTS.md` §2)'s own Methods section (fetched live,
open-access via PMC public-access policy, **PMC4886716**) states:

> "The linear stiffness of each ligament were estimated as the product of the average ligament
> cross-sectional area and an assumed ligament elastic modulus of **125 MPa**.\[13\]"

Reference [13] is **Chandrashekar N, Mansouri H, Slauterbeck J, Hashemi J. "Sex-based differences
in the tensile properties of the human anterior cruciate ligament." J Biomech. 2006;39(16):2943-50.
PMID [16387307](https://pubmed.ncbi.nlm.nih.gov/16387307/)** (citation match verified live —
author/topic/year all correct; the abstract itself only gives relative %-differences, so the exact
"125 MPa" figure is taken on Lenhart2015's own reported word, not independently re-derived digit-
for-digit). **125 MPa sits squarely inside the task's own cited ~100–150 MPa collagen-modulus
range** — the material input is reasonable in isolation.

The problem is the word **"each"**: ONE modulus value was applied to ACL, PCL, MCL, *and* LCL
alike, with only MRI-derived cross-sectional area varying between them. But Wilson et al. 2012's
own abstract reports MCL is **~680% wider** (much larger cross-section) than LCL, while the two
ligaments' *measured* structural stiffness is nearly identical (63 vs 59 N/mm) — which requires
MCL's **true** modulus to be substantially *lower* than LCL's to compensate
(`stiffness = E·A/L`). A model that applies one modulus scaled only by CSA mechanically predicts
MCL should come out stiffer than LCL by roughly the CSA ratio — **exactly** the miscalibration
measured independently here. Back-deriving implied bundle CSA from the model's own `k`
(`CSA = k / 125 MPa`) gives a **3.33x** MCL:LCL cross-section ratio (64.0 mm² vs 19.2 mm², side r,
all bundles summed) — same direction as, though smaller magnitude than, Wilson's own ~5.5–6.3x
implied by their reported width×thickness — a supportive, not tightly quantitative, secondary
check.

**Relevance to mis-distributed passive load-sharing (the task's own framing)**: a ~3.7–4.4x
medial-heavy MCL:LCL stiffness imbalance (vs ~1:1 in cadaveric literature) means any passive
coronal-plane (varus/valgus) moment in this model is reacted disproportionately by the medial
side relative to real knees — directly relevant to any free-body-cut analysis that partitions
joint moment across these 4 structures.

---

## 6. Honest gaps

1. **Woo 1991's exact old-group / alternate-orientation numbers not extracted live.** The primary
   abstract (efetch) gave only the young/anatomical figure explicitly (242±28 N/mm). The
   middle-/old-group range and the "tibial/90-90" orientation numbers came from the Peters2018
   secondary table only, cross-checked at the young-group endpoint (exact match) but not
   independently re-verified beyond that — same honest-gap class as this repo's pre-existing
   Baldwin2009 gap in `docs/MECHANISM_LIGAMENT_EPSR_SENSITIVITY.md`.
2. **Race & Amis 1994 and Harner 1995's per-bundle N/mm numbers came from the Peters2018 secondary
   review, not the primary papers' own full text** (both non-open-access; their abstracts report
   strength/kN, not stiffness/N-per-mm). Cross-validated via the failure-load-vs-abstract check
   (§3) and the 2-of-3/2-of-2 re-extraction agreement (§4), but not from the primary source
   documents directly.
3. **Robinson 2005's MCL "secant stiffness" (~52 N/mm) is my own derived approximation**
   (load-to-failure ÷ extension-to-failure), not a directly-reported linear-region tangent
   stiffness — flagged as likely an underestimate, used only as a corroborating check, never as
   the primary MCL anchor.
4. **The parallel-bundle-summing approximation** (needed to compare an 84-bundle discretized model
   against a single-structure cadaveric pull test) assumes bundles act as ideal co-directional
   parallel springs. Validated *indirectly*: applying the identical treatment to ACL and PCL (also
   multi-bundle discretizations from the same source model) reproduces literature closely, giving
   confidence the method itself is not the source of the MCL/LCL miss — but this is not an
   independently-proven-exact geometric fact for every bundle orientation.
5. **MCL at 90° flexion is not fully recruited even at 16mm of virtual pull** (13/16 bundles
   engaged-linear, R²=0.9992) — the reported 134.8 N/mm at that specific angle slightly
   underestimates the group's asymptotic ~172 N/mm capacity; not chased further with a larger pull
   (diminishing returns for the core finding, which is angle-insensitive elsewhere).
6. **WebSearch was unavailable this session** (session-wide quota already exhausted before this
   task's first call) — all literature discovery used WebFetch directly against NCBI eutils,
   Europe PMC, and PMC full text instead, the same method this repo's prior ligament docs already
   used and cite as their own convention.
7. **Single donor/source model (Lenhart2015, "a young healthy adult female" MRI subject)** — this
   task does not re-litigate that already-disclosed limitation (`docs/MECHANISM_KNEE_LIGAMENTS.md`
   §8, `docs/MECHANISM_LIGAMENT_EPSR_SENSITIVITY.md` §9); it only adds the STIFFNESS-axis finding on
   top. Consistent with this: the ACL comparison (model 283 vs literature young/healthy 218–286)
   lines up with a "young, healthy donor" regime specifically, not the full-population range.
8. **Side `r` only** — the config's own documented mirroring convention (Z-negation) means side
   `l` would differ only through small already-disclosed anthropometric asymmetry
   (`docs/MECHANISM_KNEE_LIGAMENTS.md` §2), not re-verified here (out of this task's scope).

---

## 7. Confidence tier

**Cadaveric tensile-test literature** (Woo1991, Trent1976, Noyes&Grood1976, Butler1992,
Race&Amis1994, Harner1995, Wilson2012, Robinson2005 — 8 independent studies), cross-tabulated via a
peer-reviewed systematic review (Peters et al. 2018, PMC5787350) with per-row re-verification where
a first extraction pass looked suspect (§4). All PMIDs verified live this session via NCBI eutils
(esearch/esummary/efetch), Europe PMC, and PMC full text — per the task's explicit instruction, not
recalled from training data.

---

## 8. Files

- `scripts/msk/ligament_stiffness_extract.py` — the single runnable source of every number in this
  doc: config-uniformity audit, the two independent structural-stiffness derivations (closed-form
  upper bound + geometric finite-difference virtual tensile test, with the self-caught sign-fix
  documented inline), and the full literature ledger with machine-computed ratios (nothing
  hand-typed).
- `data/msk_smoketest/ligament_stiffness/ligament_stiffness_evidence.json` — raw model-extraction
  output: config uniformity audit, per-bundle naive structural stiffness, and the full
  finite-difference sweep (5 angles × 4 groups, with engagement checkpoints and R² per fit).
- `data/msk_smoketest/ligament_stiffness/ligament_stiffness_ledger.json` — the consolidated
  deliverable: model numbers + all 8 literature citations with PMIDs and verification notes +
  Lenhart2015's stiffness-derivation trace + the falsifier verdict + all honesty flags.
- Read (unmodified, external anchors): `data/msk_smoketest/knee_ligaments/model_with_knee_ligaments.osim`,
  `scripts/msk/knee_ligament_config.json`, `scripts/msk/add_knee_ligaments.py`,
  `docs/MECHANISM_KNEE_LIGAMENTS.md`, `docs/MECHANISM_LIGAMENT_EPSR_SENSITIVITY.md`.
