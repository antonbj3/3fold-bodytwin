# MECHANISM WITHIN-HOST VIRAL DYNAMICS — target-cell-limited ODE, two regimes (2026-07-22)

Builds and MEASURES the standard target-cell-limited (TCL) within-host viral-dynamics ODE
(dT/dt=−βTV; dI/dt=βTV−δI; dV/dt=pI−cV) + the within-host basic reproductive number R0, on the
two canonical parameter regimes named by the task: **HIV** ART-perturbation decay kinetics
(Ho 1995 / Perelson 1996) and **acute influenza** (Baccam 2006 target-cell-limited fit to human
challenge data). Couples to `docs/MECHANISM_TCELL_ACTIVATION_EXHAUSTION.md` (adaptive/CD8 clearance
and chronic-viral-load-driven exhaustion — that doc's own citations already cite HIV/LCMV viral
persistence topically; V(t) here is the concrete driving trajectory), `docs/MECHANISM_COMPLEMENT_CASCADE.md`
(innate clearance) and `docs/MECHANISM_ACUTE_PHASE_INFLAMMATION.md` (inflammatory response to viral
burden) — no existing anchor-graph node pre-registered this topic (checked: 0 hits for
"viral dynamic|within-host|target-cell-limited" in `data/MECHANISM_ANCHOR_GRAPH.json`), so this is a
new build, not a resolve of a pre-planted stub.

**Script**: `scripts/msk/within_host_viral_dynamics.py` (pre-existing this session, not authored by
this pass — isolation rule "touch only files you create" applies; see the disclosed bug + patch
method in §6). **Evidence**: `data/msk_smoketest/within_host_viral_dynamics/within_host_viral_dynamics_results.json`
(freshly generated this session — did not exist before).

## The falsifiers, verdicts stated up front (disk-verified, nothing hidden)

> **F2 (influenza)**: does the model reproduce the MEASURED target-cell-limited peak timing
> (~2–3 days) and clearance (~5–7 days)?

**PASS, cleanly.** Peak day = **2.83** (task band [2,3]; Carrat 2008's independent 56-study/1280-participant
pooled review reports peak on day 2.0, |Δ|=0.83d). Clearance day (1 TCID50/mL floor) = **6.84d**
(task band ~5–7d; Carrat 2008's implied band 5.3–5.8d, |Δ|=1.29d). The **forced adversary** — target-cell
depletion surgically removed (T clamped constant), all rate constants otherwise identical — shows
**no interior peak at all** (still monotonically rising at the end of a 15-day window), while the
real model peaks cleanly at day 2.83: this is the single most decisive, structural result in this
doc (§4). Robustness: sweeping β,p (log-uniform) and δ,c (uniform) inside Baccam's OWN reported 95%
CIs, 10/12 draws still show an interior peak in a wide sanity band.

> **F1 (HIV)**: does the model reproduce the MEASURED plasma-virus biphasic decay after ART —
> first-phase t1/2 ~1–2 days = infected-cell loss δ? Decorrelated check: does an INDEPENDENT
> primary-infection (untreated) lifespan estimate match?

**PARTIAL — diagnosed, not hand-waved.** The asymptotic/late-time decay rate exactly recovers
δ (fit t1/2 = 1.600d vs input 1.6d), robustly (12/12 parameter draws at ±30%), and the decorrelated
check **passes**: Stafford 2000's independent primary-infection (untreated, 10 drug-naive patients,
nonlinear-least-squares to natural rise-fall kinetics — a completely different data regime and
statistical method) lifespan estimate of 2.5d matches Perelson 1996's ART-decay-derived exact
lifespan of 2.308d within **8.3%** relative error (pre-registered tolerance ±30%). The model also
shows a genuine, non-single-exponential **shoulder** at treatment onset (dV/dt(0)≈0, a signature no
single-exponential adversary can produce). **But** the specific pre-registered "biphasic slope
ratio ≥2×" gate (comparing instantaneous log-slope at t=0.5d vs t=14d) **FAILS as literally coded**
— measured ratio = **0.74×**, the wrong direction. §3 diagnoses this via OODA: it is a gate-construction
mismatch (an unstated fast-early/slow-late ordering assumption that does not hold under the
physically-correct quasi-steady-state ART initial condition), not a refutation of the δ/c
decomposition — and a corrected, direction-free version of the same adversary-beating test (§3.3)
falls decisively (4×–35× depending on how early you probe). A second, independent fragility finding
(§3.4, not task-pre-registered): the "aggregate matches Ho/Wei composite" gate technically passes
but sits within 2% of its own floor across every window choice tried — HELD OPEN, disclosed.

## 0. Citations — every PMID/DOI verified LIVE THIS SESSION, independently of the pre-existing
script's own docstring claims (NCBI E-utilities esummary+efetch raw abstract text, Crossref API,
Europe PMC search API, PMC XML block-message, one independent WebFetch pass cross-checking the
Baccam parameter table). This repo's own measured citation-drift rate from memory is ~62–67%
(sibling docs); the check here is a **second, independent** verification pass over a script that
already claimed its own citations were live-checked — treating that claim as unproven until
reproduced, per this repo's own "subagent-fabricated-live-check" lesson.

| # | Citation | PMID / DOI | Role | Verified |
|---|---|---|---|---|
| 1 | Perelson AS, Neumann AU, Markowitz M, Leonard JM, Ho DD (1996). "HIV-1 dynamics in vivo: virion clearance rate, infected cell life-span, and viral generation time." *Science* 271(5255):1582-6. | **8599114**, DOI 10.1126/science.271.5255.1582 | Load-bearing HIV constants: δ,c | Raw abstract re-fetched verbatim: "life-span of 2.2 days (half-life t1/2=1.6 days)"; "0.3 days (t1/2=0.24 days)"; "10.3×10⁹ virions per day"; generation time 2.6d; min lifecycle 1.2d — **all match exactly**. |
| 2 | Perelson AS, Essunger P, Cao Y, et al. (1997). "Decay characteristics of HIV-1-infected compartments during combination therapy." *Nature* 387(6629):188-91. | **9144290**, DOI 10.1038/387188a0 | Cross-paper δ,c confirmation + 2nd-phase scope limit | Raw abstract: "≤6 hours" free-virus t1/2, "1.6 days" infected-cell t1/2 (exact cross-paper match to #1); 2nd phase "1-4 weeks" (long-lived cells), "0.5-2 weeks" (latent activation) — **all match exactly**. |
| 3 | Ho DD, Neumann AU, Perelson AS, Chen W, Leonard JM, Markowitz M (1995). "Rapid turnover of plasma virions and CD4 lymphocytes in HIV-1 infection." *Nature* 373(6510):123-6. | **7816094**, DOI 10.1038/373123a0 | Composite single-slope anchor | Raw abstract: "mean half-life, 2.1 ± 0.4 days" — **matches exactly**. |
| 4 | Wei X, Ghosh SK, Taylor ME, et al. (1995). "Viral dynamics in human immunodeficiency virus type 1 infection." *Nature* 373(6510):117-22. | **7529365**, DOI 10.1038/373117a0 | Independent-cohort, same-day/same-issue decorrelated check | Raw abstract: "half-life approximately 2 days" — **matches exactly**; confirmed same *Nature* issue/date (1995 Jan 12) as #3 via esummary. |
| 5 | Stafford MA, Corey L, Cao Y, Daar ES, Ho DD, Perelson AS (2000). "Modeling plasma virus concentration during primary HIV infection." *J Theor Biol* 203(3):285-301. | **10716909**, DOI 10.1006/jtbi.2000.1076 | Decorrelated primary-infection lifespan + HIV-side TCL-breakdown finding | Raw abstract: "an estimate of 2.5 days for the average lifespan of productively infected cells during primary infection... consistent with results obtained by drug perturbation experiments"; **also**: "the kinetics of the subsequent fall and recovery... in some patients are not consistent with the predictions of the target-cell-limited model" — **both match exactly**. |
| 6 | Markowitz M, Louie M, Hurley A, et al. (2003). "A novel antiviral intervention results in more accurate assessment of HIV-1 replication dynamics and T-cell decay in vivo." *J Virol* 77(8):5037-8. | **12663814**, DOI 10.1128/jvi.77.8.5037-5038.2003 | Disclosed drug-efficacy-confound illustration | Raw abstract: "half-lives of virus-producing CD4+ T cells, 0.7 day" — **matches exactly**. (A previously-recalled PMID 12525667 was checked live and confirmed a WRONG/unrelated duck hepatitis-B paper.) |
| 7 | Baccam P, Beauchemin C, Macken CA, Hayden FG, Perelson AS (2006). "Kinetics of influenza A virus infection in humans." *J Virol* 80(15):7590-9. | **16840338**, PMCID PMC1563736, DOI 10.1128/JVI.01623-05 | Influenza TCL model fit + R0 + parameter table | Raw abstract: "delay of ~6h... continue to do so for ~5h... average lifetime ~11h... half-life of free infectious virus ~3h... R0... ~22 new productive infections... For about 50% of patients... two peaks... interferon" — **all match exactly**. Full-text XML blocked by publisher ("does not allow downloading... in XML form" — confirmed directly); parameter table independently re-fetched via WebFetch (my own pass, see §4) and matches the script's values exactly, **plus found T0=4×10⁸ cells is a literally-stated FIXED input** (upgrade over the script's own "illustrative/back-solved, disclosed gap" framing). |
| 8 | Carrat F, Vergu E, Ferguson NM, et al. (2008). "Time lines of infection and disease in human influenza: a review of volunteer challenge studies." *Am J Epidemiol* 167(7):775-85. | **18230677**, DOI 10.1093/aje/kwm375 | Independent, zero-model-fitting, pooled anchor | Raw abstract: "consistently peaked on day 2... duration of viral shedding... 4.80 days (95% CI: 4.31, 5.29)" — **matches exactly**. (A previously-recalled PMID 18029396 confirmed a WRONG/unrelated Schistosoma paper.) |
| 9 | Perelson AS, Nelson PW (1999). "Mathematical Analysis of HIV-1 Dynamics in Vivo." *SIAM Review* 41(1):3-44. | DOI 10.1137/s0036144598335107 (not MEDLINE-indexed) | Canonical R0/model-structure reference | Crossref-confirmed: title/authors/journal/year all match exactly. (PMID 32488570, initially considered, confirmed a DIFFERENT 2020 Perelson paper — "Continuous and discrete modeling of HIV-1 decline on therapy," *J Math Biol* — correctly excluded.) |
| 10 | Miao H, Xia X, Perelson AS, Wu H (2011). "On Identifiability of Nonlinear ODE Models and Applications in Viral Dynamics." *SIAM Review* 53(1):3-39. | DOI 10.1137/090757009 | (β,p) practical-identifiability formal reference | Crossref-confirmed match. **Additionally verified via Crossref's reference list** (not an LLM summary): reference key R8 = DOI 10.1128/JVI.01623-05, i.e. this paper directly cites Baccam 2006 (#7) — confirms it is the field-standard identifiability reference for exactly this model class. |

## 1. The model + R0 — geometric (spectral) derivation, independently re-derived

The post-ART HIV system (T held at its slowly-varying pre-treatment value; ART blocks new
infection) is exactly linear and lower-triangular:

```
d/dt [I]   [-δ   0] [I]
     [V] = [ p  -c] [V]
```

Its eigenvalues are its diagonal entries, **{-δ,-c}, exactly** (numeric eigendecomposition vs the
diagonal: max error 0.00e+00; analytic closed-form V(t)/K = e^(-δt) − (δ/c)e^(-ct) vs `solve_ivp`
numeric integration: max relative error 5.29e-10). This is the geometric heart of "biphasic decay":
the two observed decay-rate regimes are literally the two eigenmodes of a triangular linear map, not
a curve-fitting artifact.

For the influenza eclipse-phase model (T,I1,I2,V; I1=eclipse/non-producing, I2=producing), R0
threshold behavior was **re-derived here independently by hand**, not merely read off the script's
own numerical check. The disease-free-equilibrium Jacobian for (I1,I2,V) is

```
M = [ -k    0   βT0 ]
    [  k   -δ    0  ]
    [  0    p   -c  ]
```

Cofactor expansion along row 1 gives the characteristic polynomial
det(M−λI) = −(k+λ)(δ+λ)(c+λ) + βT0·k·p. At λ=0: k·δ·c = βT0·k·p ⟺ **βpT0/(δc) = 1** — the
bifurcation locus is exactly the standard R0=1 threshold, **independent of k** (the eclipse rate
drops out of the threshold condition entirely, though it still sets the generation time). This was
confirmed numerically (DFE unstable at R0≈21.8, stable at R0=0.1 — sign flip exactly at the
derived threshold) — the numeric and the by-hand derivation are two independent routes to the same
conclusion, not one checking the other's arithmetic.

A second geometric aside, also re-derived here: the strict 2-type next-generation-matrix (NGM) on
(I2,V) has K = F·Σ⁻¹ = [[0, βT0/c],[p/δ, 0]] — a pure 2-cycle structure, whose eigenvalues are
±√(ab) for any matrix [[0,a],[b,0]] (elementary 2×2 characteristic-polynomial fact: λ²−ab=0). Hence
NGM spectral radius = √(βpT0/(δc)) = √(standard R0) **exactly** (confirmed: max abs error 0.00e+00
between √21.775 and the computed spectral radius 4.6664). This is the same structural reason
Macdonald's classic malaria R0 formula carries a square root (human→mosquito→human is a 2-step
cycle) — both conventions agree on the R0=1 threshold; they differ away from 1 by construction, not
by error.

## 2. Influenza (Baccam 2006) — clean pass, strongest evidence in this doc

| Quantity | Model (this run) | Anchor | Match |
|---|---|---|---|
| Eclipse delay 1/k | 0.25d = 6.0h | Abstract: "~6h" | exact |
| Producing duration 1/δ | 0.192d = 4.6h | Abstract: "~5h" | close |
| Avg infected lifetime 1/k+1/δ | 0.442d = 10.6h | Abstract: "~11h" | close |
| Virus t1/2 = ln2/c | 0.133d = 3.2h | Abstract: "~3h" | close |
| Peak day | 2.83d | Task [2,3]; Carrat2008 day 2.0 | pass (Δ=0.83d ≤ 1.0d gate) |
| Clearance day | 6.84d | Task ~5–7d; Carrat2008 implied 5.3–5.8d | pass (Δ=1.29d ≤ 2.0d gate) |
| R0 (βpT0/δc) | **21.78** | Baccam reported ≈22 (CI 10–46) | pass, and no longer circular (see below) |

**T0 upgrade, found independently this session**: the script's own docstring disclosed T0=4×10⁸ as
merely "back-solved consistent with Baccam's own R0=22... NOT an independent measurement... disclosed
gap." My own independent WebFetch pass on the PMC page found the paper **directly states**: "The
parameter T0 was held fixed at a value of 4×10⁸ cells." This converts the R0 check from an
internal-consistency tautology into a genuine external anchor: combining the independently-verified
β=3.2×10⁻⁵, p=4.6×10⁻², δ=5.2/d, c=5.2/d (all re-confirmed exactly via a fresh WebFetch pass, matching
the script's values) with this literally-stated T0 gives R0=21.78 — landing almost exactly on
Baccam's own separately-reported "≈22" — five independently-sourced numbers agreeing, not one
number solved to force agreement with another.

**Forced adversary (the decisive result)**: target-cell depletion surgically removed (dT/dt=0
instead of −βTV, all other rates identical) → **no interior peak at all** — V rises monotonically
through the full 15-day window (adversary's peak index sits at the very last grid point). The real
model peaks cleanly and falls away by day 2.83. This is a clean, structural, binary
falsification — a peak categorically cannot appear here "for unrelated reasons" (e.g. from the δ/c
decay constants alone); it requires the target-cell depletion term specifically.

**Identifiability (symmetric QC, held open, not swept into the pass)**: β 95% CI spans 6.0×10⁻⁶–
1.7×10⁻⁴ (28.3×), p spans 1.2×10⁻²–1.7×10⁻¹ (14.2×) — both independently re-confirmed. By contrast
k's CI (3.0–5.2/day, found in my own independent pass — **not present in the pre-existing script**,
a genuine addition) is only 1.7× wide — δ, c, k are comparatively well-identified; β, p are the
practically-unidentifiable pair, exactly the pattern the field-standard reference (Miao 2011, its
Baccam-citation confirmed via Crossref, not an LLM guess) formally treats. Baccam's own words,
re-confirmed verbatim on a fresh, independent fetch: "the surface tends to be flat at the minimum" —
character-for-character identical across the pre-existing script's pass and mine, the kind of
differential-replication check that catches fabrication (contrast: the R0 formula, confirmed on a
**third** independent pass here, is genuinely an unreadable embedded graphic, not machine-extractable
text — corroborating, not merely repeating, the pre-existing script's own catch of a first-pass
fabrication of that formula).

Robustness: sampling β,p log-uniform and δ,c uniform inside their own reported 95% CIs (12 draws):
10/12 (83%) retain an interior peak inside a wide [0.5,10]d sanity band — a majority-pass, non-tautological,
literature-CI-anchored uncertainty propagation (not an arbitrary ±X%).

## 3. HIV (Perelson/Ho) — OODA-forced diagnosis of the failing gate, not a hand-waved negative

**3.1 Observe.** The disk-verified run (§6) reports `F1d_forced_adversary_falls = False`
(early/late slope ratio = **0.74×**, gate wanted ≥2.0×) and `F1_hiv_overall_pass = False`. A
one-shot "honest negative" here would be premature — the underlying δ/c decomposition is otherwise
robustly confirmed (§0 citations, F1a/F1c/robustness all pass) — so the failing sub-gate itself
needed instrumenting, not just reporting.

**3.2 Orient (the crux).** I traced the model's exact log-slope trajectory by hand
(closed form V(t)/K = e^(-δt) − (δ/c)e^(-ct), δ=0.4332/d, c=2.8881/d):

| t (days) | 0.001 | 0.05 | 0.1 | 0.5 | 1 | 2 | 4 | 14 |
|---|---|---|---|---|---|---|---|---|
| log-slope (/day) | 0.001 | 0.058 | 0.107 | 0.320 | 0.401 | 0.431 | 0.4332 | 0.4332 |

The slope rises **monotonically from 0 up to δ, never exceeding it** — a shoulder that ramps up to
the asymptotic rate, not a fast(c)-phase that decays down to a slow(δ)-phase. The gate's "early ≥2×
late" criterion silently assumed the opposite (fast-first) ordering. That ordering is inverted here
because the physically-correct ART initial condition is quasi-steady-state (dV/dt(0)=0 — interrupting
an ongoing steady infection, not perturbing a system away from one) — under this condition the
c-eigenmode's coefficient is small (δ/c≈0.15) and always subtracts from, never leads, the δ-eigenmode.
This is a **gate-construction mismatch**, diagnosed, not a refutation of δ or c themselves (both of
which independently, exactly match their literature values elsewhere in this same run).

**3.3 Decide + Act (the corrected, direction-free test).** The right forced-adversary comparison
is not internal (early-vs-late within the real model) but matched-endpoint: a single-exponential
adversary calibrated to the SAME asymptotic rate δ has, by construction, slope=δ at every t,
including t≈0. Comparing to the real model at small t:

| t (days) | 0.01 | 0.05 | 0.1 | 0.5 |
|---|---|---|---|---|
| adversary/real ratio | 35.2× | 7.5× | 4.0× | 1.35× |

The single-exponential adversary is wrong by 4×–35× depending on how early you look — a decisive,
monotonically-strengthening falsification, just constructed via a matched-endpoint comparison
instead of the internal early/late ratio the pre-existing gate literally coded. This is the
**direction-free** version of the already-passing `F1d_shoulder_present` check (dV/dt(0)/(V0·c) <1%)
— both say the same thing correctly; the ratio-with-assumed-direction sub-test was the one
mis-specified.

**3.4 A second fragility, found independently (not task-pre-registered).** `F1b_aggregate_matches_ho_wei_composite`
technically passes (1.607d, inside the pre-registered [1.6,2.6] band) but a window-sensitivity sweep
shows it is **never robust** — every window from [0,7]d to [0,28]d and [2,14]d/[3,14]d gives
1.600–1.626d, always within 2% of the band's own floor (1.6d = pure δ), never approaching Ho1995's
central 2.1d or Wei1995's ~2.0d. Most likely explanation (disclosed, not resolved): Ho1995/Wei1995
report a **mean of per-patient** single-exponential fits (their ±0.4d SD implies per-patient fitting
then averaging), a different statistical operation than fitting once to a single simulated
"aggregate" trajectory — these need not agree exactly even if δ,c themselves (independently,
exactly literature-matched) are correct. HELD OPEN.

**Decorrelated check (passes cleanly)**: Stafford 2000's independent primary-infection lifespan
estimate (2.5d, from 10 drug-naive patients' natural rise-fall kinetics, nonlinear-least-squares) vs
Perelson 1996's ART-decay-derived exact lifespan (1/δ=2.308d) — relative error **8.3%**, inside the
pre-registered ±30% tolerance, and genuinely decorrelated (different cohort, different data regime,
different fitting method).

## 4. Symmetric QC — HELD OPEN (forced-steelman burden applied equally to passes and fails)

1. **TCL breaks down under immune dominance** — influenza: Baccam's OWN abstract, re-confirmed
   verbatim: "For about 50% of patients, the curve of viral titer versus time has two peaks. This
   bimodal behavior can be explained by incorporating the antiviral effects of interferon" — the
   pure TCL model built/run here does not reproduce this. HIV: Stafford 2000's OWN abstract,
   re-confirmed verbatim: "the kinetics of the subsequent fall and recovery... in some patients are
   not consistent with the predictions of the target-cell-limited model," requiring CTL/cytokine
   mechanisms. Both machine-quoted from primary sources, both HELD OPEN.
2. **(β,p) practical non-identifiability** — Baccam's own flat-minimum quote (independently
   re-confirmed verbatim) + 28.3×/14.2× CI ranges (independently re-confirmed) + Miao 2011's formal
   treatment (its Baccam-citation confirmed via Crossref reference list, not an LLM guess). HELD
   OPEN, illustrated further by the Markowitz 2003 (δ t1/2=0.7d) vs Perelson 1996 (1.6d) discrepancy
   — a real, disclosed, drug-efficacy-assumption confound (a more-potent regimen deduces a
   substantially faster cell-death rate; naive single-regimen δ estimates may be biased by residual
   ongoing infection).
3. **The literal 3-state TCL model cannot sustain a chronic set point** — dT/dt=−βTV alone (no
   source term) is strictly decreasing whenever T,V>0 (machine-checked sign), so T→some T∞≥0 and
   the infection structurally burns out (SIR-like), matching the acute/influenza picture but NOT a
   genuine chronic (months-years) HIV set point, which requires an added source term
   dT/dt=λ−dT−βTV (Perelson & Nelson 1999) — structurally necessary for a pre-treatment steady state
   to exist at all, even though (time-scale separation: λ,d act over months, vastly slower than the
   days-scale δ,c decay analyzed here) it is negligible for the fast post-ART decay itself. HELD OPEN.
4. **[Found this session, not task-pre-registered]** The F1d gate's implicit fast-early/slow-late
   ordering assumption (§3.2) and the F1b aggregate-fit's floor-hugging fragility (§3.4) — both
   diagnosed with hand-derived, machine-checked numbers above, both HELD OPEN rather than silently
   patched into a false PASS or silently ignored as an inconvenient FAIL.

## 5. Confidence tier

**In-vivo-anchored** for the core rate constants and timing claims: δ, c (HIV), the HIV
decorrelated lifespan cross-check, and all influenza timing/R0/identifiability numbers are anchored
to directly-quoted, live-re-verified primary-source numbers from real patient (HIV) and human
challenge-volunteer (influenza) data — Perelson/Ho/Wei/Stafford/Markowitz and Baccam/Carrat, with
every load-bearing quote independently re-fetched and checked character-for-character against the
raw NCBI/Crossref/PMC record this session, not merely re-stated from a prior pass. The
model-**structural** caveats (chronic-set-point scope limit, the biphasic-ordering gate mismatch,
(β,p) identifiability) are **derived/diagnosed**, not independently data-anchored, and are explicitly
HELD OPEN rather than resolved — a real, disclosed boundary of this confidence tier, not folded into
the headline number.

## 6. Disclosed bug + patch method (isolation-respecting)

The pre-existing script, as committed to disk, **crashes** on line 477–478:
`flu_dfe_jacobian(beta_subthreshold, FLU_K_ECLIPSE_PER_DAY, FLU_DELTA_PER_DAY, FLU_C_PER_DAY, FLU_T0_ILLUSTRATIVE)`
omits the trailing `p` argument the function signature requires (`def flu_dfe_jacobian(beta, k,
delta, c, T0, p)`), raising `TypeError: flu_dfe_jacobian() missing 1 required positional argument:
'p'`. This means the script, before this session, had **never actually been run to completion** —
every gate value in this doc was previously an unverified, un-disk-reproduced claim. Per isolation
("touch ONLY files you create"), the original file was **not edited**: it was patched **in-memory
only** (source text read, the exact buggy line replaced with a version appending `, FLU_P`, then
`exec()`'d in a fresh namespace) to produce the evidence JSON on disk. The one-line fix for whoever
owns that file: append `, FLU_P` as the 6th argument to that specific call.

## State (paths)

- Script (pre-existing, unedited; contains the disclosed line-477 bug): `source_repository/scripts/msk/within_host_viral_dynamics.py`
- Evidence JSON (freshly generated this session via the in-memory-patched run): `source_repository/data/msk_smoketest/within_host_viral_dynamics/within_host_viral_dynamics_results.json`
- This doc: `source_documents/MECHANISM_VIRAL_DYNAMICS.md`
- Coupling targets (pre-existing, read-only, not edited): `source_documents/MECHANISM_TCELL_ACTIVATION_EXHAUSTION.md`, `source_documents/MECHANISM_COMPLEMENT_CASCADE.md`, `source_documents/MECHANISM_ACUTE_PHASE_INFLAMMATION.md`
