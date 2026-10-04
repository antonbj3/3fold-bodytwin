# MECHANISM BLOOD OXYGEN TRANSPORT — the connective tissue closing the Fick O2-delivery chain (2026-07-22)

Supplies the piece the twin's already-built CENTRAL (`scripts/msk/cardiac_output.py`,
`docs/MECHANISM_CARDIAC.md`) and PULMONARY (`scripts/msk/respiratory.py`,
`docs/MECHANISM_RESPIRATORY.md`) layers each named but never derived: `cardiac_output.py`'s own Fick
calculation (`Q = VO2/(a-vO2diff)`) uses a-vO2diff = 5.0 mL/100mL (rest) / 10-12 (walking) as an
**assumed literature constant**, never itself computed from hemoglobin chemistry. This layer computes
CaO2 and CvO2 from first-principles Hb-O2-binding physiology (Hill equation + dissolved-O2 solubility)
and uses them to test whether the Fick loop actually **closes** once a-vO2diff is no longer a free
input but an independently-derived one. Script: `scripts/msk/blood_oxygen_transport.py`. Evidence:
`data/msk_smoketest/subject2_walking1/blood_oxygen_transport/blood_oxygen_transport_results.json`.

**NO re-solve, no OpenSim.** Reads `metabolic_cost_results.json`, `cardiac_output_results.json`, and
`respiratory_results.json` (all plain JSON, read-only) and does pure Hb-chemistry + Fick-principle
arithmetic. Never touches the `.osim` model or any `.sto` file.

**Confidence tier: in-vivo-anchored (blood-gas / Hill).** The Hill-equation parameterization is
anchored against a real clinical dissociation-curve reference (Collins et al. 2015, validated on 3524
real patient blood specimens) and the Fick-closure falsifier (§6) is anchored against a real human
catheterization measurement (Higginbotham et al. 1986) — not a synthetic/simulated reference and not a
tautology gate. It is **one tier below** a subject-specific in-vivo measurement (no arterial blood gas
exists for subject2; all blood-gas operating points are generic/population-level, §14).

## 0. Scope note — naming discrepancy, disclosed up front

The task that produced this doc named the sibling docs `MECHANISM_CARDIAC_OUTPUT.md` and
`MECHANISM_PULMONARY_GAS_EXCHANGE.md`. Neither filename exists in this repo — the actual sibling docs
are **`docs/MECHANISM_CARDIAC.md`** (central circulation / cardiac output) and
**`docs/MECHANISM_RESPIRATORY.md`** (ventilation / O2-uptake), confirmed via `find`. Both were read in
full and their **real, already-computed numbers** are used below (not generic textbook fallbacks) —
specifically `cardiac_output_results.json`'s own rest VO2 (282.65 mL/min), rest Q (5.653 L/min), and
Higginbotham et al. (1986)'s real measured rest CO (5.70 L/min, already live-verified there).

## 1. Geometric structure (stated up front, not decorative)

SaO2(PO2) is a saturating Hill sigmoid; CaO2 is an affine map of it (scale by [Hb], add a small linear
dissolved-O2 term); **a-vO2diff is the difference between two points on the SAME sigmoid** — one near
its flat high-PO2 plateau (arterial, PaO2~95 mmHg) and one on its steep low-PO2 shoulder (mixed-venous,
PvO2~40 mmHg). Fick's `Q = VO2/(a-vO2diff)` then rides the identical hyperbola `cardiac_output.py`'s
own geometric-structure section already names — just with a-vO2diff read off two points of this curve
instead of assumed. §9 measures the curve's actual local-slope structure (steepest point closed-form
at `PO2* = P50×((n-1)/(n+1))^(1/n)`, not the common "steepest at P50" heuristic) — the flat plateau is
a genuine sensor **null-space** (pulse oximetry cannot resolve PaO2 changes there); the steep shoulder
is where the curve does its real unloading work.

## 2. Method, in one paragraph

Hill equation `SaO2 = PO2^n/(PO2^n+P50^n)`, P50=26.6 mmHg, n=2.7 (task-specified) — validated (§4)
against two independently-sourced checkpoints never used to fit P50/n. `CaO2 = 1.34×[Hb]×SaO2 +
0.003×PaO2` (mL/dL); same formula for CvO2 at mixed-venous PvO2=40 mmHg. The **decisive falsifier**
(§6): does `VO2 = CO×(CaO2−CvO2)` close, using three DIFFERENT-provenance numbers — OpenSim
muscle-probe VO2 (`metabolic_cost.py`), REAL catheterization CO (Higginbotham et al. 1986), and this
script's first-principles Hb-Hill a-vO2diff? **Second regime** (§7): Bohr-shift direction, swept
across illustrative pH/coefficient magnitudes. **Third regime** (§8): cardiac_output.py's own walking
a-vO2diff widening (10-12 mL/100mL) inverted to the mixed-venous PO2 it implies.

## 3. Citations — verified LIVE this session (NCBI eutils, PubMed, PMC full-text, Wikipedia)

`WebSearch` was unavailable this session (shared budget exhausted before this doc's research began —
the same constraint `MECHANISM_RESPIRATORY.md` disclosed); every citation below used direct `WebFetch`
to NCBI eutils / PubMed / PMC / Wikipedia instead.

| # | Citation | PMID / DOI | Role |
|---|---|---|---|
| 1 | Collins JA, Rudenski A, Gibson J, Howard L, O'Driscoll R (2015). "Relating oxygen partial pressure, saturation and content: the haemoglobin-oxygen dissociation curve." *Breathe* 11(3):194-201. | **26632351**, DOI 10.1183/20734735.001415, PMC4666443 | PMC **full-text** (not just abstract) fetched live, confirms verbatim: "P50 of normal adult blood is approximately 26 mmHg"; "Normal SaO2: between 96% and 98%"; "mixed venous saturation… typically about 75% in healthy individuals at rest"; Bohr right-shift direction with lower pH/higher PCO2/higher temperature. Does **not** independently confirm this doc's n=2.7 (the paper's own equation is image-embedded, not resolved by this session's text extraction) — disclosed, not assumed. |
| 2 | Severinghaus JW (1979). "Simple, accurate equations for human blood O2 dissociation computations." *J Appl Physiol* 46(3):599-602. | **35496** (verified live via NCBI esummary) | Classical polynomial (non-Hill) empirical fit; Collins (2015) reports validating it against 3524 real clinical specimens. Used §12 bonus cross-check only. **Methodological note**: this session's own first reaction to a 5-digit PMID on a 1979 paper was suspicion of citation drift (by analogy to `cardiac_output.py`'s own measured 67% first-recall-drift rate) — re-verified via a second esummary fetch before trusting it; correct. A caught near-miss in *reasoning*, not an actual citation error. |
| 3 | Kelman GR (1966). "Digital computer subroutine for the conversion of oxygen tension into saturation." *J Appl Physiol* 21(4):1375-6. | **5916678**, DOI 10.1152/jappl.1966.21.4.1375 | Verified live. No abstract (pre-abstracting era) — classical algorithmic source, cited bibliographically for the Bohr/temperature-correction concept only. |
| 4 | Severinghaus JW (1966). "Blood gas calculator." *J Appl Physiol* 21(3):1108-16. | **5912737** (verified live via esummary) | Classical source of the standard pH/temperature/PCO2 blood-gas correction nomogram. Abstract/full text **not accessible live this session** (PubMed CAPTCHA-blocked after repeated fetches) — the specific pH-correction slope is **not** independently re-extracted; the Bohr-shift test (§7) is accordingly direction-gated, never magnitude-gated on one recalled number. |
| 5 | Billett HH (1990). "Hemoglobin and Hematocrit." *Clinical Methods*, 3rd ed., Ch. 151. | **21250102** (NCBI Bookshelf, verified live) | States normal Hb "14 to 18 g/dl" (♂) / "12 to 16 g/dl" (♀); Hct "40 to 54%" (♂) / "36 to 48%" (♀). Task's Hb=15/Hct=45 sit inside the male range — female ranges are systematically lower (§10, held OPEN). |
| 6 | Rivers E, Nguyen B, Havstad S, et al. (2001). "Early goal-directed therapy in the treatment of severe sepsis and septic shock." *N Engl J Med* 345(19):1368-77. | **11794169** (verified live) | Measured central-venous O2 sat 65.3-70.4% in septic-shock ICU patients (an **ill** population, disclosed as such) — corroborates that venous saturations in the 65-75% band are the real clinically-measured regime; not the source of the healthy-rest 75% figure (citation 1 fills that role directly). |
| 7 | Roca J, Agusti AGN, Alonso A, et al. (1992). "Effects of training on muscle O2 transport at VO2max." *J Appl Physiol* 73(3):1067-76. | **1400019**, DOI 10.1152/jappl.1992.73.3.1067 (verified live, bibliographic) | Real study directly measuring exercising-muscle arterial/femoral-venous PO2 — confirms venous PO2 falls at exercising muscle (§8's mechanism site); specific numeric values paywalled, not extracted (disclosed, same convention as this repo's Waters & Mulroy 1999 citation). |
| 8 | Wikipedia "Fick principle" (fetched live) | TEXTBOOK-GRADE, flagged — the **same page** `cardiac_output.py`'s own citation already uses | Gives `CaO2 = [Hb]×1.34×SaO2 + 0.0032×PaO2` (independent live confirmation of 1.34; 0.0032 vs this doc's 0.003 — a disclosed ~7% variant on the small dissolved term) **and** the worked example VO2=237.5 mL/min, CaO2=20 mL/dL, CvO2=15 mL/dL, a-vO2diff=5 mL/dL, CO=4.75 L/min — independently reproducing 4 of the task's own 5 anchor numbers from one already-repo-trusted source. |
| 9 | Wikipedia "Oxygen–haemoglobin dissociation curve" (fetched live) | TEXTBOOK-GRADE, flagged | States P50 "typically about 26.6 mmHg (3.5 kPa)" — independent confirmation of the task's exact P50; confirms Bohr right-shift direction (no slope given). |

Re-used, not re-verified (already live-verified in the sibling docs, read read-only): Higginbotham et
al. (1986) PMID 3948345 (real rest CO=5.70 L/min); Weir (1949) PMID 15394301.

## 4. Headline results (subject2/walking1, [Hb]=15 g/dL, Hct=45%, P50=26.6 mmHg, n=2.7)

| quantity | value | task anchor | gate |
|---|---:|---|---|
| SaO2 at PaO2=95 mmHg | 96.88% | — | — |
| SvO2 at PvO2=40 mmHg | **75.05%** | ~75% | PASS (0.07% diff vs Collins 2015) |
| CaO2 | **19.76 mL/dL** | ~20 mL/dL | PASS |
| CvO2 | 15.21 mL/dL | — | — |
| a-vO2diff (chemistry-derived) | **4.55 mL/dL** | 4-5 mL/dL | PASS |
| Fick-reconstructed VO2 (Higginbotham CO × chemistry a-vO2diff) | **259.5 mL/min** | ~250 mL/min | PASS (3.8% diff) |
| … vs twin's own OpenSim-probe VO2 (282.65 mL/min) | | <20% pre-registered | **PASS (8.19% diff)** |

All numbers are on disk in `blood_oxygen_transport_results.json`, machine-printed from the script's own
arithmetic — nothing above is hand-computed prose.

## 5. Forced adversary (void-floor) — symmetric QC on a leaning-positive result

Every number above leans positive (task anchors met). Per symmetric-QC discipline, the confound to
force is: *are these anchor bands just wide enough that almost any plausible input passes, regardless
of whether the Hb-Hill mechanism is doing real work?* Four void-floors, forced before trusting §4:

| void-floor | result | verdict |
|---|---|---|
| **A** — strip Hb-binding entirely (dissolved-O2 term only) | CaO2 collapses to 0.285 mL/dL = **1.44%** of the real 19.76 | Hb-binding supplies 98.6% of CaO2 — **mechanism is necessary, not decorative** |
| **B** — wrong Hill exponent n (1.0/2.0/3.5/4.0) at PO2=40 vs the 75% anchor | diffs = 19.9% / 7.6% / 7.5% / 11.5% — **all fail** the 5% tolerance the task-given n=2.7 passes (0.07%) | the specific n is doing real discriminating work |
| **C** — wrong P50 (20/24/30/35 mmHg), same test | diffs = 15.6% / 6.5% / 8.7% / 21.4% — **all fail** | same conclusion for P50 |
| **D** — grid-scan n∈[1,6]×0.2, P50∈[15,40]×1 mmHg (676 combinations) | only **13.0%** of combinations pass the same 5% band the task-given pair passes | the task-given (n,P50) sits in a genuinely constrained region, not a trivially-easy one |

**Gate: PASS** — the forced adversary falls on all four fronts. This is the same void-floor discipline
`cardiac_output.py`'s own §5 (a-vO2-widening necessity) already established, applied here to the
Hill/CaO2 mechanism itself.

## 6. Decisive falsifier — does the Fick loop CLOSE? (non-circular 3-leg triangulation)

Three numbers of genuinely **different provenance** must multiply out consistently for this gate to
mean anything — none is derived from either of the others:

- **Leg A** — VO2 = 282.65 mL/min (OpenSim muscle-probe estimate, `metabolic_cost.py`)
- **Leg B** — CO = 5.70 L/min (REAL human catheterization, Higginbotham et al. 1986, PMID 3948345)
- **Leg C** — a-vO2diff = 4.553 mL/dL (first-principles Hb-Hill chemistry, **this doc**, never fed from
  any Fick-derived number)

`B × C` = 5.70 L/min × 45.53 mL/L = **259.51 mL/min**, vs Leg A's independently-measured 282.65 mL/min
→ **8.19% difference**. Pre-registered threshold: <20% (looser than this repo's own established 25%
band for an equivalent cross-check — i.e. not a loosened-to-pass threshold). **Gate: PASS — the loop
closes.** Also within 3.80% of the task's own classic textbook anchor (VO2≈250 mL/min).

**Honest, non-cherry-picked secondary result**: inverting the other direction (Q_chemistry =
twin's-own-VO2 / chemistry-a-vO2diff = 282.65/45.53 = **6.208 L/min**) sits *just outside*
`cardiac_output.py`'s own tight rest-Q anchor [4,6] L/min (0.21 L/min over) — reported as a marginal
miss, not smoothed over, though it clears that same doc's wider 25% band (ceiling 7.12). Both
thresholds are shown; nothing was cherry-picked to force a clean pass.

## 7. Second regime — Bohr shift, direction-forced across illustrative coefficient magnitudes

`log10(P50/P50_ref) = coef×(pH_ref−pH)`. No single coefficient magnitude was live-verifiable this
session (§3, Honest gaps) — so the direction claim is swept across coef ∈ {0.20, 0.30, 0.40, 0.48,
0.60, 0.70} at tissue PO2=30 mmHg, pH ∈ {7.40, 7.35, 7.30, 7.20, 7.10}. **For every single swept
coefficient, SaO2 strictly, monotonically decreases as pH falls** (e.g. at coef=0.48: 58.05% → 54.38%
→ 50.66% → 43.24% → 36.11%) — matching Collins (2015) + Wikipedia's live-verified right-shift
direction. **Gate: PASS at all 6 swept magnitudes** — the qualitative claim is robust to the
coefficient-magnitude uncertainty, a stronger test than pinning one possibly-fragile number.

## 8. Third regime — walking a-vO2diff widening, inverted to implied mixed-venous PO2

`cardiac_output.py`'s own walking sweep (a-vO2diff = 10/11/12 mL/100mL) is inverted (root-find,
arterial SaO2 held fixed) through this doc's Hb-Hill chemistry:

| a-vO2diff | implied PvO2 | implied SvO2 |
|---:|---:|---:|
| 10 mL/100mL | 25.89 mmHg | 48.16% |
| 11 mL/100mL | 24.04 mmHg | 43.22% |
| 12 mL/100mL | 22.28 mmHg | 38.27% |

Strictly monotonic (PvO2 falls further as a-vO2diff widens further) and physiologically plausible
(falls from the resting 40 mmHg but stays well positive/bounded — submaximal-walking territory, not
the ~20 mmHg near-maximal-exercise regime). **Gates: PASS** (monotonic + plausible). Disclosed
assumption: arterial SaO2 held fixed through submaximal walking (standard; arterial desaturation during
exercise is a near-maximal-effort/elite-athlete phenomenon) — not independently re-verified this
session. Roca et al. (1992) topically confirms venous PO2 falls at exercising muscle (the mechanism's
site) without a live-accessible number for this exact intensity.

## 9. Geometric structure — measured, not eyeballed

Closed-form location of the Hill curve's steepest slope: `PO2* = P50×((n−1)/(n+1))^(1/n)` = **19.94
mmHg** — matches the numerical argmax (19.94 mmHg, 0.014% diff) to machine precision. This is **not**
at P50 itself (a common imprecise heuristic) — a genuine, falsifiable, closed-form geometric fact.
Sensitivity ratio (slope at this point vs. slope at PO2=100 plateau): **40.8×** — the SaO2 observable
is ~41× more informative about a PO2 change in the tissue-extraction shoulder than in the arterial
plateau. Independently corroborated by the discrete PaO2 sweep (§10): shoulder slope (40-60 mmHg)
0.747 %/mmHg vs. plateau slope (top of sweep) 0.022 %/mmHg → **34.6× steeper**, measured two different
ways. **Gates: PASS** (analytical-numerical derivative match, max-slope-location match).

## 10. Non-degeneracy / void-floor sweeps

- **[Hb] swept 6→22 g/dL** (severe anemia → polycythemia), SaO2/SvO2 held fixed: CaO2 ranges
  [8.07, 28.85] mL/dL, strictly monotonic increasing; required CO for the twin's fixed VO2 ranges
  [14.72, 4.28] L/min, strictly monotonic **decreasing** with rising Hb — reproduces the real, known
  clinical phenomenon of compensatory high-output states in anemia. Analytical-vs-numerical derivative
  match: `d(CaO2)/d(Hb) = 1.34×SaO2 = 1.298` mL/dL per g/dL, confirmed to <1% (rtol). **Gates: PASS.**
- **PaO2 swept 30→150 mmHg**: SaO2 strictly monotonic [58.05%, 99.07%] — the plateau/shoulder shape
  from §9, now traced continuously rather than at 2 points. **Gate: PASS.**
- **Hb/Hct reference** (Billett 1990): task's Hb=15/Hct=45 sit inside the adult-male range; female
  ranges are systematically lower (12-16 g/dL, 36-48%) — **held explicitly OPEN**, not resolved.

## 11. Bonus/secondary (non-gating) — Severinghaus (1979) independent functional-form cross-check

| PO2 (mmHg) | Severinghaus (1979) | Hill (this doc) | diff (pct pts) |
|---:|---:|---:|---:|
| 40 | 74.95% | 75.05% | 0.11 |
| 60 | 90.58% | 89.99% | 0.59 |
| 95 | 97.39% | 96.88% | 0.50 |
| 100 | 97.75% | 97.28% | 0.47 |

Two structurally different functional forms (Hill sigmoid vs. Severinghaus's cubic-denominator
polynomial) agree to <0.6 percentage points across the physiological range. **Honest gap**: the
Severinghaus polynomial coefficients used here are **recalled**, not independently re-derived from a
live-fetched primary numeric table this session (PubMed full-text fetch for that specific paper was
CAPTCHA-blocked) — kept explicitly OFF `overall_pass`, reported as exploratory only.

## 12. couples_to — cardiac + pulmonary + metabolic (prose only; no shared-graph write this session)

- **MECHANISM_CARDIAC** (`cardiac_output.py`) — supplies the CO this doc's Fick-closure falsifier (§6)
  consumes; this doc supplies the a-vO2diff that script's own Fick calculation had to assume from
  literature.
- **MECHANISM_RESPIRATORY** (`respiratory.py`) — shares the same upstream VO2 (`metabolic_cost.py`);
  this doc adds the O2-**content** (not just ventilation-volume) side of gas exchange.
- **metabolic_cost.py / MECHANISM_METABOLIC_COST** — the ultimate VO2 source all three layers consume
  (the shared-instrument-chain caveat below applies identically here).

(Per `docs/MECHANISM_HARDENED_CONVENTIONS.md` §2/§4, promoting this coupling into a canonical graph edge
in `data/MECHANISM_ANCHOR_GRAPH.json` requires the fold pipeline — `mechanism_fold → fold_gate_v2` — a
separate step from this doc, not performed this session; stated here as prose/evidence-JSON metadata
only, consistent with `cardiac_output.py`/`respiratory.py`'s own docs, which also state their couplings
in prose without a graph-edge write.)

## 13. Pre-registered gates — 14/14 PASS (+ 1 disclosed non-gating bonus)

```
hill_equation_validated_2_independent_checkpoints:  PASS (0.07% diff @ PO2=40 vs Collins 2015's 75%)
forced_adversary_void_floor_ok:                     PASS (Hb-binding=98.6% of CaO2; wrong-(n,P50) fails; 13.0% grid-pass-rate)
cao2_in_task_anchor:                                PASS (19.76 in [18,22])
avo2diff_in_task_anchor:                             PASS (4.55 in [4,5])
svo2_in_task_anchor:                                PASS (75.05 in [70,80])
fick_loop_closes_vs_openism_vo2:                    PASS (8.19% diff, <20% pre-registered)
walking_regime_monotonic:                           PASS
walking_regime_plausible:                           PASS
bohr_direction_holds_all_swept_coefficients:        PASS (6/6 coefficients)
geometric_derivative_match:                         PASS
geometric_max_slope_location_match:                 PASS (0.014% diff)
hb_sweep_monotonic_nondegenerate:                   PASS
hb_sweep_derivative_match:                          PASS
pao2_sweep_monotonic_nondegenerate:                 PASS

BONUS (non-gating): severinghaus_routes_agree:      PASS (max 0.59 pct-pt diff)
```

**Overall: PASS** (deterministic — 2 independent runs produce byte-identical stdout; JSON-on-disk
machine-cross-checked against console-printed values, not eyeballed).

## 14. Honest gaps — symmetric QC: what this does NOT prove (disclosed, not hidden)

- **Nothing here is proven** in the strong sense. [Hb], Hct, and P50 all vary substantially
  inter-individually (Billett 1990: Hb 12-18 g/dL, Hct 36-54% across sex alone; P50 additionally
  shifts with 2,3-DPG, CO-Hb, fetal Hb, temperature — none modeled). A closed Fick loop with
  wide-plausible-band inputs is **consistency**, not proof, exactly as the task itself frames it.
- **n=2.7 is not independently re-derived from a primary numeric table this session** — it is
  task-specified and standard-secondary-sourced (Collins 2015's own P50≈26mmHg was live-confirmed;
  its explicit Hill-exponent figure was not, since the formula is image-embedded in the fetched text).
  The forced-adversary sweep (§5) shows this pairing is NOT arbitrary/interchangeable, but does not
  substitute for a primary-source re-derivation.
- **No single Bohr-effect coefficient was live-verified this session** (§3, Severinghaus 1966 —
  PubMed CAPTCHA-blocked). The direction claim (§7) is robust to this gap (swept across 6 magnitudes);
  the MAGNITUDE of the Bohr shift is not claimed at any specific precision.
- **Generic (population-level), not subject-specific**, exactly like `cardiac_output.py`/
  `respiratory.py`'s own disclosed scope: PaO2=95 mmHg, PvO2=40 mmHg, pH=7.40/7.36 are standard
  textbook operating points, not subject2's own measured blood gas (no ABG exists for subject2).
- **0-D, single-compartment, steady-state only** — no capillary transit-time/diffusion-limitation
  model, no 2,3-DPG/temperature/CO-Hb shifts, no CO2/Haldane-effect content model, no fetal-Hb or
  altitude coupling. The walking-regime inversion (§8) holds arterial SaO2 fixed — a disclosed,
  standard-but-unverified-this-session simplification.
- **SHARES THE METABOLIC-COST COMMON-MODE**, explicitly: Leg A's VO2 (§6) is the same single
  `metabolic_cost.py` OpenSim-probe output `cardiac_output.py` and `respiratory.py` already consume —
  an undiagnosed bug there propagates identically into this doc's Fick-closure number. Leg B
  (Higginbotham's real CO) and Leg C (this doc's Hb-Hill chemistry) are genuinely independent of that
  chain and of each other — the 3-leg triangulation's power comes from B and C, not from A being a
  fresh measurement.
- **Q_chemistry (§6) is a marginal, disclosed miss** against `cardiac_output.py`'s own tight rest-Q
  anchor (6.21 vs [4,6] L/min) — reported honestly, not smoothed into the wider band that does pass.
- **The Severinghaus (1979) bonus cross-check (§11) uses recalled, not live-re-verified, polynomial
  coefficients** — explicitly excluded from `overall_pass` for this reason.
- **Single subject (subject2), single trial (walking1)** — inherits, does not fix or expand, the
  upstream scripts' own scope limits.
- **No graph-edge write this session** (§12) — `couples_to` is prose/JSON metadata, not a folded
  `data/MECHANISM_ANCHOR_GRAPH.json` node/edge (would require the separate `mechanism_fold` pipeline).

## 15. Repro

```
cd ~/projects/bodytwin
source .venv-msk/bin/activate
python3 scripts/msk/blood_oxygen_transport.py
```
Requires `data/msk_smoketest/subject2_walking1/metabolic_cost/metabolic_cost_results.json` and
`data/msk_smoketest/subject2_walking1/cardiac_output/cardiac_output_results.json` to already exist
(run their producing scripts first if not); optionally reads
`data/msk_smoketest/subject2_walking1/respiratory/respiratory_results.json` (degrades gracefully if
absent — affects one printed cross-check line only, not any gate). Writes
`data/msk_smoketest/subject2_walking1/blood_oxygen_transport/blood_oxygen_transport_results.json`.
Pure Python/numpy/scipy (`scipy.optimize.brentq` for the walking-regime inversion, §8), no OpenSim
call, runs in under 2 seconds, deterministic (verified: 2 independent runs produce byte-identical
stdout). No git operations; reads all 3 upstream JSONs read-only; writes only under
`data/msk_smoketest/subject2_walking1/blood_oxygen_transport/`.
