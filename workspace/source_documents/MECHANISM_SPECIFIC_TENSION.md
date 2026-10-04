# MECHANISM Specific-Tension Decomposition — is the knee over-strength a PCSA error or a specific-tension error? (2026-07-22)

Closes the `Fmax = PCSA x specific_tension` decomposition that `docs/MECHANISM_FMAX_PCSA_VALIDATION.md`
left open. That doc found the twin's 12 knee-crossing prime movers systematically over-strong relative
to Handsfield-MRI-PCSA-implied Fmax (12/12 muscles above 1.0 at a fixed 60 N/cm^2 convention,
implied specific tension 60.5-79.3 N/cm^2) and causally confirmed the over-strength is real (correcting
it shrinks subject2's knee joint-force over-prediction 17.1%). Left open: is that over-strength a PCSA
(volume/architecture) error, or a specific-tension error? This doc backs out the model's own IMPLIED
specific tension per muscle and tests it against a THIRD, freshly-live-verified external anchor — the
primary human/mammalian specific-tension measurement literature — fully decorrelated from both PCSA
sources (Ward cadaver dissection, Handsfield MRI regression) the prior doc already compared.

**Headline (joint contribution, not a clean single-side localization):** the model's implied specific
tension for the 12 over-strong knee-crossing muscles (60.5-79.3 N/cm^2, subject2) sits **inside** the
task's own pre-registered wide falsifier band [15,100] N/cm^2 — so that test alone, taken naively,
would report "not falsified" and stop. Forcing the adversary (per this repo's watertight discipline: a
band this wide risks being a near-vacuous void-floor) turns up two sharper, decisive results the wide
band alone would miss: **(1)** the SAME [15,100] band genuinely IS violated — a real, non-noise-driven
falsification, not an artifact of tiny muscles — for 3-4 of the 12 knee-crossing muscles (gaslat_r,
gasmed_r, semimem_r, +vaslat_r for subject4) in the two smaller-bodied cross-subjects (subject3/4),
proving the band has real teeth and isn't vacuous; and **(2)** a same-muscle, same-aggregation-method
match against the two most directly relevant published studies (quadriceps femoris, pooled exactly the
way Erskine et al. 2009 and O'Brien et al. 2010 themselves pooled it) shows subject2's implied
quadriceps specific tension (70.4 N/cm^2) sits **+1.40 SD above the highest published point estimate**
for that exact muscle group (O'Brien 2010, 55+/-11 N/cm^2) — a study that is itself already 83% higher
than an equally-legitimate, equally-PMID-verified study of the identical muscles (Erskine 2009, 30+/-5
N/cm^2). **Verdict: the data are consistent with a genuine, partial specific-tension-side contribution
to the knee over-strength, layered on top of (not instead of) the PCSA-side story already established
— but this is non-falsification, not confirmation** (per this task's own symmetric-QC requirement):
the published literature's own >2x spread between equally-valid studies of the same muscles is wide
enough that "inside 1-2 SD of the most charitable published estimate" cannot be called impossible.

## 0. Non-circularity — the task's own posed puzzle, resolved explicitly

**The trap:** computing "PCSA_model" as `Fmax_model / assumed_specific_tension` and then comparing
that PCSA back to any specific-tension assumption is circular by construction — it recovers the
assumed value exactly, trivially, testing nothing.

**Why Ward's PCSA specifically cannot be used here:** the prior doc's Sec.5 established, live, that
`Fmax_model = PCSA_Ward x 61 N/cm^2` EXACTLY (Arnold, Ward, Lieber & Delp 2010, PMC2903973 — "Maximum
isometric force was calculated from measured PCSA and a specific tension of 61 N/cm^2 for all
muscles"). This script (`scripts/msk/specific_tension_validation.py`, Part 0) re-states this fact
before using it: since `PCSA_Ward (implied) = Fmax_model / 61` by pure algebra, using it as a
denominator here would trivially return `sigma_implied = 61.0` for every single muscle — a tautology,
not a measurement.

**The only non-circular path:** `sigma_implied = Fmax_model / PCSA_Handsfield`, where `PCSA_Handsfield`
comes from an INDEPENDENT in-vivo MRI-volume regression (Handsfield et al. 2014) that never entered
Fmax_model's construction at all. This is exactly the `st_implied_N_per_cm2` field the prior doc's
`fmax_pcsa_validation.py` already computed — Part 1 of this analysis re-derives it independently from
that JSON's own raw per-muscle fields (`fmax_model_N`, `pcsa_cm2`) as a machine cross-check, rather than
trusting the stored field or prior prose.

## 1. Verification gates (machine-checked, all PASS)

Run: `scripts/msk/specific_tension_validation.py` (real run, exit 0).

| gate | result |
|---|---|
| Recompute `sigma_implied = Fmax_model_N / PCSA_Handsfield_cm2` independently for all 160 (4 subjects x 40 muscles) rows vs the prior doc's stored `st_implied_N_per_cm2` field | **PASS**, max\|diff\| = 0.00e+00 (exact) |
| Literature PMIDs (8 papers) — found via NCBI eutils `esearch`, confirmed via `esummary` title/journal/author match, values pulled via `efetch` abstract text | **PASS**, all 8 confirmed live this session (not model memory) — full log `data/external/specific_tension_literature/PROVENANCE.md` |
| Maganaris 2001 kN/m^2 -> N/cm^2 unit conversion (150/155/104/658) | **machine-verified**: `1 kN/m^2 = 0.1 N/cm^2` exactly -> 15.0/15.5/10.4/65.8 N/cm^2 |
| Fukunaga 1996 "11-47 N/cm^2" and Powell 1984 "22.5 N/cm^2" cross-triangulation | **PASS** — independently re-derived via NCBI eutils author-search AND found as Arnold 2010's own references #16/#35 (re-fetched raw text live this session, `pmc.ncbi.nlm.nih.gov/articles/PMC2903973/`) — two independent retrieval paths converge on the same two PMIDs (8847297, 6511546) |
| Knee-crossing (n=12) / hip-crossing (n=25) muscle-set membership | reused, unedited, from `static_opt_knee.py`/`validate_hip_force.py` — same sets as the prior doc |

## 2. The live-verified specific-tension literature (decorrelated from both PCSA sources already used)

All PMIDs confirmed live this session (NCBI eutils `esearch`/`esummary`/`efetch`, raw `curl`, not a
summarizer guess). Full table: `data/external/specific_tension_literature/specific_tension_lit_values.json`;
retrieval log: `data/external/specific_tension_literature/PROVENANCE.md`.

| paper (PMID) | population / method | value (N/cm^2) |
|---|---|---|
| Zajac 1989 (2676342) | foundational modeling review | no number in abstract; paywalled — cited for the modeling framework only |
| Buchanan 1995 (7489126) | human, elbow flexors vs extensors, dynamometry | no number in abstract; concludes specific tension "is not a constant," flexors >> extensors |
| Powell et al. 1984 (6511546) | guinea pig hindlimb, cadaveric/architectural | **22.5** (point; = Arnold 2010's own ref #35, re-confirmed live) |
| Fukunaga et al. 1996 (8847297) | human (n=8 M) in-vivo MRI+torque, ankle plantar/dorsiflexors | **11-47** (range; = Arnold 2010's own ref #16, re-confirmed live; dorsiflexors "more than twofold higher" than plantar flexors per the abstract itself) |
| Maganaris et al. 2001 (11181594) | human (n=6 M) in-vivo, soleus+TA, electrical-stim ("realistic") vs cadaveric-param ("unrealistic", authors' own label) | **realistic: 15.0 (Sol) / 15.5 (TA)**; unrealistic: 10.4 (Sol) / 65.8 (TA) |
| Erskine et al. 2009 (19468746) | human (n=27 young M) in-vivo, whole quadriceps-femoris GROUP, voluntary MVC + activation-corrected | **30 +/- 5** |
| O'Brien et al. 2010 (19748968) | human (n=20 adults) in-vivo, same QF-group method as Erskine 2009 | **55 +/- 11 (men)**, 57.3+/-13 (women) |
| Erskine et al. 2010 (20102343) | human (n=17 M), same QF group, PRE/POST 9-wk resistance training | **+20% within-subject** (PCSA rose only +6%) |

**The wide spread IS the point (this task's own framing, now with live numbers, not just prose):**
Erskine 2009 and O'Brien 2010 measured the SAME muscle group (quadriceps femoris) with the SAME
pooling formula (tendon force / sum(PCSA_i x cos(pennation_i))), from broadly the same UK
in-vivo-muscle-physiology lab lineage (Maganaris co-authors both), and got specific-tension point
estimates **83% apart** (30 vs 55 N/cm^2, men-matched) — before any model, any PCSA-source disagreement, or any
subject anthropometry enters the picture at all. Buchanan 1995's own conclusion, independently
verified live rather than assumed, is that a single constant specific tension "in muscle models is
questionable." Maganaris 2001 demonstrated a single muscle's (TA) implied specific tension can swing
from 15.5 to 65.8 N/cm^2 (4.2x) depending on a plausible-looking but wrong methodological choice —
by the SAME six subjects, same muscle, same underlying moments.

## 3. Model's implied specific tension — per-muscle and same-method group-pooled

**Per-muscle-set (subject2, primary; distributions machine-recomputed, Part 3 of the pipeline):**

| muscle set | n | range (N/cm^2) | mean | median | n exceeding 100 |
|---|---:|---:|---:|---:|---:|
| all 40 | 40 | 38.9 - 86.7 | 64.3 | 64.0 | 0 |
| knee-crossing (the established over-strong set) | 12 | 60.5 - 79.3 | 68.8 | 67.5 | 0 |
| hip-crossing | 25 | 38.9 - 86.7 | 62.4 | 62.2 | 0 |

**Group-pooled (sum(Fmax_i)/sum(PCSA_i), matching Erskine/O'Brien's own aggregation exactly):**

| subject | quadriceps femoris (vasmed+vaslat+vasint+recfem) | triceps surae (soleus+gasmed+gaslat) |
|---|---:|---:|
| subject2 (primary) | **70.44** N/cm^2 | **69.88** N/cm^2 |
| Ward/Arnold cadaver-mean (control) | 76.34 | 75.45 |
| subject3 (cross-subject) | 94.40 | 92.24 |
| subject4 (cross-subject) | 95.93 | 93.64 |

**Matched-muscle comparison (the forced, sharper adversary — the SD-multiple above the single
highest-published, same-muscle-group, same-method point estimate found):**

| subject | quadriceps model ST | SD above Erskine 2009 (30+/-5) | SD above O'Brien 2010 men (55+/-11) | triceps-surae model ST | x Maganaris 2001 "realistic" soleus (15.0) |
|---|---:|---:|---:|---:|---:|
| subject2 (primary) | 70.4 | +8.09 SD | **+1.40 SD** | 69.9 | 4.66x |
| cadaver-control | 76.3 | +9.27 SD | +1.94 SD | 75.5 | 5.03x |
| subject3 | 94.4 | +12.88 SD | +3.58 SD | 92.2 | 6.15x |
| subject4 | 95.9 | +13.19 SD | +3.72 SD | 93.6 | 6.24x |

## 4. Cross-subject sweep — is the wide [15,100] band actually vacuous? (forced, not assumed)

The adversary this analysis is obligated to force before trusting the wide-band "not falsified" result
for subject2 alone: a band that never rejects anything isn't a test. Sweeping the SAME [15,100]
falsifier across the diverse instance-space already built by the prior doc (2 smaller-bodied
cross-subjects, 1 cadaver-mean control) shows the band DOES get violated — genuinely, not by tiny
noise-prone muscles:

| subject | knee-crossing-12 range (N/cm^2) | muscles exceeding 100 |
|---|---:|---|
| subject2 (primary) | 60.5 - 79.3 | **none** |
| Ward/Arnold cadaver-control | 66.0 - 85.7 | **none** |
| subject3 | 82.6 - 104.9 | **gaslat_r (103.6), gasmed_r (100.3), semimem_r (104.9)** |
| subject4 | 83.9 - 106.5 | **gaslat_r (105.4), gasmed_r (102.0), semimem_r (106.5), vaslat_r (100.9)** |

**These are not small-muscle regression-noise artifacts** — PCSA for the flagged knee muscles is
14.9-51.0 cm^2 (large, well-conditioned regression estimates), unlike the separately-flagged hip
outlier `glmin1_r` (PCSA 3.3-4.3 cm^2, a genuinely tiny muscle where proportional Handsfield-regression
error is largest — flagged and excluded from this "real signal" claim, not conflated with it). So the
band is confirmed non-vacuous: for smaller-bodied subjects, specific tension alone genuinely CANNOT
explain part of their over-strength (a clean, forced PCSA-side requirement for those muscles/subjects)
— while for the actual primary subject (subject2), it remains inside bounds.

**Geometric governor (same one already identified in the prior doc, re-confirmed here):** since
`sigma_implied` is just `ratio_model_over_pcsa x 60` and the ratio was already shown to be a
monotonically decreasing function of a subject's height x mass, `sigma_implied` inherits the exact
same monotonic relationship — smaller subjects mechanically push every muscle's implied specific
tension higher, independent of any new mechanism. This is why subject3/4 (H×M=107.3/105.2) cross the
100 N/cm^2 ceiling for several muscles while subject2 (H×M=153.3, closest to the cadaver template) does
not — the SAME governing relationship, not a separate phenomenon.

## 5. Symmetric QC — why "inside the wide band" is non-falsification, not confirmation

**Adversary 1 (void-floor risk in the wide band):** forced in Sec.4 — the band is not vacuous, but it
IS lenient for the primary subject specifically. Resolved by the sharper matched-muscle test (Sec.3),
which is harder to pass and still leaves subject2 above the highest published point estimate.

**Adversary 2 (measurement-family mismatch, the honest limiting factor on the matched-muscle test):**
Erskine/O'Brien's PCSA denominator is `muscle volume / fascicle length AT THE MVC JOINT ANGLE`
(dynamic, contraction-specific, from ultrasound); this repo's `PCSA_Handsfield` denominator is
`Handsfield-regression volume / the model's OWN optimal_fiber_length` (a fixed anatomical model
constant, not tied to any measured contraction state). These are not perfectly like-for-like PCSA
definitions — part of the quadriceps-group gap could reflect this definitional mismatch rather than a
"true" specific-tension difference, i.e., yet another PCSA-SIDE confound at a different level, not
necessarily proof of a specific-tension-side error. Disclosed, not resolved, exactly as the prior doc
disclosed its own analogous "shared methodological family" caveat for the Ward-vs-Handsfield
comparison.

**Adversary 3 (which comparison is more trustworthy, quadriceps vs triceps-surae — reported both, not
cherry-picked):** the quadriceps comparison has the better ANATOMICAL match (all 4 heads, exact same
pooling formula as this repo's own QUAD_GROUP) but relies on voluntary MVC (activation-corrected, but
still voluntary). The triceps-surae comparison has a weaker anatomical match (soleus-only vs this
repo's 3-muscle pooled group) but Maganaris's "realistic" value uses ELECTRICAL STIMULATION — the same
conceptual quantity as OpenSim's `max_isometric_force` (fully-activated, a=1) — arguably the cleaner
conceptual match despite the anatomical mismatch. Both point the same direction (model implied ST
clearly above the "realistic"/careful literature values); neither alone is treated as sufficient.

**Symmetric conclusion, precisely scoped:** "inside the wide [15,100] band" for the primary subject is
**non-falsification**, not confirmation, exactly as this task's own instructions required — the
measured literature's own >2x inter-study spread on the identical muscle group means "1.4 SD above the
highest of two legitimate studies" cannot be called definitively wrong. The honest reading is a
**partial, plausible, non-dispositive specific-tension-side contribution**, additive to (not a
replacement for) the PCSA-side explanation the prior doc already causally confirmed accounts for about
half of subject2's knee over-prediction gap.

## 6. Pre-registered falsifiers — resolved explicitly

- **"Model's implied specific tension for the over-strong knee-crossing muscles exceeds the measured
  15-100 N/cm^2 spread (localizing the over-strength cleanly to specific tension)"**: **NOT FALSIFIED
  for the primary subject** (max 79.3 N/cm^2, comfortably inside) — **FALSIFIED (band violated) for 3-4
  of 12 knee-crossing muscles in the two smaller-bodied cross-subjects** (subject3/4, gaslat_r/gasmed_r/
  semimem_r/+vaslat_r all >100 N/cm^2) — a genuine, muscle-set-specific, anthropometry-governed
  falsification, not a uniform yes/no across the whole instance space.
- **"Model's implied specific tension matches the measured literature (strength calibration via
  specific tension is NOT the source)"**: **NOT confirmed** — the sharper same-muscle, same-method
  quadriceps comparison shows subject2 sitting +1.40 SD above the single highest published point
  estimate (O'Brien 2010) and +8.09 SD above another equally legitimate one (Erskine 2009) for the
  identical muscle group; the triceps-surae comparison shows a 4.66x multiple over Maganaris's
  "realistic" soleus estimate. Both point toward a real, non-trivial specific-tension-side
  contribution, though the literature's own spread prevents calling this a clean proof.
- **Net localization of the knee over-strength, PCSA vs specific tension:** **BOTH, jointly** — not a
  single-side answer. This doc does not overturn or replace the prior doc's PCSA-side causal finding
  (correcting Fmax via Handsfield-PCSA at a fixed 60 N/cm^2 already closes ~50% of subject2's knee
  gap); it adds evidence that a portion of what that correction captures is itself plausibly a
  specific-tension-side effect (the "60 N/cm^2 convention" used throughout this model lineage sits
  1.4-1.9 SD above the highest careful in-vivo human estimates of the exact same muscles), meaning the
  two error sources are not cleanly separable with the data available — an honest, non-overclaiming
  answer to the task's own posed question.

## 7. Confidence tier and honest gaps

**Confidence tier: cadaveric/published-plausibility**, as specified by the task, and for the same
reasons the prior doc used that tier plus new ones specific to this literature layer:
1. **Zajac 1989 and Buchanan 1995 could not be numerically verified this session** — both are
   paywalled, no OA copy found (Unpaywall live-checked for Buchanan), so both are cited only for their
   qualitative conclusions (verified live via PubMed abstract), not a number. This is disclosed, not
   silently patched with a remembered/assumed figure.
2. **Fukunaga 1996's 11-47 N/cm^2 range is second-hand** (via Arnold 2010's own characterization,
   independently re-confirmed as that paper's actual reference #16 this session) — not this session's
   own read of Fukunaga's full-text tables.
3. **PCSA-definition mismatch between the literature anchor and this repo's Handsfield-based PCSA**
   (dynamic fascicle-length-at-MVC vs a fixed model constant) is a genuine, undecided confound on the
   matched-muscle comparison's exact magnitude (Sec.5, Adversary 2) — the DIRECTION of the finding
   (model implied ST above careful in-vivo estimates) is likely robust to it, but the precise SD
   multiples are not immune to this definitional difference.
4. **n=1 formal primary-subject conclusion + n=2 cross-subject robustness sweep + n=1 cadaver control**
   — the same scope this whole doc family has used throughout; sufficient to demonstrate the
   mechanism and its geometric governor, insufficient for a population-level magnitude claim.
5. **No independent Ward-cadaver PCSA table exists in this repo** to attempt a three-way (Ward vs
   Handsfield vs in-vivo-literature) PCSA cross-check directly — only the algebraic fact that Ward's
   implied PCSA is definitionally `Fmax_model/61` was available, which is precisely why it could not be
   used as a fresh anchor here (Sec.0).
6. **Only two anatomical groups (quadriceps, triceps surae) were matched to primary literature** — the
   other 8 of the 12 knee-crossing muscles (hamstrings, tibialis posterior, etc.) were not individually
   matched to a same-method published specific-tension study this session; the all-40/knee-12
   per-muscle-set distributions (Sec.3) stand on the wide-band test only for those muscles.

## 8. Files

- `scripts/msk/specific_tension_validation.py` — Part 0-4 pipeline (non-circularity re-statement,
  recompute cross-check, group-pooled aggregation, per-muscle-set + cross-subject falsifier sweep,
  matched-muscle SD-multiple comparison, verdict). Real run, exit 0. Reads the prior doc's own JSON;
  builds no new OpenSim model.
- `data/msk_smoketest/subject2_walking1/specific_tension_validation/specific_tension_validation_results.json`
  — every number in Sec.1-4, full per-subject/per-muscle-set/per-group detail.
- `data/external/specific_tension_literature/specific_tension_lit_values.json` — the 8-paper
  literature table (PMID, citation, verification tier, value(s)).
- `data/external/specific_tension_literature/PROVENANCE.md` — full live-retrieval log (exact curl
  commands, what was searched but not used, honest gaps).
- Reused, unedited: `scripts/msk/static_opt_knee.py` (`KNEE_ANATOMICAL_MUSCLES_EXPECTED`),
  `scripts/msk/validate_hip_force.py` (`HIP_ANATOMICAL_MUSCLES_EXPECTED`).
- Input (not modified): `data/msk_smoketest/subject2_walking1/fmax_pcsa_validation/fmax_pcsa_validation_results.json`
  (from `docs/MECHANISM_FMAX_PCSA_VALIDATION.md`).
- Prior doc this closes the loop on: `docs/MECHANISM_FMAX_PCSA_VALIDATION.md`.

Isolation respected throughout: bodytwin only; no new OpenSim model built (pure re-analysis of already
machine-verified numbers plus a fresh literature anchor); external web fetches were read-only (NCBI
eutils, PMC, Unpaywall) with no credentials and no write access. No git commit, no git push, no git add
performed.
