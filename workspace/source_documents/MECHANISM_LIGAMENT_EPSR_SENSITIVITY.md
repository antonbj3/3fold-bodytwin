# MECHANISM LIGAMENT EPSILON_R SENSITIVITY — fragility sweep of the Blankevoort1991Ligament reference strain (2026-07-21)

**Task.** Every certified ligament-engagement conclusion in `docs/MECHANISM_KNEE_LIGAMENTS.md`,
`docs/MECHANISM_HIP_ANKLE_LIGAMENTS.md`, `docs/MECHANISM_SPINE_LIGAMENTS.md`, and
`docs/MECHANISM_UNIFIED_V2.md` rests on a per-bundle `Blankevoort1991Ligament` **reference strain**
(epsilon_r) — a modeling choice, not a subject-specific measurement. This sweeps epsilon_r across a
published-plausible band and re-evaluates every headline conclusion, producing a machine-checked
survives/flips ledger. Script: `scripts/msk/epsr_sensitivity_sweep.py`. Model:
`data/msk_models/subject2_unified_v2.osim` (read-only input, never modified).
Run: `source_repository/.venv-msk/bin/python3 scripts/msk/epsr_sensitivity_sweep.py`.

**Isolation respected:** read-only against `subject2_unified_v2.osim` (every perturbed instance is a
fresh in-memory `osim.Model(...)` copy; nothing is ever `printToXML`'d back over the source); no git
add/commit/push; this doc, the harness script, and the two evidence JSONs are the only new files.

---

## 0. Headline

**Of 18 tracked conclusions: 10 are ROBUST across the full doubled exploratory range
([-0.10,+0.10]); 7 flip via a single clean sign-change; 1 (the MCL magnitude threshold) shows a
narrow, non-monotonic dip.** Of the 8 that move at all, **3 flip strictly INSIDE the task-cited
"typical" +/-0.05 band — the falsifier is triggered for exactly these 3**: knee ACL
tighter-in-extension (flips at epsilon_r ≈ baseline − 0.030), spine ISL tighter-in-flexion (− 0.038),
and — most fragile of all — spine **PLL tighter-in-flexion, which flips at only ≈ baseline − 0.012,
less than a quarter of the cited band**, independently corroborating (via a completely different
sweep/pose/model than the original build) a fragility the spine-ligament doc's own build-time OODA
process had already stumbled into and only narrowly avoided (`MECHANISM_SPINE_LIGAMENTS.md` §5). A
further 5 conclusions are only conditionally fragile — robust within the cited band but flip (or, for
the MCL case, dip) somewhere in the wider +/-0.05 to +/-0.10 exploratory margin: ATFL
tighter-in-plantarflexion (−0.063), spine ALL (−0.065), spine LF (−0.068), spine SSL (−0.071), and the
knee's **MCL "58% peak strain" headline number, which is not robust as a magnitude at all** (it ranges
49.8%-73.7% across the swept band, with an isolated dip just below 50% starting at |delta|≈0.055; the
region-depth check in §7 shows this fragility compounds with flexion depth, sharpening rather than
contradicting that doc's own pre-existing "deep-flexion magnitudes are qualitative, not validated"
caveat). By contrast, the **10 robust conclusions include every hip/ankle engagement-direction claim
the task named as headline findings** (ILFL tighter-in-extension, PTFL tighter-in-dorsiflexion, the
CFL/DELT antagonist check, PCL's reciprocal pattern) — geometrically forced, surviving the full
explored range under both adversaries. **The "0/1397 tension-only violations" count** (the task's own
name for this conclusion, from `MECHANISM_UNIFIED_V2.md`'s hip+ankle-only re-verification; this harness's
own broader aggregate spans knee+hip+ankle+spine and measures 0/4138 at every swept delta, same
qualitative result at a wider scope, not a literal re-measurement of "1397") **is also in that
robust-10, but for a reason that disqualifies it as a fragility test at all**: it is a tautological
identity of the force law's own branching (`strain<=0 -> force==0`, true for any valid positive
`slack_length`), not an emergent geometric or biomechanical fact — flagged explicitly per the task's
own "never a tautology gate" rule (§6).

---

## 1. What epsilon_r actually is here (verified live, not assumed)

`Blankevoort1991Ligament` has **no persistent `reference_strain` property** — live introspection
already recorded in `MECHANISM_KNEE_LIGAMENTS.md` §1 lists only `linear_stiffness`, `transition_strain`,
`damping_coefficient`, `slack_length` (+ `GeometryPath`). `reference_strain` is an **argument** to
`setSlackLengthFromReferenceStrain(reference_strain, state)`, called once per bundle at build time at
that fork's own reference pose (all coordinates = 0), which solves
`slack_length = L_ref / (1 + reference_strain)` where `L_ref` is the bundle's real geometric path
length at that pose. Only the resulting `slack_length` is persisted in the `.osim` file.

Sweeping epsilon_r therefore means recovering `L_ref` from each bundle's **already-published** baseline
`slack_length` and reference strain (both logged in the three fork config JSONs and present, unchanged,
in the live merged model — the graft that built `subject2_unified_v2.osim` is a pure ForceSet copy, so
`slack_length` there is bit-identical to each standalone build), then recomputing a new `slack_length`
for a shifted reference strain in closed form:

```
L_ref             = slack_length_baseline * (1 + r0)                    [recovered, per bundle]
slack_length_new  = L_ref / (1 + r0 + delta)
                  = slack_length_baseline * (1 + r0) / (1 + r0 + delta)  [applied per bundle]
```

This is **geometric**, not a heuristic substitute: it holds each bundle's real measured attachment-point
geometry (hence `L_ref`) fixed and only reinterprets how much of that fixed length is "slack" vs.
"pre-strain" — exactly isolating the one modeling choice under test, with no need to walk any model back
to its original per-fork reference pose. `r0` (baseline) differs per bundle and is read from
`scripts/msk/knee_ligament_config.json` (`reference_strain_at_laiarnold_full_extension`, 84 bundles,
range −0.14 to +0.06), `scripts/msk/hip_ankle_ligament_config.json`
(`reference_strain_at_reference_pose`, 22 bundles, uniform +0.02), and
`scripts/msk/spine_ligament_config.json` (`reference_strain_at_neutral`, 5 bundles, +0.02/+0.015/+0.06/
−0.01/−0.01) — 111 bundles total, matching `subject2_unified_v2.osim`'s own ligament count exactly
(verified live, `n_adjusted == 111` asserted in the harness).

### Citations for the swept band (live-verified this session, not recalled)

| citation | verified PMID / DOI | role |
|---|---|---|
| Blankevoort L, Huiskes R. Ligament-bone interaction in a three-dimensional model of the knee. *J Biomech Eng.* 1991;113(3):263-9. | PMID [1921352](https://pubmed.ncbi.nlm.nih.gov/1921352/), doi:10.1115/1.2894883 | source of the force law itself AND (via the vendored Lenhart2015 knee model already cited in `MECHANISM_KNEE_LIGAMENTS.md`) this repo's own per-bundle baseline reference-strain values — the "per-bundle documented ranges" the task names. |
| Baldwin MA, Laz PJ, Stowe JQ, Rullkoetter PJ. Efficient probabilistic representation of tibiofemoral soft tissue constraint. *Comput Methods Biomech Biomed Engin.* 2009;12(6):651-9. | PMID [19370459](https://pubmed.ncbi.nlm.nih.gov/19370459/), doi:10.1080/10255840902822550 | directly on-point: a probabilistic (Advanced-Mean-Value vs. Monte Carlo) sensitivity study whose own abstract names "inherent uncertainty in ligament stiffness, **reference strain**, and attachment site locations" as the studied object — confirms epsilon_r is a recognized, field-studied uncertainty axis, not an invented adversary. |
| Baldwin MA, Clary CW, Fitzpatrick CK, Deacy JS, Maletsky LP, Rullkoetter PJ. Dynamic finite element knee simulation for evaluation of knee replacement mechanics. *J Biomech.* 2012;45(3):474-83. | PMID [22209313](https://pubmed.ncbi.nlm.nih.gov/22209313/), doi:10.1016/j.jbiomech.2011.11.052 | companion evidence that soft-tissue ligament properties in this model family are routinely **calibrated per-specimen from laxity tests**, not assumed a priori — corroborates the "not tightly measured" framing. |

**Honesty flag, same class already disclosed in `MECHANISM_SPINE_LIGAMENTS.md` §6 for Pintar 1992**:
Baldwin 2009's own numeric SD/range table could **not** be extracted live this session — Europe PMC
confirms PMID 19370459 is not open access (abstract-only); the DOI redirect surfaces no visible content.
The swept magnitudes below (+/-0.05, +/-0.10) are **not** lifted from that table. +/-0.05 is the
task-specified "typical" band (and matches this repo's own already-published per-bundle spread, e.g.
knee ACLam −0.14 to −0.05); +/-0.10 is swept **additionally**, beyond the cited band, specifically to
satisfy the symmetric-QC requirement ("check you used the FULL published range, not a convenient
sub-band") — every flip point below is reported by its exact location, so a reader can see whether it
falls inside or outside the narrower cited band regardless of which band is treated as authoritative.
**Confidence tier: the sweep itself is method-only; the epsilon_r magnitude band is
cadaveric/published-plausibility (topic verified live; exact SD not extractable, honestly flagged).**

---

## 2. Two adversaries (forced, not just the convenient one)

1. **Common-mode**: every one of the 111 bundles' reference strain shifted by the **same** delta,
   21 grid values spanning [-0.10, +0.10] (fine step 0.005-0.01 near zero, coarser toward the edges).
   This is the "obvious" sweep, and — per this repo's own operating lesson on self-flattering
   common-mode-only robustness checks — is inherently friendly to directional claims (a rigid shift of
   a monotonic engagement curve rarely reorders which side is bigger).
2. **Decorrelated**: each bundle's reference strain independently jittered, IID uniform, at two
   magnitudes (+/-0.05, +/-0.10), 20 seeds each (`numpy.random.default_rng`, fixed seeds for
   reproducibility). This is the **forced adversary** against any conclusion that looks robust only
   because a rigid shift preserves inter-bundle ordering — independent per-bundle noise does not.

Both were run; results agree with each other everywhere checked (§4), which is itself a cross-check
that the common-mode sweep is not silently self-flattering here.

**Explicitly excluded from the sweep, and why (checked, not silently skipped):** the knee's
valgus/varus MCL/LCL opposite-sign-slope check, the ATFL/PTFL structural-invariance-to-subtalar-
rotation check, and the structural length-plausibility check are pure **path-geometry** facts
(attachment-point positions / lengths) that never reference `slack_length` or strain — they are
epsilon_r-**invariant by construction**, not empirically robust. Re-sweeping them would spend compute
on a foregone conclusion, so they are marked N/A rather than padding the PASS column.

---

## 3. External anchor — baseline (delta=0) reproduction, checked BEFORE trusting any delta!=0 result

The closed-form recomputation at delta=0 is a pure numerical no-op and **must** exactly reproduce the
already-published numbers from the 4 certified docs (produced by different code — the original
`add_*_ligaments.py` build scripts). This is the external, non-tautological anchor: a genuine
cross-implementation check, not circular.

| metric | published | measured (delta=0) | abs err | tol | verdict |
|---|---:|---:|---:|---:|---|
| knee ACL ext band (N) | 4.318 | 4.318 | 0.0004 | 0.05 | PASS |
| knee ACL flex band (N) | 2.143 | 2.143 | 0.0001 | 0.05 | PASS |
| knee PCL ext band (N) | 6.180 | 6.180 | 0.0004 | 0.08 | PASS |
| knee PCL flex band (N) | 355.735 | 355.735 | 0.0003 | 4.0 | PASS |
| knee MCL peak strain @140deg (%) | 58.2 | 58.207 | 0.007 | 0.6 | PASS |
| hip ILFL ext band (N) | 300.324 | 300.324 | 0.0002 | 3.5 | PASS |
| hip ILFL flex band (N) | 17.055 | 17.055 | 0.0003 | 0.3 | PASS |
| ankle ATFL plantar band (N) | 10.468 | 10.468 | 0.0002 | 0.15 | PASS |
| ankle ATFL dorsi band (N) | 0.000 | 0.000 | 0.0000 | 0.05 | PASS |
| ankle PTFL dorsi band (N) | 174.85 | 174.846 | 0.004 | 2.0 | PASS |
| ankle PTFL plantar band (N) | 88.05 | 88.049 | 0.001 | 1.2 | PASS |
| spine ALL ext band (N) | 11.718 | 11.718 | 0.0000 | 0.2 | PASS |
| spine PLL flex band (N) | 0.108 | 0.108 | 0.0003 | 0.02 | PASS |

**15/15 PASS.** Only after this gate passed did the harness proceed to report fragility (per
`epsr_sensitivity_sweep.py`'s `main()`: `if not baseline_ok: raise SystemExit(...)`).

**One real forensic finding surfaced by this gate, fixed at the source (kills/fixes are auditable, not
smoothed over):** a first version of the harness passed a literal `-35.0`/`+35.0` degree grid for the
ankle subtalar sweep, while the original `add_hip_ankle_ligaments.py` derives its grid from
`coord.getRangeMin()/getRangeMax()` (radians -> degrees), which — because the model's own stored
coordinate range in radians does not round-trip to an exact integer degree value — actually sweeps
`[-34.99999999999998, +35.00000000000002]`. That extra `2e-14` deg pushes the true edge sample just
*above* 35.0, so the original script's own `<=35.0` band filter silently **excludes** it. My literal
`-35.0/+35.0` grid does not have that overshoot and so *included* the true edge sample, giving
CFL-inversion-band-mean = 77.09 N instead of the published 74.324 N (CFL rises steeply near that edge:
91.06 N at 33deg to 96.46 N at ~35deg, so a single-sample inclusion/exclusion swings the 7-point mean by
~3.7%). This is a pure floating-point sweep-boundary artifact, **not** an epsilon_r effect (found at
delta=0, before any perturbation) — fixed by deriving the grid live from the coordinate's own range in
every sweep, exactly matching the original convention (`epsr_sensitivity_sweep.py`'s `sweep_1d`,
docstring on the fix). After the fix, CFL-inversion reproduces 74.32388... N, matching to 1e-4 N.

---

## 4. Pre-registration (set before running the sweep)

| expected ROBUST (directional / topological "tighter in X") | expected FRAGILE (absolute magnitude / near-zero count) |
|---|---|
| ACL tighter-extension, PCL tighter-flexion, ILFL tighter-extension, ISFL internal-both-levels, PFL external-in-extension (both metrics), ATFL tighter-plantarflexion, PTFL tighter-dorsiflexion, CFL tighter-dorsiflexion, CFL/DELT antagonist, spine ALL/PLL/LF/ISL/SSL engagement direction | MCL "58%" peak-strain magnitude, "0/1397" tension-only-violation count |

---

## 5. Results — survives / flips ledger

Legend: **ROBUST** = 0 sign changes across all 21 common-mode deltas in [-0.10,+0.10]. **FLIPS@d** =
single monotonic flip at common-mode delta `d` (bisected to +/-0.0005). Decorrelated columns = fraction
of 20 IID-jittered seeds that still PASS at that magnitude.

| conclusion | pre-registered | baseline | common-mode result | decorr. pass ±0.05 | decorr. pass ±0.10 | in cited ±0.05 band? |
|---|---|:-:|---|:-:|:-:|:-:|
| ACL tighter extension | robust | PASS | **FLIPS @ −0.0298** | 65% | 35% | **YES — fragile inside cited band** |
| PCL tighter flexion | robust | PASS | ROBUST | 100% | 100% | robust |
| ILFL tighter extension | robust | PASS | ROBUST | 100% | 100% | robust |
| ISFL internal both levels | robust | PASS | ROBUST | 100% | 100% | robust |
| PFL external in extension (raw gap>0) | robust | PASS | ROBUST | 100% | 100% | robust |
| PFL "amplified in ext.", raw-N gap | robust | **FAIL at baseline** | ROBUSTLY FAILS | 0% | 0% | robust (robustly wrong, not epsilon_r-fragile — already metric-fragile at baseline, `MECHANISM_HIP_ANKLE_LIGAMENTS.md` §5) |
| PFL "amplified in ext.", normalized gap | robust | PASS | ROBUST | 100% | 100% | robust |
| ATFL tighter plantarflexion | robust | PASS | **FLIPS @ −0.0633** | 100% | 85% | just outside — robust *within* cited band, not unconditionally |
| PTFL tighter dorsiflexion | robust | PASS | ROBUST | 100% | 100% | robust |
| CFL tighter dorsiflexion (informational) | robust | PASS | ROBUST | 100% | 100% | robust |
| CFL/DELT antagonist cross-check | robust | PASS | ROBUST | 100% | 100% | robust |
| spine ALL tighter extension | robust | PASS | **FLIPS @ −0.0645** | 100% | 95% | just outside cited band |
| **spine PLL tighter flexion** | robust | PASS | **FLIPS @ −0.0120** | **45%** | 50% | **YES — MOST fragile, <1/4 of cited band** |
| spine LF tighter flexion | robust | PASS | **FLIPS @ −0.0683** | 100% | 85% | just outside cited band |
| spine ISL tighter flexion | robust | PASS | **FLIPS @ −0.0380** | 95% | 75% | **YES — fragile inside cited band** |
| spine SSL tighter flexion | robust | PASS | **FLIPS @ −0.0714** | 100% | 85% | just outside cited band |
| MCL peak strain > 50% ("non-physiological") | fragile | PASS | **NON-MONOTONIC** — isolated dip below 50% only for delta in (−0.065,−0.05), min 49.67% | 100% | 95% | no — closest approach |delta|≈0.055, just outside the cited band (but only barely, and the underlying magnitude is fragile everywhere, §7) |
| 0/1397 tension-only violations | fragile | PASS | ROBUST (but TAUTOLOGICAL — §6) | 100% | 100% | robust, but not a real test (§6) |

**Falsifier verdict, exact count (18 conclusions total): 3 flip strictly INSIDE the task-cited +/-0.05
band** (knee ACL tighter-extension @ −0.030, spine ISL @ −0.038, spine PLL @ −0.012) — **the falsifier
is triggered for these 3**. **5 more are only conditionally fragile** — robust within the cited band
but flip/dip somewhere between +/-0.05 and +/-0.10 (ankle ATFL @ −0.063, spine ALL @ −0.065, spine LF
@ −0.068, spine SSL @ −0.071, and the MCL >50% boolean's narrow non-monotonic dip, whose closest
approach to zero is |delta|≈0.055 — just outside the cited band, not inside it). **10 conclusions are
robust across the entire explored +/-0.10 range**: PCL, ILFL, ISFL, both PFL rotation checks that
passed at baseline (raw-gap-sign and the normalized-gap), PTFL, CFL (informational), the CFL/DELT
antagonist cross-check, the already-known-baseline-FAIL PFL "amplified, raw-N" metric (robustly wrong
regardless of epsilon_r, not epsilon_r-fragile), and the tension-only violation count (robust, but
tautologically so — §6). 3 + 5 + 10 = 18.

---

## 6. Two things that look alike but are not: geometric robustness vs. definitional robustness

The task pre-registered the tension-only near-zero count as *expected fragile*. It measured the
**opposite** — 100% robust at every common-mode delta and every decorrelated seed at both magnitudes.
Forced diagnosis (an honest-negative-avoidance check, not accepting the flattering result at face
value): this is **not** geometric robustness like the engagement-direction claims above. It is a
**tautological identity of the force law's own branching** — `Blankevoort1991Ligament` returns
`force = 0` for `strain <= 0` **by construction**, true for *any* valid positive `slack_length`,
independent of both geometry and epsilon_r (the only failure mode would be `1 + reference_strain <= 0`,
which the harness's own guard rejects before it can happen). It verifies the **software is not buggy**,
not that the **biomechanics is real** — exactly the "never a tautology gate" anchor the task warns
against conflating with a real external check. Its trust tier should stay what it already was (an
implementation-correctness check); this sweep does not — and should not be read to — upgrade it.

By contrast, the engagement-direction conclusions that stayed robust (ILFL, ISFL, PFL-normalized,
PTFL, CFL, the antagonist cross-checks, PCL) are robust for a genuinely **geometric** reason: their
baseline force margins between the "tighter" and "looser" bands are large relative to how much a
+/-0.10 reference-strain shift can move the crossing point (e.g. ILFL 300N vs 17N — a >17x margin no
plausible epsilon_r perturbation closes). That is the "geometrically forced" case the task's
symmetric-QC section asks to distinguish from an artifact of too narrow a sweep — confirmed here
because the range explored (+/-0.10) is double the task's own cited band, and these held at 100%
decorrelated pass-rate even at the wider magnitude.

---

## 7. The MCL magnitude finding sharpens, not contradicts, an already-disclosed gap

`MECHANISM_KNEE_LIGAMENTS.md` §5/§8 already disclosed "deep-flexion magnitudes are qualitative, not
validated... quantitative trust bound ≲100deg flexion" as a kinematic-coupling limitation (LaiArnold's
prescribed 1-DOF knee vs. Lenhart2015's free 6-DOF source). This sweep adds a second, independent
reason the same region is untrustworthy, and shows the two compound:

| flexion angle | MCL peak \|strain\| range over epsilon_r in [-0.10,+0.10] | span |
|---:|---|---:|
| 90 deg (doc's own "trusted" region) | 26.8% - 38.0% | 11.2 pts |
| 100 deg (doc's own trust-bound edge) | 32.9% - 45.9% | 13.0 pts |
| 140 deg (doc's own already-flagged-untrusted endpoint) | 49.8% - 73.7% | 23.3 pts |

Geometric derivation (not a fit): per-bundle `d(strain)/d(delta) = length_at_pose / L_ref`, which grows
the further the queried pose is from the reference pose — so epsilon_r-sensitivity and
kinematic-coupling-mismatch sensitivity are not independent uncertainty axes, they multiply with
flexion depth. **No pose in this sweep is magnitude-robust** (even the "trusted" 90 deg region swings
11 points); only the qualitative engagement-DIRECTION claims are defensible as magnitude-independent.

---

## 8. Trust-ledger re-tiering — recommendation (not applied to `MECHANISM_TRUST_LEDGER.md` in this pass)

Per isolation rules for this task (repo written concurrently by another instance; touch only files
created by this task), the shared `docs/MECHANISM_TRUST_LEDGER.md` is **not** edited here. The
recommended re-tiering, for whichever process next updates that ledger:

| conclusion | current framing | recommended tier after this sweep |
|---|---|---|
| Knee MCL "58% peak strain, non-physiological" | stated as a specific number | **re-tier to diagnosed-gap**: report as "49.8-73.7% depending on the un-measured reference-strain assumption, compounding with an already-disclosed kinematic-coupling gap" — never re-quote "58%" as if it were a precise, validated figure |
| Spine PLL "tighter in flexion" | PASS, reported with its own two-round OODA fix already disclosed | **re-tier to method-only**: the sweep confirms (independently) that this PASS sits on a knife-edge less than 1/4 of the cited plausible band from FAIL — the original doc's own §5 "R=0.015 was chosen for a comfortable margin" already hinted at this; this sweep quantifies it |
| Spine ISL "tighter in flexion" | PASS | **re-tier to method-only**: flips inside the cited band (−0.038) |
| Knee ACL "tighter in extension" | PASS | **re-tier to method-only**: flips inside the cited band (−0.030) |
| ATFL / spine ALL / LF / SSL engagement direction | PASS | **keep tier, add a footnote**: robust within the task-cited +/-0.05 band, not unconditionally (flips −0.06 to −0.07) |
| ILFL / ISFL / PTFL / CFL / PFL-normalized / antagonist checks / PCL | PASS | **keep tier as-is** — geometrically forced, large margins, survives the FULL doubled exploratory range under both adversaries |
| "0/1397 tension-only violations" | PASS, cited as a headline count | **keep tier, add explicit tautology footnote**: do not let a future reader treat its epsilon_r-robustness as biomechanical evidence — it is a force-law code identity |

---

## 9. Honest gaps

1. **Baldwin 2009's own numeric SD/range table was not extractable** (paywalled, abstract-only) — the
   swept magnitudes are the task-cited band plus a doubled exploratory band, not a literature-measured
   SD. Flagged, not hidden, same class as this repo's own pre-existing Pintar-1992 gap.
2. **Common-mode and decorrelated are not the only possible perturbation models** — e.g. a
   per-*ligament-group* (not per-bundle) common shift, or a correlated-within-group/decorrelated-
   across-group hybrid, were not tested; the two adversaries used bracket the extremes (fully
   correlated vs. fully independent) but do not cover every intermediate correlation structure.
3. **20 seeds per decorrelated magnitude** is enough to resolve pass-rates to roughly +/-10-15%
   resolution (binomial), sufficient to distinguish "clearly fragile" (45-85%) from "clearly robust"
   (100%) as done here, but not to pin exact fragile-fraction confidence intervals tightly.
4. **The MCL non-monotonic dip (§5) was fine-grained only in [-0.10,-0.04]**, the region it was
   observed in; the rest of the range was only grid-sampled at the original 21-point resolution
   (adequate for the monotonic cases, whose single flip was independently bisected to +/-0.0005 in all
   7 instances).
5. **This sweep operates entirely on `subject2_unified_v2.osim`**, inheriting every honesty flag
   already disclosed in the 4 source docs (single donor subject, LaiArnold's prescribed 1-DOF knee,
   generic hip/ankle geometric constructions, generic spine stiffness values, etc.) — it tests
   sensitivity to epsilon_r specifically, not any of those other, already-disclosed gaps.

---

## 10. Files

- `scripts/msk/epsr_sensitivity_sweep.py` — the harness: baseline-reference-strain map from the 3 fork
  configs, closed-form perturbed-model builder, generalized sweep/band-mean helpers, all 18 headline
  checks, baseline-reproduction external-anchor gate, common-mode + decorrelated sweeps, survives/flips
  ledger construction.
- `data/msk_smoketest/epsr_sensitivity/epsr_sensitivity_evidence.json` — raw sweep output: baseline
  reproduction report, all 21 common-mode runs' full conclusion sets, all 40 decorrelated runs' full
  conclusion sets, the per-conclusion ledger (first-flip detection at 21-point grid resolution), MCL
  peak-strain continuous series.
- `data/msk_smoketest/epsr_sensitivity/epsr_sensitivity_ledger.json` — **the consolidated
  survives/flips deliverable**: methodology, closed-form, citations (with the honesty flag on Baldwin
  2009's unextractable table), pre-registration, baseline-reproduction anchor, the bisected
  (+/-0.0005-precision) flip points for all 7 monotonic-flip conclusions, decorrelated pass-rates, the
  MCL non-monotonic dip characterization, and the flexion-depth-vs-fragility supplementary check (§7).
- Read (unmodified, external anchors): `docs/MECHANISM_KNEE_LIGAMENTS.md`,
  `docs/MECHANISM_HIP_ANKLE_LIGAMENTS.md`, `docs/MECHANISM_SPINE_LIGAMENTS.md`,
  `docs/MECHANISM_UNIFIED_V2.md`, `scripts/msk/knee_ligament_config.json`,
  `scripts/msk/hip_ankle_ligament_config.json`, `scripts/msk/spine_ligament_config.json`,
  `data/msk_models/subject2_unified_v2.osim`.
- **Not modified** (isolation / concurrent-write rule — recommendation only, §8):
  `docs/MECHANISM_TRUST_LEDGER.md`.
