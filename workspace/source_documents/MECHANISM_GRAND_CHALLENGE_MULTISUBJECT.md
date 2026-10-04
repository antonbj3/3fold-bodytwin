# MECHANISM — Grand Challenge COHORT: multi-subject in-vivo knee validation (+ a units-bug correction to the flagship)

**Task:** scale the flagship same-subject knee validation (`MECHANISM_GRAND_CHALLENGE_DM.md`) to a cohort, using
the already-extracted First/Second/Third Grand Challenge competitions (`data/external/grand_challenge_{first,
second,third}/`). Measured anchors first (guaranteed value), Tier-1 twin only if tractable.

**Pre-registration (stated before any twin-side comparison):** falsifier = "do the distinct subjects' twin
Tier-1 lower bounds over-predict their measured in-vivo knee force ~1.5–1.8x like DM (confirming the flagship
as a cohort)?" Threshold for the measured-BW cross-check: documented vs static-standing-GRF-derived body
weight must agree within 15%. Threshold for trial identity: knee-force-csv vs trajectories-csv row counts must
match exactly.

**Headline result:** the measured cohort is delivered in full (3 distinct new subjects, §1–2). The twin side
is a diagnosed-gap this session (§5). While building the pipeline, DM's own pre-existing measured-force number
(the flagship's central anchor) was used as a calibration check — and was found to be **wrong by almost
exactly the lbf→N conversion factor** in the pre-existing flagship doc/JSON (§3), which **reverses** the
flagship's headline conclusion once corrected. This is reported with full forcing/verification detail below,
flagged for the doc owner to fix at the source (not edited here — isolation rule: touch only files created
this session).

---

## 1. Dedup: TRUE distinct-subject count = 3 (First/Second/Third), DM(6th) status unresolved

Read directly from each competition's own **"Competition Data Description.pdf"** (Fregly et al.), Section II
"Subject Information" — not inferred from filenames alone:

| competition | code | height | weight | sex | knee side | implant gen | year |
|---|---|---:|---:|---|---|---|---:|
| First | **JW** | 166 cm | 64.6 kg | unstated | right | 1st-gen EKnee (4 uniaxial load cells) | 2009 |
| Second | **DM** | 172 cm | 67.0 kg | unstated | right | 2nd-gen eTibia (6-axis load cell) | 2010 |
| Third | **SC** | 167 cm | 78.4 kg | **female** | **left** | 2nd-gen eTibia (6-axis load cell) | 2011 |

No two of {JW, DM, SC} share height+weight+sex+knee-side — SC is unambiguously a third, distinct person
(different sex and instrumented side alone settle it). **TRUE distinct-subject count across First+Second+Third
= 3.** Each competition's own PDF explicitly frames its subject as "a new patient" relative to the prior
competition (Second's PDF re: JW, Third's PDF re: DM) — direct textual confirmation this is a reuse-free trio.

**DM(second) vs the pre-existing flagship "DM"(Sixth Competition, already used in `MECHANISM_GRAND_CHALLENGE_DM.md`):**
same 2-letter code, both right-knee, weight 67.0 kg (2010) vs 69.98 kg (flagship's scaled `.osim` mass, +4.4%,
a plausible multi-year drift). New this session: DM(second)'s own independently-collected normal-gait peak
(272.6–296.2 %BW, median 294.3, n=5) closely matches the flagship DM(sixth)'s **corrected** normal-gait peak
(250.7–312.5 %BW, median 286.8, n=7 — see §3) — two independently-collected, independently-processed datasets
landing in the same range is suggestive of the same person's gait mechanics recurring across competition years,
but **not decisive** (TKA literature broadly clusters in 2–3×BW regardless of identity). No Sixth-Competition
"Data Description" PDF was available in this session's scope to settle it directly. **Verdict: flagged, not
merged** — treated as a distinct cohort entry pending that document.

## 2. MEASURED in-vivo peak knee force — the guaranteed-value deliverable

For each subject's **"ngait"/"ngait_og" trials** (each competition's own PDF confirms this label = "the
subject's normal gait pattern," i.e. level walking — not a treadmill/pathological variant), peak resultant
axial tibiofemoral contact force, own documented body weight:

| subject | n trials | min %BW | median %BW | max %BW | JCF/GRF ratio (physical sanity) | static-BW cross-check |
|---|---:|---:|---:|---:|---:|---|
| JW (1st) | 5 | 263.4 | **271.1** | 284.8 | 1.84× | +1.33% vs documented (CV 0.0011, clean) |
| DM (2nd) | 5 | 272.6 | **294.3** | 296.2 | 1.96× | −14.3% vs documented (CV 0.45, noisy — see caveat) |
| SC (3rd) | 5 | 256.9 | **262.5** | 265.5 | 2.14× | +0.91% vs documented (CV 0.0024, clean) |
| DM (6th, flagship, **corrected** — see §3) | 7 (ngait_og only) | 250.7 | **286.8** | 312.5 | 2.57× | n/a (not re-checked this session) |

**Extraction method** (per-subject, unit-verified from each competition's own README, not assumed): JW's
1st-gen EKnee reports 4 uniaxial load cells (PM/AM/AL/PL); their sum equals total axial contact force exactly
(medial+lateral regression coefficients are complementary per channel, verified algebraically). DM(2nd)/SC's
2nd-gen eTibia reports a 6-axis load cell at a known tray point; `+Fz` = inferior = compression (each
competition's own `eTibia Coordinate System.pdf`), and the medial/lateral regression coefficients on Fz sum to
~1.000 (Meyer et al. 2011 for DM; an analogous fit for SC), i.e. **Fz alone ≈ total axial force**. All raw
values are lbf/in-lbf per every competition's own README; converted via the exact CODATA factor
(4.4482216152605 N/lbf).

**Independent verification performed this session** (not just trusted from the concurrent cohort-extraction
script that had already run — `scripts/msk/grand_challenge_multisubject.py` /
`data/msk_smoketest/grand_challenge_cohort/multisubject_measured_results.json`, a prior/concurrent process's
own output, not edited here per isolation): fresh code
(`scripts/msk/grand_challenge_multisubject_verify.py`) re-parsed every raw CSV independently and reproduced
every peak to full float precision (`cross_check_vs_concurrent_script.pass_exact_match = true` for all 3
subjects) — a genuine two-pipeline agreement, not a restatement. **Trial identity**: knee-force-csv vs
trajectories-csv row counts match exactly for all 15 gait trials (0 mismatches, 0 unmatched). **Physical
sanity anchor** (new this session, decorrelated from the extraction pipeline itself): peak vertical GRF (from
the independently-recorded force-plate `*_grf.csv`) vs peak knee joint-contact-force — JCF must exceed GRF
during single-limb stance (muscle+ligament crossing force only ever *adds* to the pure reaction). Measured
ratios 1.84–2.57× — physically sane and consistent with published TKA literature (D'Lima/Fregly, typically
2–3×BW peaks).

**Honest caveat (symmetric QC):** DM(second)'s static-standing body-weight cross-check is noisy — tried 4
different static trials (`staticfor`, `staticin`, `staticout`, `staticmaxiso`); all show high CV (0.42–0.47,
vs JW/SC's 0.001–0.002) and undershoot the documented 67.0 kg by 5–20 kg. This looks like a genuine
data-collection quirk specific to this subject's Second-Competition static trials (not a bug in the
window-search algorithm — re-tried with fresh independent code, same result), not disqualifying (the
documented weight is a peer-reviewed clinical value, the more authoritative source), but disclosed rather than
hidden.

## 3. MAJOR FINDING: a units bug in the pre-existing flagship's measured-force anchor

While establishing DM(second) as an independent calibration point, the flagship's own already-published DM
number (56–83 %BW, `MECHANISM_GRAND_CHALLENGE_DM.md`) looked anomalous: **less than the subject's own peak GRF**
(110.5 %BW, from the flagship's own `grf_sanity_gate`) — biomechanically impossible during single-limb stance,
since muscle/ligament force only adds to (never subtracts from) the GRF-driven reaction, a principle the
flagship doc itself invokes elsewhere. Forced this adversary rather than accepting it as a quirky "honest"
number:

**Bit-exact reproduction of the bug** (`scripts/msk/grand_challenge_multisubject_verify.py`, function
`dm_sixth_flagship_units_bug_forensic`): the flagship's "resultant_N" values are `sqrt(Fx²+Fy²+Fz²)` computed
directly on the **raw lbf-scale** `Fx,Fy,Fz` columns of `DM_ngait_og1_knee_forces.csv`, with the lbf→N
conversion (×4.4482216152605) never applied, then mislabeled `_N`:

| flagship-claimed value | independently reproduced (raw, no conversion) | match |
|---|---:|---|
| `global_peak_resultant_N` = 437.9118301667586 | 437.9118301667586 | exact |
| `localwindow_2p30_3p10_peak_resultant_N` = 399.28207936244775 | 399.28207936244775 | exact |
| `value_at_twin_peak_instant_N` = 244.14239902155464 | 244.14239902155464 | exact |
| `DM_ngait_og3` cross-trial peak = 56.360222724197094 %BW | 56.360222724197094 %BW | exact |
| `DM_smooth1` cross-trial peak = 79.17771346972773 %BW | 79.17771346972773 %BW | exact |
| `DM_ngait_tm_ss1` cross-trial peak = 83.42406464901092 %BW | 83.42406464901092 %BW | exact |

Five independent points match to full float precision — not a coincidence. Every competition's own README in
this dataset family (First/Second/Third, read directly this session) states force columns are in lbf; the
6th-competition file shares the identical column layout (`Time(sec),Fx,Fy,Fz,Tx,Ty,Tz,GON,GRFz`) and raw
magnitude scale.

**External physical anchor (not a tautology — independently-sourced GRF, unaffected by this bug):**

| | uncorrected (as published) | corrected (×4.4482216152605) |
|---|---:|---:|
| DM(6th) global peak %BW | 63.81 | **283.84** |
| JCF / GRF ratio (GRF=110.47%BW, independently `.sto`-sourced) | **0.578× — physically impossible** | **2.57× — physically sane** |
| same-stance-phase window %BW | 58.18 | **258.80** |
| 17-trial cross-check median %BW | 68.08 | **302.83** |
| ngait_og-only (n=7) median %BW | 64.47 | **286.77** |

The corrected value is also **externally consistent** with this session's 3 freshly-derived subjects
(JCF/GRF 1.84–2.14×) and with DM(second)'s own independently-collected normal-gait cohort (median 294.3%BW,
§1) — two decorrelated anchors (physical GRF-ratio necessity + cross-subject/cross-year consistency), neither
of which is circular with the bug-detection itself.

**Revised twin-vs-measured verdict for DM:** twin Tier-1 = 104.4 %BW (unchanged, twin-side math untouched by
this bug). Against the corrected measured anchor, the ratio **reverses**:

| framing | twin | measured (corrected) | ratio | direction |
|---|---:|---:|---:|---|
| global peak | 104.4 | 283.84 | 0.368× | twin **under**-predicts by 2.7× |
| same-stance-phase window (primary) | 104.4 | 258.80 | 0.403× | twin **under**-predicts by 2.5× |
| ngait_og-only median | 104.4 | 286.77 | 0.364× | twin **under**-predicts by 2.7× |

This is the **opposite** of the flagship's published headline ("twin over-predicts by 1.5–1.8×"). **This
finding is NOT acted on here** — `MECHANISM_GRAND_CHALLENGE_DM.md` and `MECHANISM_GRAND_CHALLENGE_DM_TIER2.md`
are pre-existing files, not created this session (isolation rule: touch only files created here). Flagged for
the doc owner/orchestrator to independently re-verify and correct at the source. **Falsifier for this claim,
so it can be checked independently in under a minute:** open
`data/external/simtk_kneeloads_DM_measured_force/DM_ngait_og1_knee_forces.csv`, take row at `t≈1.675`
(`Fx=4.14, Fy=13.95, Fz=437.67`), compute `sqrt(4.14²+13.95²+437.67²) = 437.91` — compare to the flagship
JSON's `global_peak_resultant_N`. If they match to 6 figures with no ×4.4482 anywhere in between, the bug is
confirmed; if a version of the flagship pipeline is found that DOES apply the conversion and still gets
437.91, this finding is refuted.

## 4. Falsifier verdict (partial — measured side only)

The task's pre-registered falsifier ("do the cohort's twin Tier-1 lower bounds over-predict ~1.5–1.8× like
DM?") **cannot be evaluated this session** for JW/DM(second)/SC — Tier-1 twin is a diagnosed-gap (§5), so
there is no twin number to compare for the new subjects. What IS delivered: a verified, cross-checked,
externally-anchored **measured cohort** (§2) — real value on its own, and the necessary precondition for
anyone to run the twin-side comparison later. What is also delivered: the discovery that the falsifier's own
reference point (DM's measured anchor) needs correction before it can be meaningfully evaluated **on any
subject**, cohort or flagship.

## 5. Tier-1 twin: diagnosed-gap (not attempted this session)

**Reconnaissance performed** (not a zero-effort skip): compared marker-label conventions in all 3 new
subjects' own `trajectories.csv` headers against every marker-set/IK asset already in this repo.

- JW/DM(second) (1st/2nd competition) use an **older ~31/43-marker convention**
  (`R.Heel, R.Midfoot.Superior, R.Shank.Superior, R.Thigh.Superior, R.ASIS, Sacral, ...`) with no existing
  marker-set XML in this repo built for it.
- SC (3rd competition) uses a **richer ~62/70-marker convention** (`Sternum, Xiphoid, Thoracic, R.Asis, R.Psis,
  Lumbar, R.Elbow/Wrist/Radius/Ulna`, computed joint-centers `RHJC/RKJC/RAJC`) that closely resembles the
  pre-existing flagship's own `DM_markers.xml` — the most promising reuse path, but that marker set was built
  for DM's implant-specific `.osim`, not a generic model, and would still need: retargeting onto a **generic**
  LaiArnold/gait2392 model per the task's own spec, a fresh Scale step for SC's own anthropometry, and a
  GRF/motion axis-convention remap (SC's raw `*_grf.csv` lab frame ships explicit `Xglobal/Yglobal/Zglobal`
  marker triads — direct evidence the native axes are not OpenSim's Y-up convention).
- The repo's existing IK infrastructure (`scripts/msk/smoke_test_ik.py`, `subject_specific_scaling.py`) is
  built entirely around a different, already-scaled subject (`subject2`, an OpenCap `LabValidation_withVideos`
  dataset) with its own unrelated marker set — not reusable as-is.
- None of the 3 new subjects have ready OpenSim-native external-loads (`.mot`/`.sto`) or scale-setup XMLs
  (unlike the pre-existing DM-sixth flagship data, which arrived JAM-team-prepared).

**Conclusion:** a compatible marker-set + scale + IK + joint-reaction pipeline is multi-hour, per-subject
infrastructure work across two distinct marker-set generations — beyond this task's "only if tractable" bar
and this session's lean/time budget on a CPU-contended shared box. Per the task's own explicit permission,
this is a diagnosed-gap, reported honestly rather than silently skipped or forced into a low-quality attempt.

## 6. Files

- New this session (mine): `scripts/msk/grand_challenge_multisubject_verify.py` (independent re-derivation +
  units-bug forensic + final JSON), `data/msk_smoketest/grand_challenge_cohort/multisubject_cohort_VERIFIED_FINAL.json`,
  this doc.
- Pre-existing / concurrent (read-only, not edited): `scripts/msk/grand_challenge_multisubject.py` +
  `data/msk_smoketest/grand_challenge_cohort/multisubject_measured_results.json` (a concurrent process's own
  cohort-extraction, independently cross-checked and confirmed exact-match by this session's script);
  `data/external/grand_challenge_{first,second,third}/` (raw competition data, extracted before this session);
  `data/external/simtk_kneeloads_DM_measured_force/`, `data/external/opensim_jam_grand_challenge_DM/`,
  `data/msk_smoketest/DM_ngait_og1/final/dm_grand_challenge_results.json` (flagship's own data/output, site of
  the units bug); `docs/MECHANISM_GRAND_CHALLENGE_DM.md`, `docs/MECHANISM_GRAND_CHALLENGE_DM_TIER2.md` (flagged
  for correction, not edited).

**Confidence tier:** measured cohort (§2) = in-vivo-anchored, B-gradeable (cross-pipeline exact match + trial
identity + physical GRF/JCF sanity, all machine-verified). Units-bug finding (§3) = watertight (5 bit-exact
reproductions + 2 decorrelated external anchors: physical necessity + cross-subject/cross-year consistency).
Tier-1 twin (§5) = diagnosed-gap, sanctioned outcome per task spec.
