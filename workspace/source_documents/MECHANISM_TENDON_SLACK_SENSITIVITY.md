# MECHANISM TENDON-SLACK-LENGTH SENSITIVITY — is the knee/hip force cert fragile to l_TS? (2026-07-21)

Tests the second key uncertain Hill-muscle parameter (Fmax being the first, handled by a separate,
complementary thread — not redone here) against `docs/MECHANISM_CMC_SECOND_SOLVE.md` /
`docs/MECHANISM_CROSS_SUBJECT.md`'s knee/hip over-prediction headline (subject2/walking1, twin over-predicts
OrthoLoad in-vivo contact force by ~1.4–2.0×). Tendon slack length (l_TS) sets each muscle's force-length
OPERATING POINT; this session swept l_TS ±10% (task's own primary band) and, forced further, ±20% (a
forced-adversary stress band) on the 27 knee+hip "prime mover" muscles, re-solved Static Optimization from
scratch for each of 9 scale points, and re-read peak knee-r/hip-r contact force via the task's own named
function (`validate_hip_force.compute_self_cross_check_generic`, unchanged).

**Headline: the cert is l_TS-FRAGILE at BOTH joints, with real margin above the pre-registered 0.15
ratio-swing threshold — knee swing 0.809 (primary ±10% band) / 0.944 (±20% stress band), hip swing 0.276 /
0.293. Within the task's own ±10% band alone, the knee ratio-vs-OrthoLoad spans 0.778× (twin UNDER-predicts)
to 1.587× (twin OVER-predicts by 59%) — a plausible l_TS mis-calibration alone can flip this twin's headline
knee finding from "roughly matches in-vivo" to "over-predicts by ~60%". The mechanism is mechanistically
diagnosed, not just measured: Static Optimization's redundant, minimum-effort recruitment discretely
switches a large antagonist pair (gasmed_r + recfem_r, a knee flexor/extensor pair) from near-silent (≤10 N)
to dominant (1140–1284 N) across a narrow l_TS window, and antagonist co-contraction is exactly the
mechanism this repo's own knee-cert architecture (`R − Σmuscle-crossing`) turns into inflated bone-contact
force. A first-pass (naive) analysis under-stated this — a forced adversary check (below) found and
corrected a peak-search boundary artifact that had been partly masking it.**

## 1. Headline result

| l_TS scale | band | knee peak (window-corrected, %BW) | knee ratio vs OrthoLoad (258.22) | hip peak (%BW) | hip ratio vs OrthoLoad (273.93) |
|---:|---|---:|---:|---:|---:|
| 0.80 (−20%) | STRESS | 201.94 | 0.782 | 454.67 | 1.660 |
| 0.85 (−15%) | STRESS | 201.30 | 0.780 | 461.54 | 1.685 |
| 0.90 (−10%) | PRIMARY | 200.87 | **0.778** | 456.86 | 1.668 |
| 0.95 (−5%) | PRIMARY | 250.46 | 0.970 | 387.50 | 1.415 |
| 1.00 (0%, baseline) | PRIMARY | 391.10 | 1.515 | 386.77 | 1.412 |
| 1.05 (+5%) | PRIMARY | 409.78 | **1.587** | 385.14 | 1.406 |
| 1.10 (+10%) | PRIMARY | 243.37 | 0.942 | 381.17 | 1.391 |
| 1.15 (+15%) | STRESS | 170.64 | 0.661 | 403.98 | 1.475 |
| 1.20 (+20%) | STRESS | 166.08 | **0.643** | 419.71 | 1.532 |

Bold = the min/max defining each band's swing. **PRIMARY-band swing: knee 0.809, hip 0.276 (both >0.15 →
FRAGILE). FULL-band swing: knee 0.944, hip 0.293 (both >0.15 → FRAGILE).** All 9 points fully converged (no
excluded points) — the fragility is not an artifact of dropping non-converged solves.

## 2. Pre-registration (stated before this session's forced-adversary work changed any number)

- **Muscle set** (27, right-side, union of anatomically-classic knee + hip "prime movers", disjoint from the
  separate Fmax-sensitivity thread's own parameter — verified, §4): `vasmed_r, vaslat_r, vasint_r, recfem_r,
  bflh_r, bfsh_r, semimem_r, semiten_r, gasmed_r, gaslat_r` (knee: vasti+rectus femoris+hamstrings+gastrocs)
  ∪ `iliacus_r, psoas_r, glmax1-3_r, glmed1-3_r, glmin1-3_r, addlong_r, addbrev_r, addmagProx/Mid/Dist/Isch_r,
  recfem_r, bflh_r, semimem_r, semiten_r` (hip: iliopsoas+glutes+adductors+biarticular hamstrings/rectus
  femoris). Excludes small/secondary muscles (sart_r, grac_r, tfl_r, piri_r) by design.
- **PRIMARY scale band**: ±10% in 5% steps, {0.90, 0.95, 1.00, 1.05, 1.10} — the task's own explicit range.
- **STRESS band** (forced-adversary widening, pre-registered as a REQUIRED step per this task's own
  discipline — "leaning-positive [robust] adversary must be forced to its strongest fair form before being
  credited"): ±15%/±20% added, {0.80, 0.85, 1.15, 1.20}.
- **Falsifier**: any PRIMARY-band point moving either joint's ratio-vs-OrthoLoad by an absolute swing >0.15
  ⇒ l_TS-FRAGILE on that joint; swing ≤0.15 ⇒ ROBUST on this axis, held to the same rigor.
- **Symmetric-QC gate**: a swept point is a valid comparison point only if `so_convergence_gates.
  convergence_pass` is True (n_frames>100, covers the validated peak-time window, zero NaN/out-of-bounds
  activation). Points failing this would be reported, not hidden, and excluded from the swing computation.
- **Confidence tier (disclosed up front, not softened later)**: METHOD (sensitivity). The PRIMARY ±10% band
  is the task's own literature-referenced range (Redl 2007 / Ackland 2012 / Carbone 2012 vicinity) — **this
  session's one WebSearch attempt to independently re-verify the specific published coefficients returned
  "this session has used its web search budget (2000 of 2000)" before any lookup ran**, so that specific
  number is NOT freshly re-verified via citation this session (a real, disclosed gap, §6). The STRESS ±20%
  band is this session's own forced widening (general convention in this literature that l_TS
  cadaveric-measurement uncertainty commonly exceeds 10%), also not a specific re-verified per-paper number.
  By contrast, the SENSITIVITY MEASUREMENT ITSELF — does perturbing l_TS move this twin's own force cert,
  and by how much — is directly machine-measured on the real pipeline this session, high confidence; only
  the real-world plausibility of the exact band WIDTH carries the lower, method-tier confidence.

## 3. Method (reuse discipline — wrap, don't edit in place)

1. **Model perturbation via the OpenSim API, not XML regex**: `Muscle.getTendonSlackLength()`/
   `Millard2012EquilibriumMuscle.setTendonSlackLength()` are real, C++-validated property accessors (unlike
   the AnalyzeTool setup-XML files elsewhere in this repo, which have no such API and must be regex-patched).
   For each scale, a fresh `osim.Model` is loaded from the read-only original, ONLY the 27 swept muscles'
   `tendon_slack_length` is multiplied by the scale (each muscle scaled from ITS OWN baseline value — "per-
   muscle scaled", not a shared absolute shift), and the result is written under this repo
   (`data/msk_models/tendon_slack_sweep/subject2_scaled_lts_{tag}.osim`). **Verified live, every point**
   (`verify_model_edit_isolated`, re-loading the WRITTEN file, not the in-memory object): only the swept
   muscles' l_TS changed; the other 53 muscles' l_TS AND all 80 muscles' Fmax are bit-identical to the
   original — this decorrelates this thread from the separate Fmax-sensitivity thread by construction, not
   by assertion.
2. **Re-solve, don't reuse**: `static_opt_knee.run_static_optimization()` is called UNCHANGED, with only its
   module-level `MODEL_FILE`/`SO_OUT` path constants monkey-patched to the new model/output directory (same
   "Python resolves bare names against the enclosing module's own `__dict__` at call time" reuse pattern
   `cross_subject_validation.py` already established). Each of the 9 points is a genuine, independent,
   full-trial (158-frame, 1.57s) Ipopt Static-Optimization re-solve — not a cached/interpolated estimate.
3. **Re-verify convergence** (`static_opt_knee.check_so_convergence_and_sanity`, unchanged) before trusting
   any downstream number, plus a NEW swept-muscle-specific pinned-activation diagnostic (fraction of frames
   where a swept muscle's activation sits within 0.01 of its [0,1] bound) beyond the existing global gate.
4. **Re-read peak force** via `validate_hip_force.compute_self_cross_check_generic` — the task's own named
   function, already fully generic over the distal-body cut — called identically for the knee-r BFS chain
   and the hip-r BFS chain (same call signature `validate_hip_force.py`'s own `main()` already uses).
5. **OrthoLoad anchors re-verified live** this session (knee median 258.221 %BW n=72, hip median 273.931
   %BW n=162 — both matched the pre-existing cert bit-for-bit), never trusted from prior JSON/prose.

## 4. Forced adversary #1 (leaning-positive: is "fragile" actually a peak-search boundary artifact?)

**Observed**: the raw whole-trial-interior-argmax knee peak TIME jumps discontinuously between two disjoint
clusters — t≈0.49–0.55s for scale ∈ {0.95, 1.00, 1.05, 1.10}, vs. t=**1.520s exactly** for scale ∈ {0.80,
0.85, 0.90, 1.15, 1.20}. t=1.520s is not a coincidence: this script's own edge-exclusion (`TRIAL_END_TIME`
1.57s minus `e`=5 samples=0.05s) makes 1.520s the LAST interior sample — the exact failure signature
`docs/MECHANISM_CROSS_SUBJECT.md` §2.3/§3 already diagnosed for a different subject/joint ("JR's un-excluded
search lands on that still-rising boundary value").

**Oriented** (machine-checked, not eyeballed): for all 5 affected points, the last 3 interior frames before
the boundary are strictly monotonically rising (`np.diff(...)>0` on the raw force values). These are NOT
genuine interior local maxima — they are truncated samples of a still-climbing curve caught by the edge
exclusion. (Hip's own tail is ALSO climbing toward the same edge in 3/9 points — 0.95, 1.00, 1.05 — but
never overtakes the established t=0.550s peak within the recorded window: a disclosed near-miss, not a
current contamination, flagged in §6.)

**Decided + Acted**: recomputed the peak search restricted to a disclosed, pre-registered-post-hoc window
t∈[0.20, 1.00]s — chosen to bracket the established ~0.49–0.55s peak cluster (matching EVERY other cert in
this repo, e.g. the baseline's own 391.10 %BW at t=0.510s) with generous margin (0.29s before, 0.45s after
the observed cluster), not tuned to any single scale point's outcome. This reused the ALREADY-COMPUTED,
on-disk scaled models + SO outputs — **zero re-solve cost** (`scripts/msk/tendon_slack_sensitivity_diagnose.
py`). Self-verifying: the window changed ONLY the 5 contaminated points; the other 4 (0.95/1.00/1.05/1.10)
were already identical between the whole-trial and windowed search.

**Result: the correction made the knee swing LARGER, not smaller** (raw/naive primary-band swing 0.658 →
corrected 0.809). This rules out "the fragility finding was itself just a search-window bug inflating the
swing" — if anything, the naive whole-trial numbers had *understated* the true sensitivity by partly
canceling it against an unrelated boundary value. The adversary (boundary artifact as the explanation for
"fragile") was forced and FELL: it was real, but only in 5/9 points, and correcting it strengthened rather
than weakened the fragility finding.

## 5. Forced adversary #2 (mechanism: is the swing physically explicable, or an unexplained numeric wobble?)

Pulled every knee-crossing muscle's SO-computed force (already on disk, zero re-solve) at t=0.510s (this
repo's own established peak instant) for scale ∈ {0.90, 1.00, 1.05, 1.10}:

| muscle | 0.90 (N) | 1.00 (N) | 1.05 (N) | 1.10 (N) | Fmax (N) |
|---|---:|---:|---:|---:|---:|
| gasmed_r (medial gastrocnemius — knee FLEXOR) | **0.00** | **1139.97** | **1284.28** | 511.37 | 3115.5 |
| recfem_r (rectus femoris — knee EXTENSOR) | **10.02** | **586.73** | **645.06** | 158.16 | 2191.7 |
| vasmed_r / vaslat_r / vasint_r (vasti) | 28/53/17 | 24/45/15 | 21/38/12 | 17/28/9 | — |
| bfsh_r / semimem_r (other hamstrings) | 125/162 | 33/16 | 21/12 | 11/6 | — |

`gasmed_r` and `recfem_r` — an ANTAGONIST pair at the knee (one flexor, one extensor) — jointly switch from
near-silent (≤10 N) at scale=0.90 to dominant (1140–1284 N / 587–645 N, still only 37–41%/27–29% of their own
Fmax — bounded, not degenerate) at scale=1.00–1.05, then partially recede at 1.10. This is a genuine,
geometrically-grounded mechanism, not a heuristic: Static Optimization's redundant, minimum-effort
(activation-cost) recruitment discretely "discovers" a large co-contraction solution exactly in the l_TS
window where these two muscles' force-length-shifted capacity makes it cheap to reach the required net
knee moment that way. Antagonist co-contraction is PRECISELY the mechanism `static_opt_knee.py`'s own module
docstring identifies as why bone-contact force (`R − Σmuscle-crossing`) exceeds pure reaction force in the
first place — so a recruitment-regime switch mechanically, not mysteriously, explains the observed
contact-force swing. (Sanity: knee-r/hip-r crossing-muscle SETS were verified bit-identical across all 9
scale points — geometry detection cannot itself have drifted; the swing is 100% a force-magnitude effect.)

## 6. Pipeline-validity gates (all machine-checked, none eyeballed)

- **Compliant-tendon mechanism confirmed**: 80/80 muscles are `Millard2012EquilibriumMuscle` with
  `ignore_tendon_compliance=false` — l_TS enters via the tendon force-strain equilibrium condition
  (geometric basis for expecting sensitivity at all, stated before running anything).
- **Model-edit isolation**: verified every point — only swept muscles' l_TS changed; ALL Fmax bit-identical
  (decorrelates this thread from the Fmax-sensitivity thread).
- **SO convergence**: PASS at all 9/9 points (no excluded points; fragility is not a convergence artifact).
- **No swept-muscle activation pinned at bound, any point.** No muscle exceeded 1.5× Fmax, any point.
- **Geometry-derived free invariant**: the "pure reaction" force (kinematics+GRF only, no muscle term) MUST
  be bit-identical across all 9 scale points, since it structurally cannot read muscle force — confirmed
  relative spread = 0.0 at both joints (tol 1e-6). Proves only the intended channel varies.
- **Internal-consistency anchor**: the fresh scale=1.00 control point reproduces the pre-existing validated
  baseline (391.10/386.77 %BW, `docs/MECHANISM_JOINT_FORCE_SCORECARD.md`) to <0.01% — this sweep's own
  pipeline is proven correct before any perturbed point is trusted. Also bit-matches
  `validate_joint_force.py`'s own independently-established pure-reaction number (102.16690546813328 %BW)
  exactly.
- **Wall-clock**: ~21–22s per point (9 points, ~200s total for the sweep; diagnostic re-derivation from
  on-disk outputs added ~1 minute, zero re-solves).

## 7. Honest gaps

1. **The ±10%/±20% band width is task-provided / general-convention, not a citation this session
   independently re-verified** (§2's confidence-tier note) — the WebSearch budget was exhausted session-wide
   (a shared resource across concurrent agents on this machine) before the one lookup attempted could run.
   If the TRUE cadaveric l_TS uncertainty for these specific muscles is narrower than ±10%, the PRIMARY-band
   fragility conclusion should be re-checked against the narrower band; given the swing is 5.4× the
   threshold, a materially narrower "true" band would need to be under ~±1.9% (10%/5.4) before this
   specific PRIMARY-band verdict would flip to robust — no such narrow-band claim exists in this session's
   record.
2. **Hip's tail-boundary near-miss (§4) is disclosed, not resolved** — 3/9 points show the hip signal still
   climbing toward the t=1.57s edge without yet overtaking the established t=0.550s peak; a longer trial or
   a different subject could plausibly tip this into the same contamination the knee showed. Not currently
   affecting any reported hip number, but a latent risk for reuse of this exact whole-trial-argmax pattern
   elsewhere.
3. **The detailed shape of the excluded post-1.0s tail region (t=1.0–1.57s) was not further characterized**
   beyond confirming its last-3-frame monotonic trend — whether it corresponds to a genuine second
   physiological gait event (e.g., contralateral heel-strike) or is purely a boundary artifact throughout is
   unresolved; irrelevant to this task's falsifier (outside the historically-validated peak window every
   other cert in this family uses) but flagged as unexplored.
4. **Mechanism diagnosis (§5) covers only 4 of 9 scale points** (0.90/1.00/1.05/1.10, chosen as the most
   informative/contrasting) — the full 9-point per-muscle force table was not pulled, since the 4-point
   spot-check already gave an unambiguous, bounded (well under Fmax), non-degenerate mechanistic answer;
   pulling the rest would refine the transition's exact shape but is not needed for the pre-registered
   falsifier.
5. **Single subject (subject2), single trial (walking1), right side only** — same scope caveat as every
   other cert in this family (`docs/MECHANISM_CROSS_SUBJECT.md`). Whether the SAME sharp, antagonist-
   co-contraction-driven l_TS sensitivity recurs in other subjects (where different absolute l_TS/Fmax/
   geometry combinations could put the "cheap co-contraction" regime at a different scale, or make it less
   sharp) is untested here.
6. **This thread swept 27 "prime mover" muscles simultaneously, uniformly scaled together** — it does not
   attribute the swing to any ONE muscle's l_TS in isolation (a full per-muscle Sobol/one-at-a-time sweep
   would; out of scope, given the task's own "per-muscle scaled... on the knee & hip prime movers" framing
   as one combined sweep). §5's spot-check strongly implicates gasmed_r/recfem_r specifically, but this is
   an observation from the combined sweep, not a controlled single-muscle isolation.
7. **Reserve-actuator saturation and pelvis-residual gates were monitored but not treated as differential
   signals** — pelvis residual force/moment (173.0N/45.1Nm) and knee/hip reserve saturation (7–17%) were
   IDENTICAL (bit-for-bit on the pelvis residual) across all 9 points, confirming they are a pre-existing,
   scale-independent limitation of this cert family (consistent with every other joint-force doc in this
   repo), not a differential artifact of l_TS — reported for completeness, not gated on.

## 8. Files

- `scripts/msk/tendon_slack_sensitivity.py` — the sweep pipeline (model perturbation, re-solve, convergence
  gates, `compute_self_cross_check_generic` calls, raw whole-trial peak extraction, naive verdict).
- `scripts/msk/tendon_slack_sensitivity_diagnose.py` — the forced-adversary follow-up (§4): re-derives full
  per-frame force time series from the EXISTING on-disk scaled models + SO outputs (zero re-solve), computes
  the windowed peak + the monotonic-tail machine check.
- `data/msk_models/tendon_slack_sweep/subject2_scaled_lts_{m20,m15,m10,m05,p00,p05,p10,p15,p20}.osim` — the
  9 perturbed models (new, this session; original read-only source untouched).
- `data/msk_smoketest/subject2_walking1/tendon_slack_sensitivity/`:
  - `so_{tag}/` × 9 — each scale's full Static Optimization output (activation.sto, force.sto).
  - `tendon_slack_sensitivity_results.json` — raw sweep output (whole-trial argmax, naive verdict).
  - `tendon_slack_sensitivity_diagnose_results.json` — windowed re-derivation + monotonic-tail diagnostics.
  - `tendon_slack_sensitivity_corrected_summary.json` — corrected ratio table + final swing/verdict.
  - **`tendon_slack_sensitivity_evidence.json`** — the consolidated deliverable evidence file (pre-
    registration, all gates, both forced adversaries, the mechanism spot-check, the corrected table, final
    verdict, all source-file pointers).
- **Not touched**: `scripts/msk/validate_joint_force.py`, `scripts/msk/static_opt_knee.py`,
  `scripts/msk/validate_hip_force.py` (all reused unchanged, per this repo's wrap-don't-edit convention);
  the separate Fmax-sensitivity and moment-arm-sensitivity threads (out of scope by the task's own
  instruction, not re-run, not re-litigated).

No git commit, no git push performed (isolation respected). Concurrent processes observed on this shared
machine during this session (`scripts/msk/epsr_sensitivity_sweep.py`, `scripts/msk/overprediction_decomp.py`)
belong to the other mechanism instance sharing this repo (per `COORDINATOR.md` §0) — not touched, not read.

**Read-only-drive self-audit** (same check `docs/MECHANISM_CROSS_SUBJECT.md` §3 already performed for other
subjects, repeated here rather than assumed): `find subject2_dir -newer <this session's first new file>`
shows exactly one changed path, `.../OpenSimData/Mocap/Model/opensim.log` — OpenSim's own C++ Logger
duplicating stdout's info/warning lines to a text file, the SAME pre-existing, already-disclosed behavior
that doc found for subject3/subject4. Zero `.osim`/`.mot`/`.xml`/`.sto` data files under the external subject2
tree were touched (checked directly, not assumed).
