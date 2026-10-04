# MECHANISM CONTACT-FORCE MUSCLE DECOMPOSITION — where does the over-prediction live? (2026-07-21)

Localizes the twin's ~1.4-2.0x knee/hip contact-force over-prediction (`docs/MECHANISM_CROSS_SUBJECT.md`,
`docs/MECHANISM_CMC_SECOND_SOLVE.md`) to specific muscles by decomposing the already-validated
subject2/`walking1` Static-Optimization contact force, at its own established peak instants, into
per-muscle contributions along the joint-contact axis. Uses `static_opt_knee.py`'s /
`validate_hip_force.py`'s committed, **UNCHANGED** crossing-muscle geometry
(`knee_crossing_muscles_and_forces`, `compute_self_cross_check`, `compute_self_cross_check_generic`) —
zero new muscle-force solve, zero edits to either file (md5-verified identical to
`docs/MECHANISM_CROSS_SUBJECT.md`'s own recorded hashes, Sec.5 below).

**Headline: the quadriceps+hamstrings "pathological co-contraction" story named in the pre-registered
falsifier is CLEANLY REFUTED at the knee — all 4 hamstring muscles sit at 0.010-0.066 activation
(essentially the SO numerical floor) at the peak, never close to the 0.5 "maximally firing" threshold.
But the over-prediction is NOT a clean "genuine model property" story either: it is dominated by a REAL
antagonist pair (gastrocnemius, a knee flexor, 48.7% of the muscle-driven force; quadriceps, mostly
carried by rectus femoris, 22.2%), geometrically verified (live OpenSim moment arms, not a hardcoded
table) as true opposing-moment-arm antagonists — gastrocnemius' timing is a clean literature match,
rectus femoris' is a borderline 5.9-percentage-point-early case, AND rectus femoris is carrying quad
duty that published EMG assigns to vastus lateralis instead (a real, specific, redundancy-resolution
substitution, not classic antagonist co-contraction). The hip shows a THIRD, distinct pattern: no 1-2
muscle dominance at all (top-2 groups = 47.5% of the muscle sum, vs the knee's 71.0%) — a genuinely
diffuse spread across gluteus medius, iliopsoas, rectus femoris, TFL and gluteus minimus — with its own
independent anomaly (gluteus MAXIMUS, the classic hip extensor, contributes just 0.7%, essentially
silent, while gluteus MEDIUS is active well past its own published early-stance window, corroborating
docs/MECHANISM_EMG_TIMING.md's already-published FAIL for that exact muscle). **A second, decorrelated
data point sharpens this further: this subject's own REAL surface EMG (`glmed1_r` is one of only 7
electrodes in this dataset, reused from a prior session's CEINMS-style forward simulation, not
recomputed) implies a force at this exact instant of only ~100 N — vs SO's own 585 N for the same
muscle at the same instant (hybrid/method-controlled ratio 0.17-0.18) — i.e. the hip's single largest
contributor may itself be substantially over-recruited relative to the subject's actual measured muscle
activity, a second, independent redundancy-resolution-artifact candidate alongside the knee's
rectus-femoris/vastus-lateralis substitution (Sec.7).** Every per-muscle contribution is verified, not
assumed, to sum EXACTLY back to the already-published %BW total (residual 1.4e-10 relative at the knee,
7.9e-12 at the hip — floating-point noise, not an approximation).

## 0. Pre-registration (stated before any per-muscle number was computed)

**C / ¬C, thresholds fixed in advance** (`scripts/msk/contact_muscle_decomp.py` module docstring, written
before Step 9's numbers existed):

- **H_pathological** (the task's own framing): the classic quadriceps group AND hamstring group are
  BOTH "maximally firing" simultaneously at the peak — operationalized as **mean-group-activation > 0.5**
  for BOTH groups. If true: localizes the over-prediction to SO's redundancy resolution (an artifact).
- **H_genuine_property** (the task's own alternative): the **top-2 anatomical groups** account for
  **>60%** of the total positive muscle-driven contribution, AND their activation timing falls inside (or
  within a **10-percentage-point buffer** of) a published EMG ON-window for this gait phase. If true: the
  over-prediction is a property of the muscle-driven model at this instant, not a co-contraction artifact.
- **Anything else** is reported as its own, honestly-described outcome — not forced into either bucket.

**Forced adversary (this task leans toward finding "a story," so the adversary to force is a decomposition
that looks clean but doesn't actually reconstitute the number it claims to explain)**: every per-muscle
contribution is projected onto the SAME contact-force axis (`n_hat`) that defines the already-published
%BW headline, and the full set (muscle contributions + the muscle-free Newton's-law term) is checked to
sum EXACTLY back to that published total — not eyeballed, not approximated (Sec.5). A second, independent
adversary was forced for the antagonist LABELS themselves: rather than trusting a hand-written
flexor/extensor lookup table, every muscle's sagittal-plane function is read from a LIVE OpenSim moment
arm at the actual peak pose, anchored against undisputed anatomy (Sec.6) — this caught and corrected one
mislabel a naive table would have kept silently wrong (`tfl_r`, Sec.6).

**Confidence tier**: the reconstitution/geometry claims (Sec.5-6) are **method-only** (exact linear
algebra / live OpenSim kinematics — a mathematical identity, not an empirical claim, so "confidence tier"
in the empirical sense does not apply; it either holds to float precision or it doesn't, and it does).
The physiological interpretation (is the dominant pattern "normal" or "inflated"?) is
**published-plausibility**, anchored to Rajagopal et al. 2016 (PMID 27392337, direct quotes, already
verified live in this repo's `docs/MECHANISM_EMG_TIMING.md`) and Falconer & Winter 1985 (PMID 3987606, CCI
formula, already verified live in `docs/MECHANISM_EMG_DRIVEN.md`) — both reused verbatim, not
re-fetched. One narrow sub-claim (Sec.9) additionally cross-checks against this subject's own REAL
surface EMG (reused from `docs/MECHANISM_EMG_DRIVEN.md`'s prior session, not recomputed) — that specific
cross-check is **in-vivo-anchored** but covers only 4/13 knee-crossing muscles (`gasmed_r`, `vasmed_r`,
`vaslat_r`, plus the hamstring pair already used in the classic-CCI subset) and 1/25 hip-crossing muscles
(`glmed1_r`), disclosed per-claim below (Sec.9, Sec.11 item 7), not blanket-applied to the whole document.

## 1. Method (geometric, exact — not a heuristic)

`compute_self_cross_check` / `compute_self_cross_check_generic` already compute, at every one of 158
frames, the **full vector** (never previously reported past its own norm)

```
F_bone_contact(t) = R_old_vec(t) - sum_m [ tension_m(t) * unit_dir_m(t) ]
```

where `unit_dir_m` is muscle *m*'s unit pulling direction on the distal sub-chain (from
`static_opt_knee.py`'s own sign-fixed `knee_crossing_muscles_and_forces`) and `R_old_vec` is the
muscle-free Newton's-law term (gravity + inertia + GRF only). Let `n_hat = F_bone_contact(t*) /
|F_bone_contact(t*)|` at the established peak `t*` — **this is "the joint-contact direction"**, i.e. the
direction of the already-published, already-validated %BW headline. Because a dot product is linear,
projecting every additive term onto `n_hat` and summing **exactly** reconstitutes `|F_bone_contact(t*)|`
— an algebraic identity, verified in Sec.5, not assumed. Per-muscle contribution to the reported
contact-force magnitude:

```
c_m = -tension_m(t*) * ( unit_dir_m(t*) . n_hat )
```

`unit_dir_m` points from the muscle's distal (inside) crossing point toward its proximal (outside) one —
the direction tension pulls the sub-chain. A compressive bone-on-bone contact force pushes the sub-chain
the OTHER way, so a muscle whose pull is anti-aligned with `n_hat` **adds** to the compressive contact
force (`c_m>0`) regardless of whether it is a flexor or extensor. Measured, not assumed: **0/13 knee
crossers and 0/25 hip crossers had a negative contribution** — every single crossing muscle, active or
not, adds to the compressive load. This is the "moments partially cancel, forces add" mechanism already
qualitatively noted (10/13 knee muscles, raw tension only, no axis projection, no reconstitution check,
no hip table) in this module's own docstring and `docs/MECHANISM_STATIC_OPT.md` Sec.7 — made exact and
exhaustive here.

## 2. Headline numbers

| | KNEE (t=0.51s, 44.1%GC) | HIP (t=0.55s, 47.1%GC) |
|---|---:|---:|
| Total (re-derived, tie-back to published) | 391.10 %BW (2999.26 N) | 386.77 %BW (2966.05 N) |
| R_pure (muscle-free Newton's-law term, projected) | 87.67 %BW (22.4%) | 84.13 %BW (21.8%) |
| Top-1 contributor | gasmed_r: 139.94 %BW (35.8%) | glmed1_r: 75.90 %BW (19.6%) |
| Top-2 contributor | recfem_r: 75.96 %BW (19.4%) | recfem_r: 57.40 %BW (14.8%) |
| Top-3 contributor | gaslat_r: 50.72 %BW (13.0%) | iliacus_r: 55.62 %BW (14.4%) |
| Top-2 **anatomical group** share of muscle-driven force | **71.0%** (gastrocnemius+quadriceps) | **47.5%** (iliopsoas+glute_med) |
| n muscles needed for 80% of muscle-driven force | 3 (gasmed, recfem, gaslat = 87.9%) | 6 (glmed1, recfem, iliacus, psoas, tfl, glmin1 = 85.0%) |

The knee is **concentrated** (3 muscles carry 87.9% of the muscle-driven share); the hip is **diffuse** (6
muscles needed to clear 80%, and the top single muscle carries only a quarter of the muscle-driven
share). This is exactly the distinction the task asked to resolve — and it resolves DIFFERENTLY at the
two joints, not uniformly.

## 3. KNEE — full ranked table (all 13 crossing muscles, exhaustive)

| muscle | anatomical group | sagittal function (live moment arm) | activation | tension (N) | %Fmax | contribution (%BW) | % of total |
|---|---|---|---:|---:|---:|---:|---:|
| gasmed_r | gastrocnemius | flexor (ma=+0.0264) | 0.441 | 1140.0 | 36.6% | 139.94 | 35.8% |
| recfem_r | quadriceps | extensor (ma=-0.0545) | 0.227 | 586.7 | 26.8% | 75.96 | 19.4% |
| gaslat_r | gastrocnemius | flexor (ma=+0.0307) | 0.282 | 398.0 | 25.3% | 50.72 | 13.0% |
| tfl_r | tensor fasciae latae | **extensor** (ma=-0.0182) | 0.308 | 124.8 | 30.3% | 14.92 | 3.8% |
| vaslat_r | quadriceps | extensor (ma=-0.0501) | 0.010 (floor) | 45.2 | 0.9% | 5.83 | 1.5% |
| bfsh_r | hamstrings | flexor (ma=+0.0262) | 0.066 | 32.5 | 5.8% | 4.13 | 1.1% |
| vasmed_r | quadriceps | extensor (ma=-0.0511) | 0.010 (floor) | 24.2 | — | 3.16 | 0.8% |
| sart_r | sartorius | flexor (ma=+0.0199) | 0.098 | 22.0 | — | 2.87 | 0.7% |
| semimem_r | hamstrings | flexor (ma=+0.0403) | 0.010 (floor) | 14.6 | — | 1.91 | 0.5% |
| vasint_r | quadriceps | extensor (ma=-0.0523) | 0.010 (floor) | 14.6 | — | 1.90 | 0.5% |
| bflh_r | hamstrings | flexor (ma=+0.0354) | 0.010 (floor) | 8.7 | — | 1.13 | 0.3% |
| semiten_r | hamstrings | flexor (ma=+0.0493) | 0.010 (floor) | 4.9 | — | 0.64 | 0.2% |
| grac_r | gracilis | flexor (ma=+0.0312) | 0.010 (floor) | 2.4 | — | 0.32 | 0.1% |

**All 4 hamstring muscles (bflh_r, bfsh_r, semimem_r, semiten_r) sit at 0.010-0.066 activation** — the SO
numerical floor, not a substantial fraction of it. `H_pathological` (quad>0.5 AND hamstring>0.5) is
**decisively FALSE**: the classic quad/hamstring pair (vasmed_r+vaslat_r vs semiten_r+bflh_r, matching
`docs/MECHANISM_EMG_DRIVEN.md`'s own subset) reads **agonist_mean=0.0100, antagonist_mean=0.0100** — both
at floor. (Their raw CCI ratio computes to 0.999, which LOOKS like near-perfect balance — this is the
SAME scale-invariant floor-ratio artifact `docs/MECHANISM_EMG_DRIVEN.md` Sec.6 already diagnosed and
warned about: a ratio of two tiny, near-equal numbers reads "high" regardless of whether either signal is
real. The `H_pathological` test is deliberately built on the absolute 0.5 threshold, not the ratio, for
exactly this reason.)

**The real dominant pair is gastrocnemius (flexor) + rectus-femoris-carried quadriceps (extensor) — a
genuine, geometrically-confirmed antagonist pair, just not the illustrative one.** Anatomical-group
roll-up: gastrocnemius 48.7% of muscle-driven force, quadriceps 22.2%, TFL 3.8%, hamstrings 2.0%,
sartorius 0.7%, gracilis 0.1%. Broad extensor-vs-flexor CCI (moment-arm-derived groups): activation-based
0.988 (diluted by averaging several floor-level secondary muscles into each group — disclosed, not
hidden, Sec.8), **contribution-based 0.671** (the more mechanistically meaningful number: 67% of what a
perfectly-balanced antagonist pair would produce, weighted by how much force each side actually
contributes to the joint, not just how many muscles are in each bucket).

## 4. HIP — full ranked table (all 25 crossing muscles, exhaustive; top 10 shown individually, remainder <1% each)

| muscle | anatomical group | sagittal function (live moment arm) | activation | tension (N) | %Fmax | contribution (%BW) | % of total |
|---|---|---|---:|---:|---:|---:|---:|
| glmed1_r | glute_med | extensor* (ma=-0.0183) | 0.574 | 585.3 | 53.5% | 75.90 | 19.6% |
| recfem_r | quadriceps | flexor (ma=+0.0288) | 0.174 | 458.6 | 20.9% | 57.40 | 14.8% |
| iliacus_r | iliopsoas | flexor (ma=+0.0382) | 0.476 | 487.6 | 47.7% | 55.62 | 14.4% |
| psoas_r | iliopsoas | flexor (ma=+0.0267) | 0.251 | 344.2 | 24.1% | 36.48 | 9.4% |
| tfl_r | tensor fasciae latae | flexor (ma=+0.0251) | 0.337 | 133.2 | 32.4% | 17.02 | 4.4% |
| glmin1_r | glute_min | flexor* (ma=+0.0054) | 0.301 | 113.4 | 30.3% | 14.75 | 3.8% |
| glmed2_r | glute_med | extensor* (ma=-0.0305) | 0.175 | 110.5 | — | 14.41 | 3.7% |
| glmin2_r | glute_min | extensor* (ma=-0.0027) | 0.213 | 75.0 | — | 9.78 | 2.5% |
| glmin3_r | glute_min | extensor* (ma=-0.0081) | 0.133 | 41.3 | — | 5.38 | 1.4% |
| sart_r | sartorius | flexor (ma=+0.0212) | 0.113 | 23.7 | — | 3.09 | 0.8% |
| **glmax1/2/3_r (all 3 parts, summed)** | glute_max | extensor* (ma≈-0.06 to -0.08) | **0.010 each (floor)** | 20.5 total | — | 2.67 | 0.7% |
| everything else (12 muscles: glmed3_r, adductors x6, hamstrings x3 [bflh/semimem/semiten], piriformis, gracilis) | — | mixed | 0.010-0.015 each (floor) | 77.8 total | — | 10.14 | 2.6% |

\* `glmed`/`glmin`/`glmax`'s PRIMARY mechanical action is frontal-plane abduction, not sagittal
flexion/extension — the "extensor"/"flexor" tag here is their measured SECONDARY sagittal-plane lean at
this exact pose (small moment-arm magnitudes, e.g. glmin2_r's -0.0027 m vs glmax's -0.06 to -0.08 m,
correctly reflecting that abduction, not extension, is their main job). Reported because the task asks
specifically about the sagittal antagonist question, not because these muscles are being reclassified as
"really" extensors.

**No 1-2-muscle dominance**: top-2 anatomical groups (iliopsoas 23.8%, glute_med 23.7% — statistically
tied, 4.4 N apart out of 706 N) sum to only 47.5% of the muscle-driven force, well under the 60%
pre-registered dominance threshold. **Gluteus maximus — the textbook "prime mover" hip extensor — is
essentially silent** (all 3 parts at the 0.010 activation floor, 0.7% of the total contact force
combined) while gluteus medius alone (0.574 activation, 53.5% of its own Fmax) carries more contact-force
weight than any other single muscle at either joint tested this session. The classic hip
flexor/extensor pair (iliopsoas vs glute max) reads **agonist_mean=0.364, antagonist_mean=0.010, CCI=0.054**
— a real, un-diluted, clearly LOW co-contraction index (the opposite of "balanced antagonist
co-contraction": iliopsoas dominant, glute max essentially off).

## 5. Symmetric-QC: does the decomposition reconstitute the published total? (machine-checked, exact)

| check | knee | hip |
|---|---:|---:|
| Re-derived peak vs `docs/MECHANISM_STATIC_OPT.md`/`MECHANISM_HIP_FORCE.md`'s published %BW | 391.0995 vs 391.10 (0.00012% diff) | 386.7679 vs 386.77 (0.00056% diff) |
| Peak time vs published | t=0.510s vs 0.51s (exact) | t=0.550s vs 0.55s (exact) |
| Per-muscle-vector-sum recompute (single peak frame) vs full-158-frame trajectory pass, same frame | rel diff **1.05e-09** | rel diff **5.23e-11** |
| **R_pure + Σc_m vs published `\|F_bone_contact\|`** (the decomposition-sums-to-total gate) | **rel residual 1.44e-10 — PASS** | **rel residual 7.91e-12 — PASS** |
| Live moment-arm sign anchor (known extensors/flexors load the expected sign at the actual peak pose) | PASS | PASS |

All four gates PASS by many orders of magnitude of margin (residuals are floating-point noise, not
approximation error) — this is an exact linear-algebra identity, verified rather than assumed, per the
task's own symmetric-QC clause ("a decomposition that doesn't reconstitute the total is wrong").

## 6. Geometric antagonist classification — corrects one mislabel a hardcoded table would keep

Rather than a hand-written flexor/extensor lookup, every crossing muscle's sagittal function is read from
`Muscle.computeMomentArm` **at the model's actual peak pose** (not a generic/default pose), anchored
empirically: known-undisputed knee extensors (recfem_r, vasmed_r) and hip flexors (iliacus_r, psoas_r)
verified to load one sign; known knee flexors (bflh_r, semiten_r, gasmed_r, gaslat_r) and hip extensors
(glmax2_r, bflh_r, semiten_r) verified to load the opposite sign (`verify_sign_anchor`, PASS at both
joints, Sec.5). Every other muscle is then classified by that SAME verified convention.

**This caught a real correction**: `tfl_r` (tensor fasciae latae) is commonly treated as a neutral "hip
abductor/knee stabilizer" (its own comment in this repo's `KNEE_ANATOMICAL_MUSCLES_EXPECTED` doesn't even
assign it a sagittal role). Live geometry shows it loads a clear, non-borderline **negative** moment arm
about `knee_angle_r` (-0.0182 m) at the peak pose — i.e. it acts as a **weak knee extensor** at this
specific pose, not a neutral stabilizer, and its 3.8% contribution is reported as such rather than
silently mislabeled or omitted.

**A second, independent geometric check — literally "do antagonists pull in opposite directions?"**: for
every pair of measurably-active muscles (activation ≥ 0.05) at each joint, the raw Cartesian `unit_dir`
vectors' pairwise cosine similarity was computed (`pairwise_most_antiparallel`). Result, measured not
assumed: **every single pair at both joints has a POSITIVE cosine** (knee: range +0.729 to +0.996 across
15 pairs; hip: range +0.382 to +0.999 across 45 pairs) — **zero pairs pull in literally opposing Cartesian
directions.** This refines the folk "antagonists pull the tibia opposite ways" framing: every crossing
muscle, agonist or antagonist by MOMENT, pulls the sub-chain in a roughly SIMILAR Cartesian direction
(proximally, toward the joint) — it is this shared-direction pull, opposing the compressive push of bone
contact, that makes co-contracting antagonists ADD to joint compression rather than cancel (Sec.1). The
lay language of "pulling opposite ways" describes the ROTATIONAL moment (which does partially cancel),
not the Cartesian force vector on the sub-chain (which does not cancel — it accumulates). The hip's wider
range (down to +0.382, vs the knee's +0.729 floor) reflects its geometrically richer muscle arrangement
(flexors/extensors/abductors/rotators crossing at many circumferential positions, vs the knee's more
planar flexor/extensor layout).

## 7. Published-EMG literature-timing check (Rajagopal et al. 2016, PMID 27392337 — reused verbatim from `docs/MECHANISM_EMG_TIMING.md`, not re-fetched)

| joint | muscle group | %GC at peak | literature ON-window | inside? | gap to nearest edge | read |
|---|---|---:|---|---|---:|---|
| knee | gastrocnemius | 44.1 | late-stance (30-60) | **YES** | 0 pts | **clean match — physiological** |
| knee | rectus femoris | 44.1 | early-stance (0-30) + pre-swing (50-73) | no | 5.9 pts (before pre-swing opens) | **borderline / early-onset**, weakest-anchored muscle in the literature source (§3 of `MECHANISM_EMG_TIMING.md`) |
| knee | vasti (quad) | 44.1 | early-stance (0-30) | no | 14.1 pts | literature timing alone reads "consistent with OFF" — but real EMG (Sec.9) says this is likely wrong; see the substitution finding below |
| knee | hamstrings | 44.1 | early-stance + terminal-swing | no | 14.1 pts | **literature ALSO predicts OFF here — twin's floor activation MATCHES, not contradicts** |
| hip | gluteus maximus | 47.1 | broad, "relatively constant" (terminal-swing + 0-50) | yes (broad) | 0 pts | literature predicts SOME tone; twin shows **floor (0.010)** — a real mismatch, but on a low-confidence/broad literature window |
| hip | gluteus medius | 47.1 | early-stance only (0-30) | no | 17.1 pts | **twin's 0.574 activation is well past its own literature window** — corroborates `docs/MECHANISM_EMG_TIMING.md`'s already-published FAIL for glmed (Jaccard=0.439, twin ON 5-57%GC vs literature 0-30%GC), now pinned to this exact contact-force-peak instant |
| hip | rectus femoris | 47.1 | early-stance (0-30) + pre-swing (50-73) | no | 2.9 pts | borderline / early-onset, same caveat as the knee |
| hip | iliopsoas (iliacus+psoas) | 47.1 | **not in this repo's literature dict** (no channel/window exists for iliopsoas in `validate_emg_timing.py`'s `EMG_LITERATURE`) | — | — | **data gap, disclosed** — the hip's #1-tied contributor cannot be timing-checked against this repo's existing literature anchor at all |

A genuinely important, specific, real finding surfaces here that is **not** captured by the pre-registered
`H_pathological`/`H_genuine_property` binary: **the twin's SO solution substitutes `recfem_r` for `vasti`
as its knee-extensor channel.** Rajagopal et al. 2016 name **vastus lateralis** specifically as the
measured early-stance knee extensor ("knee extensor (vastus lateralis) activity in early stance"); this
twin's vasti sit at the 0.010 activation floor throughout while `recfem_r` carries essentially all of the
quadriceps-side load (0.227 activation, 19.4% of the knee's total contact force). Rectus femoris and
vasti both extend the knee (a genuine mechanical redundancy — SO's objective can substitute one for the
other without changing the net knee moment it must satisfy), so this looks like a **specific,
identifiable SO-redundancy-resolution artifact** — a real candidate mechanism for part of the
over-prediction, distinct from (and more precise than) generic "co-contraction." **This is not just a
literature-timing inference**: this subject's own real EMG (Sec.9) shows `vasmed_r`/`vaslat_r` tension
23x higher than SO's near-floor value at this exact instant — a decorrelated, in-vivo-anchored (for this
2-muscle subset, with an already-disclosed normalization-uncertainty caveat) confirmation that vasti are
not really silent, strengthening this from a plausible hypothesis toward a real, evidenced finding. Not
proven causally this session (recfem's dual hip-flexor role means some of its recruitment could genuinely reflect an
emerging hip-flexion-moment demand near mid-stance, not pure knee-extensor substitution — an untested,
disclosed alternative explanation, Sec.10).

## 8. Co-contraction indices (Falconer & Winter 1985, PMID 3987606 — CCI = 2·min(a,b)/(a+b), reused verbatim)

| pair | joint | agonist mean act. | antagonist mean act. | CCI | reading |
|---|---|---:|---:|---:|---|
| classic quad (vasmed_r+vaslat_r) vs hamstring (semiten_r+bflh_r) | knee | 0.0100 | 0.0100 | 0.999 (**floor-ratio artifact, not real balance** — same diagnosed issue as `MECHANISM_EMG_DRIVEN.md` §6) | H_pathological test uses the 0.5 absolute threshold, not this ratio, for exactly this reason — **FALSE** |
| broad extensor vs flexor (moment-arm-derived groups) | knee | 0.113 (diluted by floor muscles) | 0.116 (diluted) | 0.988 (activation, diluted) / **0.671 (contribution-weighted, more meaningful)** | real, substantial — driven by gastroc+recfem, not quad+hamstring |
| classic iliopsoas vs glute-max | hip | 0.364 | 0.010 | **0.054** | clearly LOW — iliopsoas dominant, glute max essentially silent, opposite of balanced co-contraction |
| broad extensor vs flexor (moment-arm-derived groups) | hip | 0.086 (diluted) | 0.155 (diluted) | 0.713 (activation, diluted) / **0.759 (contribution-weighted)** | real, substantial — but reflects glute_med+glute_min's secondary sagittal lean plus iliopsoas, not a classic sagittal antagonist story |

The activation-based "broad" CCI numbers are diluted by averaging several genuinely-floor-level secondary
muscles into each functional bucket alongside 1-2 dominant ones — disclosed, not hidden; the
contribution-weighted variant (which reflects how much actual force each side contributes, not how many
muscles are in each bucket) is the more mechanistically meaningful number for "how much of the
compressive load comes from balanced opposing contributions," and is reported as the lead number for that
specific question.

## 9. Bonus corroboration: this subject's own real surface EMG (reused, not recomputed)

This subject's dataset has real surface EMG on only 7/80 muscles (`docs/MECHANISM_EMG_DRIVEN.md`):
`soleus_r, gasmed_r, vasmed_r, vaslat_r, semiten_r, bflh_r, glmed1_r`. Of the muscles that dominate the
per-muscle tables above, **4 of them are covered** — `gasmed_r` (knee #1), `vasmed_r`/`vaslat_r` (the
knee's "silenced" quad channel), and `glmed1_r` (hip #1) — a materially better real-data check than a
first pass found (an earlier draft of this document under-counted this coverage; corrected here, not
silently). `recfem_r`, `iliacus_r`, `psoas_r`, `tfl_r`, and `gaslat_r` have **no** right-side channel in
this subject's file (disclosed gap, inherited from `docs/MECHANISM_EMG_DRIVEN.md` Sec.9 item 2) and
cannot be checked this way. All values below are reused verbatim from a prior session's CEINMS-style,
real-EMG-driven forward simulation (`emg_driven_results.json` / `emg_hybrid_force.sto`,
method-controlled against `naive_fwd_only_force.sto` to isolate the pure EMG-content effect from the
forward-integration method-switch effect) — not recomputed this session, and not free of that prior
session's own disclosed EMG-to-force-normalization uncertainty (Sec.9 gap #1 there, repeated below
per-muscle where it applies).

| muscle | joint role | SO tension (N) | real-EMG-hybrid tension (N) | hybrid/SO ratio | reading |
|---|---|---:|---:|---:|---|
| `gasmed_r` | knee #1 (35.8% of total) | 1140.0 | 917.7 | **0.805** | **CORROBORATES** SO — inside the pre-registered ±30-43% agreement band (`docs/MECHANISM_EMG_DRIVEN.md` §5); this muscle's own EMG-normalization ceiling was already at 100% of its 3-trial max (no known under-ceiling bias) |
| `vasmed_r` | knee, SO-silenced (0.8% of total) | 24.3 | 563.5 | **23.2** | real EMG implies **dramatically more** force than SO's near-floor value — reinforces Sec.7's recfem/vasti-substitution finding with real, decorrelated data, not just literature timing |
| `vaslat_r` | knee, SO-silenced (1.5% of total) | 45.2 | 1041.6 | **23.1** | same direction, same magnitude of disagreement |
| `glmed1_r` | hip #1 (19.6% of total) | 585.3 | 99.8 | **0.17** | real EMG implies **much LESS** force than SO chose — the hip's single largest contributor may itself be over-recruited (Sec.10) |

**Reading the vasti number carefully, not just quoting it**: `docs/MECHANISM_EMG_DRIVEN.md`'s own Sec.9
already disclosed that `vasmed_r`'s per-channel EMG normalization ceiling (this trial's own peak, used
as "100%") is only **66.2%** of the 3-trial maximum across `walking1`/`walking2`/`walking3` — i.e. even
this real-EMG number likely **understates** true excitation, not overstates it, so the 23x
disagreement's DIRECTION (vasti are not silent; SO's floor-level activation for them does not match
measured muscle activity) is credible, while its exact MAGNITUDE (23x specifically) is not
precisely calibrated — a real signal, reported with the uncertainty that already exists on it, not a
clean, exact number.

**The knee and hip point in OPPOSITE directions on this specific check**: the knee's #1 contributor
(`gasmed_r`) is corroborated by real EMG; the hip's #1 contributor (`glmed1_r`) is contradicted by it.
This is a genuinely different, joint-specific finding, not a uniform "SO always over- or
under-recruits" story — consistent with this document's broader finding that the knee and hip have
DIFFERENT mechanisms (Sec.2-4), not one shared artifact.

Hamstring real-EMG values (already reused in Sec.0/headline: 0.006 at the knee peak, 0.005 at the hip
peak, both well below SO's own 0.010-0.066) remain the strongest, cleanest real-data confirmation in this
document — an **in-vivo-anchored** (for this 2-muscle subset) confirmation that hamstring silence is not
an SO-specific artifact.

## 10. Pre-registered decision

| joint | H_pathological (classic antagonist pair both >0.5: quad+hamstring at knee, iliopsoas+glute-max at hip) | top-2 groups >60% AND literature-consistent | verdict |
|---|---|---|---|
| knee | **FALSE** (hamstring group mean = 0.024, all 4 muscles ≤0.066) | top-2 groups = 71.0% (>60%, dominance TRUE) but literature match is MIXED (gastroc clean, recfem 5.9pts early) | **MIXED**: dominance is real, mechanism is a genuine (geometrically-verified) antagonist pair — just not quad-vs-hamstring, and not cleanly literature-timed either. Refutes the task's pathological-co-contraction hypothesis; does not cleanly confirm the "genuine property, no caveats" hypothesis — lands on a third, specific, actionable finding (Sec.7's recfem/vasti substitution) that neither pre-registered bucket anticipated |
| hip | **FALSE** (iliopsoas mean=0.364, glute-max mean=0.010 — CCI=0.054, badly unbalanced, not "both maximally firing") | top-2 groups = 47.5% (<60%, dominance FALSE) | **NEITHER**: genuinely diffuse (6 muscles needed for 80%), not concentrated, not the classic antagonist-pair story either — its own specific anomaly is glute_max's near-total silence contradicted against a broad literature window, and glute_med's activity extending well past its own narrow one (independently corroborating an already-published FAIL, Sec.7) |

**Bottom line for the operator's trust question**: the twin's over-prediction is **not** explained by the
folk "antagonist pair both maximally firing" mechanism at either joint — that specific hypothesis is
refuted with large margin, not a knife-edge miss: the knee's hamstring group (SO mean 0.024, individual
muscles 0.010-0.066; real EMG 0.005-0.006 where available) sits 8-90x below the 0.5 threshold depending
on which number is used, and the hip's gluteus maximus (SO mean 0.010, all 3 parts at floor) sits 50x
below it — in every reading, by every available measure (this session's SO activations AND the subject's
own real EMG), the "both sides near-maximal" condition is not remotely close to being met at either
joint. But it is also not a clean,
caveat-free "genuine model property" either: the knee's real mechanism (gastrocnemius+rectus-femoris
antagonism, mostly-but-not-entirely literature-timed, with a specific rectus-femoris/vasti substitution
that looks like an SO redundancy-resolution choice) and the hip's real mechanism (diffuse
multi-muscle loading with an independently-corroborated gluteus-medius timing mismatch and a
gluteus-maximus silence that undershoots even a generous literature band) are two distinct, specific,
partially-actionable findings — not evidence of a single unifying artifact, and not proof the twin's
muscle recruitment is simply "correct." Both readings the task offered were partial; the decomposition's
value is in making the actual, specific, per-muscle picture visible and re-verifiable rather than in
forcing a binary verdict the data does not cleanly support.

## 11. Honest gaps

1. **Single trial, right leg, subject2 `walking1`, one instant per joint** — same scope caveat as every
   cert in this family; no claim of generality across subjects/trials/gait speeds/the left leg (though
   `docs/MECHANISM_CROSS_SUBJECT.md` already shows the over-prediction's DIRECTION and MAGNITUDE-RANGE
   generalize; this session's per-muscle MECHANISM was not re-tested on subject3/subject4).
2. **The literature-timing windows are the same population-average (Rancho Los Amigos/Perry 1992)
   convention `docs/MECHANISM_EMG_TIMING.md` already disclosed as not fit to this subject** — the 10-point
   buffer used for "borderline" calls (Sec.7) is a convention choice, not derived from a subject-specific
   EMG onset distribution.
3. **Rectus femoris's literature anchor is the weakest of those checked** (already flagged in
   `docs/MECHANISM_EMG_TIMING.md` §3 — its window is a standard convention extrapolated from a
   running-condition quote, not a walking-specific verified source) — the "5.9/2.9-point early, borderline"
   reads for recfem at both joints should be weighted accordingly, lower-confidence than the
   gastrocnemius/glmed reads.
4. **The rectus-femoris/vasti "substitution" explanation (Sec.7) is a plausible, disclosed HYPOTHESIS, not
   a proven mechanism** — this session did not decompose recfem's dual knee-extensor/hip-flexor moment
   contributions to test whether its hip role (not a pure knee-extensor substitution for vasti) explains
   some or all of its recruitment. A targeted follow-up (recfem's own moment-arm-weighted contribution to
   the hip-flexion moment at this instant, vs to the knee-extension moment) could resolve this directly.
5. **The activation-based "broad" extensor/flexor CCI is diluted by group-averaging muscles of very
   different sizes/floors** (Sec.8) — disclosed explicitly; the contribution-weighted CCI is offered as the
   more meaningful number for the specific "how much force comes from balanced opposition" question, but
   readers wanting a literature-comparable *activation*-based CCI should use the classic 2-muscle pairs
   (Sec.3/4's `classic_antagonist_pair`, not diluted) rather than the broad group version.
6. **Iliopsoas has no literature ON/OFF window in this repo's existing `EMG_LITERATURE` anchor** (Sec.7) —
   the hip's #1-tied contributor's timing plausibility could not be checked against literature at all this
   session, only against the (unbalanced) classic-pair CCI and the real-EMG cross-check's absence (no
   iliopsoas channel exists in the subject's real EMG file either, per `docs/MECHANISM_EMG_DRIVEN.md`).
7. **Real-EMG corroboration (Sec.9) covers 4/13 knee-crossing muscles (gasmed_r, vasmed_r, vaslat_r,
   plus hamstrings bflh_r/semiten_r) and only 1/25 hip-crossing muscles (glmed1_r)** — `gaslat_r` (the
   knee's #3 contributor), `recfem_r` (major at both joints), `iliacus_r`, `psoas_r`, and `tfl_r` have no
   real-EMG channel in this subject's dataset and are literature-verified only (published-plausibility
   tier), not in-vivo-anchored. The hip's dominant set beyond `glmed1_r` (iliopsoas, TFL, the rest of
   glute_med/min) is entirely in this uncovered category.
8. **`glmed`'s literature-mismatch finding here (Sec.7) reproduces, rather than independently discovers,
   an already-published result** (`docs/MECHANISM_EMG_TIMING.md`'s FAIL verdict for gluteus medius, Jaccard
   0.439) — reported as corroboration/triangulation, not as a new discovery.
9. **The pre-registered `DOMINANCE_TOP2_FRAC_THRESH=0.60` and `PATHOLOGICAL_CCI_ACT_THRESH=0.5` are
   this session's own choices**, reasonable but not themselves drawn from a published source — a
   different threshold could shift the knee's "MIXED" read toward either bucket (it clears 60% comfortably
   at 71.0%, so is not knife-edge on this specific threshold; the hip's 47.5% is further from 60% than from
   0%, so also not knife-edge, but neither number was cross-checked against a literature-derived
   concentration norm).
10. **Pelvis residual force still fails its own `<75N` band** (172.96 N here, inherited unchanged from
    every cert in this family, `docs/MECHANISM_JOINT_FORCE_SCORECARD.md` §3) — no RRA was run; joint-level
    reserve usage stays low (`joint_reserve_leaning_hard: {}`, re-verified this session) so the specific
    per-muscle numbers above are not obviously contaminated by it, but this is not a new resolution of that
    pre-existing gap.

## 12. Files

- `scripts/msk/contact_muscle_decomp.py` (new, this session) — the full, re-runnable pipeline. Imports
  `validate_joint_force.py` / `static_opt_knee.py` / `validate_hip_force.py` / `validate_emg_timing.py`
  **unchanged** (md5 recorded at run time: `e2bb15958cfe8d72f44a3e0ac3d95950` /
  `d00f962c1a8873606874fffd65e828c8` / `a6514b9f081e0b4767fb4b0f94e60200` / `94f6fa5a3ecaa14846e63c8d059f88b0`
  — the first three match `docs/MECHANISM_CROSS_SUBJECT.md` Sec.5's own recorded hashes exactly). Reuses
  `sok.knee_crossing_muscles_and_forces`, `sok.compute_self_cross_check`,
  `vhf.compute_self_cross_check_generic`, `sok.check_so_convergence_and_sanity`,
  `vet.detect_gait_cycle`/`vet.pct_gc`/`vet.EMG_LITERATURE`/`vet.MUSCLE_GROUPS` verbatim. New code: the
  peak-frame per-muscle re-derivation + axis projection (`decompose_at_peak`), the live moment-arm sign
  classification (`moment_arms_at_current_pose`, `verify_sign_anchor`), the CCI/roll-up/literature-check
  logic (`analyze_joint`), and the geometric pairwise-alignment check (`pairwise_most_antiparallel`) — all
  new orchestration around the existing, unedited geometry, per this repo's wrap-don't-edit-in-place
  convention.
- `data/msk_smoketest/subject2_walking1/contact_muscle_decomp/contact_muscle_decomp_results.json` — every
  number in this document except Sec.9, machine-written, including the full 13-row knee and 25-row hip
  per-muscle tables (this doc shows the top rows/summarized tail for the hip for brevity).
- Reused, not re-solved: `data/msk_smoketest/subject2_walking1/static_optimization/so/
  walking1_StaticOptimization_{activation,force}.sto`. Reused, not recomputed (Sec.9's real-EMG
  cross-check, via a small standalone read-only lookup, not `contact_muscle_decomp.py` itself):
  `data/msk_smoketest/subject2_walking1/emg_driven/emg_driven_results.json`,
  `.../emg_driven/emg_hybrid_force.sto`, `.../emg_driven/naive_fwd_only_force.sto` — all pre-existing
  outputs of a prior session's `scripts/msk/emg_driven.py`, read in place, not modified.
- **Not touched** (verified via md5, matches `docs/MECHANISM_CROSS_SUBJECT.md`'s own record for 3/4):
  `scripts/msk/validate_joint_force.py`, `scripts/msk/static_opt_knee.py`,
  `scripts/msk/validate_hip_force.py`, `scripts/msk/validate_emg_timing.py`.

No git commit, no git push performed (isolation respected). All new files are untracked, for the
coordinator to commit.
