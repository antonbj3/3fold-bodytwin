# MECHANISM NOCICEPTION — first mechanical-threshold nociceptor sensory layer (2026-07-21)

Adds **NOCICEPTION** as a sensory modality distinct from proprioception (which is being
built separately). `WAVE_PLAN.md` already names "(b) pain/nociception + thermoregulation"
as a planned wave item (line 257) alongside an "OpenPain brain + multimodal pain-signals"
pointer (line 245) — that CNS/brain-side pain-processing item is a **different, later**
piece of scope (out of bounds here); this build is the **peripheral mechanical-threshold
nociceptor** layer only: a sensory READOUT over loads this repo's MSK build already
computes and cert-grades, not a new physics model.

Script: `scripts/msk/nociception.py` (reads existing evidence + one live OpenSim
re-derivation, see §1). Data: `data/msk_smoketest/nociception/nociception_results.json`.
**Isolation respected:** read-only against the OpenSim install, the existing model/config
files, and every evidence JSON consumed; nothing pushed/committed; all new output is a
NEW file. Run: `source_repository/.venv-msk/bin/python3 scripts/msk/nociception.py`.

---

## 0. CRITICAL HONESTY FLAG — read this before trusting any >100° ligament result

**The single most attention-grabbing number this layer can produce — MCL right-knee
strain reaching 58.2% at 140° of passive knee flexion, an 11.3× overshoot of the
nociceptive threshold used below — is a KNOWN KINEMATIC ARTIFACT, not a physiological
finding.** This was already diagnosed, root-caused, and documented **before** this build,
in `docs/MECHANISM_KNEE_LIGAMENTS.md` §5/§8 and re-confirmed **fresh** (not cited) on the
merged whole-body model in `docs/MECHANISM_UNIFIED_MODEL.md` §3c: LaiArnold's
`walker_knee_r` is a **prescribed 1-DOF** `CustomJoint` (only `knee_angle_r` is a free
coordinate — ab-adduction/internal-external-rotation/all 3 translations are fixed
polynomial functions of flexion angle), while the Lenhart2015 model the ligament
attachment geometry was donor-sourced from uses a **free 6-DOF** knee whose secondary
kinematics *emerge* from ligament+contact force balance. Grafting real ligament geometry
onto a *different*, prescribed coupling has no guarantee of staying consistent at the
extremes of range of motion — it visibly doesn't, beyond roughly 100-110°.

This script **reproduces that number as a live cross-check** that the ligament build is
unchanged (58.21% at 140°, tension-only-corrected — see §1 — vs the source doc's
previously-published 58.2%: match) and uses it **only to illustrate this nociceptor
layer's threshold-crossing mechanism** — i.e., what the model looks like when fed a
genuinely enormous strain value. It is never reported as a real biological prediction.
Every row of the ligament channel's output carries a machine-readable
`kinematic_trust_zone` ∈ `{trusted, artifact_dominated}` (split at the source doc's own
stated 100° bound) and every event in the flat `events` list carries a boolean
`artifact_not_physiology`. The 80 events flagged `artifact_not_physiology=True` (knee
angle >100°, PCL+MCL) are illustrative-only; see §4 for the more nuanced honest treatment
of the *trusted*-zone crossings, which is not a clean "everything ≤100° is safe" story
either.

---

## 1. The model — a high-threshold mechanoreceptor (HTM) abstraction

**Geometric/physiological reasoning, not a rote lookup table.** A biological
mechanonociceptor is, to first order, a *high-threshold mechanoreceptor*: unlike a
proprioceptive low-threshold unit (fires throughout the normal working range), an HTM is
largely **silent** under ordinary physiological loading and only begins firing once the
local mechanical variable (strain, pressure, stress) crosses into a damage-relevant
range. Implemented here as the simplest faithful version of that structure, per tissue
channel:

- **threshold τ** — literature-anchored (§2), per tissue.
- **fires = (x > τ)** — the binary high-threshold readout.
- **drive = max(0, (x−τ)/τ)** — a graded relative-severity index for x>τ. This is an
  explicit **first-step simplification**: a LINEAR overshoot ratio, not a fitted
  spike-rate/recruitment curve. Real HTM units show a monotonically-increasing firing
  rate above threshold (general nociceptor physiology), but no quantitative
  firing-rate-vs-stimulus curve specific to ligament/cartilage/bone mechanonociceptors
  was found and verified within this session's scope (§6 gap).

**Four channels**, each a readout over a load this repo's MSK build already computed and
cert-graded — **no new forward-dynamics**:

| # | Channel | Tissue variable | Source |
|---|---|---|---|
| 1 | Ligament strain | per-bundle strain, ACL/PCL/MCL/LCL | live re-sweep of `model_with_knee_ligaments.osim` (see honesty flag below) |
| 2 | Cartilage contact pressure | TF/PF mean + peak (max-triangle) pressure | `data/msk_smoketest/cartilage_contact/results/joint-mechanics/passive_flexion.h5` |
| 3 | Bone (femoral cortical) stress | peak bending+axial stress | `bone_stress_results.json` |
| 4 | Joint contact force | knee/hip/ankle/spine peak %BW | `joint_force_validation`/`hip_force_validation`/`ankle_force_validation`/`spine_force_validation` JSONs |

### Channel-1 honesty flag — why a live re-sweep, not a JSON/doc-table read (a forced-adversary catch)

`knee_ligament_results.json`'s own persisted `flexion_sweep_right_knee` table stores
**group-summed force**, not strain. `docs/MECHANISM_KNEE_LIGAMENTS.md` §5's own printed
"max|strain|" table uses **abs(strain)** — and this was tested live, before trusting it:
at `knee_angle_r=0°`, that table's own "ACL max|strain|=0.140" is **not** the most-tensioned
bundle at all — it's the magnitude of a **slack** (negative-strain, zero-force) bundle;
the actual most-tensioned ACL bundle at 0° carries only +3.0% strain. A slack ligament
carries zero force (`Blankevoort1991Ligament`'s own tension-only law, already verified
0-violations-out-of-2982 in the source cert) and **cannot mechanically stimulate
anything** — naively reusing that abs-value table as a nociceptor input would have
manufactured a false signal from unloaded tissue. This script therefore **re-derives
signed, per-bundle strain directly** (reusing `add_knee_ligaments.load_config` /
`build_model_with_ligaments` for the geometry/scaling/slack-length math — zero
re-implementation — only the sweep-and-select loop, keyed on the SIGNED, max-**positive**
strain per group, is new) and applies the threshold only to `strain > 0` bundles.

Cross-check: at `knee_angle_r=90°`, this script's independent re-derivation gives
PCL=15.98%, MCL=25.70% — internally consistent with the source doc's own abs(strain)
table at the same angle (PCL 0.133 shown there is close because PCL happens to be
genuinely tensioned there; the two methods diverge exactly where a bundle is slack, as
designed).

---

## 2. Literature anchors — all live-verified this session (NCBI eutils), not recalled

Every PMID/DOI below was checked live this session via NCBI `esearch → esummary →
efetch` (abstract text), **never trusted from memory**. This mattered in practice: this
session's own **first-guess PMID for a joint-afferent paper (6827299) resolved live to a
completely unrelated 1983 visual-neuroscience monocular-deprivation paper** — caught
before it entered this document, the same failure mode `docs/MECHANISM_KNEE_LIGAMENTS.md`
§0 already found once with the task-given citations for the ligament build. Anchors are
graded a **tier**: PRIMARY = the specific number was directly quoted from live-fetched
PubMed abstract text this session; WEAKER = the citation/PMID is confirmed real and
correctly attributed, but the number is a widely-corroborated literature figure **not**
independently re-extracted from a live abstract this session (old paper, no abstract
indexed) — disclosed, not hidden, same convention as this repo's own anchor-tier legend
(`docs/MECHANISM_JOINT_FORCE_SCORECARD.md` §1).

| Anchor | Citation | PMID / DOI | Value used | Tier |
|---|---|---|---|---|
| Ligament damage onset | Provenzano PP, Heisey D, Hayashi K, Lakes R. "Subfailure damage in ligament: a structural and cellular evaluation." *J Appl Physiol.* 2002;92(1):362-71. | PMID 11744679, doi:10.1152/jappl.2002.92.1.362 | **5.14% strain** — "the onset of structural damage occurs at 5.14% strain (referenced from a preloaded length)"; "cellular damage was induced at ligament strains significantly below the structural-damage threshold" | **PRIMARY** |
| Cartilage damage stress | Torzilli PA, Grigiene R, Borrelli J Jr, Helfet DL. "Effect of impact load on articular cartilage: cell metabolism and viability, and matrix water content." *J Biomech Eng.* 1999;121(5):433-41. | PMID 10529909, doi:10.1115/1.2835070 | **15 MPa** (low end of the quoted 15-20 MPa band) — "a critical threshold stress (15-20 MPa) that caused cell death and apparent rupture of the collagen fiber matrix" | **PRIMARY** |
| Bone damage stress | Reilly DT, Burstein AH. "The elastic and ultimate properties of compact bone tissue." *J Biomech.* 1975;8(6):393-405. | PMID 1206042, doi:10.1016/0021-9290(75)90075-5 | **190 MPa** (commonly-cited ultimate compressive stress, human femoral cortical bone) | **WEAKER** — citation confirmed real live; no abstract indexed in PubMed for this pre-1976 paper (confirmed via `efetch`), so the number itself was not independently re-extracted from primary text this session (see §6) |
| Joint-afferent high-threshold property | Grigg P, Schaible HG, Schmidt RF. "Mechanical sensitivity of group III and IV afferents from posterior articular nerve in normal and inflamed cat knee." *J Neurophysiol.* 1986;55(4):635-43. | PMID 3701397, doi:10.1152/jn.1986.55.4.635 | *qualitative only* — in the **normal** (non-inflamed) knee, only 3/41 group IV and 12/18 group III units responded to local mechanical stimulation (vs. 14/36 and 24/36 respectively once inflamed) | **PRIMARY (qualitative)** — confirms the general HTM design principle; gives no absolute N/kPa number |

**What did NOT verify, tried and disclosed rather than papered over:** the task's own
suggested "~15-20% catastrophic ligament rupture strain" figure could not be
independently confirmed from a live-fetched abstract. Noyes FR & Grood ES 1976 ("The
strength of the anterior cruciate ligament in humans and Rhesus monkeys," *J Bone Joint
Surg Am*, PMID 1002748 — confirmed real, no DOI on record) gives only a **relative**
age comparison (young-adult specimens ~2-3× stronger than sixth-decade+ specimens), no
absolute strain percentage. Woo SL, Hollis JM, Adams DJ, Lyon RM, Takai S 1991
("Tensile properties of the human femur-ACL-tibia complex," *Am J Sports Med*, PMID
1867330, doi:10.1177/036354659101900303 — confirmed real) gives ultimate **load**
(2160±157 N) and stiffness (242±28 N/mm) for young specimens in anatomical orientation —
a real, force-domain failure anchor, but not a strain percentage, so not used as the
operative threshold. The 5.14% Provenzano figure is used instead: it is **more
conservative** (fires earlier, protective-direction) than the task's suggested range, not
a weaker substitute for it.

---

## 3. Results by channel

### 3a. Ligament strain (the channel that actually fires)

Full 0-140° sweep, 1° resolution, right knee, tension-only-corrected signed strain.
**First-crossing angle per ligament group** (τ=5.14%):

| Group | First fires (full range) | First fires (trusted zone, ≤100°) | Peak (140°) |
|---|---:|---:|---:|
| ACL | never | never | 2.94% (never exceeds τ anywhere in the sweep) |
| PCL | **54°** | **54°** | 23.02% |
| MCL | **56°** | **56°** | **58.21%** (the flagged artifact case, §0) |
| LCL | **0°** (fires 0-3°, then goes slack) | **0°** | n/a (slack from ~4° onward) |

Total ligament-channel firing events across the full sweep: **176** (one per
angle×group crossing) — **96 in the trusted zone (≤100°), 80 in the artifact-dominated
zone (>100°, §0)**. Representative rows (all machine-computed, cross-checked against an
independent manual probe at 90° before this script was trusted):

| angle | ACL | PCL | MCL | LCL |
|---:|---:|---:|---:|---:|
| 0° | +3.00% | +3.00% | +5.00% | **+6.00% (fires)** |
| 40° | +0.13% | +0.72% | +2.05% | −0.16% (slack) |
| 54° | −1.30% | **+5.23% (fires)** | +4.29% | −1.90% (slack) |
| 56° | −1.28% | +5.88% | **+5.16% (fires)** | −2.13% (slack) |
| 90° | −2.27% | +15.98% | +25.70% | −6.21% (slack) |
| 100° | −2.10% | +18.23% | +32.89% | −7.28% (slack) |
| 140° | +2.94% | +23.02% | **+58.21% (§0 artifact)** | −6.15% (slack) |

**Honest nuance on the "trusted-zone" crossings themselves (PCL 54°, MCL 56°, and
everything through 100°):** these are **not** re-flagged as artifact (they sit inside
`docs/MECHANISM_KNEE_LIGAMENTS.md` §5's own stated ≤100° bound), but that bound was
calibrated for a **qualitative** check (which ligament dominates at which flexion angle —
the doc's own pre-registered PASS/PASS engagement-pattern thresholds), not for
absolute-magnitude injury-threshold comparison. The same 1-DOF-vs-6-DOF coupling
mismatch that is *dominant* beyond 100° plausibly onsets **gradually**, not as a sharp
switch exactly at 100°. This build has **no decorrelated way** to separate genuine
passive-ligament mechanics from a smaller, unquantified dose of the same coupling
artifact in the 60-100° band — this is reported as an **open, undetermined question**,
not resolved in either direction. The **one crossing furthest from any artifact
concern** is LCL at 0-3° (the model's own reference/neutral standing pose): all 4 LCL
bundles are pinned at the donor model's own published +6% reference strain (Lenhart2015,
MRI-validated), a real, non-kinematically-contaminated value that happens to sit just
above (6.00% vs 5.14%) the conservative damage-onset threshold — a small, unalarming
overshoot at rest, not a warning sign, but the cleanest non-artifact data point this
build has.

### 3b. Cartilage contact pressure — silent (a non-degenerate negative)

Passive flexion trial (0→91.04°, no GRF/muscle activity), settled window, τ=15 MPa:

| Contact | Peak mean pressure | Peak max (single-triangle) pressure | Fires? |
|---|---:|---:|---|
| tibiofemoral | 2.10 MPa @ 91.0° flexion | 4.76 MPa | **No** (7.1× / 3.2× below τ) |
| patellofemoral | 2.29 MPa @ 86.2° flexion | 6.13 MPa | **No** (6.5× / 2.4× below τ) |

Correctly silent under a passive, non-injurious trial — the cartilage-damage threshold
(15-20 MPa, impact-loading literature) sits well above even the peak single-triangle
pressure this passive trial ever produces.

### 3c. Bone (femoral cortical) stress — silent

Subject2 walking1 gait trial, Euler-Bernoulli beam peak stress (already cert-graded in
`docs/MECHANISM_JOINT_FORCE_SCORECARD.md`'s sibling bone-stress cert): **54.97 MPa**
(compressive, medial subtrochanteric/proximal shaft, t=0.56s = 35.7% gait cycle) vs
τ=190 MPa → **fires=False**, ratio 0.289. Silent during ordinary gait, as expected.

### 3d. Joint contact force — diagnostic only, not receptor-anchored

No literature mechanonociceptor threshold in absolute Newtons/kPa/Nm was found and
verified this session (§2's qualitative afferent anchor gives no such number). This
channel is reported as **context only**, using a CONSTRUCTED (not receptor-anchored) "2×
published in-vivo peak" implausibility ceiling — the same style of ceiling this repo's
own scorecard already uses to flag numbers "far outside any published literature" as a
pipeline red flag (`docs/MECHANISM_JOINT_FORCE_SCORECARD.md` §3 item 2), **not** a
nociception threshold:

| Joint | Twin peak | In-vivo anchor | Constructed 2× ceiling | Exceeds ceiling? |
|---|---:|---:|---:|---|
| Knee | 391.1 %BW (t=0.51s) | 258.2 %BW (OrthoLoad, n=72) | 516.4 %BW | No |
| Hip | 387.0 %BW (t=0.55s) | 273.9 %BW (OrthoLoad, n=162) | 547.9 %BW | No |
| Ankle | 484.5 %BW (t=0.58s) | 476.7 %BW (cadaveric/FE, not in-vivo) | 953.3 %BW | No |
| Spine | 389.2 %BW (with box, static pose) | 170.9 %BW (OrthoLoad, partial-load-sharing anchor) | 341.9 %BW | **Yes** |

Spine crossing its own constructed ceiling is **not a new finding** — it is the same
joint the scorecard already flags as having the family's highest over-prediction ratio
(1.74-2.28×, vs 1.4-1.5× at knee/hip), for reasons already disclosed there (partial-
load-sharing anchor, possibly stacked with the same Moissenet-mechanism SO
over-prediction the scorecard names). Reported here as an internal-consistency check
(this channel correctly flags the one joint the repo's own independent scorecard also
flags), not as an independent nociception claim.

---

## 4. Where/when — the flat, machine-checkable events answer

`nociception_results.json["events"]` is the direct answer to "report where/when it
predicts nociceptive signaling": one entry per firing, with `channel`, `location`,
`when`, `magnitude`, `threshold`, `ratio_to_threshold`, `kinematic_trust_zone`, and
`artifact_not_physiology`. Summary: **176 ligament-strain events** (96 non-artifact /
80 artifact-flagged, §0/§3a), **0 cartilage events, 0 bone events, 0 receptor-anchored
joint-force events** (channel 4 has no receptor threshold to fire against by design,
§3d). Under every load this twin has certified for **normal, non-injurious movement**
(gait, passive flexion), the cartilage and bone channels stay silent — the nociceptor
layer does not cry wolf under benign conditions, a real (if unexciting) specificity
check. The only channel that fires is ligament strain, starting well before the
already-known 140° artifact case, with the honest ambiguity of §3a disclosed rather than
resolved.

---

## 5. Reproduction

```
source_repository/.venv-msk/bin/python3 scripts/msk/nociception.py
```

Rebuilds the ligament model in-memory (read-only against the source `.osim`, no file
written except this run's own output), re-sweeps signed per-bundle strain 0-140° at 1°
resolution, reads the persisted cartilage-contact h5 and the bone-stress/joint-force
JSONs directly, applies the four thresholds, and writes
`data/msk_smoketest/nociception/nociception_results.json` + a printed summary. No
network access, no git operations.

---

## 6. Honest gaps (full list)

1. **Channel 4 (joint contact force) has no receptor-threshold anchor** — diagnostic
   context only (§3d); an absolute-force mechanonociceptor threshold was searched for
   and not found/verified live this session.
2. **Bone-damage threshold is WEAKER-tier** (§2) — the citation is real
   (Reilly & Burstein 1975) but its ~190 MPa figure could not be independently
   re-extracted from a live PubMed abstract (none indexed for this pre-1976 paper); two
   follow-up searches (Carter & Hayes 1977, PMID 561786; Turner & Burr 1993, PMID
   8274302 — both confirmed real) also did not yield a live-quotable absolute number.
3. **The task's suggested ~15-20% "catastrophic" ligament-rupture strain was not
   independently verified live** (§2) — 5.14% (subfailure-damage-onset, Provenzano 2002)
   is used instead: more conservative, not a substitute for the same claim.
4. **The 60-100° ligament-channel crossings carry an open, undetermined artifact-vs-
   physiology question** (§3a) — not resolved either way, disclosed as ambiguous rather
   than silently classed as either "safe" or "artifact."
5. **Ligament channel inherits every honesty flag already logged against the source
   build** (`docs/MECHANISM_KNEE_LIGAMENTS.md` §8): literature-scaled not
   subject-specific attachment geometry, `linear_stiffness` left unscaled, single donor
   subject/model family.
6. **Cartilage channel uses a passive (no-GRF, no active muscle) trial only** — a
   loaded-gait cartilage-pressure trial was never built in this repo
   (`docs/MECHANISM_CARTILAGE_CONTACT.md` §7 gap 1); this script cannot and does not
   report a loaded-gait cartilage-nociception result. Loaded gait would very plausibly
   raise the cartilage channel's peak pressure closer to (though, per that cert's own
   literature anchor, still likely below) the 15 MPa damage threshold.
7. **Graded "drive" is a linear overshoot ratio**, an explicit first-step simplification
   — no fitted firing-rate/recruitment curve for connective-tissue mechanonociceptors
   was found and verified within this session's scope.
8. **No sensitization/inflammation model** — the Grigg/Schaible/Schmidt 1986 anchor
   itself shows thresholds drop sharply once inflamed (14/36 vs 3/41 group IV units
   responsive); this build models only the naive/uninflamed threshold.
9. **Static thresholds only** — no habituation, adaptation, or loading-rate dependence
   (Torzilli's own cartilage finding, and connective-tissue viscoelasticity generally,
   both suggest rate-dependence; not modeled here).
10. **Proprioception (the sibling modality) is out of scope by design** — built
    separately per the task; no shared-afferent or cross-modality interaction is
    modeled here.
11. **Single subject, single trial family** — inherits subject2/walking1's own scope
    limits (`docs/MECHANISM_JOINT_FORCE_SCORECARD.md` §7): one healthy young-ish subject,
    right side primary, no claim of generality across subjects/trials/populations.

---

## 7. Files

- `scripts/msk/nociception.py` — the model + all 4 channels + literature anchors +
  machine cross-checks + JSON evidence writer. Reuses `add_knee_ligaments` (ligament
  geometry/scaling, unmodified) and `cartilage_contact` (h5 path/loader, unmodified) —
  zero re-implementation of either's proven physics.
- `data/msk_smoketest/nociception/nociception_results.json` — full machine-readable
  evidence: anchors, thresholds, all 4 channels' raw series/summaries, the flat `events`
  list, and `honest_gaps`.
- Read for parent/dependency verification (not re-litigated):
  `docs/MECHANISM_KNEE_LIGAMENTS.md` (ligament build + the artifact diagnosis, §5/§8),
  `docs/MECHANISM_UNIFIED_MODEL.md` (§3c, fresh reconfirmation of the same artifact),
  `docs/MECHANISM_CARTILAGE_CONTACT.md` (contact-pressure build + its own literature
  anchor), the bone-stress cert behind `bone_stress_results.json`,
  `docs/MECHANISM_JOINT_FORCE_SCORECARD.md` (joint-force numbers + its own anchor-tier
  legend, reused here), `WAVE_PLAN.md` (lines 245/257, the pre-existing pain/nociception
  wave-plan pointer this build fulfills).
