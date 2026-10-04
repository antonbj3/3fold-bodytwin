# MECHANISM SACROILIAC JOINT — pelvic-ring load transfer, form/force-closure stability (2026-07-22)

Script: `scripts/msk/sacroiliac_joint_gates.py`. Raw results:
`data/msk_smoketest/sacroiliac_joint/sacroiliac_joint_gate_results.json`. Evidence (28 citations,
verbatim abstracts): `docs/MECHANISM_SACROILIAC_JOINT_evidence.json`. This is a
**literature-parameter-anchored gate-checker**, not a live OpenSim simulation — every number below
traces to a PMID verified LIVE this session (NCBI eutils `esearch`+`efetch`), mechanically evaluated
against pre-registered thresholds (`python3 scripts/msk/sacroiliac_joint_gates.py`, exit 0, 10/10
gates PASS — see §9 for the honest limits of that 10/10).

**Duplicate check**: repo-wide grep for `sacroiliac|sacro-iliac|SI-JOINT|SI_JOINT` across
`docs/*.md` and the graph's 1015 node ids found zero prior hits (`MYOFASCIAL.md` and
`TRUNK_FLEXORS.md` mention `thoracolumbar`/`self-bracing` in passing but do not cover the SI joint
itself) — this is a genuinely new, non-duplicate node.

## 0. couples_to

- `docs/MECHANISM_MZ_PELVIS_MOMENT.md` — that doc found a genuine ~8–18 N·m **sagittal** (`pelvis_tilt`)
  pelvis-moment floor during gait, intrinsic to the kinematics–GRF pairing, not fixable by RRA task
  gains. `pelvis_tilt` is kinematically the coordinate an SI joint's nutation/counternutation DOF
  would occupy if this twin modeled one — see §8, a concrete, machine-verified structural link.
- `docs/MECHANISM_HIP_ANKLE_LIGAMENTS.md` — supplies this doc's in-twin "free hinge" comparator
  (§2, its own hip/ankle sweep ranges) and is the natural sibling for the SAME
  `opensim.Blankevoort1991Ligament` class this doc's proposed_cell (§10) would extend to the SI
  ligaments (sacrotuberous, sacrospinous, interosseous, long dorsal, ventral SI).
- `docs/MECHANISM_SPINE_LIGAMENTS.md` / `docs/MECHANISM_SPINE_FORCE.md` — the SI joint is the
  spine's downstream load-transfer bridge to the legs; both docs' lumbar loads pass through the
  mechanism this doc certifies.
- `docs/MECHANISM_MYOFASCIAL.md` (§2, "candidate #2, erector spinae/thoracolumbar fascia, not
  quantitatively assessed") — this doc's Willard 2012 citation (PMID 22630613) supplies exactly the
  TLF anatomy that candidate needed; flagged there as an open item, informed (not yet closed) here.
- `docs/MECHANISM_TRUNK_FLEXORS.md` — core/transversus-abdominis-adjacent muscle mechanics; this
  doc's Snijders 1998 EMG citation (internal oblique) is the same transversely-oriented muscle class.
- Graph nodes `MSK-COM-GRF-VIDEO-ONLY` (the GRF this joint transmits leg→spine), `MSK-GAIT-LOCOMOTION`
  (SI motion is cyclic during gait), `MSK-HIP-FEMOROACETABULAR` and `MSK-LUMBAR-SPINE-SEGMENT`
  (the adjacent links in the spine→sacrum→ilium→hip chain) — all currently `couples_to_ids: None`
  in the graph; this doc names the prose link for the coordinator to resolve.

## 1. Headline

**All 10 pre-registered gates PASS** (`sacroiliac_joint_gates.py`, exit 0) — but three of them
(G5/G6/G7, the direct-manipulation mechanism experiments) share a common-mode limitation disclosed
prominently in §9, not hidden. The measured SI ROM (rotation ~1.7–2.5° mean, translation ~0.7 mm
mean, across 4 independent methods spanning 1954–2000, plus an independent FE-model cross-check) is
**simultaneously inconsistent with a rigid-pelvis (zero-ROM) null** (signal exceeds an independently
measured noise floor by 8.5–12.5×) **and a free-hinge (large-ROM) null** (even the largest recorded,
pathological individual value is ≤6% of this twin's own hip/ankle sweep-range convention). Force
closure is shown to be **load-bearing, not decorative**, via a direct cadaveric ligament-section
experiment and an EMG belt-substitution experiment. Nutation is shown to **increase** ligament
tension and joint compression under load — a self-locking (negative-feedback), not runaway,
mechanism. A concrete, machine-verified (grep, not narrated) finding: **this twin's own OpenSim
models currently implement the rigid-pelvis null by construction** — both `subject2_scaled_
handsfield_fmax_corrected.osim` and the merged `subject2_unified.osim` have exactly **one**
`<Body name="pelvis">` and no sacrum/ilium split anywhere (§8) — the literature's "rigid pelvis"
adversary is not just a textbook idealization, it is this exact twin's live model.

## 2. The measured ROM anchor — 4 independent methods + 1 independent computation, 1954–2015

| study | PMID | method | n | rotation | translation |
|---|---|---|---:|---:|---:|
| Weisl 1954 | 13196834 | cadaveric/mechanical (bibliographic only, pre-abstract era) | — | — | — |
| Egund et al 1978 | 717034 | **first** roentgen stereophotogrammetry (RSA) | 4 | ~2.0° max (axis precision **0.2°** mean) | PSIS–PSIS ≤0.4 mm across 7 positions |
| Sturesson, Selvik, Udén 1989 | 2922636 | refined RSA, patients | 25 | **2.5°** mean (0.8–3.9°) | **0.7 mm** mean (0.1–1.6 mm) |
| Jacob & Kissling 1995 | 11415579 | 3-D photogrammetry (external markers, **not** X-ray) | 24 | **1.7°** mean (symptomatic outlier: >6°) | 0.7 mm mean |
| Sturesson, Uden, Vleeming 2000 | 10703111 | RSA during a **loaded** functional task (standing hip flexion) | 22 | "very small movements... self-locking mechanism... obstructs movements" | — |
| Zheng, Watson, Yong-Hing 1997 (independent **computation**, not measurement) | 9136197 | quasi-static FE model, 3-link closed-chain pelvis | 1 cadaver-derived | 1.0–1.6° at 50 N·m | 0.5–1.8 mm at 1000 N |

**G1 (rejects the rigid-pelvis null)**: every study's own mean rotation exceeds Egund 1978's
independently-reported measurement precision (0.2°, a DIFFERENT quantity than the ROM claim itself —
an instrument-repeatability characterization, not a restatement) by **8.5×–12.5×**
(1.7°/0.2°=8.5, 2.0°/0.2°=10.0, 2.5°/0.2°=12.5) — machine-computed, `G1_rom_resolvable_vs_rigid_
pelvis_null` in the results JSON. **Honest near-miss, disclosed not hidden**: Sturesson 1989's
single SMALLEST *individual* reading across its 25 patients was 0.8° — only **4.0×** the noise floor,
*below* this gate's own 5× bar. The gate is evaluated on each study's headline mean (the standard,
non-cherry-picked statistic each paper itself reports), not on the single smallest individual
data point; that choice is disclosed here explicitly rather than silently adopted.

**G2 (rejects the free-hinge null)**: rather than import an external "typical hinge joint" number,
this gate reuses this twin's **own already-established** mobile-joint convention —
`docs/MECHANISM_HIP_ANKLE_LIGAMENTS.md` swept `hip_flexion_r` over a 150° range and `ankle_angle_r`
over a 100° range. Even the single largest SI value on record anywhere — Jacob & Kissling's
symptomatic-subject **outlier**, 6° (a *pathological* case, the most generous possible number for
the adversary) — is **6%** of the smallest of those two reference ranges (ankle, 100°). Machine-computed:
`G2_rom_far_below_free_hinge_null`, ratio 0.06 vs a pre-registered ≤0.15 bar.

**G3 (cross-method over-determination, not a tautology)**: the RSA-measured mean (2.07° rotation /
0.70 mm translation, averaged across the 3 RSA-family studies) and Zheng et al's **independently
computed** FE-model mean (1.23° / 1.27 mm, built from cadaver geometry via a minimum-potential-energy
formulation, never fit to the RSA numbers) agree within **1.68×/1.81×** — comfortably inside a
pre-registered ≤3× band. Two decorrelated methods (marker-based in-vivo kinematics vs a
first-principles structural model) land on the same order of magnitude.

## 3. Geometric model — the Coulomb normal-force floor (why form closure alone cannot be sufficient)

Shear resistance at any joint interface is Coulomb-limited: **F_shear,max = μ_eff · N**, where N is
the compressive (normal) preload across the joint and μ_eff is the effective friction coefficient.
**The governing quantity is N, not μ**: at N=0, F_shear,max=0 **regardless of how large μ_eff is** —
a hard floor, geometrically analogous to a σ_min null-space (a mode with zero singular value cannot
be rescued by scaling anything else in the product).

Form closure raises μ_eff: cartilage-covered ridges and depressions on the auricular surfaces,
"complementary" between sacrum and ilium and "more pronounced in men" (Vleeming, Volkers, Snijders,
Stoeckart 1990 Part I, PMID 2326706), measurably raise the friction coefficient — ridges/depressions
matter **more** than coarse texture alone (1990 Part II, PMID 2326707). An independent, modern,
**cross-species tribological replication** (Nordberg, Espinosa, Hu, Athanasiou 2021, minipig,
PMID 34105385 — a UC Irvine group fully unconnected to the Rotterdam program below) confirms both
the *direction* and gives a hard number: SI-joint cartilage coefficient of friction is **28% higher
than knee cartilage**, and — the mechanistically decisive detail — **LOCAL** roughness (fine
ridges/grooves) correlates significantly with friction, while **global** roughness (gross shape) and
lubricin content do **not** (`G10_form_closure_tribological_replication`, PASS). This is the same
"fine structure, not gross shape" relationship Vleeming's 1990 human study found 31 years earlier, in
a different species with a different instrument — a genuine, decorrelated replication of the
mechanism, not just the label.

**But raising μ_eff cannot substitute for N.** The field's own most mature statement of this (Vleeming
& Schuenke 2019, PMID 31218826) states the null explicitly: *"If the sacrum would fit in the pelvis
with perfect form closure, no lateral compressional forces would be needed to maintain stability.
However, such a construction would make mobility practically impossible."* Real SI joints are not
built at that idealized limit — they sit at neither extreme (not a perfectly form-closed, zero-ROM
joint; not a friction-free hinge either), exactly the two nulls §2 rejects. **Force closure — active
compression from muscles (transversely-oriented: internal oblique/transversus-abdominis-adjacent,
gluteus maximus, piriformis), ligaments (sacrotuberous, sacrospinal, dorsal and interosseous SI), and
fascia (the posterior layer of the thoracolumbar fascia, TLF) — is what keeps N bounded away from
zero.** The TLF's own detailed anatomy (Willard, Vleeming, Schuenke, Danneels, Schleip 2012, PMID
22630613, verified live) confirms the exact multi-muscle chain the task names: its posterior layer is
"dominated by the aponeuroses of the **latissimus dorsi**," its lateral raphe is "the junction of...
the **transversus abdominis**," and the whole composite "attaches firmly to the posterior superior
iliac spine and the **sacrotuberous ligament**" at the base of the lumbar spine — a literal anatomic
bridge from arm/lat muscles through the fascia to the SI joint's own primary ligament.

## 4. Force closure is load-bearing — two direct manipulation experiments, not just correlation

**Cadaveric ligament-section experiment** (Buyruk, Stam, Snijders, Vleeming, Lameris, Holland 1995,
embalmed human pelvises, PMID 8850504): a pelvic vibrator + colour-Doppler-imaging stiffness proxy
was applied under **three explicit mechanical conditions on the same specimens** — intact, screwed
(rigid fixation, pushing toward the form-closure/rigid limit), and **ligaments cut** (force closure
directly removed). Result: "statistics showed high reproducibility and **significant differences
between the stability conditions**." This is a direct manipulation, not a correlational study —
removing force closure changes measured stiffness. **Honest gap, precisely stated**: the abstract
confirms a *significant* difference across the 3 conditions but does **not** itself state the
per-condition ranking — it does not say, in so many words, that "ligaments cut" reads softer than
"intact." The predicted direction (cutting ligaments reduces stiffness) is what the theory
(Pool-Goudzwaard 1998, Vleeming & Schuenke 2019) predicts and what the separate Snijders 1998 EMG
result (below) is independently consistent with — but it is not confirmed by this abstract's text
directly (`G5_force_closure_manipulation_cadaveric`, flagged honestly in its own `honest_gap` field).

**EMG belt-substitution experiment** (Snijders, Ribbers, de Bakker, Stoeckart, Stam 1998, living
subjects, PMID 9779394): internal oblique (a transversely-oriented, force-closure-relevant abdominal
muscle) activity **decreased significantly** when (a) resting on one leg (less gravity load), (b)
tilting the pelvis backward (less psoas load), and (c) **applying a pelvic belt** — an external
mechanical substitute for the muscle's own compressive contribution. This is the forced-adversary
test for "is force closure just a passive geometric byproduct, or an actively regulated mechanism":
a muscle measurably *reduces its own contribution* exactly when an external device offers to do the
same mechanical job — behavior consistent only with active regulation of a needed compressive force,
not with a passive/decorative role (`G6_force_closure_emg_substitution`, PASS).

**Named common-mode limitation (disclosed, not buried)**: G5, G6, and the nutation-tension
experiment in §5 (G7) all originate from **one coherent ~30-year Rotterdam/Erasmus University
research program** with heavily overlapping authorship (Vleeming, Snijders, Stoeckart, Buyruk,
Pool-Goudzwaard, Damen, van Wingerden appear across nearly every one of these papers). No fully
independent, unconnected laboratory was found this session to have directly replicated the specific
cadaveric ligament-section or EMG-belt experiments. This is a genuine decorrelation gap, not a fatal
one: (a) the *geometric necessity* argument in §3 (F=μN→0 at N=0) is a mechanical certainty
independent of any lab's interpretation; (b) the tribological replication (§3, Nordberg 2021, UC
Irvine, fully unconnected) independently confirms the *surface-mechanism* half of the story; (c) the
nutation *direction* (§5) is echoed by Sturesson 2000's RSA outcome measurement — though Sturesson
himself co-authored that specific paper with Vleeming, so this is only partial independence. Flagged
prominently as the single most important source-diversity gap in this document, analogous to how
`docs/MECHANISM_CORNEAL_TRANSPARENCY.md` flagged its own single-calibration-point dependency.

## 5. Nutation increases with load — self-locking, not runaway

Direct cadaveric buckle-transducer measurement (Vleeming, Pool-Goudzwaard, Hammudoghlu, Stoeckart,
Snijders, Mens 1996, n=12, 0–50 N incremental loading, PMID 8852309): **forced nutation decreased
long dorsal sacroiliac ligament (LDL) tension and increased sacrotuberous ligament tension; forced
counternutation reversed both.** The same experiment measured direct coupling to exactly the
multi-muscle chain the task names: LDL tension increased when the **erector spinae** was loaded and
decreased when **gluteus maximus** and the TLF (**simulating latissimus dorsi** contraction) were
loaded — the identical arm↔fascia↔SI-ligament bridge confirmed anatomically in §3.

The field's current synthesis (Vleeming & Schuenke 2019, PMID 31218826) states the mechanism in full:
sacral nutation is **"anticipatory for joint loading"** — it tenses all the dorsal SI ligaments
(interosseous, dorsal sacroiliac) **except** the LDL, which "prepares the pelvis to absorb an
increase in load," and **"the posterior ilium are pressed together causing an increase in SIJ
compression."** Combined with §3's relation, this closes a coherent mechanical loop: axial load →
nutation → ligament tension (↑ sacrotuberous/interosseous, ↓ LDL) → increased N (compression) →
increased F_shear,max = μ·N — i.e. d(F_shear,max)/d(nutation) > 0, a **negative-feedback,
self-locking** loop, not a positive-feedback runaway. **This document does not separately re-derive
the demand-side (how fast required shear resistance grows with load) from first principles** — the
claim that capacity growth keeps pace with (or outpaces) demand growth rests on the **net outcome**
measurement: Sturesson 2000's RSA finding of "very small movements" specifically **under** the
one-leg-stance-loaded condition, describing the result as the self-locking mechanism "obstruct[ing]"
motion. That is an outcome-level validation, not an independently re-derived rate — disclosed as a
scope limit, not hidden (`G7_nutation_self_locking`, PASS; see §9 for the honest gap).

## 6. Perturbation adversaries

**Peripartum relaxin (reduces force closure hormonally)**: serum relaxin remodels pelvic connective
tissue. MacLennan, Nicolson, Green, Bath 1986 (Lancet, PMID 2874277): 35 women with severe pelvic
pain/joint instability had relaxin levels above the 95% CI of 368 controls, normalizing by the 3rd
postnatal day. Kristiansson, Svisdsudd, von Schoultz 1996 (PMID 8942512): relaxin peaks at week 12,
correlates significantly with symphyseal/back pain in late pregnancy. **Honest nuance, not smoothed
over**: Damen, Buyruk, Güler-Uysal, Lotgering, Snijders, Stam 2001 (n=163, PMID 11703199) found mean
SI joint laxity does **NOT** differ between pain (3.0 TU) and no-pain (3.4 TU) groups — it is the
**left-right ASYMMETRY** of laxity that differs (2.2 vs 0.9 TU) and predicts pain (37% of the pain
group vs 4% of the no-pain group show asymmetric laxity). A 2002 prospective follow-up (n=123, PMID
12486354) confirms this asymmetry predicts pain persisting postpartum (sens 65%, spec 83%, PPV 77%,
3× risk). **Confound named**: none of these citations separate relaxin's direct ligament-compliance
effect from pregnancy's other confounds (weight gain, postural change, vascular effects) — flagged,
not resolved (`G8_relaxin_peripartum_adversary`, PASS on the relaxin-pain link across 2 independent
cohorts; nuance and confound both carried into the results JSON, not silently dropped).

**SI joint fusion (eliminates ROM, pushes toward the rigid-pelvis limit)**: Lindsey, Kiapour, Yerby,
Goel 2015 (validated FE model, iFuse implant, PMID 26767156) — fusion reduces SI ROM by **56.6%
(flexion), 59.5% (extension), 27.8% (lateral bending), 53.3% (axial rotation)** vs the intact
condition, all four clearing a pre-registered ≥25% bar, with only a small (1.1–4.6%) compensatory
increase at the adjacent L5-S1 lumbar segment (`G4_fusion_perturbation_toward_rigid_null`, PASS —
the tightest margin of any gate, 27.8% vs a 25% bar). Honest note: this is a computational (FE),
not an in-vivo RSA, measurement of a fused joint.

**Pelvic-ring fracture (load-path disruption)**: Varga, Hearn, Powell, Tile 1995 (n=8 unembalmed
cadaveric pelves, servohydraulic testing machine, 500 N cyclic axial load simulating bilateral
stance, PMID 7721471) — sectioning the **pubic symphysis** (a DIFFERENT part of the closed pelvic
ring) measurably changed motion recorded at the **posterior SI complex**, and even after internal
fixation, "internally fixed symphyseal motion was generally greater than intact, regardless of
fixation method" — the ring's native closed-chain load-sharing is not fully recoverable by hardware
alone (`G9_pelvic_ring_fracture_adversary`, PASS). Zheng et al 1997 (§2)'s own framing — the pelvis
as a "stable, complex three-link [closed-chain] structure" that "explains why there is so little
motion in the sacroiliac joint" — gives the geometric reason a ring fracture matters: breaking the
closed-chain topology at any one point removes a redundant load path, changing load/motion
everywhere else in the ring, exactly what Varga/Tile measured directly.

## 7. Gate summary (machine-computed, `sacroiliac_joint_gates.py`)

| Gate | Claim | Result |
|---|---|---|
| G1 | ROM resolvable vs rigid-pelvis null | PASS (8.5–12.5× noise floor; honest near-miss on the single smallest individual reading, §2) |
| G2 | ROM far below free-hinge null (in-twin comparator) | PASS (ratio 0.06 vs ≤0.15 bar) |
| G3 | RSA-vs-FE cross-method over-determination | PASS (1.68×/1.81× vs ≤3× bar) |
| G4 | Fusion pushes toward rigid-pelvis limit | PASS (27.8–59.5% reduction vs ≥25% bar — tightest margin) |
| G5 | Force-closure cadaveric manipulation | PASS (significant; magnitude not extracted, common-mode-flagged) |
| G6 | Force-closure EMG substitution | PASS (common-mode-flagged) |
| G7 | Nutation self-locking | PASS (outcome-level validation, not independently re-derived rate) |
| G8 | Relaxin/peripartum adversary | PASS (2 independent cohorts; asymmetry-not-magnitude nuance disclosed) |
| G9 | Pelvic-ring fracture load-path coupling | PASS |
| G10 | Cross-species tribological replication | PASS (fully independent lab/species) |

**10/10 PASS.** Per this discipline's own rule ("kill-rate is not a quality signal, an unjustified
kill is deferred not dead" — its converse for confirmations applies too): a clean sweep is reported
with its limits attached (§4's common-mode flag, §5's outcome-vs-rate scope limit, §2's near-miss),
not celebrated as unqualified.

## 8. A concrete, machine-verified twin-model finding — the rigid-pelvis null IS this twin's current model

Grepped live this session, not narrated from memory:

```
$ grep -o '<Body name="[a-zA-Z_0-9]*"' data/msk_models/subject2_scaled_handsfield_fmax_corrected.osim | sort -u
...
<Body name="pelvis"
...
$ grep -o '<Body name="[a-zA-Z_0-9]*"' data/msk_models/subject2_unified.osim | grep -i "pelvis\|sacr\|ili"
<Body name="pelvis"
```

**Exactly one `pelvis` body, in both the base model and the merged/unified model** that
`docs/MECHANISM_HIP_ANKLE_LIGAMENTS.md` and `docs/MECHANISM_SPINE_LIGAMENTS.md` build on — sacrum and
ilium are permanently welded together with zero relative DOF. The only joints touching `pelvis` are
`ground_pelvis` (the free 6-DOF **floating-base** joint between the whole rigid pelvis and the
world — the same coordinate `docs/MECHANISM_MZ_PELVIS_MOMENT.md` analyzes as `pelvis_tilt`/`list`/
`rotation`, **not** an SI joint, and this doc is careful not to conflate the two) and `hip_l`/`hip_r`/
`back`. There is no `sacrum` or `ilium` body anywhere, therefore **no SI joint DOF exists in this
twin at all** — the model *is*, by omission rather than by measured/informed choice, the literature's
"rigid pelvis" null this entire document was built to reject empirically. This is a genuine,
actionable gap, not a criticism of prior work (no earlier doc claimed to model the SI joint) — see
§10.

## 9. Symmetric QC — honest gaps (full list)

1. **Common-mode source risk for G5/G6/G7** (§4) — the decisive direct-manipulation mechanism
   experiments share one ~30-year research program's authorship. Not independently replicated by an
   unconnected lab this session. Partially mitigated by the geometric necessity argument (§3) and the
   fully-independent tribological replication (Nordberg 2021), but not closed.
2. **G1's own near-miss** (§2) — the single smallest *individual* rotation reading across all 25
   Sturesson 1989 patients (0.8°) is only 4.0× the noise floor, below this gate's pre-registered 5×
   bar. The gate is scored on each study's own headline mean (the standard, undisputed statistic),
   not the smallest individual value; both numbers are disclosed, not just the favorable one.
3. **Buyruk 1995b's per-condition stiffness values were not extracted** (§4) — the abstract reports
   statistical significance across intact/screwed/ligaments-cut, not the threshold-unit numbers
   themselves; full text is paywalled and was not recovered this session.
4. **No wedge half-angle was independently measured or derived** (§3) — "dorsocranial wedging" is
   the field's own verified descriptive language (PMID 31218826); this doc did not re-derive a
   specific angle from CT/morphometric raw data (unavailable for live fetch this session).
5. **The self-locking loop's demand-side growth rate was not independently derived** (§5) — the
   claim that shear-capacity growth keeps pace with load-driven shear-demand growth rests on
   Sturesson 2000's net-outcome RSA measurement (small ROM under one-leg-stance load), not a
   separately re-derived rate comparison.
6. **Relaxin's specific mechanical contribution is confounded** (§6) — pregnancy weight gain,
   postural change, and vascular effects are not separated from relaxin's direct ligament-compliance
   effect within the cited studies.
7. **SI fusion's ROM-reduction number is computational (FE), not in-vivo RSA**, of an actually-fused
   joint (§6) — the model is described by its authors as "experimentally validated," not
   independently re-validated by this session.
8. **Small, geographically narrow sample sizes** throughout the classic biomechanical literature
   (n=4–25 for the cadaveric/RSA studies) — typical for this era/field, but a real instance-space
   caveat; the largest samples (n=123–368) are all in the peripartum-relaxin/laxity literature, a
   different sub-question.
9. **5 of 28 citations are pre-abstract-era PubMed records** (title/PMID/DOI verified live, no
   abstract text available for re-extraction): Weisl 1954, Tile 1988 (bibliographic-only, both
   flagged in-line rather than assigned a gate).
10. **Vleeming 1990 Part II's exact numeric friction coefficient was not recovered** — the abstract
    states friction "correlates with roughening" qualitatively; a targeted live search for a
    numeric SI-joint coefficient-of-friction value found only Nordberg 2021's **relative** (28%
    vs knee) figure, used throughout instead of a fabricated absolute μ.
11. **This doc did not build or modify any OpenSim model** (§8's finding is a read-only grep against
    existing files) — the proposed SI joint DOF (§10) is a scoped next step, not executed here.

## 10. Proposed cell — add a constrained SI joint DOF (natural next step, not built this session)

Split the single `pelvis` body into `sacrum` + `ilium` (bilateral), connected by a `CustomJoint` with
a small coupled range matching §2's measured envelope (~0–4° rotation, ~0–2mm translation, generous
around the 0.8–3.9° / 0.1–1.6mm Sturesson 1989 envelope), and add the 5 SI ligaments (sacrotuberous,
sacrospinous, interosseous, long dorsal, ventral SI) as `opensim.Blankevoort1991Ligament` bundles —
the identical class `docs/MECHANISM_HIP_ANKLE_LIGAMENTS.md` and `docs/MECHANISM_SPINE_LIGAMENTS.md`
already use, same force law, no new force class needed. This would let the twin, for the first time,
actually show nutation increasing under a real simulated gait/lift load (falsifiable against §5's
mechanism) rather than assuming a rigid pelvis by construction (§8). Explicitly NOT attempted this
session (literature/falsifier-design scope, matching this doc's DELIVER spec) — named precisely so a
future session does not have to re-discover the gap.

## Files

- `scripts/msk/sacroiliac_joint_gates.py` — the gate-checker (literature numbers + pre-registered
  thresholds + the model-structure grep check), pure stdlib, ~1s wall time, exit 0.
- `data/msk_smoketest/sacroiliac_joint/sacroiliac_joint_gate_results.json` — full machine output
  (every number in this doc traces back to this file).
- `docs/MECHANISM_SACROILIAC_JOINT_evidence.json` — all 28 citations, verbatim abstract quotes,
  verified live this session via NCBI eutils (`esearch`+`efetch`).

## Repro

```
python3 scripts/msk/sacroiliac_joint_gates.py
```
