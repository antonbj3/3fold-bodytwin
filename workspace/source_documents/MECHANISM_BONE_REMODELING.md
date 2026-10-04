# MECHANISM BONE REMODELING — Wolff's law / Frost mechanostat, first step (2026-07-21)

Adds a **BONE REMODELING** layer above the existing BONE STRESS layer: converts the twin's
predicted femur cortical stress during gait into strain (µε) and classifies it against **Frost's
mechanostat** — the strain-magnitude threshold model that gates whether bone tissue nets
**resorption** (disuse), **maintenance** (adapted homeostasis), or **formation** (mild-overload
modeling). Reports which window this twin's own gait prediction lands in.

Script: `scripts/msk/bone_remodeling.py` (read-only consumer, run: `python3
scripts/msk/bone_remodeling.py`). Data: `data/msk_smoketest/subject2_walking1/bone_remodeling/
bone_remodeling_results.json`.

**Isolation respected:** bodytwin only; no git operations; reads
`data/msk_smoketest/subject2_walking1/bone_stress/bone_stress_results.json` **read-only** and
writes only under `data/msk_smoketest/subject2_walking1/bone_remodeling/`. New files only.

---

## 0. Coupling — what this reuses, and a naming correction

Per task instruction, this layer couples to the bone-stress layer's own **results JSON**, reusing
it rather than recomputing: `data/msk_smoketest/subject2_walking1/bone_stress/
bone_stress_results.json` (`scripts/msk/bone_stress.py`, mtime 2026-07-21 19:14) — Euler-Bernoulli
beam free-body of the femur under hip + knee joint reaction + all 26 live femur-attaching muscles,
subject2/walking1.

**Naming/timing note, stated plainly (not silently glided past):** `bone_stress.py`'s own
docstring references `docs/MECHANISM_BONE_STRESS.md` for further honesty detail — that file **did
not exist on disk** when this build started (confirmed via `find`, start of this session) and
**this layer was built by reading `bone_stress_results.json` directly**, exactly as the task's own
fallback logic anticipates ("if its results JSON exists reuse it" — it does, so that branch was
taken; the corrected-hip/knee-contact-force fallback was not needed). `docs/MECHANISM_BONE_STRESS.md`
subsequently **appeared mid-session** (mtime 19:16:11, between this script's and this doc's own
timestamps) — this repo is explicitly shared between two concurrent mechanism instances
(`COORDINATOR.md` §0), so a sibling instance most likely wrote it concurrently. Spot-checked, not
blindly trusted: its headline table (−54.97 MPa, medial cortex, subtrochanteric, s=0.156,
t=0.56s/36%) matches exactly what this layer independently pulled from the same raw JSON — a
reassuring over-determination, not a dependency (this layer never read that doc; both arrived at
the same numbers from the same underlying JSON independently).

**This script recomputes nothing upstream.** Every stress number below is copied verbatim from
that JSON; this layer only adds the σ→ε conversion and the mechanostat classification on top.

---

## 1. Method — why strain (not stress), why peak (not mean), derived not assumed

Frost's mechanostat is a **feedback-control** model: a single sensed scalar (peak strain per
loading cycle) gates which of three disjoint effector regimes a bone's remodeling machinery sits
in, exactly like a thermostat's sensed temperature gates heating/cooling/neither. Two consequences
follow directly from that structure, not from convention:

- **Strain, not stress, is the sensed variable.** Osteocytes are embedded in the mineralized
  matrix and respond to the matrix's own deformation, not to load per unit area directly. The
  conversion used here is the standard linear-elastic one, ε = σ/E, with E a property of the
  tissue (assumed uniform), so it is geometry-invariant given a fixed E — the same conversion
  bone_stress.py's own literature-anchor comparison already uses (§2).
- **Peak (not mean/RMS) is the correct sample.** Frost's own thresholds are themselves defined as
  peak-strain-per-cycle criteria (§2 quote) — so evaluating at each configuration's own worst-case
  section/frame (which the upstream beam-theory sweep already computed) is the mechanostat-
  consistent choice, not a convenience.

**Scope, stated up front:** this is a **static** strain→remodeling-**signal** map — an
instantaneous classification of predicted peak strain against literature thresholds. It is **not**
a time-integrated remodeling simulation: no bone-mineral-density/cortical-area state variable is
advanced through simulated days or loading cycles, and no remodeling **rate** is modeled. It
inherits the upstream layer's generic (literature-range, not this subject's CT/DXA-measured) femur
geometry.

---

## 2. Literature anchors — live-verified this session via NCBI E-utilities + Crossref REST

Every citation below was queried **directly against the primary bibliographic record** this
session (`eutils.ncbi.nlm.nih.gov/entrez/eutils/{esearch,esummary,efetch}` for PubMed metadata +
abstract text, `api.crossref.org/works/<doi>` as an independent cross-check of title/journal/
volume/issue/pages/year/author) — not recalled, and not taken from a WebSearch results summary
(WebSearch's own budget was exhausted this session; the eutils/Crossref REST APIs are a strictly
stronger anchor anyway: raw primary-source text, machine-queried, not an LLM's paraphrase of a
search snippet).

| citation | PMID | DOI | verified content |
|---|---|---|---|
| Frost HM. "Bone 'mass' and the 'mechanostat': a proposal." *Anat Rec.* 1987 Sep;219(1):1-9. | 3688455 | 10.1002/ar.1092190104 | Abstract (fetched verbatim): *"bone strains in or above the **1500-3000 microstrain** range cause bone modelling to increase cortical bone mass, while strains below the **100-300 microstrain** range release BMU-based remodeling which then removes existing cortical-endosteal and trabecular bone."* |
| Frost HM. "Bone's mechanostat: a 2003 update." *Anat Rec A.* 2003 Dec;275(2):1081-101. | 14613308 | 10.1002/ar.a.10119 | Conceptual update of the same two-threshold negative-feedback picture; abstract has no restated µε numbers — cited for conceptual corroboration only. |
| Frost HM. "A 2003 update of bone physiology and Wolff's Law for clinicians." *Angle Orthod.* 2004 Feb;74(1):3-15. | 15038485 | 10.1043/0003-3219(2004)074<0003:AUOBPA>2.0.CO;2 | Existence/context verified; abstract has no restated µε numbers (full text not accessed — paywalled). |
| Duda GN et al. "Influence of muscle forces on femoral strain distribution." *J Biomech.* 1998 Sep;31(9):841-6. | 9802785 | 10.1016/s0021-9290(98)00080-3 | Re-fetched, full abstract: *"the model with all thigh muscles showed peak surface strains **below 2000 mu epsilon** (45% gait cycle). Under simplified load regimes surface strains reached values **close to 3000 mu epsilon**."* Carried from `bone_stress_results.json`. |
| Taylor ME et al. "Stress and strain distribution within the intact femur: compression or bending?" *Med Eng Phys.* 1996 Mar;18(2):122-31. | 8673318 | 10.1016/1350-4533(95)00031-3 | Re-fetched, full abstract: *"The results of this investigation strongly support the hypothesis that the femur is loaded **primarily in compression**, and not bending as previously thought."* Used here as the §5 consistency-check anchor. |
| Reilly DT, Burstein AH. "The elastic and ultimate properties of compact bone tissue." *J Biomech.* 1975;8(6):393-405. | 1206042 | 10.1016/0021-9290(75)90075-5 | Existence/title/journal/year verified; **no abstract on file at PubMed** — the E=17 GPa figure attributed to it (reused verbatim from `bone_stress.py`'s own `CORTICAL_E_GPA`) is not independently re-extracted from primary text. Same honest gap as upstream, not resolved here. |
| Turner CH. "Three rules for bone adaptation to mechanical stimuli." *Bone.* 1998 Nov;23(5):399-407. | 9823445 | 10.1016/s8756-3282(98)00118-5 | Re-fetched, full abstract: adaptation is driven by **dynamic**, not static, loading; needs only **short duration**; bone cells **accommodate** to customary/routine loading (desensitize). Used in §6 honest gaps. |

**Secondary/tertiary triangulation (explicitly lower-confidence, not a peer-reviewed primary
source):** the Wikipedia "Mechanostat" article was fetched live (raw wikitext) this session and
gives a four-window table with its own numbers: disuse <~800 µε, adapted 800–1500 µε, overload
>1500 µε, fracture >~15000 µε. Used **only** to triangulate a second threshold set, never as the
anchor. Critically, that **same article** states: *"a typical bone (e.g., the tibia) has a
security margin of about 5 to 7 between typical load (**2000 to 3000 microstrain**) and fracture
load"* — i.e., by its **own** table, habitual/typical peak long-bone strain sits inside its
**overload** band, not its own "adapted" (maintenance) band. This is the direct evidence behind
§4's finding on the task's anchor claim.

**On the task's own specified thresholds (1000 / 1500 / 3000 µε):** implemented exactly as
specified (§3), but flagged plainly against the verified Frost 1987 primary text: Frost's own
disuse/resorption threshold is **100-300 µε**, three-to-ten-fold lower than the task's 1000 µε.
This is a widely-repeated **simplified/secondary-literature rounding**, not Frost 1987 verbatim.
Both are reported side by side (`threshold_variants` in the JSON), never silently reconciled.

---

## 3. Results — the map

### 3.1 Primary configuration (the upstream layer's own central-estimate geometry)

Subtrochanteric/proximal shaft, s=0.156 (normalized hip→knee position), **35.7% of the gait
cycle**, literature neck-shaft offset=40mm, cross-section OD=27.0mm/thickness=6.5mm:

| quantity | stress | microstrain | classification (all 3 variants agree) |
|---|---:|---:|---|
| medial cortex (compression) | −54.97 MPa | **3233 µε** | formation-or-higher |
| lateral cortex (tension) | +48.71 MPa | **2865 µε** | formation-or-higher |
| axial-only hypothetical floor (bending zeroed) | −3.13 MPa | **184 µε** | resorption |

### 3.2 Offset sweep (0mm = null/no hip eccentricity → 60mm literature upper bound)

| offset | stress | microstrain | classification | axial-dominated (Taylor-1996-consistent)? |
|---:|---:|---:|---|:---:|
| 0 mm (null) | −26.75 MPa | 1573 µε | formation-or-higher | No |
| 30 mm | −47.19 MPa | 2776 µε | formation-or-higher | No |
| 40 mm (primary) | −54.97 MPa | 3233 µε | formation-or-higher | No |
| 50 mm | −62.82 MPa | 3695 µε | formation-or-higher | No |
| 60 mm | −70.72 MPa | 4160 µε | formation-or-higher | No |

### 3.3 Cross-section (geometry) sweep, at primary offset

| OD × thickness | stress | microstrain | classification |
|---|---:|---:|---|
| 24.0 × 5.0 mm | −81.83 MPa | 4814 µε | formation-or-higher |
| 24.0 × 6.5 mm | −75.30 MPa | 4429 µε | formation-or-higher |
| 24.0 × 8.0 mm | −72.59 MPa | 4270 µε | formation-or-higher |
| 27.0 × 5.0 mm | −60.85 MPa | 3579 µε | formation-or-higher |
| 27.0 × 6.5 mm (primary) | −54.97 MPa | 3233 µε | formation-or-higher |
| 27.0 × 8.0 mm | −52.20 MPa | 3070 µε | formation-or-higher |
| 30.0 × 5.0 mm | −47.02 MPa | 2766 µε | formation-or-higher |
| 30.0 × 6.5 mm | −41.82 MPa | 2460 µε | formation-or-higher |
| 30.0 × 8.0 mm | −39.17 MPa | 2304 µε | formation-or-higher |

**16 raw sweep entries, 14 distinct configurations** (offset=40mm and OD27.0×th6.5 both
re-express the primary point — deduped by exact stress value, not hand-counted).

### 3.4 Headline — which window does gait sit in

**FORMATION-OR-HIGHER (mild overload / modeling), not maintenance and not resorption** — true for
**all 14 distinct swept configurations**, under **all 3 independently-sourced threshold variants**
(task-specified, Frost-1987-verified, Wikipedia-secondary). The single most conservative real
configuration in the entire sweep (offset=0mm, i.e. zero hip eccentricity) still reads 1573 µε —
**~5% above** the 1500 µε formation threshold, not below it. Zero of the 14 configurations land in
maintenance (1000/800/100-300–1500 µε, depending on variant) or resorption.

**Honest caveat on "robust across 3 variants":** all 3 sourced threshold sets share the **identical
1500 µε** formation-onset boundary — they were found to disagree only on the *lower* (disuse)
threshold (100-300 vs 800 vs 1000 µε), which turned out to be irrelevant here since no
configuration is anywhere near it. "Robust across variants" should be read as *robust to that
disuse-threshold disagreement*, not as three independently-conflicting formation-threshold
hypotheses all happening to agree.

---

## 4. The critical caveat: this verdict is likely biased HIGH, not settled

This is the single most important honest finding of this build, so it is its own section rather
than buried in a gap list.

At **every** configuration where the axial/bending split is available (primary + all 5 offsets),
this coarse beam model is **bending-dominated**: |σ_axial| ≪ σ_bend_amp (3.13 vs 51.84 MPa at the
primary point; even at the offset=0 null, 2.06 vs 24.69 MPa). This directly **contradicts** the
verified Taylor et al. 1996 anchor (§2): *"the femur is loaded primarily in compression, and not
bending as previously thought."* Notably, the contradiction **persists even at zero hip
eccentricity** — so it is not solely an artifact of the (swept, literature-range, not
subject-measured) neck-shaft offset; some of the excess bending is coming from elsewhere in this
coarse model (muscle-geometry moment arms, the circular-annulus idealization, or beam theory
itself vs a subject-specific CT-FE mesh).

Since the mechanostat classification above is driven almost entirely by that same bending term,
the "formation, not maintenance" verdict inherits this reliability gap directly:

- If Taylor 1996's compression-dominated regime held **here** instead of this model's own
  bending-dominated one, the pure-axial component alone (184 µε at the primary section) would
  classify as **resorption** under every threshold variant checked.
- Duda et al. 1998's own subject-specific finite-element model — which, like this one, includes
  **all** thigh muscles — found peak strain **below 2000 µε** (comfortably lower than this
  model's 3233 µε at a similar, though not identical, gait-cycle phase: 45% vs 35.7%).

**Conclusion:** the true, subject-specific answer plausibly sits somewhere between this layer's
"formation" estimate and a lower "maintenance"-adjacent estimate closer to Duda 1998's own FE
result. This first-step layer **cannot adjudicate between them** — it inherits, rather than
resolves, the upstream bone-stress layer's own already-flagged reliance on generic (not
subject-measured) geometry. The "formation" headline should be read as this specific coarse
model's honest prediction, not asserted as a settled biological fact.

---

## 5. Machine cross-checks / gates (all PASS, `bone_remodeling_results.json.gates`)

| gate | result |
|---|---|
| self-test (unit conversion + bin-boundary logic, exact round-number construction: 17 MPa @ 17 GPa = 1000 µε exactly) | PASS |
| medial-cortex reconstruction (axial − bend_amp) exactly reproduces upstream's stored peak | PASS (exact) |
| ALL 14 distinct swept configs ≥ 1500 µε (task formation threshold) | True |
| ZERO configs classify maintenance/resorption under task thresholds | True |
| verdict robust across all 3 threshold variants | True (see §3.4 caveat) |
| primary config bending-dominated (contradicts Taylor 1996) | True |
| axial-only hypothetical floor would classify resorption | True |

The self-test is a correctness check on this script's **arithmetic and binning logic only** — it
does not and cannot validate the upstream physics; that burden is carried by §4's caveat and by
`bone_stress_results.json`'s own gates (headline cross-check against hip/knee %BW, synthetic
self-tests, stress sanity band — all already PASS, re-read not re-run here).

---

## 6. Honest gaps

1. **Static map, not a time-integrated simulation** — no bone-mass/BMD/cortical-area state is
   advanced through simulated time; only an instantaneous strain→regime classification at each
   upstream configuration's own peak. No remodeling *rate* is modeled.
2. **Generic femur geometry** (inherited, not new) — literature-range neck-shaft offset and
   cortical cross-section, not this subject's CT/DXA-measured geometry. This is the parameter the
   verdict is most sensitive to: 1573→4814 µε (3.1×) from geometry/offset assumptions alone, all
   of it still within one classification bin (formation-or-higher) only because the low end
   (1573) barely clears 1500.
3. **Task thresholds (1000/1500/3000 µε) are a secondary-literature rounding**, not verbatim Frost
   1987 (100-300/1500-3000 µε) — reported side by side, not reconciled (§2).
4. **No independently-verified 4th "pathologic overload" bin.** Some clinical literature splits
   the formation-or-higher zone further at ~3000 µε; neither fetched Frost abstract (1987, 2003)
   states this split, so it is reported only as a hedge annotation in the JSON
   (`secondary_severity_hedge`), never as a verified boundary.
5. **"Habitual gait sits in maintenance" (the task's own anchor claim) could not be confirmed
   live.** If anything, the Wikipedia secondary source's own numbers place typical/habitual
   long-bone strain in its *overload* band, consistent with the mechanostat's own
   cost-minimizing-controller logic (tune habitual peak strain near the modeling threshold, not
   deep in a quiescent zone) — flagged, not silently assumed (§2, §3.4).
6. **§4's bending-vs-Taylor-1996 contradiction is the single largest reliability caveat** on the
   headline and is explicitly unresolved by this layer.
7. **Rate/cycle-count/accommodation effects omitted** — Turner 1998 (verified, §2) shows dynamic
   loading, short duration, and habituation to customary loading all matter; this map only
   classifies peak magnitude.
8. **Single subject, single trial, right leg only** (subject2/walking1) — no inter-subject or
   inter-trial variability characterized.
9. **E=17 GPa reused, not re-verified** — Reilly & Burstein 1975 citation exists; the specific GPa
   figure has no retrievable PubMed abstract this session, so it is not independently
   re-extracted from primary text (same gap as upstream, carried forward unchanged).
10. **Geometry-sweep points have no independently-stored gait-cycle location** — inherited from
    the primary offset per `bone_stress.py`'s own source comment, not re-verified per-point.
11. **`docs/MECHANISM_BONE_STRESS.md` did not exist when this build started** (§0) — coupled
    directly to the results JSON instead, per the task's own reuse-if-present instruction; the doc
    appeared mid-session from a concurrent sibling instance and was spot-checked (not depended
    on) after the fact.

---

## 7. Files

- `scripts/msk/bone_remodeling.py` — read-only consumer of `bone_stress_results.json`; unit
  conversion, 3-variant threshold classification, self-test, gates, JSON writer.
- `data/msk_smoketest/subject2_walking1/bone_remodeling/bone_remodeling_results.json` — full
  machine-readable evidence: every swept configuration's stress/microstrain/classification under
  all 3 threshold variants, gates, literature anchors (with live-verification method noted
  per-citation), headline text, honest gaps.
- `docs/MECHANISM_BONE_REMODELING.md` — this file.

---
---

# PART B — BMU turnover KINETICS (rate) + a second, decorrelated site (2026-07-22)

Part A (above) is explicitly a **static** strain→regime **classifier** — its own honest_gaps[0]
states plainly: *"no remodeling RATE is modeled."* Part B fills exactly that gap: a first
**quantitative BMU (basic multicellular unit) kinetics** model, cross-checked against **live-verified
in-vivo histomorphometry** (tetracycline double-labeling, Parfitt/Eriksen/Clarke lineage), plus a
**second, decorrelated anatomical site** — the femoral **neck** (`docs/MECHANISM_FEMORAL_NECK_STRESS.md`)
— run through Part A's own mechanostat classifier, reused verbatim, which Part A itself never touched
(Part A only ever read the shaft's `bone_stress_results.json`).

**Script**: `scripts/msk/bmu_turnover_kinetics.py` (imports `bone_remodeling.py` unedited; md5
`3481c9b029f40d7a4cfa54b743a5f2dc`, confirmed identical to the value already recorded in
`docs/MECHANISM_BONE_STRESS_CONSEQUENCE.md` §10 — bit-for-bit reuse, not a fork). **Data**:
`data/msk_smoketest/subject2_walking1/bone_remodeling/bmu_turnover_kinetics_results.json` (a NEW
file; `bone_remodeling_results.json` from Part A is untouched — `git diff --stat` on it is empty).
**Isolation respected**: no git add/commit; only new files written; this repo is confirmed, live,
shared with a concurrent sibling instance this session (`git status` shows unrelated in-flight edits
to other docs/scripts neither read nor touched here).

---

## 8. Two independent falsifiers, run separately, never merged

**FALSIFIER 1** (the task's own framing): *does a BMU-kinetics model built from live-verified
histomorphometry reproduce the independently-reported annual bone-formation/turnover rate (~10%/yr
cortical, ~25%/yr trabecular) and the resorption/formation BALANCE?*

**FALSIFIER 2** (the task's own "decorrelated strain-vs-turnover cross-check"): *does Frost's
mechanostat put the femoral-**neck** gait strain in the maintenance window, using a SECOND
anatomical site/free-body as an independent check on Part A's shaft-only finding?*

---

## 9. Falsifier 1 — turnover rate vs. live-verified histomorphometry

### 9.1 The geometric derivation (not a looked-up ratio)

A point on bone surface is modeled as a two-state renewal process: **idle** (quiescent, lining
cells) or **active** (currently inside a BMU's resorption/formation envelope). By the
renewal-reward theorem — a standard ergodic-theory result, not an assumed formula — the long-run
**fraction of time a point is active** equals its own mean active-sojourn duration (`T_bmu`) divided
by the mean time between successive activations at that same point (`T_recur`):

```
occupancy_fraction = T_bmu / T_recur
```

This is algebraically identical to Parfitt's own `Ac.f × Σ` formalism (`Ac.f` = activation
frequency = `1/T_recur`; `Σ` = BMU total lifespan = `T_bmu`) — derived here from first-principles
renewal theory rather than asserted as a black-box histomorphometric formula (GEOMETRIC THINKING:
the derivation is a time-occupancy argument, the same structure that underlies stereological
surface↔volume conversion). If resorption and formation are volumetrically **coupled** (balanced,
non-pathological) and spatially homogeneous, this occupancy fraction is a first-pass proxy for the
**annual volumetric turnover rate**.

### 9.2 Trabecular — PASSES within the pre-registered factor-of-2 band

Using **two independently live-verified** Eriksen 2010 numbers (PMID 21188536, DOI
10.1007/s11154-010-9153-1, fetched via NCBI efetch + Crossref + PMC full-text this session):
cancellous BMU cycle **200 days** (*"remodeling occurs on the surface of trabeculae and lasts about
200 days in normal bone"*) and full-surface renewal **2 years** (*"the total surface of cancellous
bone is completely remodeled over a period of 2 years"*):

```
occupancy_fraction = 200 / 730 = 27.40 %/yr
```

vs. the task-stated **25%/yr** — ratio 1.10×, **PASS** (pre-registered tolerance: factor of 2).

### 9.3 Cortical — FAILS; reported plainly, not reconciled

Clarke 2008 (PMID 18988698, DOI 10.2215/CJN.04151206, PMC3152283, full text fetched live) states
directly, verbatim: *"The relatively low adult cortical bone turnover rate of **2 to 3%/yr** is
adequate to maintain biomechanical strength of bone."* Midpoint **2.5%/yr** vs. the task-stated
**10%/yr** — ratio **4.0×**, **FAILS** the same factor-of-2 band. This is reported as a genuine,
live-verified discrepancy, not smoothed over: the specific "10%/yr cortical" figure could not be
independently reproduced from what this session verified.

### 9.4 The cross-compartment RATIO itself also fails, even though the trabecular number alone passes

Task's own two numbers imply trabecular turns over **2.5×** faster than cortical. The
live-verified/derived numbers imply **11.0×** (27.40/2.5). **FAILS** the factor-of-2 band on the
ratio, even though §9.2's trabecular absolute number passed. Reported separately — a genuine
non-triviality, not a contradiction (a number and a ratio are different claims).

### 9.5 Geometric explanation for the compartment gap (informational, not independently 3rd-source-checked)

Cortical and trabecular BMU cycle **lengths** differ by <2× (120 days median cortical vs. 200 days
trabecular, both Eriksen 2010, verified). Inverting §9.1's relation with Clarke's own directly-verified
2.5%/yr gives an **implied** cortical mean recurrence interval: `120 / 0.025 = 4800 days ≈ 13.15
years` — vs. trabecular's own **verified** 2-year interval. The compartment turnover-rate gap is
therefore driven mostly by **activation frequency** (how often a given site is revisited), not BMU
cycle speed — the standard literature mechanism for this is trabecular bone's much higher
surface-to-volume ratio (BS/BV) vs. cortical bone (qualitatively corroborated, not numerically
re-verified this session, by Parfitt et al. 1983, PMID 6630513, live-fetched abstract: age-related
bone loss proceeds via resorption-cavity-driven trabecular-plate perforation — a surface-accessible
process). This 13.15-year figure is a **derivation**, not an independently cross-checked 3rd number.

### 9.6 Resorption/formation BALANCE (coupling) check — the task's own "osteoclast vs osteoblast balance" question

```
wall_thickness_derived = MAR × formation_phase_duration
```

Using Eriksen 2010's own verified cancellous formation-phase duration (**150 days**, from *"The
resorption period has a median duration of 30-40 days and is followed by bone formation over a
period of 150 days"*) against the task-stated MAR range (0.5–1.0 µm/day), compared to Eriksen
2010's own verified resorption-lacuna depth (**40–60 µm**, *"depth ... varies between 60 in young
individuals and 40 µm in older individuals"* — the amount a BALANCED BMU must refill):

| MAR (µm/day) | derived wall thickness | vs. 50 µm anchor (mid) | within factor-2? |
|---:|---:|---:|:---:|
| 0.5 (low end, task) | 75 µm | 1.5× | **PASS** |
| 1.0 (high end, task) | 150 µm | 3.0× | **FAIL** |

The LOW end of the task's own stated MAR range reproduces the coupling anchor; the HIGH end
over-predicts by 3×. Reported per-endpoint, not averaged into one verdict.

### 9.7 BMU lifespan — task's 6–9 months is compartment-dependent, not uniform

| source | value | vs. task's 182.6–274.0 days (6–9 mo) |
|---|---:|---|
| Eriksen 2010, cortical median | 120 d (3.94 mo) | **below** the task's range |
| Eriksen 2010, trabecular | 200 d (6.57 mo) | **within** the task's range |
| Clarke 2008, cortical (resorption 2–4wk + formation 4–6mo) | 135.8–210.6 d (4.46–6.92 mo) | straddles the **low edge** |
| Parfitt 1994 (PMID 7962158) | "many months" (qualitative only) | consistent, not quantitative |

"~6–9 months" is best supported for the **trabecular** compartment and Clarke's own upper cortical
estimate — not Eriksen's cortical median point estimate. Reported per-compartment.

### 9.8 A note on a different 8%: osteoclast nuclei ≠ bone volume

Parfitt 1994 (PMID 7962158, verified) also states individual **osteoclast nuclei** turn over "about
**8% per day**" — a cell-biology quantity (nuclear replacement inside one multinucleated osteoclast),
**not** tissue-level bone-volume turnover. Reported for completeness, explicitly **not** used in any
turnover-rate cross-check above, to avoid conflating two different kinetic quantities that happen to
share a "%" unit — a real trap this analysis deliberately avoided, not merely disclaimed.

---

## 10. Falsifier 2 — femoral-NECK mechanostat, a second decorrelated site

Part A classified only the femoral **shaft** (`bone_stress_results.json`). The task's own falsifier
asks about a **decorrelated cross-check** — this repo family's own `docs/MECHANISM_FEMORAL_NECK_STRESS.md`
independently derived a **neck**-axis free body (different geometry, same underlying per-frame
loads) and is the anatomically higher-value site for hip-fracture mechanostat reasoning. Part A's
own classifier (`THRESHOLD_VARIANTS`, `classify_all_variants`, `stress_mpa_to_microstrain`) is
imported **unedited** and applied to **18 real (bending-included) neck configurations**
(primary AS-IS, null control, 6-point NSA sweep, 9-point cross-section sweep, de-inflation
counterfactual) plus **2 axial-only hypothetical floors** (primary, null), none of which Part A ever
computed. A `bit_for_bit_consistency` gate re-derives microstrain from each config's own stored
`stress_MPa` and checks it against the neck JSON's own independently-stored `microstrain_magnitude`
— **PASS**, exact, for all 18 (catches any transcription slip before trusting the classification).

**Result: 18/18 real configurations classify `formation_or_higher_overload`; zero classify
`maintenance_adapted`.** The task's own anchor claim ("habitual gait strain sits in the maintenance
window") is **not supported** at the neck — the same verdict Part A already reached at the shaft,
now independently reinforced at a second site with its own separately-derived geometry (NSA/offset/
anteversion, Carmona et al. 2019) and its own separately-derived valid-domain boundary (§3,
`docs/MECHANISM_FEMORAL_NECK_STRESS.md`).

**A structural feature, not previously stated this explicitly**: both axial-only hypothetical floors
(what the verdict would be if Taylor et al. 1996's compression-dominated regime — PMID 8673318,
verified — held here instead of this coarse model's own bending-dominated one) classify
`resorption_disuse` (487 µε primary, 659 µε null; both < the 1000 µε task threshold), **not**
maintenance either. Across **both** anatomical sites this repo family has built (shaft: Part A;
neck: here), across **every** tested variant (real sweeps AND the axial-only hypothetical), this
class of coarse beam model **never lands a configuration inside the "maintenance" bin** — it is
either resorption (axial-only-hypothetical) or formation-or-higher (real, bending-included), with
"maintenance" apparently jumped over rather than landed in. This is reported as an honest structural
observation about this model class, not oversold as a biological finding.

---

## 11. Symmetric QC — held OPEN, not resolved by this build

1. **Mechanostat thresholds have wide, partly phenomenological spread** (Part A §2: Frost 1987's own
   100–300 µε disuse threshold vs. the task's 1000 µε vs. Wikipedia's 800 µε — only the ~1500 µε
   formation-onset boundary is shared across variants). Falsifier 2's "zero maintenance" finding is
   therefore robust to the *disuse*-threshold disagreement specifically (nothing tested is near it)
   but the maintenance/formation BOUNDARY itself (1000–1500 µε) is exactly the zone this coarse
   model's own bending-vs-axial split brackets over without landing in — genuinely unresolved by this
   analysis, not adjudicated.
2. **Histomorphometry is biopsy-site-specific.** Every histomorphometric number in §9 (Clarke 2008,
   Eriksen 2010, Parfitt 1994/1983, Hauge 2001) derives from **iliac crest** biopsy — the
   conventional histomorphometry site, not the femoral neck/shaft this twin models mechanically.
   Compartment-level (cortical vs. trabecular) kinetics are assumed transferable; **regional**
   (iliac vs. femoral) kinetics are not independently verified to match — an open, unresolved gap,
   consistent with the cross-bone caveat `docs/MECHANISM_BONE_STRESS_CONSEQUENCE.md` §7 already
   disclosed for its own tibial-vs-femoral strain-anchor comparison.
3. **The occupancy-fraction derivation (§9.1) assumes coupled, spatially homogeneous remodeling** —
   true at steady state in normal bone, not in osteoporosis/high-turnover disease (Hauge 2001's own
   verified finding: primary hyperparathyroidism patients show a 50% BRC-surface change after
   parathyroidectomy — turnover state is NOT a fixed constant even within one skeleton over time).
4. **This session's literature search missed on 5 of ~14 attempted PMIDs** recalled from memory
   before live verification (Manolagas 2000, a Recker-et-al histomorphometry-methods guess, a Frost
   1969 tetracycline guess, a Robling 2006 guess, all resolved by title-mismatch — full list:
   `bmu_turnover_kinetics_results.json.literature_search_misses_this_session`). Each was caught, not
   silently used — but this is stated plainly as a reminder that the surviving citations are the
   ones that **passed** verification, not an a-priori-guaranteed-correct set.
5. **Nothing here is proven.** Falsifier 1 PASSES for trabecular turnover-rate (within 2×) and the
   low-MAR coupling check; FAILS for cortical turnover-rate, the cross-compartment ratio, and the
   high-MAR coupling check. Falsifier 2's task-anchor claim FAILS (is not supported) at both tested
   sites. A mixed, honestly-reported scorecard — not a clean win, not a clean loss.

---

## 12. Confidence tier

**In-vivo-anchored (histomorphometry)** for the literature turnover-rate/phase-duration figures
themselves (Clarke 2008, Eriksen 2010, Parfitt 1994/1983, Hauge 2001 — all live-verified
primary-source abstracts or PMC full text this session, PMID+DOI cited per number, §9). **Method-only**
for (a) the renewal-theory occupancy DERIVATION (a first-pass geometric model, not independently
cross-validated against a genuinely distinct 3rd method this session) and (b) the neck mechanostat
classification, which inherits `femoral_neck_stress.py`'s own "partially in-vivo-anchored crossed
with method-only" tier unchanged (this script performs no new mechanical modeling for Falsifier 2,
only re-classifies already-gated numbers).

---

## 13. Honest gaps (Part B additions — Part A's own gap list, §6 above, still stands unchanged)

1. The task's exact "10%/yr cortical, 25%/yr trabecular" figures were not found verbatim in any
   primary source reachable this session — retained as a labeled task-stated variant, cross-checked
   (not replaced) against Clarke 2008 (cortical, direct) and an Eriksen-2010-derived renewal
   calculation (trabecular) — §9.2–9.4.
2. Cortical turnover FAILS its own falsifier check (Clarke 2008's 2–3%/yr is ~4× lower than the
   task's stated 10%/yr) — reported plainly, §9.3.
3. The cross-compartment RATIO fails even though the trabecular absolute number passes — §9.4.
4. The coupling/balance check passes only at the LOW end of the task's stated MAR range — §9.6.
5. "BMU lifespan ~6–9 months" is compartment-dependent: supported for trabecular and Clarke's upper
   cortical estimate, not for Eriksen's cortical median point estimate — §9.7.
6. The neck mechanostat cross-check inherits every honest gap `femoral_neck_stress.py` and
   `bone_stress.py`/Part A already disclosed unchanged (generic geometry, bending-vs-Taylor-1996
   tension, boundary-value neck peak, quasi-static/translational-inertia-only, population-level
   n=628 Carmona geometry not this subject's own) — no new mechanical modeling was performed here.
7. Iliac-crest histomorphometry vs. femoral mechanics is a real, unresolved cross-region gap (§11.2).
8. Single subject (subject2), single trial (walking1), right leg only — same scope boundary as every
   MSK cert in this repo family.

---

## 14. Files (Part B)

- `scripts/msk/bmu_turnover_kinetics.py` — imports `bone_remodeling.py` UNEDITED (md5-verified);
  renewal-theory occupancy derivation, coupling/balance check, neck-mechanostat re-classification,
  self-test, gates, JSON writer. Reads `femoral_neck_stress_results.json` read-only.
- `data/msk_smoketest/subject2_walking1/bone_remodeling/bmu_turnover_kinetics_results.json` — full
  machine-readable evidence: both falsifiers, all literature anchors (PMID+DOI+verified quotes),
  the 5 literature-search misses, gates, honest gaps.
- `docs/MECHANISM_BONE_REMODELING.md` — this file (Part B appended; Part A, §0–§7 above, untouched —
  `git diff --stat` on both `bone_remodeling_results.json` and this doc's Part A confirmed empty
  before this append).
- No git add/commit/push performed (isolation respected). Only new files created; no existing
  tracked file was modified.
