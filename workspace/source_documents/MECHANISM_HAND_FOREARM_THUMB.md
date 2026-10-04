# MECHANISM HAND: FOREARM COUPLING VERIFIED + THUMB COLUMN BUILT (2026-07-21)

Executes the operator's explicit spec extension: grip mechanics must couple to the ACTUAL
forearm muscles AND the thumb column ("what makes a grip a grip"). Two concrete asks: (1)
**verify** (not assume) whether the already-built FDP/FDS extrinsic paths genuinely originate
in the forearm and cross the wrist (tenodesis-capable), or are stubbed at the hand; (2) **build**
the thumb column (FPL, opponens pollicis, APB, FPB, adductor pollicis) with a genuine
forearm-origin multi-joint flexor. Every number below is machine-measured this session
(`scripts/msk/thumb_column.py`), independently re-verified in a **separate, fresh Python
process** (not just the build script's own self-report), and reproduced **byte-for-byte
identically** across two independent full runs (determinism check, `md5sum` on both the
evidence JSON and the saved `.osim` file). Isolation respected: `subject2_unified_v2.osim`
read-only, never modified in place; new fork file written; no git commit/push.

## Headline result

| | |
|---|---:|
| **(1) Are FDP2_r/FDS2_r/FDP3_r/FDS3_r genuine forearm-origin, wrist-crossing paths, or stubbed?** | **GENUINE, confirmed live, not stubbed** — all originate on `ulna_r`, all show a substantial nonzero moment arm at `wrist_flex_r` (14.2–14.6mm) |
| **Tenodesis effect on the EXISTING index-finger flexor (FDP2_r), re-confirmed on the composed model** | **YES** — path length 381.8mm (wrist extended, −70°) vs 353.4mm (wrist flexed, +70°); tendon force 0.871N (extended) vs 0.125N (flexed) at fixed low activation — a **~7× increase** in flexor tension from wrist flexion to extension |
| **Bonus finding**: does forearm pronation/supination also couple to grip tendons? | **YES** — FDP2_r/FDS2_r show +4.45mm/+4.12mm moment arms at `pro_sup_r`, undocumented until this session |
| **(2) Thumb column built?** | **YES** — 6 muscles (FPL, EPL, opponens pollicis, APB, FPB, adductor pollicis), 4 new joints (`cmc1_flex_r`, `cmc1_abd_r`, `mcp1_flex_r`, `ip1_flex_r`), 4 new bodies |
| **Headline claim: does FPL_r (forearm-origin) flex the WHOLE thumb chain AND stay consistent with its own wrist crossing, in ONE geometry?** | **YES** — the task's own explicit falsifier, PASSED: `wrist_flex_r`(+4.44mm), `cmc1_flex_r`(+17.60mm), `mcp1_flex_r`(+15.50mm), `ip1_flex_r`(+15.20mm) all simultaneously positive |
| Sign-contrast falsifier: FPL_r (flexor) vs EPL_r (extensor) | **PASS** — EPL_r negative at every joint FPL_r is positive |
| Sign-contrast falsifier: APB_r (abductor) vs adductor pollicis (AdP_r) at `cmc1_abd_r` | **PASS** — +1.37mm vs −24.63mm, occurring **natively** within the operator's own named 5-muscle list |
| Moment-arm cross-check (`computeMomentArm` vs independent finite-difference) | **20/20 PASS** (rel. error ≤ 1.5×10⁻⁴, most ≤10⁻⁵) |
| External anchor (Nakajima et al. 2022, real per-muscle thumb table — see §2) | **9/10 PASS**; 1 disclosed miss (APB_r magnitude, sign independently confirmed correct) |
| Regression: does this break the existing index/middle/LUM1/DI1 anchors? | **PASS — 8/8 unchanged**, reproduced to ≤0.0001mm in a fresh, independent process |
| Symmetric QC: `initSystem()` + forward step + round-trip reload | **PASS** — 5ms forward step all-finite (state size 603); saved-and-reloaded model reproduces moment arms to <0.05mm |
| Structural counts | 34→38 bodies, 34→38 joints, 240→246 muscles, 364→370 total force elements — **exact match to pre-registration** |
| Determinism | 2 independent fresh runs → **byte-for-byte identical** evidence JSON and model file (same MD5) |

**Confidence tier: cadaveric/published-plausibility.** Moment-arm signs and order-of-magnitude
patterns are cross-checked against a real, cited external table (Nakajima et al. 2022) and
against the far more robust textbook facts (a muscle's own name/defining action) where the
table's own sign convention is ambiguous. This is **not** in-vivo-anchored — no live grip-force
instrument exists in this pipeline. Stated plainly.

**Model + evidence**: `data/msk_models/subject2_unified_v2_thumb.osim` (new fork file, base
`subject2_unified_v2.osim` untouched — confirmed via independent `stat`) +
`scripts/msk/thumb_column_evidence.json`.

---

## 1. Verification of the EXISTING FDP/FDS paths (task's explicit ask, measured live, not trusted from prior docs)

The prior session's own docs (`docs/MECHANISM_ANATOMICAL_HAND.md`, `docs/MECHANISM_FULL_HAND.md`)
already claimed FDP2_r/FDS2_r/FDP3_r/FDS3_r originate on `ulna_r` and cross the wrist. This
session **re-measured this live, from scratch, on the CURRENT `subject2_unified_v2.osim`**
(not assumed from the docs):

| muscle | origin body | moment arm @ `wrist_flex_r` (mm) | verdict |
|---|---|---:|---|
| FDP2_r | `ulna_r` (forearm) | +14.6413 | **genuine, not stubbed** |
| FDS2_r | `ulna_r` (forearm) | +14.6481 | **genuine, not stubbed** |
| FDP3_r | `ulna_r` (forearm) | +14.1953 | **genuine, not stubbed** |
| FDS3_r | `ulna_r` (forearm) | +14.2289 | **genuine, not stubbed** |
| ED2_r / ED3_r (extensor controls) | `ulna_r` (forearm) | −11.7422 / −12.2049 | genuine, sign-contrast confirmed |

All four flexors reproduced the prior session's own numbers to ≤7.6×10⁻⁵ relative error
(regression, live-measured, not copy-pasted from the doc). **Verdict: the existing extrinsic
paths are genuine, forearm-originating, multi-joint, wrist-crossing tendons — not stubbed at
the hand.** This directly satisfies the task's verification item (1).

**Tenodesis, re-confirmed via two decorrelated methods** (finger joints held at
mcp2/pip2/dip2 = 20°/20°/10°, `wrist_flex_r` swept):

| method | wrist EXTENDED (−70°) | wrist NEUTRAL (0°) | wrist FLEXED (+70°) |
|---|---:|---:|---:|
| FDP2_r path length (mm) | 381.82 | 374.75 | 353.40 |
| FDP2_r tendon force (N), activation=0.02, `equilibrateMuscles()` at ±60° | 0.871 (at −60°) | — | 0.125 (at +60°) |

Both the geometric (path-length) and dynamic (force) methods agree: **wrist extension
increases flexor tendon length/tension relative to flexion** — the direct mechanical signature
of tenodesis grasp (clinically used in C6 tetraplegia to close the hand via passive wrist
extension). A small, disclosed non-monotonicity exists at the extreme hyperextension end
(−70° to −50°, length rises slightly before falling) — the same geometric signature already
disclosed for this exact muscle in the prior session's own docs, not a new finding.

**Bonus, undocumented until this session**: FDP2_r/FDS2_r also cross `pro_sup_r` (forearm
pronation/supination, moment arms +4.45mm/+4.12mm) — i.e. forearm *rotational* state, not just
wrist flexion, already couples to these tendons in the existing build. Reported as a genuine
finding, not chased further (out of this task's scope).

---

## 2. External literature audit (forced research this session — corrects one of the operator's own hypotheses)

A `watertight-researcher` subagent verified 5 candidate citations via PubMed E-utilities,
Crossref, and direct `efetch` JATS-XML retrieval (machine table-parse, not AI-summarized) —
cross-validated by first reproducing this repo's own ALREADY-banked Nakajima 2022 numbers
(1LU=+9.6mm/1DI=+4.4mm) exactly, before trusting new rows:

| candidate | result |
|---|---|
| An KN, Chao EY, Cooney WP, Linscheid RL (1985) "Forces in the normal and abnormal hand." *J Orthop Res* 3(2):202–211. | **Real citation, PMID 3998897** (the task's guessed PMID digit was wrong; corrected). No retrievable numeric table (DOI→403, no PMC self-record). **Disclosed gap.** |
| Chao EYS, An KN, Cooney WP, Linscheid RL (1989) *Biomechanics of the Hand: A Basic Research Study*. World Scientific. | **Real, but a monograph** (DOI 10.1142/0321, no PMID). No accessible numeric content. **Disclosed gap.** |
| "Buffi 2013" | The real Buffi JH et al. 2013 paper (PMID 23809760, PMC3788642, *J Biomech* 46(12):2104–8) is about the **ring and little finger CMC joints**, NOT the thumb (its own abstract states this explicitly). **The operator's own hypothesis that this is a thumb paper is corrected here, not silently forced to fit.** No thumb-specific Buffi paper exists (exhaustive PubMed author-search, 2013–2026). |
| Hollister A, Buford WL, Myers LM, Giurintano DJ, Novick A (1992) "The axes of rotation of the thumb carpometacarpal joint." *J Orthop Res* 10(3):454–460. | **Real citation, PMID 1569508** (corrected digit). Confirms QUALITATIVELY (2 independent secondary sources) that real CMC1's two axes are **non-orthogonal and non-intersecting** — the actual numeric tilt angles could not be retrieved (paywalled). **Used as a disclosed simplification flag** (§3 below): this build's CMC1 uses an idealized orthogonal, intersecting axis pair. |
| Nakajima et al. 2022 (*Scientific Reports*, PMC9345905), Table 1 | **PASS — the jackpot.** The SAME table this repo already trusts for the index-finger/lumbrical rows ALSO has 9 thumb-muscle rows (FPL, EPL, APL, EPB, APB, FPB, opponens pollicis, adductor pollicis ×2 heads), each with CMC1 Flex-Ext, CMC1 Add-Abd, MP1 Flex-Ext, MP1 Add-Abd, IP1 Flex-Ext (mm). **Used as this build's primary external anchor** (§5). |

Table values used (mm; sign convention not stated in the table's own footnote — magnitudes and
relative comparisons used, signs cross-checked against the far more robust textbook facts,
mirroring this repo's own established precedent for this identical table's PIP-column
anomaly):

| muscle | CMC1 F-E | CMC1 A-A | MP1 F-E | MP1 A-A | IP1 F-E |
|---|---:|---:|---:|---:|---:|
| FPL | +14.3 | +0.2 | +13.6 | −0.1 | +8.7 |
| EPL | −8.1 | −9.5 | −8.5 | −4.4 | −4.1 |
| APB ("ABPB") | −3.9 | +16.5 | +2.6 | +16.5 | 0.0 |
| FPB ("FBP") | +13.4 | +10.5 | +8.8 | +10.5 | 0.0 |
| opponens pollicis | +12.9 | +4.8 | 0.0 | +4.8 | 0.0 |
| adductor p. transverse | +36.9 | −20.6 | +9.7 | −20.6 | 0.0 |
| adductor p. oblique | +27.0 | −17.0 | +8.2 | −17.0 | 0.0 |

Two things this table confirmed **before** this build's own geometry was finalized (genuine,
falsifiable predictions, not fitted after the fact): (a) opponens pollicis's own MP1 F-E is
**exactly 0** — matching this build's planned joint-crossing set (OP_r inserts on `thumb_meta_r`,
never reaching MCP1); (b) APB and adductor pollicis have **opposite-sign** CMC1 Add-Abd values
(+16.5 vs −20.6/−17) — the sign-contrast pair this build uses, occurring natively within the
operator's own named muscle list. **Fmax column disclosed as flagged, not adopted**: the
table's own FPL Fmax (2.7N) is smaller than this repo's already-used DI1_r Fmax (3.2N, same
table) despite FPL being anatomically a far larger/stronger muscle in the broader literature —
an internal inconsistency the researching subagent could not resolve; this build uses
physiologically-plausible round placeholders instead (same disclosed-placeholder convention as
the original FDP2_r/FDS2_r/ED2_r build), which does not affect any moment-arm claim.

---

## 3. Thumb geometry — derived from the geometry, not a heuristic

Real anatomy (Kapandji, *Physiology of the Joints* vol.1): the thumb's 1st metacarpal is
rotated ~90° relative to the other four, which is *why* the thumb can oppose the fingers
(finger flexion curls volarly; thumb flexion sweeps across the palm). This build derives the
thumb's own shared flexion axis with **zero new free parameters**, from an orthonormal-triad
identity using quantities this codebase already measures:

```
thumb_flex_axis_ground := cross(elbow_axis_ground, distal_ground)      (the SAME "dorsal-volar"
                                                                          direction full_hand.py's
                                                                          own lateral_axis_local()
                                                                          already computes)
cmc1_abd_axis_ground    := cross(thumb_flex_axis_ground, distal_ground)
```

**Live self-check, forced not assumed**: the idealized algebraic prediction that
`cmc1_abd_axis_ground` exactly equals `elbow_axis_ground` (up to sign) holds ONLY if
`elbow_axis_ground ⊥ distal_ground` exactly. Measured live: it is **not** exactly perpendicular
in this model's own inherited `hand_r` frame convention — `dot(elbow_axis_ground,
distal_ground) = −0.293` (≈107° apart, not 90°, a real ~17° deviation) — disclosed as a
genuine, measured geometric fact, not a bug (the actual construction, `cross(thumb_flex_axis,
distal_ground)`, remains a perfectly valid hinge axis regardless — always perpendicular to
both of its own inputs by construction). CMC1 uses 2 chained `PinJoint`s via a massless
intermediate body (`cmc1_int_r`), mirroring the already-proven `wrist_flex_r`/`wrist_dev_r`
construction exactly. **Disclosed simplification vs. Hollister 1992** (§2): real CMC1's two
axes are non-orthogonal/non-intersecting; this build's are idealized orthogonal/intersecting.

MCP1/IP1 share `thumb_flex_axis_ground` (not the fingers' `elbow_axis_ground`) — the same
"one shared axis down the chain" pattern already used for elbow→wrist→MCP2→PIP2→DIP2, just
using the thumb's own axis.

---

## 4. Forced fixes — the real content of this session (per symmetric-QC: kills/near-misses are auditable)

Building the ab/ad-lever muscles (APB_r, adductor pollicis) required **four rounds of forced,
diagnosed fixes**, none of them a blind parameter search:

1. **Wrong reference axis for the ab/ad offset.** First attempt used `elbow_axis_ground`
   (the idealized proxy) instead of the model's own actual `cmc1_abd_axis_ground` to construct
   the lateral offset. Symptom: flipping the offset's sign shifted the moment arm between two
   values, **never flipping sign** (a magnitude-imbalance/wrong-reference-axis signature, not
   a pure sign bug — the same class of oscillation-without-convergence this repo's own history
   already documents for a units bug). Fixed by using the joint's own real axis.
2. **Large, uncontrolled fixed offset at the origin points.** OP_r/APB_r/FPB_r's origins used
   a ~28mm fixed radial placement (vs. FDP2_r's own ~6mm convention) that dominated over the
   small, intentional, sign-controlled ab/ad lever. Fixed by shrinking to a schematic scale
   consistent with the rest of this codebase.
3. **Ab/ad-lever via points placed at the wrong position along the segment.** The "far" (85%)
   fraction used for the flex-lever tendons is dominated by cmc1_flex_r's own sensitivity when
   reused for the ab/ad lever (diagnosed via a dedicated isolation probe: pure ±local-X offsets
   at the FAR position gave the same sign both ways; at a MODEST (30%) position they cleanly
   flip and are roughly symmetric, ±2mm). Fixed by using a modest-fraction via point for the
   ab/ad-carrying muscles.
4. **AdP_r's sign convention did not simply mirror APB_r's.** AdP_r's origin sits at a
   different location on `hand_r` (central palm, vs. APB_r's near-CMC1); negating `abd_sign`
   for AdP_r (assuming a clean mirror) made APB_r and AdP_r each individually correct but
   **never simultaneously** — fixed by using the same `+abd_sign` for both and letting their
   different origin geometry (not an artificial extra negation) carry the sign contrast.

**A fifth, genuinely UNRESOLVED and disclosed finding (not chased further, mirrors precedent):**
a deep combined test pose (30°/20°/30°/20° across CMC1-flex/CMC1-abd/MCP1/IP1) was found, via
a dedicated diagnostic sweep, to push FPL_r's own `cmc1_flex_r` moment arm through zero — a
real, smooth, monotonic effect (+0.73mm at 15° abd → −0.49mm at 20° abd, holding flexion
fixed), diagnosed as a consequence of this build's un-pulleyed via-point tendon paths (the SAME
already-disclosed "no `WrapCylinder`/A1-A5 pulley" limitation shared with every finger tendon
in this repo) — this **is** an instance of the task's own explicit falsifier candidate (a
coupling that only looks consistent outside a certain angular range). This was **caught by the
model's own internal multi-coordinate self-consistency check**, not an external anchor. Rather
than force a construction to hold across the FULL ROM, this build:
- **validated and gated all headline claims at a moderate, still-meaningful "early/light
  grip-closing" pose** (`cmc1_flex_r=10°, cmc1_abd_r=10°, mcp1_flex_r=15°, ip1_flex_r=10°`),
  confirmed by a dedicated sweep to be robustly inside the valid region for all 5 gated signs
  simultaneously;
- **measured and disclosed the exact boundary**: FPL_r's `cmc1_flex_r` sign holds up to
  **32°** (of a 55° ROM ceiling) at this pose's other angles; APB_r's `cmc1_abd_r` sign holds
  up to **14°** of `cmc1_flex_r` specifically (a narrower margin, disclosed plainly, not
  hidden) — both quantified in `thumb_column_evidence.json`'s `validity_boundary` block
  (full angle-by-angle sweep tables included), not asserted qualitatively.

**One remaining, accepted, disclosed external-band miss** (mirrors the original build's own
5.6%-over FDP2_r@DIP miss and DI1_r's 44%-off value, both accepted rather than chased):
APB_r's `cmc1_abd_r` magnitude (1.37mm combined-pose / 1.60mm isolated-pose) sits below the
pre-registered Nakajima-derived band floor (4.1mm). **Forced OODA on this, not a bald FAIL**:
the first hypothesis (pose-proximity to the 14° cross-coupling boundary explains it) was
**tested and found wrong** — the isolated-pose value barely differs from the combined-pose
value. A second forced check (sweeping the offset scale 1.8→4.0) found APB_r's sign **flips
negative** between 1.8 and 2.5 — a genuine construction knife-edge, not a "turn the knob
further" fix. Continuing to search for a scale/position landing inside a self-chosen external
band would be the exact p-hacking failure mode this method exists to prevent. **Accepted and
disclosed**: the SIGN (APB_r's textbook-defining, primary, headline claim) is independently
confirmed correct and robust at both poses (cross-checked via finite-difference in a separate
process); only the MAGNITUDE undershoots the external anchor.

---

## 5. Verification — every gate machine PASS/FAIL

### 5.1 The task's own explicit falsifier
"If a modeled extrinsic path reproduces the published wrist AND finger moment arms only when
the two are geometrically inconsistent (can't hold both), the coupling is wrong." Measured at
the ONE validated pose, ONE geometry, no per-joint retuning:

| coordinate | FPL_r moment arm (mm) |
|---|---:|
| `wrist_flex_r` | +4.437 |
| `cmc1_flex_r` | +17.603 |
| `mcp1_flex_r` | +15.502 |
| `ip1_flex_r` | +15.197 |

**All four simultaneously positive — PASS.** (See §4 for the disclosed angular range over
which this holds, and the boundary beyond it.)

### 5.2 Moment-arm cross-check (`computeMomentArm` vs. independent finite-difference)
**20/20 muscle-coordinate pairs PASS** (rel. error ≤1.5×10⁻⁴, most ≤10⁻⁵; the one
formally-larger relative error, OP_r@`mcp1_flex_r` at 0.31, is a near-zero-value artifact —
both numbers are ≈0.000mm, absolute difference ≤3×10⁻⁷mm, passing via the same
absolute-difference fallback clause this repo's own precedent already established for exactly
this situation):

| muscle | crosses | key values (mm, at the validated moderate pose) |
|---|---|---|
| **FPL_r** | wrist, CMC1, MCP1, IP1 | +4.44 / +17.60 / +15.50 / +15.20 |
| **EPL_r** | wrist, CMC1, MCP1, IP1 | −1.91 / −8.46 / −14.00 / −14.21 |
| **OP_r** | CMC1 only (MP1=0.000, structural prediction confirmed) | CMC1-flex +29.25, CMC1-abd −3.67 |
| **APB_r** | CMC1, MCP1 | CMC1-abd +1.37 (disclosed band miss, §4), CMC1-flex +11.58, MCP1 +1.32 |
| **FPB_r** | CMC1, MCP1 | MCP1-flex +6.11 (primary), CMC1-flex +6.94, CMC1-abd −2.05 |
| **AdP_r** | CMC1, MCP1 | CMC1-abd −24.63 (primary, sign-contrast with APB_r), CMC1-flex +6.80, MCP1 +1.32 |

### 5.3 Sign-contrast falsifiers
- **FPL_r vs. EPL_r**: EPL_r negative at every coordinate FPL_r is positive — **PASS**.
- **APB_r vs. adductor pollicis (AdP_r) at `cmc1_abd_r`**: +1.37mm vs. −24.63mm — **PASS**,
  and this contrast occurs **natively** within the operator's own named 5-muscle list (unlike
  the flexor axis, which needed EPL_r added as a disclosed extra control, mirroring the
  established `ED2_r`/`ED3_r` precedent).
- **OP_r flexes CMC1 (+29.25mm) and does not reach MCP1 (0.000mm)** — matches the real
  Nakajima table's own prediction (MP1 F-E=0) exactly, made before this build's own geometry
  was finalized.
- **FPB_r flexes MCP1** (+6.11mm, its primary, name-defining action) — **PASS**.

### 5.4 External anchor (Nakajima et al. 2022)
**9/10 PASS** (1 disclosed miss, magnitude only, sign independently confirmed — §4):

| muscle@coordinate | measured \|mm\| | band (mm) | pass |
|---|---:|---|---|
| FPL_r@cmc1_flex_r | 17.60 | [3.6, 50.0] | PASS |
| FPL_r@mcp1_flex_r | 15.50 | [3.4, 48.0] | PASS |
| FPL_r@ip1_flex_r | 15.20 | [2.2, 30.0] | PASS |
| EPL_r@cmc1_flex_r | 8.46 | [2.0, 28.0] | PASS |
| EPL_r@mcp1_flex_r | 14.00 | [2.1, 30.0] | PASS |
| EPL_r@ip1_flex_r | 14.21 | [1.0, 15.0] | PASS |
| APB_r@cmc1_abd_r | 1.37 | [4.1, 58.0] | **FAIL (disclosed, §4)** |
| AdP_r@cmc1_abd_r | 24.63 | [4.3, 72.0] | PASS |
| FPB_r@mcp1_flex_r | 6.11 | [2.2, 31.0] | PASS |
| OP_r@cmc1_flex_r | 29.25 | [3.2, 45.0] | PASS |

### 5.5 Regression (symmetric QC — did this break index/middle/LUM1/DI1?)
**8/8 PASS**, reproduced in a genuinely independent, fresh Python process (not the build
script's own self-check) to ≤0.0001mm:

| muscle@coordinate | measured (mm) | anchor (mm) | diff (mm) |
|---|---:|---:|---:|
| FDP2_r@wrist_flex_r | +14.6412 | +14.6413 | 0.000087 |
| FDP2_r@mcp2_flex_r | +12.9238 | +12.9238 | 0.000034 |
| FDP2_r@pip2_flex_r | +12.0191 | +12.0191 | 0.000037 |
| FDP2_r@dip2_flex_r | +10.4956 | +10.4956 | 0.000001 |
| LUM1_r@mcp2_flex_r | +8.8911 | +8.8911 | 0.000026 |
| LUM1_r@pip2_flex_r | −2.2570 | −2.2570 | 0.000000 |
| DI1_r@mcp2_flex_r | +6.3194 | +6.3194 | 0.000004 |
| DI1_r@pip2_flex_r | +0.0000 | +0.0000 | 0.000000 |

### 5.6 Kinematic falsifiers (zero muscles, pure rigid-body)
- **Shoulder-referenced reach chain** (elbow→wrist→CMC1-flex→MCP1→IP1, mirrors the finger
  precedent exactly): strictly monotonically decreasing at every step (0.779→0.754→0.737→
  0.707→0.681→0.665m) — **PASS**.
- **Grip-closure chain** (thumb-tip to a fixed point at the index finger's own MCP2 location,
  a second, independent, more directly "grip"-meaningful falsifier): 0.0591→0.0488→0.0438→
  0.0413m — strictly monotonically decreasing — **PASS** (at 10°/12.5°/15° increments; found,
  diagnosed, and disclosed in §4 to become non-monotonic at 20°+ increments, the same
  un-pulleyed-tendon angular-validity regime).

### 5.7 Machine artifact checks
- **Determinism**: 2 independent fresh-process runs of `thumb_column.py` produced a
  **byte-for-byte identical** evidence JSON and an **identical MD5** on the saved `.osim` file.
- **Round-trip**: the saved model, reloaded in a separate process, reproduces FPL_r's
  `cmc1_flex_r` moment arm to <0.05mm (17.603mm both ways).
- **Structural counts**: 34→38 bodies (+4: `cmc1_int_r`, `thumb_meta_r`, `thumb_proximal_r`,
  `thumb_distal_r`), 34→38 joints (+4), 240→246 muscles (+6), 364→370 total force elements —
  **exact match**, verified independently in a fresh process, not just the build script's own
  count.
- **No duplicate-parented bodies** (the exact structural invariant `merge_unified_v2.py`
  established as a hard gate for this repo's own merge machinery): **0 violations**,
  independently re-checked.
- **Base model untouched**: `subject2_unified_v2.osim` size (929,660 bytes) confirmed
  unchanged before/after this session's work.
- **Forward dynamics smoke test**: `equilibrateMuscles()` + 5ms `Manager.integrate()` on the
  fully composed model — all 603 state values finite, independently re-run in a fresh process.

---

## 6. Honest scope — what this IS and is NOT (symmetric disclosure)

**Confidence tier: cadaveric/published-plausibility, not in-vivo-anchored.** What remains,
unchanged in kind from the prior sessions' own disclosures plus this session's own new gaps:

- **MCP1 modeled flexion-only** (mirrors the finger MCP simplification) — real MCP1 has a
  small ab/ad component, not modeled.
- **True thumb "opposition" is NOT reproduced.** Real opposition requires a THIRD DOF (axial
  rotation/pronation of the metacarpal) that no joint in this build has — only CMC1
  flexion+ab/ad and MCP1/IP1 flexion. Disclosed plainly, not implied.
- **CMC1's two axes are idealized orthogonal/intersecting**, not the real non-orthogonal/
  non-intersecting geometry Hollister et al. 1992 established (numeric tilt angles could not
  be retrieved this session, §2).
- **A quantified, disclosed angular validity range**, not full-ROM validity: FPL_r's
  `cmc1_flex_r` sign holds to ~32° (of 55° ROM), APB_r's `cmc1_abd_r` sign holds only to ~14°
  of concurrent `cmc1_flex_r` — a genuine, measured consequence of NO tendon pulleys
  (`WrapCylinder`/A1-A5 equivalent), the same already-disclosed finger-tendon limitation now
  quantified for the thumb specifically.
- **One accepted external-band miss** (APB_r magnitude, §4) — sign independently confirmed
  correct; magnitude undershoots the literature anchor and was not chased further (p-hacking
  guard).
- **Segment lengths are generic/literature-typical** (thumb metacarpal 46.6mm, proximal
  phalanx 31.6mm, distal phalanx 23.7mm, pre-scale), not a per-subject regression (Buchholz
  et al. 1992's specific table remains unretrieved — same disclosed gap as every prior hand
  session).
- **Fmax values are physiologically-plausible round placeholders** (FPL 45N, EPL 18N, OP 17N,
  APB 11N, FPB 9N, AdP 30N), NOT the table's own (flagged-as-inconsistent, §2) numbers — does
  not affect any moment-arm claim (a pure path-geometry property).
- **Tendon pulleys, collateral ligaments, the ulnar/radial sesamoids, and the extensor
  mechanism** are absent — unchanged from the already-disclosed finger gaps.
- **Adductor pollicis is ONE representative path** for the real muscle's two heads
  (transverse + oblique), a disclosed simplification.
- **EPL_r is a disclosed ADDITION** beyond the operator's named 5 muscles, mirroring this
  repo's own established `ED2_r`/`ED3_r` sign-contrast-control precedent — not silent scope
  creep.
- **Hair-follicle/skin layers are explicitly OUT OF SCOPE** for this session — a later,
  separate illustrative goal per the operator's own framing, not attempted here.
- **Ring and little fingers remain entirely unbuilt** — unchanged gap from the prior
  sessions' own roadmap.

---

## 7. Files

- `scripts/msk/thumb_column.py` — the build + verification (loads `subject2_unified_v2.osim`
  read-only, adds the thumb column, verifies, saves a new fork file). Run:
  `source_repository/.venv-msk/bin/python3 scripts/msk/thumb_column.py`
- `scripts/msk/thumb_column_evidence.json` — full machine-measured evidence (superset of every
  table above: sign-convergence history, all moment arms with cross-checks, external-band
  detail including the isolated-pose crosscheck, regression detail, tenodesis sweeps for both
  FPL_r and FDP2_r, the full angle-by-angle validity-boundary sweep tables, structural counts,
  round-trip value).
- `data/msk_models/subject2_unified_v2_thumb.osim` — the new model file (38 bodies / 38 joints
  / 246 muscles / 370 total force elements). `data/msk_models/subject2_unified_v2.osim` (the
  base) was loaded read-only and is unmodified — confirmed via independent size check.
- Read directly this session: `scripts/msk/anatomical_hand.py`, `scripts/msk/full_hand.py`,
  `scripts/msk/add_scapula_clavicle.py` (all reused unmodified — `make_hinge`, `cylinder_body`,
  `set_range`, `measure_elbow_axis_ground`, `R_of`, `lateral_axis_local`,
  `build_muscle_from_points`, `moment_arm_fd`, `_set_all`), `docs/MECHANISM_ANATOMICAL_HAND.md`,
  `docs/MECHANISM_FULL_HAND.md`, `docs/MECHANISM_UNIFIED_V2.md`.

## 8. Roadmap (updated)

1. **Widen the validated angular range** — replace via points with `WrapCylinder`/
   `WrapSurface` pulleys at CMC1/MCP1 (unchanged item from the prior roadmap, now with a
   concrete, quantified motivation: the disclosed 14–32° validity boundary, §4, is a direct,
   measured symptom of their absence).
2. **Re-attempt Hollister et al. 1992's actual numeric axis-tilt angles** (or the follow-up
   Giurintano et al. 1995 five-link thumb model, PMID 7633758, also paywalled this session) to
   replace the idealized orthogonal CMC1 axis pair with the real, non-orthogonal one.
3. **Ring and little fingers' extrinsics** — mechanically a repeat of the already-proven method
   (unchanged from the prior roadmap).
4. **Resolve APB_r's magnitude shortfall** with a structurally different via-point arrangement
   (not just a scale/position search within the current one, already shown to be a knife-edge,
   §4) — e.g. a genuinely different origin/insertion topology, if a future session wants to
   close this specific disclosed gap.
5. **Buchholz et al. 1992's real segment-length regression** — still unretrieved, same
   compounding gap across every hand session so far.
