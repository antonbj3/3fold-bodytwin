# MECHANISM HAND: COMPLETE — one model, five digits, one pulley method (2026-07-21)

Executes the operator's explicit ask: stop having scattered hand forks. Three partial,
sibling forks existed (all built off the same common ancestor,
`subject2_unified_v2_thumb.osim`, 38 bodies/246 muscles: index+middle+thumb):

| Fork | Has | Missing |
|---|---|---|
| `subject2_unified_v2_ringlittle.osim` | ALL FIVE digits (44 bodies/257 muscles) | NO pulleys anywhere — moment-arm sign breaks at ~14-59° of flexion, disclosed per-finger in `docs/MECHANISM_HAND_FIVE_DIGIT.md` §5.7 |
| `subject2_unified_v2_pulleys.osim` | The annular-pulley via-point-tightening fix, verified 13/13 sign-hold across FULL ROM | ONLY on digits 2/3 + FPL — built before ring/little existed |
| `subject2_unified_v2_thumb.osim` | Common ancestor of both | — |

This session composes ONE model: starts from `ringlittle.osim` (all 5 digits) and extends
`hand_pulleys.py`'s own tightening method — imported and reused, not re-derived — from
digits 2/3+FPL to digits 4/5 too. **Every existing fork file is untouched** (read-only,
confirmed by size check); the composed result is a genuinely new file,
`data/msk_models/subject2_hand_complete.osim`.

**Status: HYPOTHESIS awaiting independent QC.** Confidence tier: **cadaveric/published-
plausibility** — inherits the tier from both parent docs (`docs/MECHANISM_HAND_PULLEYS.md`,
`docs/MECHANISM_HAND_FIVE_DIGIT.md`): the Landsmeer/An/Armstrong-Chaffin geometric
tendon-excursion model is real, cited, textbook-standard biomechanics; the specific radius
values are a disclosed, schematic, not-directly-cited-from-a-paper's-own-table choice
(reusing this repo's own already-established `R_RATIO`/`PULLEY_VIA_OFFSET_RATIO` numbers a
second time on digits 4/5, not a new tunable knob).

## Headline result

| Falsifier (pre-registered) | Result |
|---|---:|
| (1) Every finger's FDP simultaneously positive at wrist+MCP+PIP+DIP, ED negative, ONE geometry (all 5 digits + wrist posed together) | **PASS**, 4/4 fingers |
| (2) FPL positive at wrist+cmc1+mcp1+ip1, SAME geometry | **PASS** |
| (3) Moment-arm SIGN holds across the FULL declared ROM, all 5 digits' flexors | **PASS — 23/23** (was 13/13 on digits 2/3+FPL per `hand_pulleys.py`; **10/10 new** on digits 4/5; original 13/13 still hold unchanged) |
| (4) Lumbrical double action (LUM1_r/LUM3_r/LUM4_r: MCP+, PIP−) | **PASS**, 3/3 |
| (5) ALL moment arms, `computeMomentArm` vs. independent finite-difference | **PASS — 57/57** (pre-registered union of every crossing any prior session in this repo already established for that muscle) |
| (6a) `initSystem()` + forward-dynamics step, all-finite | **PASS** |
| (6b) IK sanity (real marker IK on `walking1.trc` vs. reference `walking1.mot`) | **PASS**, worst unaffected-coordinate RMSE 1.056° (< 2.0° established tolerance) |
| (7) No duplicate-parent-body (reused `check_no_duplicate_parent_bodies` from `merge_unified_v2.py`, independently re-derived here too) | **PASS**, 0 dupes both ways |
| Structural-count arithmetic (pre-registered sum) | **PASS** — bodies 44→44, joints 44→44, muscles 257→257, wrap objects 44→59 (+15 = 9 replicated + 6 new), all exact |
| Regression vs. `ringlittle.osim` base (untouched muscles/joints reproduce exactly) | **PASS — 20/20** |
| External anchor: composed model vs. the ALREADY-independently-verified sibling `pulleys.osim` (digit 2/3/thumb) | **PASS — 13/13**, diff 0.000000mm (byte-identical reproduction of an external, independently-built file — not a tautology) |
| Round-trip (fresh reload) | **PASS**, all 23 pairs < 0.05mm |
| Tenodesis re-confirmation (FPL_r, FDP2_r, + NEW: FDP4_r, FDP5_r) | **PASS**, 4/4 |
| Determinism (2 independent OS-level builds) | **PASS** — byte-identical `.osim` (md5 `b057a3c79198cf36e0b50b94512740ac` both times), **after a forced correction to a self-inflicted false-negative in this session's own first-attempt test** (§5) |
| **A genuine, forced, pre-adoption geometric risk found and fixed** | Naively reusing `hand_pulleys.py`'s method on ring/little's shared via point would have made the MCP sign break EARLIER than doing nothing (§2) — caught and fixed BEFORE it reached the delivered model |

**Model + evidence**: `data/msk_models/subject2_hand_complete.osim` (44 bodies / 44 joints /
257 muscles / 61 coordinates / 59 wrap objects) + `scripts/msk/hand_complete_evidence.json`.
Build script: `scripts/msk/hand_complete.py`.

---

## 1. What was composed, precisely

**Base**: `data/msk_models/subject2_unified_v2_ringlittle.osim` (read-only; confirmed
byte-size-unchanged before/after, 2,149,577 bytes). **Method reused unmodified from**
`scripts/msk/hand_pulleys.py`: `rescale_lateral()`, `add_wrap_cylinder()`, the exact
`PULLEY_VIA_OFFSET_RATIO=0.5`/`THUMB_VIA_OFFSET_M=0.0045`/`PULLEY_WRAP_RADIUS_RATIO=0.4`
constants, and the finger_plans/thumb_plan construction pattern.

**Touched**: the joint-flanking `PathPoint` locations of `FDP2_r`/`FDS2_r`/`FDP3_r`/`FDS3_r`
(re-tightened — reproduces `pulleys.osim` exactly, see §4), `FPL_r` (re-tightened, same),
and **newly**, `FDP4_r`/`FDS4_r`/`FDP5_r`/`FDS5_r` (ring/little's own extrinsic flexors,
tightened for the first time). Plus 15 `WrapCylinder` objects (9 replicated: mcp2/pip2/dip2,
mcp3/pip3/dip3, cmc1/mcp1/ip1; 6 new: mcp4/pip4/dip4, mcp5/pip5/dip5).

**Not touched**: extensors (`ED2_r`-`ED5_r`, `EPL_r` — real extensors use a structurally
different hood/sagittal-band mechanism, out of scope, matches every prior session's own
disclosed boundary), intrinsics (`LUM1_r`/`DI1_r`/`LUM3_r`/`DI4_r`/`LUM4_r`/`PI3_r`/`ADM_r`),
thumb's other named muscles (`OP_r`/`APB_r`/`FPB_r`/`AdP_r`) — all regression-confirmed
exactly unchanged (§4).

---

## 2. Forced OODA — a real geometric risk, found and fixed BEFORE it reached the model

**Identified before building** (geometric reasoning, then decisively measured — not
assumed either way): `hand_pulleys.py`'s tightening plan moves path-point index 2 (the
"`hand_r`, MCP-side" via point shared by every finger's extrinsic bracket) via a
single-axis decomposition against `[0,-1,0]`: project onto Y ("along"), shrink the WHOLE
(X,Z) remainder ("lateral") to a small target. For index/middle this is safe because their
OWN point-2 X-coordinate is already 0 (measured directly off the model). But
`ring_little_fingers.py`'s own build (`docs/MECHANISM_HAND_FIVE_DIGIT.md` §3) had to give
ring's/little's point-2 a **large, deliberate, non-zero X-shift** (−44.2mm / −68.6mm,
matching that finger's own MCP joint X-translation) as a forced fix for a real bowstringing
bug. That X is a **positional correction identifying which finger this tendon belongs to**,
not a lateral pulley offset — collapsing it during tightening would undo the fix that made
ring's/little's MCP moment-arm sign correct in the first place.

**Decisive test, run BEFORE adopting a design** (`probe_naive_vs_fixed()` in
`hand_complete.py`, throwaway/unsaved models):

| | naive reuse (falsified) | derived fix (adopted) |
|---|---:|---:|
| FDP4_r@mcp4_flex_r sign holds to | **60.3°** of an 84.2° ROM ceiling | **never flips** across the full ROM |
| FDP5_r@mcp5_flex_r sign holds to | **46.3°** of an 85.0° ROM ceiling | **never flips** across the full ROM |
| FDP4_r ma range (mm) | [−2.01, +4.17] (crosses zero) | [+1.51, +8.71] |
| FDP5_r ma range (mm) | [−2.55, +2.95] (crosses zero) | [+1.04, +8.71] |

Naive pulley-tightening would have made little's own sign-holding boundary **regress**
below the un-pulleyed baseline (`docs/MECHANISM_HAND_FIVE_DIGIT.md` §5.7 measured the
un-pulleyed little finger already holding sign to 59°) — i.e. the naive "fix" is actively
**worse than doing nothing**. This is exactly the false-positive/premature-negative trap the
watertight method exists to catch, caught here before it reached the delivered artifact.

**The fix** (`rescale_z_only()` in `hand_complete.py`, geometrically derived, not a new
tunable number): for point-2 only, preserve X (finger identity) and Y (along-metacarpal
bracket spacing) exactly; rescale ONLY Z (`hand_r`'s own volar↔dorsal axis — the genuinely
lateral, Landsmeer-relevant component here) to the SAME already-established target
(`r_hand_m × PULLEY_VIA_OFFSET_RATIO`, identical number `hand_pulleys.py` already uses for
index/middle's own point-2). This is the general form of the exact same operation
`hand_pulleys.py` already performs on index/middle — confirmed: index/middle's own X is
already 0, so the two code paths agree to machine precision there. All OTHER via points
(indices 3-7, on each finger's own dedicated phalanx bodies) have no shared-body collision
and use `hand_pulleys.py`'s unmodified `rescale_lateral()` — byte-identical method to
digits 2/3, just parameterized on ring's/little's own measured distal directions.

Re-measured on the REAL, adopted model: **sign never flips across the full declared ROM**
for FDP4_r@mcp4_flex_r (0-84.2°, range [1.344, 8.706]mm) and FDP5_r@mcp5_flex_r (0-85.0°,
range [1.109, 8.706]mm) — see full curves in the evidence JSON.

---

## 3. Falsifiers — every gate, machine PASS/FAIL

### 3.1 Falsifier (1)+(2): ALL FIVE DIGITS posed simultaneously, ONE geometry
Pose: `wrist_flex_r`=20°, each finger's MCP/PIP/DIP=30/30/20°, thumb `cmc1_flex_r`/
`cmc1_abd_r`/`mcp1_flex_r`/`ip1_flex_r`=20/10/20/10° — ALL SIXTEEN coordinates nonzero at
once (a strictly stronger test than any prior session's own per-finger-only combined pose).

| digit | FDP (wrist/MCP/PIP/DIP, mm) | ED (mm) | verdict |
|---|---|---|---|
| 2 (index) | +17.62/+9.12/+5.73/+4.70 | −6.48/−3.62/−9.94/−9.38 | **PASS** |
| 3 (middle) | +17.36/+5.98/+5.73/+4.70 | −7.05/−2.72/−9.94/−9.38 | **PASS** |
| 4 (ring) | +17.19/+8.56/+5.26/+4.37 | −7.26/−4.01/−9.80/−9.20 | **PASS** |
| 5 (little) | +17.02/+8.56/+5.26/+4.37 | −7.32/−4.01/−9.80/−9.20 | **PASS** |
| thumb (FPL) | wrist +8.19, cmc1 +8.66, mcp1 +5.63, ip1 +4.99 | — | **PASS** |

All FDP/FPL values simultaneously positive, all ED values simultaneously negative, across
**all four fingers + thumb, in one shared pose** — PASS.

### 3.2 Falsifier (3): moment-arm sign across the FULL declared ROM (the pulley property)
23/23 (muscle, coordinate) pairs, single-joint swept 1°-resolution across each coordinate's
own full declared range, all others neutral:

| muscle@coordinate | ROM (deg) | ma range (mm) | sign holds |
|---|---:|---:|---:|
| FDP4_r@mcp4_flex_r | [−22.3, 84.2] | [1.344, 8.706] | **PASS** |
| FDP4_r@pip4_flex_r | [−16.2, 96.0] | [1.180, 5.395] | **PASS** |
| FDP4_r@dip4_flex_r | [−8.1, 82.1] | [2.476, 4.659] | **PASS** |
| FDS4_r@mcp4_flex_r | [−22.3, 84.2] | [1.344, 8.706] | **PASS** |
| FDS4_r@pip4_flex_r | [−16.2, 96.0] | [3.186, 6.055] | **PASS** |
| FDP5_r@mcp5_flex_r | [−26.3, 85.0] | [1.109, 8.706] | **PASS** |
| FDP5_r@pip5_flex_r | [−13.2, 91.8] | [1.937, 5.395] | **PASS** |
| FDP5_r@dip5_flex_r | [−6.8, 84.6] | [2.082, 4.659] | **PASS** |
| FDS5_r@mcp5_flex_r | [−26.3, 85.0] | [1.109, 8.706] | **PASS** |
| FDS5_r@pip5_flex_r | [−13.2, 91.8] | [3.677, 6.055] | **PASS** |

The original 13 (digit2/3/thumb, from `hand_pulleys.py`) **all still hold, unchanged** —
13/13. Combined: **23/23**. Combined "closing the fist" sweep (all joints of a chain moving
together, 19 synchronized steps) also PASS on all **9** chains (FDP2/FDS2/FDP3/FDS3/FDP4/
FDS4/FDP5/FDS5/FPL).

### 3.3 Falsifier (4): lumbrical double action
| muscle | MCP (mm) | PIP (mm) | verdict |
|---|---:|---:|---|
| LUM1_r | +8.891 | −2.257 | **PASS** |
| LUM3_r | +3.101 | −2.247 | **PASS** |
| LUM4_r | +1.993 | −2.247 | **PASS** |

### 3.4 Falsifier (5): ALL moment arms, `computeMomentArm` vs. independent finite-difference
**57/57 PASS** (tolerance: rel. error <1% OR abs. diff <0.01mm, applied uniformly — guards
near-zero crossings, e.g. `DI1_r@pip2_flex_r`, a muscle that structurally doesn't reach that
joint, from a spurious relative-error blow-up near a true value of ~0). Scope: the
pre-registered union of every (muscle, coordinate) crossing any prior session in this repo
already established for that muscle (the 9 pulleyed muscles' own 23 crossings + wrist
crossings for all 5 flexor chains + all 4 fingers' + thumb's extensor sign-contrast controls
+ all 5 intrinsics' own crossings + thumb's other 4 named muscles) — principled and
auditable, not a cherry-picked subset, not an exhaustive sweep over irrelevant pairs either.

### 3.5 Falsifier (6): initSystem + forward step + IK sanity
6a: `equilibrateMuscles()` + 5ms `Manager.integrate()` — all state values finite. **PASS.**
6b: real marker IK on `walking1.trc` vs. reference `walking1.mot` (reused convention from
`merge_unified_v2.py`) — worst unaffected-coordinate RMSE = **1.056°** (`lumbar_bending`),
inside the established 2.0° tolerance. **PASS.** (`walking1.trc` has no hand markers,
already disclosed in `docs/MECHANISM_UNIFIED_V2.md` — this checks that the REST of the body
still assembles/IK-solves correctly with the larger hand subtree attached, not meaningful
hand kinematics; not a new gap.)

### 3.6 Falsifier (7): no duplicate-parent-body
`check_no_duplicate_parent_bodies()` reused unmodified from `merge_unified_v2.py`: `{}`
(zero dupes). **Independently re-derived a second time** (fresh `Counter` over
`Joint.getChildFrame().findBaseFrame()`, not just imported) directly on the delivered file
in a separate process: `{}` — confirmed both ways. **PASS.**

### 3.7 Structural-count arithmetic (pre-registered BEFORE building)
| | ringlittle (base) | hand_complete (composed) | expected | pass |
|---|---:|---:|---:|---|
| bodies | 44 | 44 | 44 (+0) | PASS |
| joints | 44 | 44 | 44 (+0) | PASS |
| muscles | 257 | 257 | 257 (+0) | PASS |
| wrap objects | 44 | 59 | 59 (+15 = 9 replicated + 6 new) | PASS |

Zero new bodies/joints/muscles — pulley-tightening only relocates existing path points and
adds wrap objects, exactly mirroring `hand_pulleys.py`'s own precedent on its own base.

### 3.8 Regression vs. `ringlittle.osim` base
**20/20 PASS** (diff <0.001mm, same-process comparison): every untouched muscle/joint
(`LUM1_r`, `DI1_r`, `ED2_r`-`ED5_r`, `EPL_r`, `OP_r`, `APB_r`, `FPB_r`, `AdP_r`, `LUM3_r`,
`DI4_r`, `LUM4_r`, `PI3_r`, `ADM_r`) reproduces the base model exactly.

### 3.9 External anchor: composed model vs. the independently-built `pulleys.osim`
The strongest available anchor for the "reused" half of this method: `pulleys.osim` is a
SEPARATE, already-independently-verified file (13/13 sign-hold PASS in its own session).
**13/13 PASS, diff 0.000000mm** — the composed model's digit-2/3/thumb moment arms
reproduce `pulleys.osim`'s own saved values to machine precision, confirming the extension
did not perturb the already-certified portion at all (not a tautology — an external,
independently-built artifact).

### 3.10 Round-trip, tenodesis, determinism
Round-trip (fresh reload): **PASS**, all 23 pairs <0.05mm. Tenodesis (path-length
extended>flexed + tendon-force extended>flexed): **PASS, 4/4** — FPL_r (325.9→317.2mm),
FDP2_r (387.0→358.6mm), **new**: FDP4_r (374.4→347.6mm), FDP5_r (369.3→342.8mm).

**Determinism — a forced self-audit, reported in full.** This session's own FIRST
determinism check FAILED (4/23 pairs passed). Not accepted on sight (a one-shot fail is
premature surrender, not honest negative): diagnosed decisively that the check itself
compared a sweep-grid sample nearest to 0° (typically ±0.1-0.5°, pose-dependent) against a
separate build's EXACT 0.0° pose — an apples-to-oranges pose mismatch, not a real
non-determinism (confirmed: the 4 "passing" pairs were EXACTLY the 4 whose sweep grid
happened to land on 0.000°). Corrected two ways: (a) explicit identical exact-zero pose on
both models, (b) md5 of the saved `.osim` after two independent OS-level process builds —
**byte-identical**, md5 `b057a3c79198cf36e0b50b94512740ac` both times. `hand_complete.py`'s
own check is fixed in-place so this will not recur silently on a future re-run.

---

## 4. Honest scope — unchanged, disclosed gaps (symmetric, not hidden)

Carried forward from every parent doc, unchanged in kind: segment lengths remain
generic/literature-typical; Fmax/optimal-fiber-length remain physiologically-plausible
placeholders (except the 5 intrinsics' Fmax, sourced from Nakajima et al. 2022); the
extensor mechanism (hood/sagittal bands), collateral ligaments, MCP ab/ad DOF, and ~10 more
real intrinsics remain absent (`docs/MECHANISM_HAND_FIVE_DIGIT.md` §6); the 9 replicated wrap
cylinders (digit 2/3/thumb) were already measured to be a passive/non-engaging safety layer
behind the tightened via points (`docs/MECHANISM_HAND_PULLEYS.md` §1.3) — the 6 new ring/
little wrap cylinders are structurally present with the same analytically-exact
axis/position construction but were **not** independently re-measured for engagement this
session (the via-point tightening is the established, already-verified load-bearing
mechanism; not re-litigated here, disclosed as an unchecked-but-precedented assumption, not
silently assumed). The secondary, stronger external-magnitude-band bar (not just sign)
partially misses near full flexion for the ORIGINAL digit2/3/thumb pairs, unchanged from
`hand_pulleys.py`'s own disclosed §5 finding (not chased further, p-hacking guard, same as
that session) — digit 4/5's own per-digit Nakajima bands were not re-swept against the
FULL ROM curve this session (only checked at the reference pose, §"external_band_check" in
the evidence JSON), a narrower check than digit 2/3's own; flagged, not chased.

A cosmetic, harmless, pre-existing OpenSim quirk: `GeometryPath` subcomponents named
"pathwrap" get auto-renamed (`pathwrap_0`, `pathwrap_1`, …) whenever a muscle carries more
than one `PathWrap` (e.g. `FDP2_r` crosses 3 pulleys) — confirmed present identically
whether or not this session's own edits are applied; does not affect any structural count,
moment arm, or falsifier (all passed); not a defect, not chased.

---

## 5. Files

- `scripts/msk/hand_complete.py` — the build + verification (loads
  `subject2_unified_v2_ringlittle.osim` read-only, imports and reuses `hand_pulleys.py`'s
  `rescale_lateral`/`add_wrap_cylinder`/constants, `merge_unified_v2.py`'s
  `check_no_duplicate_parent_bodies`/IK-smoke-test helpers, `thumb_column.py`'s
  `measure_validity_boundary`/`tenodesis_check`, `ring_little_fingers.py`'s
  `EXTERNAL_BAND_MM`; adds the new `rescale_z_only()` fix + the pre-adoption
  naive-vs-fixed probe). Run:
  `source_repository/.venv-msk/bin/python3 scripts/msk/hand_complete.py`
- `scripts/msk/hand_complete_evidence.json` — full machine-measured evidence: every
  full-ROM curve, the combined-all-digits-pose falsifier, the lumbrical check, all 57 FD
  cross-checks, the 20-pair regression, the 13-pair external anchor vs. `pulleys.osim`, the
  duplicate-parent-body result, structural counts, round-trip, IK sanity, tenodesis, the
  naive-vs-fixed diagnosed gap, and the corrected determinism self-audit.
- `data/msk_models/subject2_hand_complete.osim` — the composed model (44 bodies / 44
  joints / 257 muscles / 61 coordinates / 59 wrap objects). md5
  `b057a3c79198cf36e0b50b94512740ac` (confirmed byte-identical across 2 independent
  builds). All 3 source forks (`subject2_unified_v2_ringlittle.osim`,
  `subject2_unified_v2_pulleys.osim`, `subject2_unified_v2_thumb.osim`) untouched.
- Read directly this session (reused unmodified): `scripts/msk/hand_pulleys.py`,
  `scripts/msk/ring_little_fingers.py`, `scripts/msk/anatomical_hand.py`,
  `scripts/msk/full_hand.py`, `scripts/msk/thumb_column.py`, `scripts/msk/merge_unified_v2.py`,
  `scripts/msk/merge_unified_model.py`, `docs/MECHANISM_HAND_PULLEYS.md`,
  `docs/MECHANISM_HAND_FIVE_DIGIT.md`, `docs/MECHANISM_UNIFIED_V2.md`.

## 6. Roadmap (updated)

1. **Fold into the whole-body unified model.** This composed hand is a fork of
   `subject2_unified_v2_ringlittle.osim` (hand+arm only, per `docs/MECHANISM_HAND_FIVE_DIGIT.md`'s
   own lineage), not yet grafted onto `subject2_unified_v2.osim` (the full-body 34-body
   unified model that already carries hip/ankle/spine/knee ligaments). A future session
   could repeat `merge_unified_v2.py`'s own graft pattern to fold this hand in.
2. **Re-measure the 6 new ring/little wrap cylinders' engagement** (currently assumed, not
   independently re-checked this session, unlike the original 9 which were measured
   non-engaging in `docs/MECHANISM_HAND_PULLEYS.md` §1.3).
3. **Sweep digit 4/5's own external bands across the FULL ROM curve**, not just the
   reference pose (mirrors digit 2/3/thumb's own already-completed, stronger check).
4. Every unchanged gap from `docs/MECHANISM_HAND_FIVE_DIGIT.md` §6/§8 and
   `docs/MECHANISM_HAND_PULLEYS.md` §5/§7 (extensor mechanism, MCP ab/ad DOF, the ~10
   remaining intrinsics, Buchholz et al. 1992's segment-length regression).
