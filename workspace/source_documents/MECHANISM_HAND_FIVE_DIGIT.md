# MECHANISM HAND: ALL FIVE DIGITS — RING + LITTLE FINGERS BUILT (2026-07-21)

Executes the operator's explicit spec: "add the ring + little fingers so the twin has a
full anatomical hand." Builds digit-4 (ring) and digit-5 (little) the SAME way index/middle
were built (`scripts/msk/anatomical_hand.py` / `scripts/msk/full_hand.py`): forearm-origin
extrinsic flexors (FDP/FDS) crossing wrist→MCP→PIP→DIP, an extensor (ED) on the dorsal side,
one lumbrical + one interosseous each, plus the little finger's own hypothenar abductor
digiti minimi (ADM). Loads `data/msk_models/subject2_unified_v2_thumb.osim` (38 bodies/246
muscles: index+middle+thumb) **read-only**, writes a **new fork file**
(`data/msk_models/subject2_unified_v2_ringlittle.osim`, deliberately distinct from a
separate, concurrent pulley-thread's own output). Every number below is machine-measured
this session (`scripts/msk/ring_little_fingers.py`), **independently re-verified in a
separate, fresh Python process** (not just the build script's own self-report), and
**reproduced byte-for-byte identically** across two independent full runs (`md5sum` on both
the evidence JSON and the saved `.osim` file).

## Headline result

| | |
|---|---:|
| **Ring (digit 4) extrinsics built?** | **YES** — `FDP4_r`, `FDS4_r`, `ED4_r` (forearm-origin, `ulna_r`), mirroring the already-verified index/middle method |
| **Little (digit 5) extrinsics built?** | **YES** — `FDP5_r`, `FDS5_r`, `ED5_r` |
| **Ring intrinsics?** | **YES** — `LUM3_r` (3rd lumbrical), `DI4_r` (4th dorsal interosseous, real anatomical numbering — confirmed by this session's own fetched external table) |
| **Little intrinsics?** | **YES** — `LUM4_r` (4th lumbrical), `PI3_r` (3rd palmar interosseous — the anatomically-correct choice, since little has no dorsal interosseous of its own), **`ADM_r`** (abductor digiti minimi, the task's own explicit hypothenar ask) |
| **Task's explicit falsifier: FDP4_r/FDP5_r simultaneously positive at wrist+MCP+PIP+DIP, ONE geometry (the same multi-joint test FPL_r passed)?** | **PASS** — Ring: `+17.19/+10.38/+9.51/+9.28mm`; Little: `+17.02/+10.38/+9.51/+9.28mm`, all four simultaneously positive, both fingers |
| ED4_r/ED5_r negative at the same 4 joints (sign-contrast)? | **PASS** — Ring: `−7.26/−4.02/−9.80/−9.20mm`; Little: `−7.32/−4.02/−9.80/−9.20mm`, all negative |
| LUM3_r/LUM4_r textbook double action (flexes MCP, extends PIP)? | **PASS** — LUM3: MCP `+2.39mm`, PIP `−1.25mm`; LUM4: MCP `+1.31mm`, PIP `−1.25mm` |
| DI4_r / PI3_r / ADM_r flex their own MCP? | **PASS, all 3** — `+4.52mm` / `+3.46mm` / `+6.71mm` |
| Moment-arm cross-check (`computeMomentArm` vs. independent finite-difference), all 11 new muscles | **34/34 PASS** (rel. error ≤6.0×10⁻³, most ≤10⁻⁴) |
| Sign convergence | **1 iteration, 0 flips** — every joint axis and muscle side correct on the first attempt |
| External anchor (Nakajima et al. 2022, real **per-digit** table — a genuine upgrade over the 2-anchor precedent index/middle's own intrinsics had) | **18/23 PASS**; 5 disclosed misses, all traced to a pre-existing, already-accepted method property (§4), not a new defect |
| Regression: does this break index/middle/thumb anchors? | **PASS — 19/19 unchanged**, reproduced independently to ≤0.0005mm |
| Symmetric QC: `initSystem()` + forward step + round-trip reload | **PASS** — 5ms forward step all-finite; two independent fresh reloads agree to <0.0001mm |
| Structural counts | 38→44 bodies, 38→44 joints, 246→257 muscles, 370→381 forces, 55→61 coordinates — **exact match to pre-registration**, independently re-counted |
| Determinism | 2 independent fresh runs → **byte-for-byte identical** evidence JSON and model file (same MD5: `8f34b17d...`/`6bbc806b...`) |
| Base model untouched | `subject2_unified_v2_thumb.osim` — **2,057,010 bytes before and after**, confirmed via independent `stat` + `md5sum`, not just self-report |

**Confidence tier: cadaveric/published-plausibility.** Every new muscle's moment-arm sign
AND magnitude is cross-checked against a real, freshly-fetched, **per-digit** external table
(Nakajima et al. 2022) covering all 11 new muscles — a stronger evidentiary base than the
2-anchor precedent the index/middle intrinsics had. Not in-vivo-anchored — no live
grip-force instrument exists in this pipeline. Stated plainly.

**Model + evidence**: `data/msk_models/subject2_unified_v2_ringlittle.osim` (new fork file,
base `subject2_unified_v2_thumb.osim` untouched — confirmed via independent `stat`/`md5sum`)
+ `scripts/msk/ring_little_fingers_evidence.json`.

---

## 1. Real, freshly-fetched external data (raw HTML grep, not AI-summarized — this repo's
established method), cross-validated against 7 already-known numbers before trusting new rows

### 1.1 ROM: Mohamed Ibrahim BK et al. 2024 (PMID 39345665, PMC11436331), Table 2
The SAME already-cited source as index's/middle's own ROM (its own title covers "Index,
Middle, Ring, and Little Fingers"). Before trusting the new Ring/Little rows, the fetch was
cross-checked against the ALREADY-known Index (MCP 86.0/25.7°, PIP 97.2/13.7°, DIP
81.6/6.0°) and Middle (MCP 86.6/24.7°, PIP 96.2/15.6°, DIP 84.5/8.4°) numbers — **matched
exactly, byte-for-byte identical to the values already cited in `anatomical_hand.py`/
`full_hand.py`** — before reading off the new rows:

| finger | MCP flex/ext (°) | PIP flex/ext (°) | DIP flex/ext (°) |
|---|---:|---:|---:|
| Ring | 84.2 / 22.3 | 96.0 / 16.2 | 82.1 / 8.1 |
| Little | 85.0 / 26.3 | 91.8 / 13.2 | 84.6 / 6.8 |

### 1.2 Moment-arm + Fmax anchor: Nakajima et al. 2022 (Scientific Reports, PMC9345905), Table 1
The SAME table already trusted for FDS2/FDP2/EDC2/1LU/1DI. Column parsing re-verified
against **7** already-known numbers (FDS2 MP=11.9/PIP=13.6, FDP2 MP=10.2/PIP=−8.5/DIP=4,
EDC2 MP=−9.4/PIP=0/DIP=−1.6, 1LU MP=9.6/PIP=+2.6 [the already-disclosed
textbook-contradicting sign], 1DI MP=4.4/PIP=0) — **all matched exactly** before the new
Ring/Little rows were trusted. This table's own finger-groups **independently confirm** the
anatomical structure this build uses: Ring = FDS4/FDP4/EDC4/**3LU**/2PI/**4DI**; Little =
FDS5/FDP5/EDC5/EDM/**ABDM**/FDMB/**4LU**/**3PI**/ODM — i.e. little has **no dorsal
interosseous of its own** (confirming `PI3_r`, not a "5DI", is the anatomically-correct
choice) and little's real hypothenar set is **exactly** {ABDM, FDMB, ODM} — matching the
task's own explicit naming verbatim, found here as an independent confirmation, not assumed.

| muscle | Fmax (N) | MP Flex-Ext | MP Add-Abd | PIP Flex-Ext | DIP Flex-Ext |
|---|---:|---:|---:|---:|---:|
| FDS4 | 2.0 | +9.9 | −1.2 | +5.0 | 0 |
| FDP4 | 3.0 | +8.9 | −0.8 | +6.2 | +4.1 |
| EDC4 | 1.7 | −8.1 | −0.5 | −2.4 | −1.2 |
| 3LU | 0.1 | +6.6 | −7.5 | **−2.0** (textbook-correct, unlike 1LU) | −1.5 |
| 4DI | 1.7 | +4.7 | +7.1 | −1.2 | −0.9 |
| FDS5 | 0.9 | +8.6 | +3.2 | +4.7 | 0 |
| FDP5 | 2.8 | +8.5 | +4.0 | +5.9 | +3.2 |
| EDC5 | 0.9 | −4.9 | +0.9 | −2.6 | −1.3 |
| 4LU | 2.0 | +6.3 | +7.2 | −2.2 | −1.7 |
| 3PI | 1.0 | +2.1 | +7.7 | −2.0 | −1.5 |
| ABDM | 1.4 | +4.7 | −8.0 | −2.0 | −1.5 |

Fmax column **flagged, not adopted, for the 6 extrinsics** (SAME already-established reason:
this table's own FDP2 Fmax was already 2.7N vs. this repo's 45N placeholder — an internal
inconsistency flagged once already in `thumb_column.py`, reconfirmed here for FDP4/FDP5
too). **Adopted directly for the 5 intrinsics** (mirrors the `LUM1_FMAX_N`/`DI1_FMAX_N`
precedent) — `LUM3_r`=0.1N, `DI4_r`=1.7N, `LUM4_r`=2.0N, `PI3_r`=1.0N, `ADM_r`=1.4N.

---

## 2. Forced fix #1: a genuine OpenSim/Simbody segfault, root-caused (not geometry)

The first full run **segfaulted** (`dumped core`) partway through verification. A one-shot
"crash, give up" would be premature — forced OODA:

- **Observe**: bisected with `faulthandler` — the crash is inside
  `Muscle.computeMomentArm()`, specifically when called on a **pre-existing** muscle
  (`FDP2_r`) at a **pre-existing** coordinate (`wrist_flex_r`), on the in-memory `Model`
  object returned by `build()`.
- **Orient (the crux)**: is this a geometry bug in this build, or something else? A
  decisive test: reload the **just-saved** `.osim` file **fresh** (`osim.Model(path)`,
  one `initSystem()` call) in the SAME process, while the old, crash-prone object still
  exists — `FDP2_r@wrist_flex_r` computes cleanly (**14.6411mm**, matching the
  already-established anchor 14.6413mm) with **zero crash**. This proves the saved model is
  structurally sound — the crash is an in-memory staleness property of a `Model` object
  that has been through **multiple** `initSystem()` calls during construction (bodies+joints,
  then muscles, then tendon-slack-length solve) — a genuine, if obscure, OpenSim/Simbody
  quirk (an already-documented class of fragility in this codebase's own history —
  `add_scapula_clavicle.py`'s CustomJoint/SpatialTransform segfault precedent — now a second,
  distinct instance of "don't trust a long-lived Model object across repeated
  `initSystem()` calls").
- **Decide/Act**: all POST-convergence verification (moment-arm cross-check, falsifiers,
  band checks, regression check, validity sweep, structural counts) now runs on a **fresh
  reload** of the just-saved model file, never the in-memory object that survived multiple
  `initSystem()` calls. The sign-convergence loop itself is unaffected (confirmed safe by
  construction: each candidate rebuild only ever queries newly-added muscles/coordinates on
  that build's own single most-recent `initSystem()` state).

This is a reusable methodological finding for any future multi-step build in this codebase:
**verify moment arms on a fresh file reload, not a long-lived in-memory Model that has
called `initSystem()` more than once.**

---

## 3. Forced fix #2: a real geometric bug caught by the falsifier itself, root-caused

At the pre-registered combined falsifier pose (all 6 new coordinates nonzero at once, plus
`wrist_flex_r`=20°), `FDP4_r`/`FDP5_r`'s own MCP moment arm came out **negative**
(`−0.11mm`/`−1.41mm`) — failing the task's own explicit falsifier. A bald accept of "this is
just the disclosed un-pulleyed limitation" would have been premature — forced OODA:

- **Observe**: sweeping MCP alone (PIP/DIP fixed) showed the SAME monotonically-decreasing
  "bowstringing" signature in **all four fingers**, including the already-verified
  index/middle — ruling out a qualitatively new mechanism. Sweeping PIP alone (MCP fixed)
  showed **zero** effect on the MCP moment arm for any finger — ruling out cross-axis
  coupling as the cause (unlike the thumb's own CMC1-flex/CMC1-abd cross-coupling issue).
- **Orient (the crux)**: ring's/little's own **baseline** (0°) MCP moment arm was itself
  anomalously small (5.33mm/3.23mm vs. index's 14.87mm/middle's 10.98mm). Quantified root
  cause: the `hand_r`-side via points sat at `x=0` regardless of finger, while ring's/
  little's own MCP joint sits at `x=−47.9mm`/`−72.2mm` (2/3 adjacent-metacarpal-spacing
  units) — a mismatch that **dwarfs** the intended lever's own ~14.4mm scale, so the
  tendon's direction at the MCP crossing was dominated by this large, physiologically
  arbitrary x-travel rather than the intended volar/dorsal lever.
- **Decide**: give the `hand_r`-side via points the SAME finger-specific x-offset as that
  finger's own MCP joint translation — real flexor tendons fan out toward their own finger
  as they cross the palm, they do not all leave from one shared point. A geometric
  correction, not a fitted parameter.
- **Act**: verified numerically on an isolated single-muscle probe BEFORE adopting: `+1.82mm`
  (unshifted) → `+12.67mm` (shifted) at MCP=20° — a clean, decisive improvement, not a
  marginal tweak. Applied **only** to the 6 extrinsic muscles (the ones with this
  two-point-on-`hand_r` bracket construction) — the 5 intrinsics (`LUM3_r`/`DI4_r`/`LUM4_r`/
  `PI3_r`/`ADM_r`) use a single `hand_r`-to-finger segment and were measured to **not** show
  this sign-flip at the same deep pose, so left unchanged (mirrors `LUM1_r`/`DI1_r`'s own
  precedent — not fixing what isn't broken).

After the fix, the falsifier passes cleanly (§5.1) and the validity-sweep (§6) shows the
sign now holds to 58-59° of MCP flexion (of an 84-85° ROM ceiling), comfortably covering the
falsifier's own 30° test pose.

---

## 4. The 5 disclosed external-band misses — forced, not hand-waved

18/23 per-digit Nakajima band checks pass. The 5 misses (`ED4_r@pip4/dip4`,
`ED5_r@pip5/dip5`, `LUM4_r@mcp5`) were **not** accepted on sight — forced check:

**ED4_r/ED5_r@PIP/DIP** (measured 9.2-9.8mm vs. bands of 0.3-9.1mm): is this new, or does it
already exist in the already-verified `ED2_r`? Live cross-check on the composed model:
`ED2_r@pip2_flex_r` = **−10.63mm**, which would ALSO fail a comparably tight per-muscle band
derived from the SAME table's own `EDC2` anchor (DIP=−1.6mm → ceiling ≈5.6mm). **Confirmed,
not assumed**: this is a pre-existing, shared property of the "extensor mirrors flexor's
exact via-point construction" method used since the ORIGINAL index-finger build — the real
extensor mechanism (central slip, lateral bands) is explicitly not modeled anywhere in this
codebase (a disclosed gap since `anatomical_hand.py`'s very first session), and building it
would be new-scope, not "the same method as index/middle." Accepted and disclosed, not a
regression, not chased.

**LUM4_r@mcp5_flex_r** (measured 1.31mm vs. band floor 1.575mm, a 17% shortfall): traced to
the SAME already-anticipated mechanism disclosed in this script's own module docstring
*before* measuring — little sits 3 spacing-units from `hand_r`'s unshifted intrinsic origin
(vs. ring's 2), and the intrinsics' `hand_r`-side points were deliberately left unshifted
(§3). A small, honest, pre-anticipated shortfall — sign independently confirmed correct
(the headline double-action falsifier passes), only the magnitude undershoots a wide,
pre-registered band. Not chased further (continuing to search for a scale/offset that lands
inside a self-chosen band is the exact p-hacking failure mode this method forbids —
established precedent: `thumb_column.py`'s own accepted `APB_r` magnitude miss).

---

## 5. Verification — every gate machine PASS/FAIL

### 5.1 The task's own explicit falsifier (independently reproduced in a separate process)
At the ONE combined pose (`wrist_flex_r`=20°, all 6 new finger coordinates nonzero
simultaneously — no per-joint retuning):

| muscle | wrist_flex_r | MCP | PIP | DIP |
|---|---:|---:|---:|---:|
| **FDP4_r** (ring) | +17.186mm | +10.376mm | +9.508mm | +9.278mm |
| **ED4_r** (ring) | −7.260mm | −4.015mm | −9.801mm | −9.197mm |
| **FDP5_r** (little) | +17.024mm | +10.376mm | +9.508mm | +9.278mm |
| **ED5_r** (little) | −7.322mm | −4.015mm | −9.801mm | −9.197mm |

**All FDP values simultaneously positive, all ED values simultaneously negative — PASS**,
for both fingers, independently reproduced in a fresh process (values above, second column
onward, match the build script's own report to <0.001mm).

### 5.2 Lumbrical double action + interossei/hypothenar flexion
| muscle | MCP | PIP | verdict |
|---|---:|---:|---|
| **LUM3_r** (ring) | +2.387mm | −1.251mm | **PASS** (flexes MCP, extends PIP) |
| **LUM4_r** (little) | +1.306mm | −1.251mm | **PASS** (flexes MCP, extends PIP) |
| **DI4_r** (ring) | +4.522mm | 0.000mm | **PASS** (flexes MCP, doesn't reach PIP) |
| **PI3_r** (little) | +3.456mm | 0.000mm | **PASS** |
| **ADM_r** (little) | +6.714mm | 0.000mm | **PASS** |

### 5.3 Moment-arm cross-check (`computeMomentArm` vs. independent finite-difference)
**34/34 muscle-coordinate pairs PASS** (rel. error ≤6.0×10⁻³, most ≤10⁻⁴) — every crossing
of all 11 new muscles, at the same combined pose as §5.1.

### 5.4 Sign convergence
**1 iteration, 0 flips** — every joint axis (`mcp4_flex_r`/`pip4_flex_r`/`dip4_flex_r`/
`mcp5_flex_r`/`pip5_flex_r`/`dip5_flex_r`) and every muscle side (`FDP4`/`FDS4`/`ED4` group,
`FDP5`/`FDS5`/`ED5` group, `LUM3_r`, `DI4_r`, `LUM4_r`, `PI3_r`, `ADM_r`) correct on the
first attempt — the elbow-axis-transport method continues to get every sign right when
applied consistently, matching the established precedent's own experience.

### 5.5 External anchor (Nakajima et al. 2022) — 18/23 PASS, 5 disclosed (§4)
| muscle@coordinate | measured \|mm\| | band (mm) | pass |
|---|---:|---|---|
| FDP4_r@mcp4 | 10.376 | [2.2, 31.2] | PASS |
| FDP4_r@pip4 | 9.508 | [1.6, 21.7] | PASS |
| FDP4_r@dip4 | 9.278 | [1.0, 14.3] | PASS |
| FDS4_r@mcp4 | 10.376 | [2.5, 34.7] | PASS |
| FDS4_r@pip4 | 10.723 | [1.25, 17.5] | PASS |
| ED4_r@mcp4 | 4.015 | [2.0, 28.3] | PASS |
| ED4_r@pip4 | 9.801 | [0.6, 8.4] | **FAIL (disclosed, §4)** |
| ED4_r@dip4 | 9.197 | [0.3, 4.2] | **FAIL (disclosed, §4)** |
| LUM3_r@mcp4 | 2.387 | [1.65, 23.1] | PASS |
| LUM3_r@pip4 | 1.251 | [0.5, 7.0] | PASS |
| DI4_r@mcp4 | 4.522 | [1.175, 16.45] | PASS |
| FDP5_r@mcp5 | 10.376 | [2.125, 29.75] | PASS |
| FDP5_r@pip5 | 9.508 | [1.475, 20.65] | PASS |
| FDP5_r@dip5 | 9.278 | [0.8, 11.2] | PASS |
| FDS5_r@mcp5 | 10.376 | [2.15, 30.1] | PASS |
| FDS5_r@pip5 | 10.723 | [1.175, 16.45] | PASS |
| ED5_r@mcp5 | 4.015 | [1.225, 17.15] | PASS |
| ED5_r@pip5 | 9.801 | [0.65, 9.1] | **FAIL (disclosed, §4)** |
| ED5_r@dip5 | 9.197 | [0.325, 4.55] | **FAIL (disclosed, §4)** |
| LUM4_r@mcp5 | 1.306 | [1.575, 22.05] | **FAIL (disclosed, §4)** |
| LUM4_r@pip5 | 1.251 | [0.55, 7.7] | PASS |
| PI3_r@mcp5 | 3.456 | [0.525, 7.35] | PASS |
| ADM_r@mcp5 | 6.714 | [1.175, 16.45] | PASS |

### 5.6 Regression (symmetric QC — did this break index/middle/thumb?)
**19/19 PASS**, self-measured baseline (on the untouched loaded base model, before any
addition) vs. post-build, within the SAME script run (more watertight than a hardcoded
anchor dict — no dependency on a stale external number), reproduced again independently in
a separate process:

| muscle@coordinate | baseline (mm) | post-build (mm) | diff (mm) |
|---|---:|---:|---:|
| FDP2_r@wrist_flex_r | 14.6412 | 14.6416 | 0.000434 |
| FDP2_r@mcp2/pip2/dip2_flex_r | 12.9238/12.0191/10.4956 | (unchanged) | 0.000000 |
| ED2_r@mcp2/pip2_flex_r | −4.7557/−10.6305 | (unchanged) | 0.000000 |
| LUM1_r@mcp2/pip2_flex_r | 8.8911/−2.2570 | (unchanged) | 0.000000 |
| DI1_r@mcp2_flex_r | 6.3194 | (unchanged) | 0.000000 |
| FDP3_r@wrist_flex_r | 14.1951 | 14.1955 | 0.000430 |
| FDP3_r@mcp3/pip3/dip3_flex_r | 6.6045/11.5300/10.5650 | (unchanged) | 0.000000 |
| FPL_r@wrist/cmc1_flex/mcp1/ip1 | 4.4373/17.6031/15.5020/15.1973 | (unchanged) | 0.000000 |
| APB_r@cmc1_abd_r | 1.3706 | (unchanged) | 0.000000 |
| AdP_r@cmc1_abd_r | −24.6290 | (unchanged) | 0.000000 |

The two ~0.0004mm diffs (FDP2_r/FDP3_r@wrist_flex_r) are the SAME ordinary cross-process
floating-point non-reproducibility in OpenSim's own multibody solver already diagnosed as
benign in `docs/MECHANISM_FULL_HAND.md` Sec.4.4 — not a real regression.

### 5.7 Validity sweep (forced, disclosed — mirrors the thumb column's own precedent)
The same un-pulleyed via-point angular-range limitation already disclosed for the thumb
column applies here too (measured, not assumed):

| finger | FDP flexor sign holds to | ED extensor sign holds to | ROM ceiling |
|---|---:|---:|---:|
| Ring | 58° | 68° | 84.2° |
| Little | 59° | 64° | 85.0° |

Both comfortably cover the falsifier's own 30° test pose and the fetched real-population ROM
midpoint. A separate, concurrent pulley-thread is adding `WrapCylinder`/`WrapSurface`
geometry to widen this range for the whole hand at once — this build does **not** duplicate
that effort, only measures and discloses the same already-known limitation.

### 5.8 Machine artifact checks
- **Determinism**: 2 independent fresh-process runs → **byte-for-byte identical** evidence
  JSON (MD5 `6bbc806b8a2434241742ef06112e885f`) and `.osim` file (MD5
  `8f34b17d57cec94c7862f3b521487ba0`).
- **Round-trip**: two independent fresh reloads of the saved model agree on
  `FDP4_r@mcp4_flex_r`/`FDP5_r@mcp5_flex_r` to <0.0001mm (10.3757mm both times).
  Independently re-confirmed in a separate process too.
- **Structural counts**: 38→44 bodies (+6: `ring_proximal_r`/`ring_medial_r`/
  `ring_distal_r`/`little_proximal_r`/`little_medial_r`/`little_distal_r`), 38→44 joints
  (+6), 246→257 muscles (+11), 370→381 forces total (+11), 55→61 coordinates (+6) — exact
  match, independently re-counted in a separate process, not just the build script's own.
- **Base model untouched**: `subject2_unified_v2_thumb.osim` — 2,057,010 bytes before AND
  after (identical), confirmed via independent `stat` + `md5sum`, not just self-report.
- **Forward dynamics smoke test**: `equilibrateMuscles()` + 5ms `Manager.integrate()` on the
  fully composed model — all state values finite.
- **No segfault** on the fix (§2) — verified 3 separate full runs post-fix, all clean exits.

---

## 6. Honest scope — what this IS and is NOT (symmetric disclosure)

**Confidence tier: cadaveric/published-plausibility, not in-vivo-anchored.** What remains,
unchanged in kind from every prior session's own disclosures, plus this session's own new
gaps:

- **All five fingers + thumb now exist structurally** — this IS "a full anatomical hand" in
  the sense the task named (5 digits, each with extrinsic flexor/extensor + at least one
  intrinsic; little additionally has a hypothenar abductor). It is **not** a complete
  muscular hand: ~10 more real intrinsics are still missing (middle finger's own 2nd
  lumbrical/2nd+3rd dorsal interossei, ring's own 2nd palmar interosseous, little's own
  flexor/opponens digiti minimi beyond ADM, the remaining palmar interossei) — unchanged
  gap, disclosed in every prior session too.
- **MCP ab/ad DOF is still absent for every finger** — so `DI4_r`'s/`PI3_r`'s/`ADM_r`'s real
  PRIMARY action (abduction/adduction) is not modeled or tested; only their real, secondary
  MCP-flexion contribution is (mirrors `DI1_r`'s own already-established simplification
  exactly, now applied to 3 more interossei/hypothenar muscles).
- **Segment lengths reuse the same generic index/middle values unchanged** for both ring
  and little (disclosed simplification, mirrors `full_hand.py`'s own middle-finger
  precedent) — real ring segments are close to index/middle's; real little segments are
  known to be noticeably shorter (a robust qualitative fact) but no citable per-digit
  regression was fetched this session (Buchholz et al. 1992 remains unretrieved, the same
  persistent gap across every hand session in this repo).
- **The 5 new intrinsics' `hand_r`-side origin points do not carry the finger-specific
  lateral (x) spacing offset** (§3) — mirrors `LUM1_r`/`DI1_r`'s own established
  simplification on index exactly — more schematic ("diagonal") for the more-ulnar digits,
  a real, disclosed, quantified consequence (§4's `LUM4_r` miss), not hidden.
- **5 external-band misses** (§4), all forced-diagnosed and traced to pre-existing,
  already-accepted method properties (the un-modeled extensor mechanism; the intrinsics'
  un-shifted hand_r origin) — not new defects, not chased further (p-hacking guard).
- **Tendon pulleys (A1-A5), collateral ligaments, the true extensor mechanism** remain
  absent — unchanged from every prior session's disclosed gap list. A separate,
  **concurrent** thread is adding `WrapCylinder`/`WrapSurface` pulleys — this build does
  **not** duplicate that effort.
- **Ring's interosseous uses only the dorsal one** (`DI4_r`); **little's uses only the
  palmar one** (`PI3_r`) — matches the task's explicit "one interosseous each" ask; the
  OTHER real interosseous each finger anatomically also has (ring's own 2PI, e.g.) is not
  built.
- **Fmax for the 6 new extrinsics are physiologically-plausible round placeholders**
  (45N/40N/20N, matching FDP2/FDP3's own convention), **not** this table's own (flagged as
  inconsistent) numbers. **Fmax for the 5 new intrinsics ARE real, sourced numbers**
  (Nakajima et al. 2022) — a genuine improvement in anchoring, consistent with the
  established precedent.
- **A genuine, root-caused OpenSim/Simbody segfault** (§2) was found and fixed at the
  process level (fresh reload for verification) — a reusable lesson for any future
  multi-`initSystem()` build in this codebase.

---

## 7. Files

- `scripts/msk/ring_little_fingers.py` — the build + verification (loads
  `subject2_unified_v2_thumb.osim` read-only, adds ring + little, verifies via a fresh
  reload, saves a new fork file). Run:
  `source_repository/.venv-msk/bin/python3 scripts/msk/ring_little_fingers.py`
- `scripts/msk/ring_little_fingers_evidence.json` — full machine-measured evidence (superset
  of every table above: sign-convergence history, all 34 moment-arm cross-checks, all 23
  external-band checks including the isolated-pose-equivalent detail, all 19 regression
  pairs, both fingers' validity-boundary sweeps, structural counts, round-trip values).
- `data/msk_models/subject2_unified_v2_ringlittle.osim` — the new model file (44 bodies / 44
  joints / 257 muscles / 381 total force elements). `data/msk_models/
  subject2_unified_v2_thumb.osim` (the base) was loaded read-only and is unmodified —
  confirmed via independent `stat`/`md5sum` (2,057,010 bytes, unchanged).
- Read directly this session: `scripts/msk/anatomical_hand.py`, `scripts/msk/full_hand.py`,
  `scripts/msk/thumb_column.py`, `scripts/msk/add_scapula_clavicle.py` (all reused
  unmodified — `make_hinge`, `cylinder_body`, `set_range`, `measure_elbow_axis_ground`,
  `R_of`, `lateral_axis_local`, `build_muscle_from_points`, `moment_arm_fd`, `_set_all`,
  `R_RATIO`, `OFFSET_TO_RADIUS_RATIO`), `docs/MECHANISM_HAND_FOREARM_THUMB.md`,
  `docs/MECHANISM_FULL_HAND.md`, `docs/MECHANISM_ANATOMICAL_HAND.md`.

## 8. Roadmap (updated)

1. **Tendon pulleys** (`WrapCylinder`/`WrapSurface`) — owned by a separate, concurrent
   thread; once merged, re-run this build's own validity sweep (§5.7) to confirm the
   angular-range limitation is resolved for all five fingers at once.
2. **The remaining ~10 intrinsics** (middle's own 2LU/2DI/3DI, ring's own 2PI, little's own
   FDMB/ODM beyond ADM, the palmar interossei not yet built) — mechanically closer to a
   repeat of this session's own method now that both archetypes (single-segment
   interosseous/hypothenar, cross-to-dorsal lumbrical) are established for 3 of 4 fingers.
3. **MCP ab/ad DOF** for all fingers — needed to exercise the interossei's/ADM's real
   PRIMARY action (unchanged gap from every prior session's own roadmap item).
4. **Resolve the 5 disclosed band misses' root cause structurally** (a real extensor-hood
   sub-model for ED2-5_r; a finger-specific x-shift extended to the intrinsics too) — not
   attempted here per the p-hacking guard; a future session's explicit, scoped task.
5. **Buchholz et al. 1992's real segment-length regression** — still unretrieved, same
   compounding gap across every hand session so far, now affecting 4 of 5 fingers.
