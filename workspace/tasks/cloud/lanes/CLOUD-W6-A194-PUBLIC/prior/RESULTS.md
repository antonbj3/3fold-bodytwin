# CLOUD-W4-A194-ARMS — Moment-arm geometry independent audit

## Verdicts (PREREG.md, sha256 d6f016ea…2daa, unchanged)

| Criterion | Verdict |
|---|---|
| Empirical (public versioned geometry + independent measured moment arms, matched by definition, improvement in ≥4/5 muscle groups) | **UNKNOWN**: licensed public geometry exists (OpenSim Gait2354, CC BY 3.0, demonstration model), but no matched independent measured moment arms were tested |
| Synthetic method (virtual-work error < 0.1 % away from switches) | **PASS**: 5/5 synthetic groups, 25 configurations |

**Only the synthetic method check passed. Nothing here is empirical validation.** All numbers come from original synthetic geometry and a model-free comparison of three estimators. None of them was compared with measured anatomy.

## Why the empirical criterion is UNKNOWN
**Correction.** An earlier version of this file said no explicitly licensed public geometry was found. That was wrong. My guessed Gait2392 path returned 404, but that was my wrong guess, not evidence that such a model is missing. A licensed public model exists and I have now inspected it:

- **https://raw.githubusercontent.com/opensim-org/opensim-models/master/Models/Gait2354_Simbody/gait2354_simbody.osim**
  - Fetched 2026-09-24: HTTP 200, 387 781 bytes, sha256 `b96156dcd721777c4d1ac5fbc2706fac7ee23be64c9c1669dcca39556a416702`.
  - The repo's `master` branch was at commit `d9b05d470b1a481c222372c85b75772faf8f7792` (from `git ls-remote`).
  - `OpenSimDocument Version="40000"`.
  - **License:** the `<credits>` element says "License: Creative Commons (CCBY 3.0) … http://creativecommons.org/licenses/by/3.0/".
  - **Scope:** the same `<credits>` element says the number of muscles "was reduced by Anderson to improve simulation speed for demonstrations and is **not intended to be used in research**".
  - **Lineage** (from `<credits>`/`<publications>`): Delp et al. 1990, IEEE TBME 37:757–767; planar knee of Yamaguchi & Zajac 1989, J Biomech 22:1–10; Anderson & Pandy 1999/2001.
  - **Contents (counted by grep):** 54 `Thelen2003Muscle`, 12 `CustomJoint`, 12 `ConditionalPathPoint`, 4 `MovingPathPoint`, and no wrap geometry (13 `WrapObjectSet` elements, none holding a wrap object).

Other sources checked:
- `.../master/README.md` (200): describes the repo; the README itself has no license text.
- `.../master/LICENSE.txt` and `.../master/LICENSE` (404).
- `.../master/Models/Gait2392_Simbody/gait2392_simbody.osim` (404): wrong guessed path.
- `https://api.github.com/repos/opensim-org/opensim-models/contents/...` (403): directory listing blocked.

**Corrected reason for UNKNOWN.** Licensed, versioned public geometry exists (Gait2354, CC BY 3.0). The criterion still cannot be met, for three reasons:
1. I did not obtain or test any **independent measured moment arms** (cadaver or imaging) that are matched by anatomical definition to this model's muscles and joint coordinates.
2. I computed **no moment arms from Gait2354**. That would need its CustomJoint splines, MovingPathPoints and ConditionalPathPoints, which my synthetic code does not implement, and this follow-up was limited to no new heavy computation.
3. The PREREG does not define the baseline that "improvement" is measured against.


The inspection is still relevant to the synthetic audit: Gait2354 uses 12 ConditionalPathPoints, the same switch type as failure case 1 below. Near those switches, finite-difference moment arms from this model would be unreliable.

## Method (synthetic)
The joint is a hinge with angle q. The attachment Q is fixed to the distal body. Three estimators are compared:
- **VW (virtual work):** r = −dL/dq. The derivative is exact, computed by forward-mode automatic differentiation (dual numbers) through the path-length code, including atan2 and arccos in the wrap.
- **FL (force line):** r = axis · ((Q − j) × u), where u is the unit direction of the last path segment at Q. This uses positions only.
- **FD:** central difference with h = 1e-5 rad.

Five synthetic "muscle groups". Each has 5 attachment configurations: a base plus 4 seeded jitters of ±1 cm (seed 194). These are the "controlled attachments":
G1 straight planar; G2 via point on the distal body; G3 cylinder wrap centred on the joint (R = 25 mm); G4 cylinder wrap offset from the joint (R = 18 mm); G5 straight 3-D path with an oblique axis and an off-origin joint centre.
The grid is q = −120…120° in 0.1° steps. Settings fixed in the code before the run: 0.5° exclusion band around each switch (switches found by bisection to 1e-12 rad), and relative error |x − ref| / max(|ref|, 1 µm).

## Results (moment arms in mm)
| Group | worst rel. err VW vs FL | worst rel. err VW vs FD | worst abs VW vs FD (mm) | wrap switches (deg) |
|---|---|---|---|---|
| G1 straight | 1.8e-13 | 2.8e-7 | 4.0e-9 | none |
| G2 via point | 4.7e-13 | 8.9e-7 | 6.0e-9 | none |
| G3 centred cylinder | 1.4e-15 | 2.5e-10 | 7.5e-9 | −4.4, −9.6, −4.1, −7.7, −4.3 |
| G4 offset cylinder | 5.6e-15 | 1.4e-9 | 7.8e-9 | −3.2, −5.3, −5.3, −7.3, +0.4 |
| G5 oblique 3-D | 2.5e-13 | 1.1e-7 | 4.5e-9 | none |

Closed-form check: in the wrapped region of G3, VW = R = 25.000 mm (maximum absolute difference 0.0 mm).

**Uncertainty.** VW and FL agree to floating-point level (≤5e-13 relative). The VW–FD differences (≤9e-7 relative, ≤8e-9 mm) are consistent with O(h²) truncation plus round-off. The FD step study at G4, q = 0 shows the error falling as h² down to h = 1e-5, then rising to about 1e-7 mm at h ≤ 1e-7. The PASS margin is more than 1000× below the 0.1 % threshold. This conclusion does not depend on the seed or the grid. It does depend on the geometry being planar-cylinder or straight: spheres, ellipsoids and multi-object wraps were not tested.

## Failure cases / counterexamples
1. **Conditional via point** (OpenSim-style, active for −30° ≤ q ≤ 40°). L jumps at the switch. A central FD that straddles the switch gives 160 050 mm against VW = 31.3 mm. VW and FL are one-sided and stay finite, but the moment arm itself is discontinuous. No estimator gives a meaningful value *at* the switch.
2. **Wrap lift-off with a coarse FD step.** At 1e-4 rad past the G4 lift-off, the FD error is 0.21 mm at h = 0.01 rad and 0.017 mm at h = 1e-3 rad. The moment arm is continuous here, but its derivative is not. A coarse-step FD can therefore produce sub-mm moment-arm discrepancies that are purely numerical. That is the same order as the A194 context values.
3. **Exploratory correction (does not change the criterion).** I first coded the G3 closed-form reference as −R. It showed a 50 mm (= 2R) discrepancy, which was a sign-convention error in the reference, not in the estimators. It was corrected to +R. The PASS rule does not use the closed form.

## A194 context (not used for any verdict)
The private values 0.39 mm and 0.77 mm were not accessed. As an illustration only: in synthetic G1 at q = 30°, moving the insertion 1 mm changes r by 0.35 mm (x) or 0.72 mm (y). A sub-mm moment-arm disagreement is therefore within typical attachment-digitisation uncertainty. Coarse FD near a wrap switch can also produce it (failure case 2).

## Limitations
- Planar hinge and single-DOF paths only. No coupled knee kinematics, sphere or ellipsoid wraps, or multi-wrap paths.
- The "muscle groups" are synthetic labels, not anatomy.
- Gait2354 was inspected for license and structure only. It was not simulated.
- VW and FL agreeing is expected from theory: for a geodesic path, dL/dQ = −u. Their agreement checks the implementation, not the physiology.
- "f attachments" in BRIEF.md was read as controlled attachment perturbations. This is an interpretation.

## Reproduce
```
pip install numpy          # tested with Python 3.11.15, numpy 2.4.6
python3 arms_audit.py      # ~10 s; writes results_synthetic.json, prints summary
```
Files: `arms_audit.py` (sha256 16b9fbd6…7c7c), `results_synthetic.json` (raw per-configuration output), `results.json` (summary), `run_log.txt` (stdout).
