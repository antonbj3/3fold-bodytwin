# CX-MOTIONSCAN — leave no stone unturned in 3fold-motion-engine: what in it can solve BodyTwin's open nodes? (read-only inventory; nothing built)

Anton: "there may be more in the motion engine… I'm counting on no stones being left unturned."
Search: `../3fold-motion-engine/` (src/, docs/, docs/papers P1–P4 + reproducibility/, tests/, examples/, `_private/romi_collab/HUNT.md` and build/ — read-only). Also grep `~/projects/3fold_staging/` for other engine repos with the same relevance.

## BodyTwin's open nodes (match each module against these)
1. Knee contact: two layers/non-parabolic gap/patch coupling (A1400 KNEEMERGE2: a geometric error gate is needed), force-controlled contact, secondary kinematics (Lemma G, A1402).
2. Recruitment / the gauge: what the body chooses (CX-INVERSEOC), a stability/stiffness constraint, co-contraction.
3. Determined quantities / identifiability / OED (A1384, A1394, A1400).
4. Speed for population/what-if: certified surrogate + fallback (A1379 EMUCERT), batched GPU (A363), dual in joint space.
5. Derivatives through contact events (A369, A907; certified gradient radius).
6. Individual geometry: the CT↔marker/fluoro frame (one missing DOF, W6-FRAME-GAUGE), patella tracking/rollback (A912, A1380).
7. Uncertainty/certificates: intervals, determinism, f32 bands.
8. Soft tissue / cartilage / ligaments: materials, thin layers, MPM.
9. Sensors: IMU/video → motion/load.

## Deliverables in `results/CX-MOTIONSCAN/`
- `MOTION_MAP.md` (starting with `# CX-MOTIONSCAN`): for each relevant module or result, give the path, what it does, its status (tested with a number / partial / idea), which node 1–9 it touches, and how BodyTwin would use it (concretely: which function, which input/output). State which ones are ALREADY in use (grep bodytwin results/ and bodytwin_core/).
- A ranking of the top 10 by (expected value for BodyTwin) / (integration cost), with one line of justification each.
- Explicitly: "not relevant" for large parts, with a search log (a negative requires a search log).
- `motion_map.json`: the same data, machine-readable.
No builds, no computation over 60 s, no writes outside `results/CX-MOTIONSCAN/`. lane_runner has full permissions to read.
