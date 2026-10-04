# MECHANISM — Fiber-length correction: causal test of the muscle-architecture root cause

> **Salvage note:** the agent completed the run + wrote the evidence JSON
> (`data/msk_smoketest/subject2_walking1/fiber_length_correction_test/fiber_length_correction_test_results.json`,
> `scripts/msk/fiber_length_correction_test.py`) but was stopped (a mistaken coordinator kill during a CPU-starvation
> stall) right before writing this doc. Authored by the coordinator from the verified JSON.

**Question:** `MECHANISM_QUAD_HAM_FASCICLE.md` found the model's optimal fiber lengths are systematically too long
(+0.8 to +2.9 SD vs Ward 2009, PMID 18972175). Does *causally* correcting them (Ward2009 values, pennation +
Fmax held at baseline) reduce the knee/hip over-prediction and fix the tendon-strain/fascicle symptoms?

## Result: fiber-length is a MAJOR knee lever (over-corrects), WORSENS the hip, does NOT fix tendon strain
| joint | baseline ratio vs OrthoLoad | fiber-length corrected | Fmax-only | **fiber-length + Fmax compound** |
|---|---:|---:|---:|---:|
| **Knee** | 1.515× (391 %BW) | **0.945×** (244 %BW) — *over-corrects below in-vivo* | 1.255× | **0.921×** |
| **Hip** | 1.412× | (worsens) | 1.417× | **1.757× — WORSE** |

- **Knee:** correcting fiber-length alone drops contact force **391 → 244 %BW** (ratio 1.515 → 0.945) — it more
  than closes the gap, over-shooting *below* OrthoLoad. Compounded with the Fmax correction: **0.921×**. So the
  knee over-prediction is **fully closable (indeed over-closable) by muscle-architecture calibration** (fiber
  length + Fmax) — strong confirmation that over-long fibers + over-strong Fmax are the knee's real, coupled cause.
- **Hip:** fiber-length correction makes it **worse** (compound 1.757× vs baseline 1.412×) — opposite direction.
  Consistent with Fmax also doing nothing at the hip: the hip over-prediction is **not** a muscle-architecture
  issue; it is the glute-med recruitment issue (`MECHANISM_HIP_STRUCTURAL_MECHANISM.md`). Knee and hip have
  genuinely different root causes.
- **Tendon strain:** barely moves (Achilles 2.50 → 2.63 %, still below the 4–6 % in-vivo band) — fiber-length
  correction does **not** fix the low-tendon-strain symptom; that is a separate (tendon force/stiffness) axis.
- Convergence: the compound-corrected solve **passes** (158 frames, pelvis residual 173 N) — not a
  de-convergence artifact.

## What it means (symmetric-QC)
Using Ward2009 values *directly* (not tuned) **over-corrects** the knee (0.945×, below in-vivo) — so this is
**not** a clean "fix"; it demonstrates that fiber-length is a *high-gain* lever whose true value likely sits
between the model's too-long values and the cadaveric reference, and that it interacts strongly with Fmax
(compound 0.92×, i.e. correcting both over-shoots). The honest reading: **the knee over-prediction is dominantly
a muscle-architecture-calibration problem (fiber length + Fmax), and it is correctable — but the exact
calibration is under-determined** (over-correction with cadaveric values; do not tune to the anchor). The hip is
a different problem entirely. **Confidence tier:** cadaveric-anchored (Ward2009) correction × in-vivo-anchored
(OrthoLoad) outcome. **Gap:** over-correction means the true fiber-length is not pinned; single-subject.
