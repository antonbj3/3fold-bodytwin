# X24 CBCT HU CALIBRATION

Review: LANE_XREVIEW_BATCH3, ACCEPTED_WITH_CORRECTION. Reviewed result SHA256: `ca9e36c25374dab4cf7345d57ddc6d1308f4cf7858024dcb40b139b7d213b899`.

Scope: "All 9 explicitly requested deliveries; README/results/raw comparison, one independent headline recount each, copied run_all replay, existing fault injections, source identity/quantity, resolution and readable limits. No new design/model/experiment. X24 raw 22-volume remeasurement and X30 full345s atlas not rerun; cached raw-value guards and default scripts rerun. Prior7 imported reviews excluded from rates."

Corrections: [{"job": "LANE_X24_CBCT_HU_CALIBRATION", "target_id": "DENT-MAT-BONE-HU-MODULUS", "claim_as_written": "Python 3, NumPy, SciPy and Matplotlib needed", "finding": "Accurate equality for recalculated polyfit coefficients fails with NumPy 1.21.5 despite only rounding difference; NumPy 2.2.6 passerar.", "evidence_locator": ["@DENTAL_IMPLEMENTATIONS@/X24_CBCT_HU_CALIBRATION/run_r2.py:55", "@DENTAL_INPUT_ROOT@/artifacts/LANE_XREVIEW_BATCH3/evidence/X24_COEFFICIENT_PROBE_1.21.5.json"], "corrected_statement": "Lock the existing environment, check the anchor accurately and allow only declared coefficient rounding; then use frozen coefficients. +1 in coefficient/anchore is folded.", "class": "ACCEPTED_WITH_CORRECTION", "scope": "Executable copy/runtime contract", "conclusion_change": "Numerical headlines unchanged; copy/control guarantee corrected"}]

Physical validation: UNKNOWN. Code and source tests are retained for inspection. Full original replay requires externally supplied data and optional dependencies; see ../../docs/DATA.md. Portable profiles, where available, use declared synthetic or published aggregate inputs.
