# XBREAK HUNT 5

Review: LANE_XREVIEW_BATCH35, ACCEPTED_WITH_CORRECTION. Reviewed result SHA256: `39f4ca5e07a8226dd37799821de4e03407f38d1a4a04bf3e5d937bd168eb8b4e`.

Scope: "All11 supplied test preparations in R1 negative cone result, two explicitly selected R3/R4 pilots, stored STL witness checks and corrected continuous-gap replay. Physical scale, whole crown and any unsupplied same-summary sufficiency test are not admitted."

Corrections: "Apply verify.patch to recompute whole-source-surface lower/upper gap bounds against saved STL. Original verifier accepted upper_max=0 below measured positive lower_max; patched verifier rejects it and unchanged baseline passes. Replay currently reads producer source VTP/NPZ/STL paths outside the bundle. No standalone replay, required exact-summary sufficiency witness or completed README at the cutoff; source units stay UNKNOWN. target_id is absent."

Physical validation: UNKNOWN. Code and source tests are retained for inspection. Full original replay requires externally supplied data and optional dependencies; see ../../docs/DATA.md. Portable profiles, where available, use declared synthetic or published aggregate inputs.
