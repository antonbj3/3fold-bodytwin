# XSRC VERIFY

Review: LANE_XREVIEW_BATCH35, ACCEPTED_WITH_CORRECTION. Reviewed result SHA256: `307f8d6c2c7e3f29430b584ad7e5e09c06b4f66661b4cb8eba2a48dfd42b9345`.

Scope: "80-bundle/429-field accounting, 38/113 attrition, main source-table probes in all six source-job strata, consumer transfer exclusions, existing sufficiency witnesses and offline replay. This is not individual independent semantic certification of all 80 bundles."

Corrections: "Apply portable.patch plus scoped_source_gate.py. Use package-local recipient snapshots and raw/SOURCE_DRAFT.jsonl. The original copied run failed at the external draft directory. Number-presence control accepted C04 occlusal/axial swap; labelled primary-table gate rejects it. Printed values and consumer conclusions are unchanged."

Physical validation: UNKNOWN. Code and source tests are retained for inspection. Full original replay requires externally supplied data and optional dependencies; see ../../docs/DATA.md. Portable profiles, where available, use declared synthetic or published aggregate inputs.
