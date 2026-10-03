# CX-PROGRESSWAVE — 80 The swarm + 30 reserve_worker packets that PRODUCE PROGRESS for BodyTwin (new ability or new knowledge about the body), no reanalysis

Anton: "we're supposed to make progress here". The previous waves filled the queue with second-order questions about our own results (reanalyses of DXA, KKT gaps, certificate status, recalculations). Those are banned here.

## Every packet must do ONE of these
- **Build**: a new, small, working model component with a test against a public reference. Examples:
  - a patella tracking function;
  - a disc pressure model;
  - a tendon stiffness curve;
  - a cartilage layer;
  - a cell/tissue transport module.
- **Measure/extract**: a new quantity from real data we have that is not yet used. Examples:
  - rollback from the Grand Challenge fluoro trials (stairs, step-up, openfe, twist);
  - knee angle–moment from isokinetics;
  - muscle CSA from the VSD CT;
  - GRF features from OpenCap.
- **Test against the world**: compare a model with a published/public reference that we have not tested against. The DOI source must actually be accessible (Wilke disc pressure, Krevolin patella arm, Bergmann OrthoLoad in new activities, OpenCap EMG).

## Areas (spread)
- **Knee chain:** patella/rollback/operating range/extensor mechanism; cartilage in cooperation with Field U384 (only as a consumer).
- **Spine/L5–S1:** disc pressure, lever arms, intra-abdominal pressure.
- **Individual signals:** strength at several angles, CT muscle geometry, body shape from video/scanning.
- **Sensors:** IMU/video → joint load with public data.
- **Hyperrealism roadmap** (`results/HYPERREALISM_ROADMAP_20260924/`: HYPERREALISM_PLAN, IDEA_CATALOG, CELL_GEOMETRY_*, SPATIAL_CELL_SIMULATION): bounded cell/tissue/organ modules with a public reference, e.g. muscle energetics, tendon adaptation, bone remodelling under load, perfusion.

## Banned
- recomputing/reanalysing our own results;
- questions about certificates/registers/audits;
- the DXA estimand;
- KKT gaps;
- D1 details;
- "is X usable when Y is UNKNOWN";
- packets that can only establish that data is missing.

## Format
Same as CX-DSWAVE:
- a complete `inputs/` with DATA_SUFFICIENCY and a load test;
- the brief must NOT refer to files outside the packet;
- The swarm packets ≤ 60 s of computation;
- reserve_worker packets somewhat heavier.
- IDs: `BT-PG-001` … `BT-PG-110`.

**Do NOT queue.** Write `results/CX-PROGRESSWAVE/QUEUE_PROPOSAL.txt` (one line: `<A|B|C> <swarm|reserve_worker> BT-PG-xxx | area | one-sentence question | build/measure/test`). The coordinator reads the list and queues it.

lane runner has full permissions in the workspace; ~/projects/bodytwin is read-only.
