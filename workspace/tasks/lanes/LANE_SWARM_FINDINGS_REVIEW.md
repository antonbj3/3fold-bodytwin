# LANE_SWARM_FINDINGS_REVIEW

A stronger model reviews the night's most promising constructive swarm results (The swarm/swarm_worker, new source families and surgical families). Weaker models' findings do not count as "holding" until verified; the XSEED review showed that roughly half of the swarm's claims were PARTIAL or WRONG. Results folder `results/LANE_SWARM_FINDINGS_REVIEW/`.

## Candidates

`tasks/build_night/swarm/TOP_CONSTRUCTIVE_20261001_0600.txt` — 25 non-negative results ranked by research value (job ID | target | source keys | decision). Job directories: `external_mount<ID>/` (results.json, RESULTS.md, code). Prioritise: surgical couplings (SURG_INCISION × SURG_HEMOSTASIS: vessels per incision area, wall shear, fraction of vessels that do not self-seal), Q080×Q100, Q077×Q088 (ATP synthase → oxygen supply), Q049→Q154 (series conductance), Q009×Q154 (free/total binding with efflux), Q058×Q084.

## Per result

1. Rerun the job's code; does it reproduce the decisive number?
2. Derive the central claim yourself (dimensions, conservation, limiting cases). Check that the source family's model does not carry one of the errors found by LANE_SOURCE_CONSERVATION_AUDIT (AUDIT_TABLE.json) — Q005, Q012, Q052, Q090, Q115, Q146, Q168, MITOSTRESS have known errors.
3. Verdict: CONFIRMED / PARTIAL / WRONG / INDETERMINATE with a recalculated number.

## Delivery

REVIEW_TABLE.json (job → claim → verdict → recalculated number → file:line), RESULTS.md for readers: which swarm findings hold and can become new ports in BodyTwin, and which must be redone. Better 8 thoroughly reviewed than 25 superficially.
