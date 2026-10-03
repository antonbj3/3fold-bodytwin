# LANE_SWARM_FINDINGS_REVIEW

A stronger model reviews the night's most promising constructive swarm results (The swarm/swarm_worker, new source families and surgical families). Weaker models' findings do not count as "holding" until verified; the XSEED review showed that roughly half of the swarm's claims were PARTIAL or WRONG. Results folder `results/LANE_SWARM_FINDINGS_REVIEW/`.

## Kandidater

`tasks/build_night/swarm/TOP_CONSTRUCTIVE_20261001_0600.txt` — 25 non-negative results ranked by research value (jobb-ID | Objectives | source keys | beslut). Jobbkataloger: `/mnt/games-240/research/bunny48_20260926/bodytwin/<ID>/` (results.json, RESULTS.md, kod). Prioritera: kirurgiska kopplingar (SURG_INCISION × SURG_HEMOSTASIS: vessels per cut surface, wall shear, percentage of vessels not self-sealing), Q080×Q100, Q077×Q088 (ATP-syntas → oxygen supply), Q049→Q154 (seriekonduktans), Q009×Q154 (fri/total-bindning med efflux), Q058×Q084.

## Per resultat

1. Rerun the job's code; does it reproduce the decisive number?
2. Derive the central claim yourself (dimensions, conservation, limiting cases). Check that the source family's model does not carry one of the errors found by LANE_SOURCE_CONSERVATION_AUDIT (AUDIT_TABLE.json) — Q005, Q012, Q052, Q090, Q115, Q146, Q168, MITOSTRESS have known errors.
3. Verdict: CONFIRMED / PARTIAL / WRONG / INDETERMINATE with a recalculated number.

## Leverans

REVIEW_TABLE.json (job → claim → verdict → recalculated number → file:line), RESULTS.md for readers: which swarm findings hold and can become new ports in BodyTwin, and which must be redone. Better 8 thoroughly reviewed than 25 superficially.
