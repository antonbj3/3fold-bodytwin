# CLOUD-W2-A359-MISSING: NHANES A359 selection and missingness audit

Use only public data bundled here or downloaded from primary public sources and original code written for this lane. Do not open other repositories or use restricted model data, TLEM, the collaborator data, LHDL or internal A matrices. Prior register numbers below are unreviewed claims, not input observations. Distinguish proof, synthetic test, public empirical evidence, and UNKNOWN. Give exact commands, units, source URL/DOI, uncertainty and counterexamples. Do not invent unavailable data or references. No push, publishing, email, credentials or subagents. Final answer: print RESULTS.md, results.json and small code exactly between =====FILE <relative path>===== and =====END FILE===== markers (total <=200 KB).

## Prior evidence and strongest baseline
Wave 1 replicated a 13.1% mass-only improvement with public CDC 1999–2002 data, but the stronger adjusted baseline gained only 2.34%; original 19.1% did not reproduce. Eight CDC XPT files and hashes are bundled.

## Bounded work
Using the bundled data, quantify age>=50 inclusion, DXA and strength missingness by age/sex/cycle; inverse-probability weighted sensitivity and pessimistic pattern-mixture bounds for the adjusted DXA increment. Keep train 1999–2000 and untouched 2001–02 holdout. Report force in N and assumed torque arm separately. Compare adjusted baseline and adjusted+DXA on identical rows.

Deliver RESULTS.md, results.json and original runnable code; state exact commands, sample sizes, baselines, deviations and limitations.
