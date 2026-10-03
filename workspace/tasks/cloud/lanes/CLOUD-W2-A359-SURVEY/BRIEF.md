# CLOUD-W2-A359-SURVEY: NHANES A359 survey-design robustness

Use only public data bundled here or downloaded from primary public sources and original code written for this lane. Do not open other repositories or use restricted model data, TLEM, the collaborator data, LHDL or internal A matrices. Prior register numbers below are unreviewed claims, not input observations. Distinguish proof, synthetic test, public empirical evidence, and UNKNOWN. Give exact commands, units, source URL/DOI, uncertainty and counterexamples. Do not invent unavailable data or references. No push, publishing, email, credentials or subagents. Final answer: print RESULTS.md, results.json and small code exactly between =====FILE <relative path>===== and =====END FILE===== markers (total <=200 KB).

## Prior evidence and strongest baseline
Wave 1 adjusted DXA increment: 84.43 to 82.46 N, 2.34% holdout improvement; descriptive prediction among completers. Eight public CDC XPT files and hashes are bundled.

## Bounded work
Independently implement weight-aware training and evaluation with NHANES strata/PSU and cycle weights. Avoid combining incompatible weights without stating assumption. Compare weighted and unweighted adjusted models, paired uncertainty by PSU, and sex/age strata. Audit coding, merge, imputation and leakage. Report exact estimands.

Deliver RESULTS.md, results.json and original runnable code; state exact commands, sample sizes, baselines, deviations and limitations.
