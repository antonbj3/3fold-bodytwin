# CLOUD-W2-A363-BATCH: Independent A363 batched NNLS algorithm audit

Use only public data bundled here or downloaded from primary public sources and original code written for this lane. Do not open other repositories or use restricted model data, TLEM, the collaborator data, LHDL or internal A matrices. Prior register numbers below are unreviewed claims, not input observations. Distinguish proof, synthetic test, public empirical evidence, and UNKNOWN. Give exact commands, units, source URL/DOI, uncertainty and counterexamples. Do not invent unavailable data or references. No push, publishing, email, credentials or subagents. Final answer: print RESULTS.md, results.json and small code exactly between =====FILE <relative path>===== and =====END FILE===== markers (total <=200 KB).

## Prior evidence and strongest baseline
A363 register claims a batched GPU solver around 499 steps/s at 1000 individuals with agreement ~5.8e-9. No internal systems or matrices are allowed.

## Bounded work
Implement original batched nonnegative convex quadratic/NNLS solver and scalar reference on synthetic positive, correlated, rank-deficient and changing-active-set systems. Report primal/dual KKT residual, objective gap, timing with hardware and scaling; include warm-start and failure cases. Use public NumPy/SciPy only or explain alternatives.

Deliver RESULTS.md, results.json and original runnable code; state exact commands, sample sizes, baselines, deviations and limitations.
