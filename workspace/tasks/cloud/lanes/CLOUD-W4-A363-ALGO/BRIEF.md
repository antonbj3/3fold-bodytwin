# CLOUD-W4-A363-ALGO: Independent batched NNLS algorithm audit

Only public primary data or original synthetic code. Do not access private repositories, restricted model data, TLEM, the collaborator data, LHDL or internal A matrices. Do not use subagents, push, email, publish or enable paid usage. Cite sources actually checked, including exact URL/DOI and which value they support. Clearly distinguish synthetic method checks from empirical validation. If required public data are unavailable, report UNKNOWN and do not invent values. Deliver RESULTS.md, results.json and small runnable code, and print files in final message between =====FILE <path>===== and =====END FILE=====, total <=200 KB.

## Bounded work
Derive a batched projected active-set or block principal-pivot method for original synthetic NNLS, benchmark against SciPy NNLS on 1000 seeded systems at 4 conditioning levels, including rank defects. Report KKT, objective and CPU throughput scaling. A363 private GPU throughput is context only.

Report baseline, independent check, exact command, units, uncertainty, failure cases and limitations.
