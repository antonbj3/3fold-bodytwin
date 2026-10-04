# CLOUD-W4-A907-PARITY: Implicit KKT derivative parity

Only public primary data or original synthetic code. Do not access private repositories, restricted model data, TLEM, the collaborator data, LHDL or internal A matrices. Do not use subagents, push, email, publish or enable paid usage. Cite sources actually checked, including exact URL/DOI and which value they support. Clearly distinguish synthetic method checks from empirical validation. If required public data are unavailable, report UNKNOWN and do not invent values. Deliver RESULTS.md, results.json and small runnable code, and print files in final message between =====FILE <path>===== and =====END FILE=====, total <=200 KB.

## Bounded work
Independently derive the derivative of a strongly convex equality constrained quadratic force solver with respect to geometry. Construct a seeded 2-body synthetic case with changing mass and wrench, compare implicit derivative with converged central finite differences at 5 step sizes, and give a cancellation/error budget. A907 reports 2.6e-8 on its private model, which is context only.

Report baseline, independent check, exact command, units, uncertainty, failure cases and limitations.
