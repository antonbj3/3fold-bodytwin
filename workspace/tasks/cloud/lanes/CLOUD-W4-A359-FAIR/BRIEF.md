# CLOUD-W4-A359-FAIR: Strength prediction cohort transport

Only public primary data or original synthetic code. Do not access private repositories, restricted model data, TLEM, the collaborator data, LHDL or internal A matrices. Do not use subagents, push, email, publish or enable paid usage. Cite sources actually checked, including exact URL/DOI and which value they support. Clearly distinguish synthetic method checks from empirical validation. If required public data are unavailable, report UNKNOWN and do not invent values. Deliver RESULTS.md, results.json and small runnable code, and print files in final message between =====FILE <path>===== and =====END FILE=====, total <=200 KB.

## Bounded work
Use the bundled public NHANES CDC XPT archive, verified by MANIFEST.sha256, and the bundled baseline code a359.py. Do not change the frozen cohort or split. Compare A359 adjusted mass+height+sex+age baseline with DXA across sex and age cells; quantify eligibility and missingness. Do not alter holdout split after inspection.

Report baseline, independent check, exact command, units, uncertainty, failure cases and limitations.
