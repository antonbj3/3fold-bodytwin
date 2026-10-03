# CLOUD-W6-A363-COLSCALE: Column-scaled NNLS stopping counterexample

Read SOURCE_NOTES.md when present. Only public primary data or original code. No restricted model data, TLEM, the collaborator/LHDL data, private A matrices or private repositories. No subagents, Git push, email, publishing or paid usage. Cite exact checked source URLs/DOIs and what each supports. Distinguish synthetic method checks from empirical/anatomical validity. If required data are unavailable, say UNKNOWN without fabricating. Deliver RESULTS.md, results.json and runnable code; print small files between =====FILE <path>===== and =====END FILE=====, total <=200 KB.

## Bounded task
Wave 4 A363 audit passed 1000/1000 on normalized matrices, but only 31/100 when columns spanned 1e-6 to 1e6. Build original seeded benchmark comparing original stopping, column-normalized solver and per-column dual tolerances with independent SciPy/BVLS references. Freeze seeds and tolerance before comparing.

Report strongest baseline, source/version, exact commands, units, uncertainty, counterexamples and limitations.
