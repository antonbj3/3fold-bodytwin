# CX-INTEGRATE — gather what demonstrably works into ONE usable BodyTwin package with tests and a single end-to-end command

Plan (`~/research/PLAN_BODYTWIN_FIELD_8H_20260924.md` §3, §9): 35 % development/integration, "results introduced into the right repo", "used via the project's intended interface". That is what has been missing. Everything that works is spread across ~100 result directories.

## Package `bodytwin_core/` in the workspace (`./bodytwin_core/`)
Only components with audited or reproduced evidence. For each module, copy the code with a source reference (file + sha256); do not reimplement:
1. `geometry` — geomgr (results/N7c/geomgr: instantiation, stable ID, parameter editing ops, certificate v2), A194/A183. The fixed transformation to the reference model export.
2. `solver` — the exact full-body solver (results/CX-SOLVER2/cx_solver.py, A295) + the implicit KKT derivative (results/CX-WHATIF2/implicit_kkt.py, A369) + scaling robustness from CX-SOLVER-ROBUST if it is ready (otherwise a clear TODO).
3. `whatif` — parameter → hip force with derivative (CX-WHATIF2 + CX-D1PARITY, A907; U380 is Field's code: call it via an adapter with a source reference, do not copy it in without a reference).
4. `emulator` — per active set (BT-XF4, A303, audited HOLDS).
5. `population` — N2b batched (A363) + the band (N50/N2b).
6. `nullmodels` — B24, N1, N1g (the reference everything is measured against).

## Requirements
- A single command: `python -m bodytwin_core.demo` runs on 1–2 individuals with public data or the collaborator's internal data (locally): instantiate → solve the collaborator's box lift → what-if CCD/AV with derivative → band → compare against the null models. Timing per step.
- pytest: one test per module that reproduces the key number from the source's results.json (tolerance stated). Everything green; skipped tests are counted and reported (no hidden skips).
- `bodytwin_core/README.md` (≤ 60 lines): what each module does, its evidence (A-number, audit), known limitations (for example: knee force — N1g beats muscle models, A312/A324/A358; the DXA strength law was corrected, A359). No trust rhetoric; verb + object + number.
- restricted model data/the collaborator data must not be copied into the package. Read it from its current paths (a config file). The package must be runnable without internal data for the public parts.

Do not write in `~/projects/bodytwin` (read-only tonight). Merging the package into the original repo is a separate step.

Deliverables: the package, pytest output, the demo log with timing, and `results/CX-INTEGRATE/RESULTS.md` starting with `# CX-INTEGRATE`.
