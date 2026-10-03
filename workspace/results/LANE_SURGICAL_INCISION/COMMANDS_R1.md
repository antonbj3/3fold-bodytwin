# Executed commands and reproduction

All numerical commands were run from the project root with `PYTHONDONTWRITEBYTECODE=1` and OMP/OPENBLAS/MKL/NUMEXPR_NUM_THREADS=2. Python3, NumPy, SciPy; no GPU/bigmem/cloud job. No private code was executed as a script. HX modules are imported read-only with bytecode disabled.

```bash
python3 results/LANE_SURGICAL_INCISION/incision_r1.py --tag first --sizes 32,64,96
python3 results/LANE_SURGICAL_INCISION/nonlinear_exterior_r1.py
python3 results/LANE_SURGICAL_INCISION/incision_r1.py --tag physical_patch --sizes 32,64,96
python3 results/LANE_SURGICAL_INCISION/residual_compiler_r1.py
python3 results/LANE_SURGICAL_INCISION/centered_residual_r1.py
python3 results/LANE_SURGICAL_INCISION/hx_refinement_r1.py
python3 results/LANE_SURGICAL_INCISION/compiled_query_r1.py --length-mm 9
python3 results/LANE_SURGICAL_INCISION/compiled_query_r1.py --length-mm 8
```

Actual stdout/stderr logs are in `r1/*_run.log`. The commands created new, exclusive result directories/files and should not be run again with the same output. `incision_r1.py --tag physical_patch_replay` can be run directly without replacing earlier raw data. For other new experiments: choose a new explicit directory suffix in a lane-local copy before the run. The warm NumPy consumer only writes stdout and can be run directly against the existing package.

The first incision version, with a five-row strip and constant illustrative HX seed, is preserved as `r1/incision_r1_first_snapshot.py`. The first residual compiler before affine centering is preserved as `r1/residual_compiler_uncentered_snapshot.py`. The snapshot files are preserved source text; their relative LANE/import paths refer to the original location and they are not standalone entrypoints from r1. The latest incision version has fixed observation depth8mm and actual cohesive mean damage as the HX port. The residual compiler's new centering is optional; default remains uncentered, and `centered_residual_r1.py` selects the new operation explicitly.

Parents: Q033/Q036 were reused; Q049 and the private SURG scripts were read. The graph's `working packet`/`working rank` were run read-only; the packet was saved locally as GRAPH_CONTEXT_R1.json. No `dispatch`, `feedback`, queue write or graph mutation was run because the lane write boundary forbids their external output.

A primary article was archived with `curl -L --fail --max-time 45 https://arxiv.org/pdf/2412.02021` to `literature/goda_2024.pdf`, followed by pdftotext. Network errors at some PMC/DOI/PDF views are source acquisition states, not numerical or scientific negative results. WORLD_REFERENCES_R1.json distinguishes frozen comparisons from later context sources.

`finalize_r1.py` compiles preserved artifacts and their costs and performs integrity/numerical verifications. Direct comparison against full FE, residual identity and blind certificate rejection are the core verifications; no implementation-mirroring unit test suite has been added.
