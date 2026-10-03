# CX-FIELDSHARE — femur test cases + standalone edit operator for the field lane (shape derivatives)

Consumer: the field lane (anton-12/anton-51a9) computes d(V, centre of mass, I)/dθ with error bounds and wants real femurs with anatomical parameter editing as test cases.

## Deliverables: `/mnt/shared_data/bodytwin_work/share_field/` (+ copy of the report in `results/CX-FIELDSHARE/`)
1. 4 femurs: 3 Imperial (public zenodo 167808) + 1 VSD. For each, `femur_<id>_base.npz` (V mm, F, named landmarks, measures CCD/AV/length per geomgr measures.py MEASURE_VERSION) and edited versions `_AV+10`, `_CCD-10`, `_LEN+10` via `results/N7c/geomgr/ops.py:op_edit_osteotomy` (with y0/y1, u = twist/varus/lengthening, held_maxabs, distal_bit_identical).
2. `femur_edit.py`: standalone (numpy/scipy only) copy of `edit_frame`, `osteotomy`, `_measure` and what they need from geomgr (core/measures). Include a source reference per function (file:line). API: `edit(V, F, reg, target, delta) -> V1, report` plus `osteotomy(V, frame, u)` for arbitrary u. The latter is what the field lane differentiates: θ = u (twist, varus, lengthening), with smooth transition (smoothstep weight w).
3. `test_femur_edit.py` (pytest): the standalone version reproduces geomgr's V1 at ≤ 1e-9 mm on all 4 × 3; u = 0 gives identical V; central FD of V with respect to u is finite and smooth.
4. `README.md` (≤ 30 lines): units, frame, what θ means, how w works, the known effect (length "held" 92.9 % in A183), and the data licence (Imperial CC-BY per zenodo; VSD according to its licence; say UNKNOWN if you cannot verify).

## Rules
- Read-only in `~/projects/bodytwin` and other lanes' directories. No computation over 60 s per step. No emails or pushes.
- `results/CX-FIELDSHARE/RESULTS.md` starting with `# CX-FIELDSHARE`: list of files, sha256, test results.
