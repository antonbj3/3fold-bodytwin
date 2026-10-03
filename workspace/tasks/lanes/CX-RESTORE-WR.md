# CX-RESTORE-WR — restore the collaborator's full-body equilibrium system BEFORE Q/SVD (W and R, with named rows/columns) for Field's U401 gauge analysis

Field's U401 (`the public staging tree/3fold-motion-engine/_private/romi_collab/build/U401/RESULTS.md`, A1385) cannot analyse the graph: the stored A is dense after SVD rotation, and the links to frames/act_raw/refs/reac2 in BT-XF5/N40b are broken.

## Task
- Find the build path for the collaborator's system: `results/CX-SOLVER2/n40b_build.py` / `cx_solver.py`, `results/N40b`, `results/N14b` (n14_solver.py, n14_solve.py), `results/BT-XF5`, and `bodytwin_core/vendor/solver`.
- Rebuild, without the Q/SVD step, per time step (at least steps 1, 75 and 141, preferably all):
  - W: the physical equilibrium matrix, i.e. force columns → equations;
  - R: the external load/right-hand side;
  - named ROWS (segment × DOF) and COLUMNS: muscle element name, contact coordinate (medial/lateral TF, PF, …), joint reaction, and the type of each (muscle/contact/reaction), with the one-sidedness (≥ 0) marked.
- Check: after your own QR/SVD, the solution must equal CX-SOLVER2's forces to ≤ 1e-8.
- Repair the broken links: write a MANIFEST with the correct paths + sha256.
- Output: `/mnt/shared_data/bodytwin_work/share_field/john_WR/` (npz + json names + MANIFEST + README).
  - When done, write `results/CX-RESTORE-WR/READY_FOR_FIELD.md` with the paths.
- `results/CX-RESTORE-WR/RESULTS.md` starting with `# CX-RESTORE-WR`.
