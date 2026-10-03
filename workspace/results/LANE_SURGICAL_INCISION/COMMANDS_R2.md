# Executed R2 commands

All Python calculations use `OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2 NUMEXPR_NUM_THREADS=2`; local CPU only. Writes stayed in `results/LANE_SURGICAL_INCISION/`. New data paths are exclusive. These commands reproduce their original runs only in a fresh destination; do not rerun against existing numerical files.

```
pdftotext -layout literature/r2/pissarenko_2020.pdf literature/r2/pissarenko_2020.txt
python layered_skin_r2.py --nx 32 --out r2/layers_n32
python layered_skin_r2.py --nx 64 --out r2/layers_n64
/usr/bin/time -v python layered_skin_r2.py --nx 96 --out r2/layers_n96
python acquisition_r2.py
python layered_skin_r2.py --nx 32 --out r2/layers_v2_n32
/usr/bin/time -v python layered_skin_r2.py --nx 64 --out r2/layers_v2_n64
/usr/bin/time -v python layered_skin_r2.py --nx 96 --out r2/layers_v2_n96
/usr/bin/time -v python layer_extension_r2.py
python finalize_r2.py
```

Commands above are relative to this lane for readability; actual invocations were from the workspace root with the full lane path. Raw logs are in r2/*.log. The `layers_n*` version is preserved as `r2/layered_skin_r2_initial_snapshot.py`; final `layers_v2_n*` uses the current module. The first layer-extension energy report came from `r2/layer_extension_r2_initial_snapshot.py`. Current extension corrects the reporting sign; original raw results remain untouched.

Additional inline Python commands, saved outputs:

- Primary PDF download using urllib from `https://meyersgroup.ucsd.edu/papers/journals/Meyers%20477.pdf`; Chu2019 primary HTML from Nature. Manifest records hashes and failed acquisitions. Four prereg/targets files hashed before their calculations, then a final targeted-refinement supplement hashed before its two cases.
- Corrected energy identity recomputed from existing rows → `r2/energy_balance_check_corrected.json`, without rerunning the mechanical solver.
- `Layers(128,64,0).solve([.020,.008,0])` → `r2/targeted_refinement_n128.json` and raw NPZ.
- `Layers(128,64,angle).solve([.008]*3)` for angle0/45 only → `r2/full8_targeted_refinement_n128.json` and two raw NPZs.
- `verify(Layers(32,16,45))` → VERIFICATION_R2.json via finalize.
- Every new JSON parsed, prereg hashes compared, required handoffs/files checked, raw initial/final versions preserved. This validates reproducibility/integrity, not scientific admission.

Final report supplements the two full8mm refinement cases after preserving the failed64→96 criterion. No new material parameters were fitted. The ordinary matched FE/LP control receives identical data and equations.
