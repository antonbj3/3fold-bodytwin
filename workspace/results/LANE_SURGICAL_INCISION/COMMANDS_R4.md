# Executed R4 commands

All scripts and results are in the lane's own folder. Run environment: OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2 NUMEXPR_NUM_THREADS=2 PYTHONDONTWRITEBYTECODE=1. No large-RAM/GPU/heavy_run, no FE/LP.

1. Read all of COMMON, steerR4, the lane brief, latest checkpoints/results/next and peer artifacts. No new graphdispatch because writing only in the lane's own folder; local GRAPH_FEEDBACK_R4 is the handoff.
2. Freeze PREREG_R4.json before the new research calculation; acquire_r4.py acquires originals. Run through the Python2thread environment. Manual legitimate second acquisition uses the same `fetch` with new suffixes: owen2022_html, owen2022_bioc_retry and owen2022_mdpi_pdf. All manifests including misses are under literature/r4.
3. Read the Barnett original with pdftotext; extract original raster with PyMuPDF. Own png-format original images and pages in literature/r4. Freeze PREREG_R4_SERIES.json before extraction/test; freeze_series_r4.py writes SOURCE_TARGETS_R4.json. Correct Owen range to 0,40–0,58 before prediction in new SOURCE_TARGETS_R4_v2.json and FREEZE_RECEIPT_R4_v2.json.
4. `python results/LANE_SURGICAL_INCISION/needle_series_r4.py` writes exclusively to r4/. For reproduction: `python results/LANE_SURGICAL_INCISION/needle_series_r4.py --out results/LANE_SURGICAL_INCISION/repro_r4_new` (choose a directory that does not yet exist).
5. Freeze PREREG_R4_WORK_AREA.json before budget/nuisance calculation. `python results/LANE_SURGICAL_INCISION/work_area_r4.py` writes exclusively to r4/work_area_requirements.json.
6. `python results/LANE_SURGICAL_INCISION/handoff_r4.py` writes final PORTS, source/peer/decomposition/attempts/cost/feedback. `python results/LANE_SURGICAL_INCISION/plot_r4.py` gives standalone PDF/PNG; MPLCONFIGDIR within the lane.
7. `python results/LANE_SURGICAL_INCISION/verify_r4.py`:42/42 PASS. It validates integrity/scope and the exact bounded-ratio error budget, not biological admission or independent scientific review.
8. finalize_r4.py runs after verification. Snapshot previous WORK_STATUS/NEXT before updating, append-only new RESULTS section, final CHECKPOINT/manifest/nightr4.

Fixed-destination scripts refuse existing files. Reproducing them requires scripts/frozen inputs to be copied to a new directory within the lane with their relative dependencies; original r4 must not be overwritten. Reproduction in a separate directory may change timing and acquisition metadata. R1–R3 have been preserved, including previous misses.
