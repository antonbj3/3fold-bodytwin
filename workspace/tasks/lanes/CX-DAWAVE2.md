# CX-DAWAVE2 — 30 more ANALYSIS packets (cross-trial, held-out) over DATAMATRIX_TABLES + the datasets not yet included

Read first: `results/CX-DMCOMPUTE/RESULTS.md` + `QUEUE_PROPOSAL.txt` (BT-DA-001…008 already running ×3). Do NOT duplicate them.
Extend with:
- (1) OrthoLoad (`~/projects/bodytwin/scripts/msk/index_orthoload_forces.py`, grep the data path; hip/knee/spine per activity and patient): compute the tables deterministically first, then analyses;
- (2) OpenCap LabValidation (check the mount point; the previous attempt found it disconnected — look for an alternative path in tasks/index/DATASETS.json);
- (3) the geometries Imperial/VSD/Keast (shape parameter tables already exist in results/ → cross-dataset analyses).
The same filter as CX-DATAMATRIX (named consumer, raw data/deterministic tables, ≥ 2 consumers, a non-circular held-out check, not already done).
Numbering BT-DA-009…038. Do NOT queue: `results/CX-DAWAVE2/QUEUE_PROPOSAL.txt`, which the coordinator reviews.

