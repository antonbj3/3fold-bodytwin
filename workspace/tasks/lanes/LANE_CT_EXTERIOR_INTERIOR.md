# LANE_CT_EXTERIOR_INTERIOR

Antons seed (X Bookmarks 28–30/9, batch2 TOP5 #2 and batch3 TOP5 #3; sources in ~/research/X_BOOKMARKS_20260930/). Resultatmapp `results/LANE_CT_EXTERIOR_INTERIOR/`.

## Missing capability

- `MSK-VISCERAL-FAT-IDENTIFIABILITY` has R² ≈ 0,65–0,68 for visceral fat (VAT) from the body's exterior, but the label (DXA-VAT) is itself derived from the abdomen's external geometry — the result is circular and the ceiling unknown.
- `ORGAN-DEPTH-VISCERAL` lacks an organ-imaging cohort; spleen volume varies 112–215 cm³ between cohorts.
CT has the outside and inside in the same scan: a non-circular label.

## Decomposition and operation

Per CT from TotalSegmentator v2 (1 228 CT, 117 structures, CC BY 4.0, Zenodo 10047292, zip 23,6 GB): body mask (HU > −500), L3 from vertebral masks, VAT area at L3 (HU −190…−30 inside the abdominal wall), silhouette measures from orthographic AP and lateral projection (width, depth, circumference per level), body-size mode s = PC1 of the silhouette, organ volumes from the dataset's masks. Relations: VAT_L3 = f(silhuett); V_organ = g(s) + residual. 10×5-CV with ridge and gradient boosting.

## Strongest control and frozen target

Baseline sex + height + body volume; subset without contrast and with L3 fully in view (report exclusions); residual within waist decile; shuffled labels between patients should give R² ≈ 0; report the cohort's pathology mix. Frozen target: R²(silhouette → CT-VAT_L3) with 95 % CI. |R² − 0,65| ≤ 0,05 → the identifiability ceiling is no label artifact; R² < 0,60 → the DXA circle carried part of it. Plus fraction of organ-volume variance explained by s for liver, spleen, kidneys and muscle area at L3. Write proposals about node numbers in the result directory, not into the graph.

## Data, disk and resources — IMPORTANT

- **All data on `/mnt/games-240/research/bodytwin_ct/`** (94 GB ledigt). **Never on `/` eller `/mnt/shared_data`** — The root counter is: 98 % full.
- Download the zip file there; do **not** unpack everything. Read members streaming from zip (Python zipfile) or extract only a subset at a time and delete extracted volumes after features are calculated. Save only features (CSV/JSON) and small summaries.
- Start with a pilot of ~50 CT, measure time/memory per CT, then scale.
- RAM per CT may be a few GB: run heavy steps via `tasks/heavy_run.sh` if > 8 GB, otherwise directly with 2 threads. No GPU is needed.
- NHANES- the comparison (`/media/anton/8838D60F38D5FBDE/mechanism_data/vat_test`) may only be read locally, and only summarized numbers may be quoted; it contains no personal data that may leave the machine.
- Inga interna BodyTwin-personer, the collaborator, Grand Challenge eller restricted model data.
