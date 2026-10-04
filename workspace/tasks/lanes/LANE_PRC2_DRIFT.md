# LANE_PRC2_DRIFT

Anton's seed (X bookmarks, batch1 TOP5 #2; external_research_path and SEEDS_BODYTWIN.md B1). Results folder `results/LANE_PRC2_DRIFT/`.

## Missing ability

`bodytwin:AGING-ETA-DRIFT-CELL` stands as REFUTED on a null result against a composite clock (Horvath2, partial Spearman −0,208, p = 0,31, n = 4; proof file `~/projects/bodytwin/reports/probes/mt_aging_null.json`, read-only). We lack an outcome measure that separates slow chromatin drift from fast layers. Hypothesis (Yücel, Molière & Gladyshev, Nature 658:45–54, 2026): operation in PRC2 bound low-methylated regions follows **cell divisions**, not calendar time.

## Operation

Load GSE179847 (fibroblasts in culture with population doublings and calendar days, including contact-inhibited samples). Annotate CpG with PRC2 amount from public H3K27me3 peaks for fibroblasts. Operating measure = mean deviation from the donor's own baseline. Per donor: drift against duplicates and drift against days, separately in the PRC2 set and in the control set.

## Strongest control

Quantity-matched random set: same number of CpGs, same distribution of baseline methylation and coverage as the PRC2 set (otherwise regression to the mean is measured at low-methylated positions). Second check: contact-inhibited samples (days pass, no splits) should follow duplications and not days.

## Frozen target

PREREG with pre-run hash: the difference in slope against duplicates between the PRC2 set and the matched control (with intervals), plus the same difference against calendar time. The result is reported regardless of sign. If operation follows divisions: write a proposal for new outcome variable for AGING-ETA-DRIFT-CELL in the results folder (not in the graph).

## Data and resources

Data on `external_mount` (never `/` — root disk 98 % full). The methylation array is some GB: read chunkwise/column wise, keep only the CpGs needed, save summaries. Over 8 GB RAM → `tasks/heavy_run.sh`. Public data; no internal data.
