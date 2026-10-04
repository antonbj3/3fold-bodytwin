# LANE_PRC2_DRIFT

Anton's seed (X bookmarks, batch1 TOP5 #2; external_research_path and SEEDS_BODYTWIN.md B1). Results folder `results/LANE_PRC2_DRIFT/`.

## Missing capability

`bodytwin:AGING-ETA-DRIFT-CELL` is REFUTED on a null result against a composite clock (Horvath2, partial Spearman −0,208, p = 0,31, n = 4; evidence file `~/projects/bodytwin/reports/probes/mt_aging_null.json`, read-only). We lack an outcome measure that separates slow chromatin drift from fast layers. Hypothesis (Yücel, Molière & Gladyshev, Nature 658:45–54, 2026): drift in PRC2-bound low-methylated regions follows **cell divisions**, not calendar time.

## Operation

Load GSE179847 (fibroblasts in culture with population doublings and calendar days, including contact-inhibited samples). Annotate CpG with the PRC2 set from public H3K27me3 peaks for fibroblasts. Drift measure = mean deviation from the donor's own baseline. Per donor: drift against doublings and drift against days, separately in the PRC2 set and in the control set.

## Strongest control

Size-matched random set: same number of CpGs, same distribution of baseline methylation and coverage as the PRC2 set (otherwise regression to the mean is measured at low-methylated positions). Second control: contact-inhibited samples (days pass, no divisions) must follow doublings rather than days.

## Frozen target

PREREG with hash before the run: the difference in slope against doublings between the PRC2 set and the matched control (with intervals), plus the same difference against calendar time. The result is reported regardless of sign. If drift follows divisions: write a proposal for a new outcome variable for AGING-ETA-DRIFT-CELL in the results folder (not in the graph).

## Data and resources

Data on `external_mount` (never `/` — root disk 98 % full). The methylation matrix is a few GB: read in chunks/by column, keep only the CpGs needed, save summaries. Over 8 GB RAM → `tasks/heavy_run.sh`. Public data; no internal data.
