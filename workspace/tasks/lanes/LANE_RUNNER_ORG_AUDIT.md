# Task: verify that ALL old mechanism/3FOLD material is available to the graph and workers

You're a lane runner. `lane-model`Agent with full mandate from the coordinator. (Anton)You own the decisions. **not** att bygga nytt — it is to verify that the existing is **Available, readable and usable** of graphene and of running workers.

## Background (the coordinator’s context, verbatim claims to verify)
Anton previously worked with a program described as follows (mechanism / 3FOLD BodyTwin, `~/projects/bodytwin`, a fork of cad-to-simulation):
- **RNA-seq/GEO ingestion at scale:** 26 source-catalog waves in `projects/bodytwin/data/mechanism_catalog/sources_wave*.json` (~310+ sources, ~all BT-HOLD). Molecular legs densest: wave 2 (signal time series: NF-κB, ERK/Akt, Ca²⁺, insulin+EGF phosphoproteomics, GH pulse), wave 7 (patch-seq, CITE-seq, 10x-Multiome, SHARE-seq), wave 21 (same-sample multimodal), wave 25 (Reactome/WikiPathways/GO/STRING/BioGRID; KEGG/MSigDB flagged as license traps).
- `data/MECHANISM_ANCHOR_GRAPH.json` has 643 strings referencing transcriptome/RNA-seq/scRNA/GEO/proteome, with named accessions (GSE130977, GSE143704, GSE150482, GSE249746, GSE124395, GSE122541, PXD006182…).
- **DE/gene-set/pathway executed:** rank-based single-sample gene-set enrichment AUC on 7 428 common-coverage genes, with cross-organ held-out transfer (heart/intestine) + adversary stripping literature overlap; transcriptome→phenotype R² work including self-caught collapse of "same-cell coupling" R²=0,72 to cell-type identity; pathway/network DBs treated as "FILLED-decoders" (priors) — WAVE_PLAN.md arch-note 6.
- **Multi-omics executed:** `reports/probes/mt_multiomic_calibration.json` (TEA-seq RNA+ATAC+ADT, ARI 0,376, z≈286, h5ad intermediates 108 MB RNA + 326 MB ATAC), `reports/probes/mt_patchseq_calibration.json` (z 16,6).
- **Signaling ODEs:** `scratch_glut4sig_model/` (insulin→PI3K/Akt→GLUT4, 21-state), cardiac CICR, mitochondrial OXPHOS, Na/K-ATPase in `~/outputs/`. Program = "geometry → physical model → evidence", first closed cell = knee connected to physics-biomechanics anchors (COORDINATOR.md §5).
- **The genuine gap:** no end-to-end "disease/tissue-specific signaling map → physical model’s BC/parameter set". Omics is used as anchors/legs in a certification pipeline, not as a standalone deliverable.
- **Two warnings:** much of the catalog is BT-HOLD hypothesis, not executed (WAVE_PLAN.md); some probes self-report WEAKENED/IN-FLIGHT — the graph’s own notes are the honest source.

## Task (verification, not new construction)
1. **Inventory** the above against disk: do `projects/bodytwin/data/mechanism_catalog/sources_wave*.json`, `data/MECHANISM_ANCHOR_GRAPH.json`, `reports/probes/mt_*`, `scratch_glut4sig_model/`, `~/outputs/`, `WAVE_PLAN.md`, `COORDINATOR.md`, the `bt_memory/` probe exist? Count lines/strings/records and note exact paths + sizes.
2. **Availability to the graph:** can the graph reach and consume this? Run the `./graph` tools (working rank/packet/dispatch/feedback) and verify that relevant nodes (MOL-*, ORG-*, signaling-MODEL-*) exist and that the sources are linked. Report broken or missing links.
3. **Availability for workers:** the material is located on a path which: SEED-/field-workers actually read? (T.ex. `inputs/`, kataloger under `/mnt/games-240/research/...`, `results/<job>/inputs/`.) Verify that a worker can find and read a representative file without a special right.
4. **Executed vs planned:** produce the exact split for the transcriptomic cells (executed analysis with result file vs BT-HOLD/WAVE_PLAN plan). The graph’s own notes are ground truth.
5. **Gaps:** list exactly what is NOT available (missing, wrong path, broken link, permission problem, outdated snapshot) and what is required to make it available.
6. If something is easy to fix (symlink, manifest, copy to the right inputs directory): **do it** — but do not touch SEED results or running jobs, and document every intervention.

## Utdata
1. `~/research/ORGANIZATION_AUDIT.md` — overall picture: what exists, where, available for graph+workers or not, run-vs-planned-split, glitch.
2. Machine-readable inventory (JSON) with path, size, status (available/missing/broken), graph link, worker-readable.
3. `tasks/lanes/LANE_RUNNER_ORG_AUDIT_RESULT.md` — decisions, exact findings, what was fixed, what was NOT done.

## Scope
- This is **verification**, not new research. No invented data. The graph’s own notes are the honest source.
- No mail/pushes/credentials; no SEED results are deleted or changed.
- Write/change in the cloud as ubuntu if needed; otherwise locally in the workspace.

---

## THE ORIGINAL QUOTE (verbatim, for context)
"Yes — substantially. 'Our research' here is the mechanism / 3FOLD BodyTwin program (~/projects/bodytwin, a fork of cad-to-simulation) … [see the brief: sources, accessions, executed probes, ODEs, the gap, the two warnings] … Want me to pull the exact executed-vs-planned split on the transcriptomic cells, or draft a short 'here's what we already have / here's the genuine gap' note for that conversation?"

The coordinator’s addition: "Old stuff I have worked with is described here, I just want to ensure everything is available to the graph and our workers."
