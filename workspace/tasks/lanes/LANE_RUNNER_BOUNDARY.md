# Task: expand BodyTwin at the BOUNDARY — a maximally insightful whole that can be built on

You're a lane runner. `lane-model`Agent with full mandate from the coordinator. (Anton). You own all decisions. The coordinator only serves context. Perform without asking; report default choices.

## Core insight (the coordinator's direction)
What is missing is not more inner molecular nodes. What is missing is **the connection at the boundary itself**: an end-to-end path that takes external/own omics → DE → pathway/gene-set → tissue/disease-specific signalling map → **boundary conditions/parameter set that the physics-based body model consumes** — and the reverse port. All expansion must happen **at the boundary**, not in the interior.

## Task
Read the original message below (verbatim). Make **your own, freely chosen expansion** whose combined result becomes **maximally insightful for further building**. You may choose between:
- (a) creating **~300 new boundary mutations**, and/or
- (b) letting new mutations **build on the existing 500** (v6) and each other.
Justify your choice. Prioritise **depth and reusability** over count — a smaller, sharp set that connects the boundary is better than 300 weak ones.

## Konkreta krav
1. **Boundary focus:** every new mutation must hit an explicit **port** between the omics world and the physics-based body model (e.g. omics → BC/parameter set; simulation → testable omics prediction). Name port variables, units, time windows and consumer.
2. **Resolution:** make the difference from existing attempts **measurable**: exact observable quantity, exact baseline, exact discriminating test at matched cost.
3. **Connectivity:** every new mutation must bind to **at least two existing** nodes/mutations AND to the graph (`MOL-*`, signalling `MODEL-*`, BC/parameter ports). Use the project's `./graph working rank` and `./graph working packet --id <ID>` to ground the connections. No invented coverage.
4. **First principles:** follow `~/research/FREE_AUTONOMY_20260926/FIRST_PRINCIPLES.md` (elements, interactions, controlling relationships, leaf status, stopping argument, discriminatory test).
5. **Evidence integrity:** no fabricated data/citations; negative results preserved; no automatic evidence admission. Follow `AGENTS.md` and `START.md`.

## Indata
- **Originalmeddelandet** (koordinatorns, ordagrant) — se nedan.
- Befintliga 500 (v6): `~/research/BODYTWIN_MUTATIONS_IMPROVED_v6.md`, `_TOP40.md`; and `./tasks/free48/SEED_PROGRAM_500_IMPROVED_20260929.json`.
- Originalserien: `~/research/BODYTWIN_100_MUTATIONS_20260929*.md`.
- Grafen: `./graph working rank`, `./graph working packet --id <ID>`, `notes/GRAPH_WORKFLOW.md`.
- Run context: `tasks/free48/{STATE.json,CATALOG.json,INITIAL_JOBS.json}`.

## Utdata
1. `~/research/BODYTWIN_BOUNDARY_REPORT.md` — en **insightful** overall picture: the boundary gap, the gates, what is done now and why, as well as an explicit map of connectivity to v6/grafen.
2. The boundary mutations (count according to your choice), with numbering that does not collide with M1–M500 (e.g. `BM1…`), machine-readable file + markdown.
3. A short "buildability" section: which next agent/loop can take this further without you, and what it then needs.
4. `tasks/lanes/LANE_RUNNER_BOUNDARY_RESULT.md` — decisions, justification, exactly what was created, what was NOT done.

## Constraints
- You own the decisions; the coordinator touches nothing.
- No emails/pushes/credentials. No SEED results deleted. Do not touch running jobs (only new files + possibly BRIEF for unstarted jobs if you choose).
- Write `/opt/agents` files as ubuntu if you touch the cloud; otherwise work locally in the workspace.

---

## ORIGINALMEDDELANDET (ordagrant)

"Sorry — here is the deep dive into the subject itself (transcriptomics/molecular signalling), not repo archaeology.
What actually exists at that level
117 molecular/signalling nodes in the graph (MOL-* + signalling-MODEL-*): 105 OPEN, 10 ASSUMED, 1 DEFERRED, 1 REFUTED. The vast majority are literature/mechanistic anchor cells (monogenic diseases, one signalling pathway each), not RNA-seq analyses.
Executed molecular probes (all 2026-07-18, under MOL-* cells — reports/probes/):
Cell / file | What it did | Result
MOL-PATHWAY-NETWORK-FILLED-DECODER (mt_pathway_cert.json) | Pathway/network analysis: Reactome pairs vs DepMap CRISPR co-essentiality | 24 574 pairs, 0,04 vs essentiality-matched null 0,0057 = 7×; BioGRID 9,02×
MOL-PHARMACOLOGY-DRUGRESPONSE-TOX (mt_pharmacology_overlap.json) | L1000 transcriptome signatures + PRISM viability + ToxCast, MoA clustering | 629 comp., 138 MoA clusters, κ 0,115 vs null 0,093 (z 4,86 → 18,85 after mega-cluster exclusion)
MOL-DNA-EPIGEN-PROGRAM (mt_dnaepigen.json + _leg2) | DNAm clocks + GSE113957 fibroblast transcriptome, age prediction | held-out R² 0,712 (n=141, 1–96 years), z 7,5 above null (2 decorrelated legs)
MOL-SAMESAMPLE-MULTIOMIC-CALIBRATION | TEA-seq RNA+ATAC+ADT, same cell | ARI 0,376, z 286 (recalculated by the coordinator)
MOL-SAMESAMPLE-PATCHSEQ-CALIBRATION | transcriptome+ephys+morphology | ARI 0,177, z 16,6
mt_liver.json (GSE124395), mt_somatosensory.json (GSE249746, 1136 DRG neurons, acc 0,938) | tissue/cell-specific transcriptome classification | executed
Signalling models (scripts/msk/): insulin_pi3k_akt_signaling.py, nfkb_signaling_dynamics.py (+Hoffmann2002 BIOMD140-rebuild), mapk_erk_cascade.py, mtorc1_feedback_rebuild.py, hedgehog_gli_8state_ode.py, wnt_betacatenin_rebuild.py, angiogenesis_vegf.py, coagulation_hemostasis.py, thrombin_full_hockin_mann_falsifier.py, nephron_tgf_oscillation_rebuild.py, insulin_ultradian_sturis_hopf.py.
Data actually ingested: data/hpa_gtex_ingest/rna_tissue_consensus_subset.tsv (HPA/GTEx tissue RNA, 20 162 genes × 51 tissues) for MOL-TISSUE-TRANSPORTER-LOCALIZATION; ~30+00 GEO accessions referenced (GSE113957, GSE122541, GSE84133, GSE124395, GSE116222, GSE241132, GSE158013, GSE179848, GSE136831, GSE92742/LINCS…), plus PXD proteomics and phs/dbGaP.
Against your collaborator’s taxonomy
- RNA-seq/GEO ingestion: yes. Multi-omics integration: yes, run E (RNA+ATAC+ADT, patch-seq, DNAm+RNA). Pathway/gene-set: yes, but as individual cert probes (Reactome-co-essentiality, L1000-MoA), not as a product. Differential expression: weak/tangential — the analyses are classification/age prediction/clustering, not DE tables. Signalling maps: yes, as ODE models. Disease/tissue-specific: yes — MOL-* are almost exclusively disease-gene-specific.
The gap (what you asked about)
There is no end-to-end pipeline that takes external/own RNA-seq → DE → pathway/gene-set → tissue/disease-specific signalling map → boundary conditions/parameter set consumed by the physics-based body model. The molecular layer is used as anchors/legs to certify hidden states (one mechanism per cell, external anchor) — pathway DBs are explicitly "FILLED-decoders" (priors, not measured maps), and the connection to physics goes via geometry→physics→proof + cert/fusion, not omics→BC/parameters. In addition, 105/117 MOL-* are still OPEN.
So: overlapping components yes, and some genuinely executed analyses (pathway-co-essentiality, L1000-MoA+transcriptome, same-sample-multiomics, transcriptome-based age/tissue classifiers, signalling ODEs) — but specifically "signalling map as the model’s boundary conditions/parameter set" does not exist.
Do you want me to write this up as a short gap analysis (what we have vs. what they offer, with file references) in research/?"
