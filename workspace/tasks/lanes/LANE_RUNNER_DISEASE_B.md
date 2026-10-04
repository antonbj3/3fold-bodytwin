# Mission: diseases as high-dimensional objects (Sol B)

You're a lane runner. `lane-model`Agent with full mandate from the coordinator. (Anton)You own all decisions. Perform without asking; report the default selection. You are working **oberoende** av andra agenter — write your own interpretation, even if it clashes with others'.

## Goal
Describe **individual diseases in high dimensions** and lay the groundwork to actually **solve them** — not as labels, but as structures with interactions and edge profile. Produce a **scaled disease library** as high-level goals to input and densify the boundary work.

## Representation (start from, improve freely)
A disease \(d\) is a structural element in a product space:

  D_d = ( g, e, h, x, s, theta, tau ; C, B )

- g — genotype/molecular configuration (perturbed nodes/pathways)
- e — environment/exposure (microbiome, toxins, nutrition, mechanical load, lifestyle; TIDSBEROENDE e(tau))
- h — heredity (germline, epigenetic inheritance, family history) — separate from somatic g
- x — molecular expression (transcriptome/proteome/metabolome)
- s — physiological cell/tissue/organ state (kinetics, transport, mechanics)
- theta — boundary conditions/parameters that the physics-based body model consumes
- tau — time/scale structure (time constants, exposure window)
- C — interaction structure (graph g->x->s->theta->y, incl. g×e gene-environment interaction and g×h)
- B — edge profile (what leaks/blocks/accumulates at the edge of the system)

Key aspects: \(g\times e\) (same genotype → different phenotype depending on environment), \(h\) as prior/covariate rather than cause (P(y | g,e,theta; h)), and \(e(τ)\) time-dependent → the disease becomes a BANA/trajectory, not a point.

## "Solve" a disease = find a structure preservation map
  M: D_d -> y (predictive on hold-out, discriminating against healthy/baseline)
plus identifiability (which components of D_d are determinable from y, which are latent) and intervention u: y -> which variable in D_d moves y ("cure/stop").

## Task
1. Select a set of diseases (scale: as many as provide value — justify number). Cover different subject areas (metabolic, neuro, immune, cardiovasc, oncology, degenerative, infection...).
2. For each: describe D_d high-dimensional (all axes), interaction graph C, edge profile B, the structure that makes the disease "severe", and a discriminant test.
3. **Candidate maps:** for each disease, point out which maps (omics→theta, g→x, e(tau)→s, etc.) are required to solve it.
4. **Densification:** identify which maps are REUSED across multiple diseases — that's the highest value signal.
5. Scale: produce a machine-readable library (JSON) + an insightful markdown report.

## Inputs
- The graph: `./graph working rank`, `./graph working packet --id <ID>`, `notes/GRAPH_WORKFLOW.md`.
- Molecular/signaling layer: MOL-* and signaling MODEL-* (105 OPEN, 10 ASSUMED, 1 DEFERRED, 1 REFUTED).
- Existing probes run: pathway-co-essentiality (mt_pathway_cert.json), L1000-MoA (mt_pharmacology_overlap.json), DNAm-epigen (mt_dnaepigen.json), same-sample-multiomics (TEA-seq, patch-seq), tissue classifier (mt_liver/mt_somatosensory); signaling ODEs in scripts/msk/.
- Ingested: data/hpa_gtex_ingest/rna_tissue_consensus_subset.tsv (20 162 genes × 51 tissues).
- First principles: external_research_path

## Output
1. `external_research_path` — insightful overview.
2. Machine-readable disease library (JSON).
3. `tasks/lanes/LANE_RUNNER_DISEASE_RESULT.md` — decisions, numbers, maps, density, what was NOT done.

## Constraints
- No fabricated data/citations; negative results are preserved; no automatic admission of evidence. Follow AGENTS.md and START.md.
- No emails/pushers/credentials; no SEED results are deleted.
