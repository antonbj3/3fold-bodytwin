### Q127 — Formation and retention of myonuclei

**Question:** When do the formation and loss of fiber-bound myonuclei add information about retraining beyond the muscle's previous size?

**Resolution target:** Nuclei's placement within a fiber, fusion, cytoplasmic domains and nuclear/fiber identity over time.

**Physical basis:** Reaction kinetics, matter/energy balance and spatial cell mechanics; stochastic transitions when small counts require them.

**Mechanism:** Key nuclei to fiber and training history; retention is a nuclear counting phenomenon and should not be reported as functioning motor memory.

**Mathematical start:** `dM/dt=J_fuse−k_loss M. M [myonuclei/fiber], J_fuse [myonuclei/(fiber·day)], k_loss [day⁻¹]. If F is the number of fusion events per day, a dimensionless nuclear yield per event is used, not an extra day⁻¹ factor.`

**Detailed reference to develop:** Spatial cell and niche models with individual states and lineage tracing, interconnected with intracellular/extracellular transport.

**Laws/parameters to determine:** Regulatory network, division/death/fusion, lineage history, epigenetics and niche-dependent rates.

**Outputs:** Retention fraction, fraction 0–1; Nuclear density, count per mm²; Fusion events, count per fiber and day

**Reuse:** `src/bodytwin/cells/musculoskeletal/muscle_memory_substrate_function_dissociation.py`

**Next work/data:** Operation: longitudinal identity, section bias control and spatial mixture model for fusion.

**Comparison:** Calibrated model of nuclear formation and loss without persistent training effect, tested against the same longitudinal fiber-identified data.

**Test that can defeat the proposal:** Dismiss retention if fiber-bound nuclear amount does not differ on held individuals or sections.

**Cross-connections to test:** Q032, Q067, Q071, Q076. Mechanistic development proposal without a specific external source in this review.

