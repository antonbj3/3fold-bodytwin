# Perturbations, substances and modeling methods

These are model inputs and research mechanisms. Each needs its own target-binding, exposure, tissue-access and downstream-function model; no treatment regimen is proposed.

| Class / example | Mechanism to resolve | Required comparison or unresolved input |
|---|---|---|
| Checkpoint antibodies, e.g. anti-PD-1 / PD-L1 / CTLA-4 | Occupancy, receptor turnover, co-stimulation, cell subset | Tumor and tissue-specific antigen context; calibrated PK axis (I09–I11) |
| Tocilizumab / anakinra | IL-6 receptor versus IL-1 receptor intervention | CRS and neural/meningeal injury are separate endpoints; animal evidence I01 |
| Infliximab / vedolizumab | TNF signaling versus alpha4beta7-associated gut trafficking | Local versus systemic effect; observational confounding and exposure (I06) |
| Abatacept / ruxolitinib | CD80/86 co-stimulation versus JAK signaling | Combination, organ function and respiratory-muscle effects need separation (I07) |
| Glucocorticoids | Receptor kinetics, genomic delays, metabolic and immune effects | Link HPA and cell-specific receptor abundance; temporal/tissue heterogeneity (M01, M03) |
| Complement intervention | C3/C5 or amplification regulator perturbation | Injury suppression versus antimicrobial function; existing local complement models |
| Oxygen, glucose, lactate, arginine, itaconate, redox pools | Material and energy fluxes in specific cell states | Derive stoichiometry and distinguish substrate amount from signaling effect |
| Microbial metabolites, antigens, PAMPs, DAMPs, extracellular DNA | Production, access to sensor, transport and degradation | Alternative causal pathways and clearance; mouse-to-human transfer remains a question (I08) |
| Growth factors, antibodies and cell replacement/depletion | Production, maturation, function and persistence | Expand as underactivity mechanisms only after primary-source and native-code checks |

Method ladder: stoichiometric reaction network and receptor kinetics; CTMC/SSA when copy numbers require it; ODE/delay or population-balance models; PBPK/target-mediated PK; reaction-diffusion and chemotaxis in tissue geometry; agent-based contact/trafficking; QSP coupling; structural/practical identifiability; uncertainty-aware parameter inference; organoids and organ-on-chip perturbations; independent aggregate study validation. Inspect PopelLab QSP before rebuilding its compartments. Translate model formats only with unit/stoichiometry/observable checks. Atomistic or electronic calculations are a deeper branch where reaction rates/binding energetics dominate uncertainty, not a substitute for measuring those rates.
