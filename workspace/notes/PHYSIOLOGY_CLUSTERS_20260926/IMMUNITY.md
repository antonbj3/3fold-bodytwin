# IMMUNITY — first-principles research cluster

Status: proposed work, PENDING_INDEPENDENT_REVIEW. Equations below are model specifications requiring closure and units.

## I01: Receptor binding and antigen discrimination

- Constituents: peptide-MHC; TCR; CD28; PD-1; PD-L1; CTLA-4.
- Equation skeleton: `dB/dt=kon*C*(Rtot-B)-koff*B; proofreading chain with reverse/reset rates`.
- Observable: occupancy; activation probability; dwell time s.
- Discriminating test: Compare kinetic discrimination with existing two-signal gate at matched receptor and antigen numbers.

## I02: Antigen processing and dendritic priming

- Constituents: tumor debris; uptake; proteolysis; MHC; cDC1; lymph-node T cells.
- Equation skeleton: `dn/dt=S*v; antigen flux -> presentation -> priming with finite delay`.
- Observable: presented complexes/cell; primed cells/s.
- Discriminating test: Separate antigen abundance from presentation failure using perturbation of uptake versus MHC loading.

## I03: Tolerance, regulatory cells and autoimmunity

- Constituents: self/foreign clones; Tregs; co-stimulation; suppressive signals.
- Equation skeleton: `birth-death-transition rates by clone/antigen; explicit suppressive interaction`.
- Observable: clone counts; healthy-cell killing cells/s.
- Discriminating test: Can anti-tumor activation rise without the same change in tissue-self clones? Reject a single global immune strength parameter.

## I04: T-cell activation, exhaustion and memory

- Constituents: naive; effector; progenitor exhausted; terminal exhausted; memory.
- Equation skeleton: `dN_i/dt=division_i-death_i+sum(transitions_ji-transitions_ij)`.
- Observable: counts; killing per viable cell; recovery time.
- Discriminating test: Reuse existing T-cell code; separate proximal activation from functional killing and preserve native claim/code numerical mismatch.

## I05: Tumor killing and escape

- Constituents: cancer clones; antigen loss; MHC loss; cytotoxic cells; checkpoints.
- Equation skeleton: `dT_a/dt=growth_a-kill_a(E,contact)+mutation/selection fluxes`.
- Observable: tumor cells; antigen fractions; viable burden.
- Discriminating test: Compare equal tumor burden with altered antigen distribution; simple bulk-volume model is baseline.

## I06: Myeloid suppression and plasticity

- Constituents: macrophage state; MDSC; arginine; itaconate; oxygen; cytokines.
- Equation skeleton: `continuous state dynamics plus shared substrate balances; avoid fixed M1/M2 labels`.
- Observable: arginine flux; secreted cytokines; killing suppression.
- Discriminating test: Test whether continuum state or resource competition adds predictive value beyond existing polarization model.

## I07: Complement, amplification and host protection

- Constituents: C3; C3b; convertases; C5; regulators; terminal complex.
- Equation skeleton: `mass-action cascade with catalytic amplification, depletion and receptor/target binding`.
- Observable: C3/C5 turnover mol/s; lysis fraction.
- Discriminating test: Reuse complement and bistability code; compare target injury reduction with lost pathogen clearance, including dose-blind gate failure.

## I08: Innate sensing and inflammasomes

- Constituents: PAMP; DAMP; mtDNA; cGAS-STING; NLRP3; IL-1 processing.
- Equation skeleton: `separate priming and activation events with compartment-specific sensor access`.
- Observable: mature IL-1 mol/s; death events/s.
- Discriminating test: Challenge redundant-pathway assumptions with cell-type-specific ablations; inherited AIM2/cGAS/NLRP3 node remains open.

## I09: Neutrophils, NETs and vascular injury

- Constituents: neutrophils; adhesion; ROS; NET DNA; platelets; endothelium.
- Equation skeleton: `reaction-transport plus adhesion/detachment and NET release/clearance`.
- Observable: NET density; permeability; thrombotic occupancy.
- Discriminating test: Distinguish antimicrobial benefit from endothelial injury at equal neutrophil count; define missing native-code coverage.

## I10: B-cell immunity and immunoglobulin deficits

- Constituents: B cells; plasma cells; antibody; antigen; memory; marrow supply.
- Equation skeleton: `dN_B/dt=production-expansion-depletion; dN_Ig/dt=secretion-clearance`.
- Observable: B-cell count; Ig concentration; neutralization.
- Discriminating test: For cell-depleting immunotherapy compare target-cell removal with infection-control loss and delayed antibody recovery.

## I11: Marrow, thymus and underactivity

- Constituents: progenitors; neutrophils; lymphocytes; trafficking; exhaustion.
- Equation skeleton: `age-structured or transit-compartment birth/maturation/death model`.
- Observable: cell counts; functional response; recovery time.
- Discriminating test: Distinguish too few cells, wrong location and exhausted cells despite similar blood count; evaluate infection and vaccine response as separate outputs.

## I12: Immune trafficking in real geometry

- Constituents: blood; endothelium; chemokines; ECM; lymph; tissue-resident cells.
- Equation skeleton: `partial_t n=div(D grad n-chi*n*grad c)-div(u*n)+sources-sinks`.
- Observable: cells/m3; influx cells/s; intercell distance m.
- Discriminating test: Compare spatial model to well-mixed compartment at identical cell totals and exposure; use synthetic geometry first.

## I13: CAR-T cytokine release

- Constituents: CAR-T; leukemia; monocytes; macrophages; IL-1; IL-6; NO.
- Equation skeleton: `tumor-contact activation -> cell-specific secretion -> receptor dynamics and clearance`.
- Observable: cytokine time course; killing; endothelial response.
- Discriminating test: Reproduce a conditional IL-1/IL-6 distinction against a single cytokine pool; Norelli/Giavridis are animal mechanisms.

## I14: Neurotoxicity and blood-brain interfaces

- Constituents: endothelium; tight junctions; pericytes; cytokines; brain/meningeal compartment.
- Equation skeleton: `dN_brain/dt=PS(C_blood-C_brain/K)-clearance; permeability has explicit injury/repair law`.
- Observable: transport mol/s; permeability m/s; neural proxy.
- Discriminating test: A lower systemic IL-6 signal must not automatically imply lower CNS injury; distinguish BBB, blood-CSF and meningeal routes.

## I15: Checkpoint colitis and gut barrier

- Constituents: resident T cells; epithelium; IFN-gamma; renewal; microbial products.
- Equation skeleton: `epithelial birth-loss with immune contact; barrier mass transport coupled to cytokine exposure`.
- Observable: barrier flux mol/s; viable epithelial fraction.
- Discriminating test: Compare systemic cytokine suppression with gut trafficking blockade; preserve regeneration, spatial and microbiome confounding.

## I16: Myocarditis, myositis and electrical function

- Constituents: heart clones; muscle clones; self-antigen; injury; conduction; respiratory muscle.
- Equation skeleton: `antigen-specific tissue kill -> damage/repair -> stated organ-function closure`.
- Observable: injury marker; contractility; conduction time; muscle force.
- Discriminating test: Test shared-antigen versus distinct tissue-clone models; human 2024 non-recognition prevents universal myosin assumption.

## I17: Endocrine, hepatic, pulmonary, renal and skin injury

- Constituents: organ-specific immune cells; parenchyma; ducts/barriers; repair capacity.
- Equation skeleton: `one organ-specific injury/repair port per compartment with exposure and function mapping`.
- Observable: hormone output; clearance; gas exchange; filtration; barrier.
- Discriminating test: Build an organ matrix and implement one discriminating mechanism; never use a common damage scalar as all organ endpoints.

## I18: Target-mediated drug distribution

- Constituents: free drug; receptor; complex; internalization; tissue transport.
- Equation skeleton: `dB/dt=kon*C*R-koff*B-kint*B; drug amount balance includes bound pools`.
- Observable: free drug mol/m3; occupancy; clearance m3/s.
- Discriminating test: Compare explicit binding/turnover with Hill exposure response; detect dose-axis arbitrariness in inherited checkpoint model.

## I19: Selective versus systemic immunomodulation

- Constituents: IL-6R; IL-1R; TNF; alpha4beta7; CD80/86; JAK; glucocorticoid receptor.
- Equation skeleton: `distinct target engagements feeding cell-specific signaling and trafficking ports`.
- Observable: tumor kill; tissue injury; pathogen clearance.
- Discriminating test: Specify one matched-exposure comparison with uncertainty; effects of combined interventions cannot be assigned to one drug from observational outcome alone.

## I20: Infection and reactivation during inflammation control

- Constituents: pathogen; effector cells; barrier leakage; latent reservoir; treatment.
- Equation skeleton: `dP/dt=replication-clearance(E,Ab)-drug_effect; activation/reactivation transitions`.
- Observable: pathogen burden; tissue injury; clearance time.
- Discriminating test: Construct high inflammation plus weak clearance as a counterexample to one-axis immune activity; preserve species/pathogen-specific closures.

## I21: Allergy, mast cells and hypersensitivity

- Constituents: IgE; Fc receptors; antigen cross-linking; mast cells; histamine; complement.
- Equation skeleton: `binding/cross-linking -> degranulation -> mediator transport and clearance`.
- Observable: mediator release mol/s; vascular/airway response.
- Discriminating test: Separate immediate hypersensitivity from T-cell-mediated irAEs under the same gross symptom proxy.

## I22: Resolution, fibrosis and persistent effects

- Constituents: apoptotic cells; efferocytosis; macrophages; fibroblasts; ECM; memory.
- Equation skeleton: `injury/recovery and matrix deposition/degradation with persistent cell states`.
- Observable: collagen mass; tissue stiffness Pa; functional recovery.
- Discriminating test: Test whether acute inflammation normalization can coexist with lasting function loss; maintain separate active disease and accumulated damage.

## I23: Microbiome, metabolism and immune response

- Constituents: microbial strains; metabolites; barrier; dendritic cells; tumor immunity.
- Equation skeleton: `metabolite production/transport + host receptor response with explicit causal alternatives`.
- Observable: metabolite mol/s; barrier flux; tumor response.
- Discriminating test: Reuse Q168; distinguish microbial abundance from metabolic function and permeability; mouse commensal result is not human efficacy.

## I24: Experimental methods and causal observability

- Constituents: flow cytometry; TCR sequence; spatial omics; cytokines; organ-on-chip; functional assay.
- Equation skeleton: `y=H(x,sampling,assay)+noise; represent shared cohort and batch dependence`.
- Observable: identifiable parameters; prospective prediction error.
- Discriminating test: Choose a minimal measurement panel distinguishing tumor efficacy, tissue toxicity and infection protection; hold out study/perturbation rather than rows.

