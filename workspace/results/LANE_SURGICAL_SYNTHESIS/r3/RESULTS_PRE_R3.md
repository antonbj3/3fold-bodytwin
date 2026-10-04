# Surgical synthesis: executable capability and open measurement ports

Tools, incision path and tissue stack can now be run through a declared scenario chain to vessel-radius-dependent bleeding, persistent hemostatic memory, oxygen edge, healing inventories and directional strength proxy. **51/51 tests pass. Joint biological prediction is UNKNOWN; strongest equally informed control is TIE.** Native rupture strength remains null. [Run instructions](surgical_chain/README.md), [33-port table](PORT_TABLE_R1_FINAL.md), [combined measurement specification](MEASUREMENT_SPEC_R1.md) and [figure/PDF](CHAIN_FIGURE_R1.pdf).

## Night 1/10, round 1

**Capability before:** Four separate lanes had measurement/mechanism ports and negative tests. INCISION/BINDINGS lacked native scalpel toughness and biological/vascular injury width; GAP lacked a mode-specific bridge inventory; RESPONSE's joint chemistry/strength missed the world. The new synthesis lane lacked its own checkpoint and code. Existing models and frozen observations are reused.

**New operation:** [surgical_chain/](surgical_chain/README.md) consumes ports with dimensions and MEASURED/SYNTHETIC/UNKNOWN, source, uncertainty, assay scope and shared scenario ID. Evidence-default leaves outcomes null. An explicit scenario runs prescribed gap and **separate** biological/perfusion widths, layer work, radius histogram, state-dependent shear, frozen coagulation/platelet kinetics, O2 inventory and FV history, as well as Q_U/Q_I/Q_M and a directional rupture proxy. Full tool work is unknown when contact/release/pre-tip-work is missing. Radius/speed is registered but has no calibrated law to gap/injury.

A new persistent per-bin seal-state retains formed closure when thrombin later falls. The first construction that forgot seal is preserved in `r1/scenario_nominal`; new PREREG and v2/v3/v4 preserve the construction change. Seal changes a **synthetic atmospheric conductance**, never the perfusion of cut vessels. U/I/M are downstream observers; they do not replace the previous FV C-state that feeds back to O2. Final RESPONSE R4 is included as a separate conservative D/H/A branch with explicitly unknown chemistry→native traction. This is not an identified joint collagen law.

**What holds against the world, with scope:**

| Basis | Locked result | What it actually supports |
|---|---|---|
| [INCISION R4](sources/LANE_SURGICAL_INCISION/r4/fracture_transfer.json) | Same needle across speeds at most2,8409% error; new tools35,96/22,14/54,64% error | Tool-specific effective paired-work assay, not universal Γ_cut/Γ0 |
| [INCISION R3](sources/LANE_SURGICAL_INCISION/r3/geometry/needle_geometry_results.json) | Held-out needle geometry11,4286% central,38,3333% extreme reading | Bounded force/geometry capability; no independent skin-scalpel port |
| [BINDINGS R4](sources/LANE_SURGICAL_BINDINGS/INTERFACE_BUDGET_RESULTS_R4.json) | Synthetic microbudget ceiling6965,37J/m²; native interfacearea/law null | Refutes the specified borrowed-interface/V-over-A construction, no universal upper skin bound |
| [GAP R2](sources/LANE_SKIN_TOUGHNESS_GAP/MECHANISM_TABLE_R2_FINAL.json) | Hydraulic≤0,11047kJ/m² in chosen volume/modulus; crimp-box≈0,14028; mm bridge24,5/19,1375kJ/m² | Conditional narrowing of magnitude; mode split and native terminal work still unknown |
| [RESPONSE R3](sources/LANE_SURGICAL_RESPONSE/r3/crosslink_port/summary.json) | HP-input-proxyRMSE8,30637pp; nominal gateFAIL. Early chemistry predictsHP42≈0,020 against0,043 | 8,31pp is externally supplied HP + previous C, **not** U/I/M's joint early chemistry prediction |
| [RESPONSE R4](sources/LANE_SURGICAL_RESPONSE/r4/SCORES_R4_v2.json) | HP inventory42/10=2,22026; denominator-onlyFAIL;0/16glyco profiles pass the joint gate | Pure denominator explanation fails under fixed assay reference; nativepool and kinetic order missing |

The targets come from different species/cohorts/assays; no joint posterior has been created. **The night log's bleeding direction is corrected locally:** the raw file's count-matched control **overestimates**25,3757×; it does not underestimate. The same source's equation text has R³ for shear, while the code correctly uses R. [Corrections and preserved originals](SOURCE_CORRECTIONS_R1.json).

**Strongest control:** Separately calibrated curves per observable get the same observations and can interpolate them; this is no held-out cross-prediction. `SeparateCurve` leaves new interventions unknown without an extra law. The chain can give a **conditional** coherent answer when several scenario parameters change, but the corresponding conventional FV/cohort/reaction/Poiseuille/history control has the same state, information and operations: **TIE**, no unique gain. New native cross-prediction is still UNKNOWN. The count-only map and factorized chemistry-angle are information-loss diagnostics, not the strongest control.

**Actual outcome:** [Main run](r1/nominal_complete_v4/summary.json) has273cells. At120min, synthetic edge pressure is54,2560Torr, blood118,768µL and remaining flow-weightedseal0,617795. One vessel bin has first closure4,49min; the whole wound is not closed within120min, so hemostasis time is null/right-censored. The strength proxy on day7/14/28/42/90 is0,172/2,068/10,019/16,698/22,945%. These are model scenarios, not biological measurements.

The integrated strength proxy gives **43,7037pp RMSE: FAIL** against the preserved, unmatched rat-incision reference, compared with previous legacyHill14,91pp and the different HP-input proxy8,31pp. No parameter was rescued by a new strength fit. [All comparisons](r1/INTEGRATED_STRENGTH_DIAGNOSTIC.json). The standalone R4 branch givesHP42=0,02633845 with the same cell-free controlrate/yield and synthetic D10; native chemicalpool/rupturelaw is null.

**Tested mathematics:** [51/51 tests](VERIFICATION_R1.json) lock positive and negative source results, strict unit/null policy, r⁴-flow/r¹-shear, zero flow/zero birth/zero bridging, pathway blockade, healthy-O2, PSD/inventories, causal seal and checkpoint. U/I/M balance error≤2,89e−15; early relative O2 inventory balance8,78e−11. [Refinement](REFINEMENT_R1.json): halved time steps give≤0,01308pp strength difference and0,00175Torr; halved edge step≤0,00172pp and0,27461Torr. Restart from day21 matches the full prefix-history run within7,11e−15pp and D/H/A within1,39e−17. The prefix is not replayed; the entire spatial suffix is solved again. No certifiedlocalquery/globalstate-free cost capability follows.

A seal→atmosphere ablation changes the day90 proxy2,63117pp, with exactly equal early fields before seal. Three joint synthetic scenario directions give20,1115–23,8522%; the interval is not a CI or joint biological posterior.

**Reasonable defaults and cost:**10mm×2mm incision, dermis1,5mm/subcutis0,5mm, explicit Γ150J/m² in both synthetic layers, gap200µm and separate injury/perfusion widths250µm; radii8/25µm with density1,18e6/2e4m⁻², Δp4000Pa, μ0,0035Pa·s, vessel path1mm. Domain8mm, edge step5µm,120min/90days, dt0,1min/0,1days; Q036/FV/coagulation defaults are openly labeled closures in [CLOSURE_MANIFEST](CLOSURE_MANIFEST_R1.json). Main solverwall2,284s; six refinement/scenario cases18,307s total. All earlier constructions, reruns and restarts are in [cost accounting](COST_ACCOUNTING_R1.json). Import, development, actual acquisition and the entire work round are unknown separate costs, not zero. No10×gain is claimed; twoCPUthreads, no heavy jobs.

**Graph binding:** Four dispatch+feedback have been executed with frozen standard receipt code in [the lane's local review workspace](graph_review_workspace/bodytwin/tasks/graph_runs/), one per BT-CTX-SURG-INCISION/COLLAGEN/HEMOSTASIS/HEALING. Kindreview, all negative results remain, scientific_admission=false. The workspace's standard CLI writes `tasks/graph_runs` and is outside the task's write boundary; **global context import is therefore PENDING_COORDINATOR**, not executed or claimed. Source graph/code/services are untouched.

**Remaining obstacles:** Native tool→gap/biological/perfusion injury, radius/connectivity/flow map, shear/seal/erosion→O2, absolute dryvolume/Hb, native glyco-pool/turnover and connected chemistry-age-angle→traction/rupture are missing. The pilot ports and proxy chain give engineeringPASS, **innovation controlTIE**, integrated strength proxyFAIL and empiricaljointUNKNOWN. GRAPH_INVERSE/MULTIPHYSICS are OPEN. Everything **PENDING_INDEPENDENT_REVIEW**.

**Next construction change:** [Prioritized measurement specification](MEASUREMENT_SPEC_R1.md): start with a registered incision that gives gap, separate cell viability/perfusion and vessel radius/pressure/flow on the same clock. It opens three early transitions; a new isolated Γ value does not. Then dryvolume/O2 and chemistry×connectedbridge with separate ultimate-load. [NEXT_ROUND](NEXT_ROUND.md) gives an immediately executable start and falsifiers.


## Night 1/10, round 2

**Capability before:** R1 gave an executable scenario chain and51 green regressions but strength-proxy RMSE43,7037pp. R3's8,3064pp was a separate component diagnostic; prioritization in the measurement specification was qualitative. Native strength, full hemostasis time and edge necrosis lacked a joint observation law.

**New operation:** [Exact port attribution](ATTRIBUTION_R2.md) replaces one port at a time and runs the two-factor interaction. [Quantitative measurement specification](MEASUREMENT_SPEC_R2.md) propagates256 frozen synthetic scenarios with declared measurement precision, shared calibration error and cost; after a linear hypoxia miss, a positive nonlinear observation operator is executed on64 new scenarios. Physics/biology in surgical_chain is unchanged and verified againstR1's code hashes. [Figure/PDF](DIAGNOSIS_FIGURE_R2.pdf).

**Strongest control:** The same conventional FV/cohort, port intervention, optimal linear forecast and positive kernel operator get the same observations, prior, noise/covariance and setup: **TIE**. No new global or local cost gain is claimed. The separate R3 curve must not be treated as the same forward chemistry: `75 C_R3 x_HP` uses external HP and inherited amount calibration against older strength3% day7/20% day21. R3 made no new strength fit, but C_R3 is not strength-free measured collagen information. [Source corrections](SOURCE_CORRECTIONS_R2.json).

**Actual outcome — which port carries the strength error:** ΔRMSE is reduction from43,703696; negative numbers worsen it.

| Port change | RMSE pp | ΔRMSE pp |
|---|---:|---:|
| Amount only A→R3 C | 13,000632 | +30,703064 |
| Maturation only m→R3 HP | 39,456781 | +4,246914 |
| Both A and m | 8,306367 | +35,397329 |
| Amount A→the chain's own FV C | 29,560901 | +14,142795 |
| Early state/hazard→R3 | 45,697008 | −1,993312 |
| k_deposit→R3 | 59,827729 | −16,124033 |
| All upstream R3 ports | 60,007651 | −16,303955 |

The full table includes gap, separate widths, face-conductance, legacy maturation, D/H/A-HP, bridge and load-angle in [ATTRIBUTION_R2.md](ATTRIBUTION_R2.md). Two-factor Shapley gives **30,926739pp amount and4,470590pp maturation**; the sum35,397329pp and interaction−0,447351pp are preserved. The identity `S=75 A m` reproduces R1 within1e−11pp and both port changes reproduce8,306367pp without a new fit. This is a diagnostic reconstruction; proxy/native gates are not upgraded.

**Located obstacle:** FV C and U/I/M have different reference inventories. U/I/M formation uses `(1−C_FV)` and adds an extra immature-loss on top of the FV flows; internal U→I→M fluxes cancel but ΣUIM is not FV C. Day90 has FV C0,640170, ΣUIM0,313012 and R3 C0,998320. Importing k_deposit≈5/day fills FV C≈1 but leaves ΣUIM0,009881: the birth port closes when the wrong inventory fills. All upstream R3 ports reproduce R3's saved FV C but the U/I/M observer gives only0,620360% strength day90. [Equations and raw state](r2/attribution_v1/INVENTORY_DIAGNOSIS.json).

**Information value per measurement:** With open default costsM1/M2/M3/M4/M5/M6=8/3/5/12/16/20 relative work units, the moment order becomes **M2 oxygen/flux → M1 registered incision/vessels → M5 connected bridges/mechanics → M3 plug → M4 chemistry → M6 work ledger**. Prior-SD for strength proxy28/90 is2,8418/6,7295pp; for injury-proxy width196,680µm. M5's conditional direct observer port gives1,2791/4,2652pp, M1 gives89,356µm injury-SD and M2 gives116,851µm. Score/costM2/M1/M5 is0,06325/0,02649/0,02305. The table with **every R1 entry and its subobservations**, chosen precision, shared measurement error, SD reduction,0,5–2× noise/cost and unresolved law floor is in [MEASUREMENT_SPEC_R2.md](MEASUREMENT_SPEC_R2.md). Actual lab hours/prices are null; factor-cost intervals prove no unconditional top ranking.

M4's chemistry and M6's Γ have **zero implemented** path to the requested outcomes. Their value for a future identified law is UNKNOWN, and absolute collagen moles/connected traction are still load-bearing for native strength. M5's numerics assume that the assay can actually observe synthetic Q/bridge/reference ports; this is a conditional information scenario, not calibrated native precision. [All represented and missing subentries](MEASUREMENT_COVERAGE_R2.json).

All320 scenarios are right-censored for full hemostasis at120min. SD0 for `min(T,120)` is a constant model observable, **not** certainty about later hemostasis. Full time-SD and its information value are null. Edge-necrosis SD is also null: Q036-hazard>0,20 and O2<10Torr are declared proxies, not a viability law. No numerical biological posterior is fabricated.

**Preserved information gate and executed construction change:** The linear128→128 test worsened hypoxia-MSE68,16% forM2. The new positive kernel operator improves hypoxia-MSE39,39% forM2 and19,42% forM1 on64 new scenarios with32 measurement-noise draws per scenario. But M1's injury-RMS becomes121,984µm against linear control92,292µm; the frozen joint gate is therefore **FAIL**, not PASS. New risk/cost orderM2→M5→M1→M3 is preserved together with the original failure. M5 gives strength-RMS1,3256/3,9914pp on new scenarios against prior2,6807/5,2368pp. A shared kernel width is therefore not a free joint uncertainty operator for all targets. [Frozen test and raw residuals](r2/nonlinear_information_v1/VALIDATION.json).

**Verification and cost:** **70/70 tests PASS**, including all51 fromR1, attribution identities, a new chemistry-intervention countertest, covariance/noise, shared calibration error and preserved information misses. [VERIFICATION_R2](VERIFICATION_R2.json). Four halved-time-step tests give at most0,011149pp strength difference and0µm injury-width difference; numerical gatePASS.256+64 full scenario solutions cost851,753+253,888s solverwall; four refinement cases34,272s; attribution and test/training/evaluation costs are reported separately in [COST_ACCOUNTING_R2](COST_ACCOUNTING_R2.json). Highest reported RSS336792KiB, threads2. Full development/acquisition/setup cost is UNKNOWN, never0. Earlier misses, initial kernel training and allR1 raw files are preserved.

**Graph binding and status:** The coordinator's fourR1 imports have been verified read-only. Four newR2 dispatch+feedback exist in the lane's local review workspace; [R2 receipts](GRAPH_BINDING_STATUS_R2.json) and concrete coordinator records are delivered for import. No change to source graph, canonical code or services. EngineeringPASS, attributionPASS, integrated strength proxyFAIL43,7037pp, information pivotFAIL, strongest controlTIE and empiricaljointUNKNOWN. Everything **PENDING_INDEPENDENT_REVIEW**; GRAPH_INVERSE/MULTIPHYSICS are OPEN.

**Next construction change:** A single declared reference amount and shared formation/removal ports must carry FV C and the U/I/M partition. Acquire absolute collagen amount on fixed area/volume and matched early O2/viability, then chemistry×connected traction. The information operator needs separate target-adapted observations/bandwidths and a new test; a shared kernel must not hide M1's loss. [Immediately executable continuation](NEXT_ROUND.md).
