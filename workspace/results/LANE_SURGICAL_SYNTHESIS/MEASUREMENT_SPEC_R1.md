# Combined measurement specification: open the couplings that carry the chain

PENDING_INDEPENDENT_REVIEW. This is a prioritized acquisition specification. No experiments are ordered. Prioritization is qualitative: value means which downstream ports an observation can identify, cost means relative instrument/specimen burden. The budgets below are requirements, not achieved precision or power calculations.

**First measurement package:** registered incision with actual gap, separate cell viability and perfusion, and vessel radii/pressure/flow on the same clock and in the same skin model. It opens injury zone→bleeding→oxygen edge, the earliest broken link. An isolated new Γ value does not open this transition.

| Priority | Acquisition and unlocked port | Value per cost / prerequisite | Decisive control |
|---|---|---|---|
| 1 | Before/after incision: x,z,Ft,Fn, actual front/gap; cell viability and perfusion loss separately; vessel radius histogram, connectivity, before/after flow and pressure per layer | Three early links. Synchronous optical/perfusion imaging with calibration; moderate–high cost. Radii must be measured, not just counts | Held edge/donor. Count-matched control is information-weak; strongest control gets the entire vessel map and pressure boundary conditions |
| 2 | Oxygen profile and atmospheric flux, Hb/blood advection, wet volume and dry mass, outer skin pressure and sensor response; same early wound model | Makes absolute O2 consumption and transport identifiable. Oxygen sensor + weighing relatively cheap; spatial flow more expensive. Requires step1 | Same observed pressure/flux boundary conditions in conservative FV/Hb control. Surface Robin-law should give actually measured skin pressure |
| 3 | Radius-matched plug dynamics: adhesion/aggregation, fibrin, local r(t), Q(t), γ(t), seal and erosion under fixed pressure boundaries | Identifies shear→seal and seal→O2. Flowchamber and synchronous imaging; moderate cost. Requires known geometry/pressure | GPIIb/IIIa blockade separate from adhesion; same state-driven shear control. Formed plug should not disappear solely when thrombin falls |
| 4 | Absolute GG/G/free-DHLNL, HP, other aldehyde product, catalyticLOX, collagen moles and fixed reference area/volume at10/14/21/42/61/90; verified LOX blockade from10 | Distinguishes maturity from denominator/sampling. Mass spectrometry + cohort design: higher cost. Independent assayday14 is a cheaper first discriminant for two rate laws | D+2H+A with separate import/loss. First/second-order law control gets the same early data. No fit of nativepool to HP42 |
| 5 | Chemical age×direction×connected crossings over incision front; anchor survival, local full traction–opening and ultimate load and reference area | Opens strength prediction. 3D-SHG + local mechanical assays: high cost but the strength law’s load-bearing information. Requires step4 for chemical identification | Same conventional cohort/cohesive-history. Measure low strain and failure load separately; HP or Fourier eccentricity alone is insufficient |
| 6 | Full Ft/Fn work ledger with 3D newly created area, unload/reload/repeat/precut drift, release and pre-tip-work; active new interfaces/fullterminal-work in ModeI/III | Opens new tool→work and micro→tearing. High cost. Can be run jointly with step1 but nanocensus is separate | Held actual radius and donor; same contact/cohesive law. SingleΓ for all needles should retain its previous FAIL |

## Shared lineage and current default choices

Start from a shared species/site/age/layer definition; do not use mouse excision, rat incision, human dermis, juvenile porcine full skin and cell-free culture as a joint parameter posterior. Reuse public data where its raw curves/metadata actually fill the port. Same animal is not same specimen. Destructive chemical/mechanical assays need registered specimen IDs and animal/cohort blocks; connect them through observed blocking and keep sampling/confounding explicit. Ex vivo mechanics alone cannot give a90-day healing series.

Six independent animals, PBS22±1°C and adjacent registered specimens are earlier BINDINGS/INCISION pilot proposals, not a power or precision guarantee. Tool series deliberately differ: INCISION proposes1/5/20/50µm and BINDINGS1/10/50/100/500µm; start with shared1/50µm plus a held radius and documented wedge angle. INCISION’s1/20mm/s and BINDINGS2/20/200mm/min are different speeds; a shared20mm/min case can link mechanical assays, while the needle result’s20mm/s is not transferred. Actual radius/speed drive the analysis; no conflation of shaft diameter and edge.

## Accuracy and instrument floor

For incision work, a geometric crack-plane area is used, not the sum of both free surfaces. With Ahat>εA:

`εΓ0 ≤ (εW+εΨ+εcontact+|Γeff_hat| εA)/(Ahat−εA)+εB`.

INCISION’s design |Γeff|≤3500J/m², area0,2%, work/release/contact8/5/5J/m² and terminalB4J/m² gives **29,0501J/m²**. At planned5×1,5mm incision, work8J/m² means60µJ and tangential constant-force error only12mN. BINDINGS1%thickness+1%crackadvance is insufficient for this0,2%area budget; the different budgets apply to different observables. If the requirement cannot be measured, Γ0 remains UNKNOWN.

For tearing energy, the previous25%design budget is another contract: gamma10% and active-area10% give product error21%; with other terms it becomes23,98%. Local opening≤0,8µm at *synthetic*6MPa traction is needed for5J/m² work resolution. BINDINGS100/1000fps starting choices do not replace the requirement≤5µm observed frontadvance/frame. Imaging depth/unobserved core and unique newly-detached area must be counted; a visible fiber count is no active area census.

For the radius port, Poiseuille gives `Q∝r⁴`: a5%radius error gives up to21,55%flow error even before count/pressure/length errors. The distribution tail should have an explicit upper censoring bound, not the image field of view as maximum radius. `γ=Δp r/(2μL)` has linear radius dependence. In the swarm’s assumed map, count-matched loguniform control overestimates25,3757×; it is not a measured biological gain.

The O2 target≤25µm hypoxic width and strength’s previous±5pp-plus-declared-envelope are prospective constraints, not clinical thresholds. TwoCPU threads and5µm edge grid are numerical defaults, not a measured biological resolution requirement. The measurement threshold10Torr for the hypoxia proxy does not identify cell necrosis.

## Raw data contract and acceptance

Deliver CSV/HDF5 with timestamp, force/movement/front, r/Q/pressure, localO2 and sensor response; JSON with species/site/donor/specimen/layer, reference area/volume/mass, actual radius/speed/prestrain/hydration and unit. Images should have scale/clock calibration, registered axes and terminalevent-ID. Chemistry should state mole equivalents, internal standards, specimen reference, import/loss and assay sampling. Preserve covariation, batch and shared systematic errors; unknown covariance is null, never an automatic independent product.

Calibrate blade contact from separate work, chemistry from independent assay and traction from separate mechanics. Freeze heldoutdonor/tool and late time points before fit. Approve no joint empirical chain before the same validity domain carries all ports. Strongest control gets the same observations and restart states. Laboratory acquisition, prep, imaging, setup, fit, storage, cert, rebuild and fallback must be counted; missing cost is stated as null.

Basis: [INCISION R4](sources/LANE_SURGICAL_INCISION/MEASUREMENT_SPEC_R4.md), [BINDINGS R4](sources/LANE_SURGICAL_BINDINGS/MEASUREMENT_SPEC_R4.md), [GAP R2](sources/LANE_SKIN_TOUGHNESS_GAP/MEASUREMENT_SPEC_R2.md), [RESPONSE R4](sources/LANE_SURGICAL_RESPONSE/NEXT_ROUND_R4_FINAL.md). RESPONSE had no standalone MEASUREMENT_SPEC; its final files are the measurement specification source.
