# Measurement specification after round 4

Status: PROPOSED_MEASUREMENT / PENDING_INDEPENDENT_REVIEW. No measurement, order or data collection has started. Consumer: LANE_SURGICAL_INCISION; parents GRAPH_INVERSE, MULTIPHYSICS and BT-CTX-SURG-COLLAGEN. SOURCES_R4.json separates measured analogues from missing skin parameters.

**Decision the measurement should enable:** distinguish interface work from fiber rupture, bulk/leg work and blade contact, and test whether the same measured state gives both Mode I/III tearing work and an unfitted cutting prediction. First test whether delamination has sufficient *newly created area*. A photograph of fiber counts does not identify this.

## Specimens and blocking

For an initial pilot: fresh, unfixed juvenile porcine full skin from **six independent animals**, same anatomical region, documented age, slaughter time and longitudinal/transverse orientation. Six is a chosen pilot count, not a statistical power calculation. Transport on ice; no freeze cycles; specimen time after slaughter recorded. PBS and 22±1 °C are proposed defaults. Document water content, thickness, prestrain and bath time; weigh before/after. Measure actual thickness in each work segment, including epidermis and dermis separately. 1,5 mm is only a calculation example, not a prescribed specimen thickness.

Take matched adjacent specimens for (a) full-skin tearing Mode I and III, (b) isolated dermis, (c) DEJ-peel, (d) local fiber bridge/interface, (e) blade series and (f) heat treatment. A destroyed specimen cannot be reused as virgin control. Same animal is not same specimen: retain position and shared animal blocking and measure differences between adjacent specimens. No species/donor marginals may be treated as a joint independent posterior.

Primary mechanical speeds: 2, 20 and 200 mm/min, a chosen series to expose rate dependence, not skin parameters transferred from Tang2025. Record actual crack speed separately from grip/blade speed. Primary scenario is unheated tissue; heat is a separate intervention.

## 1. Matched work ledger and crack area

Use published Mode I/III geometry from Pissarenko2020 when reproducing their targets; document grip dimensions, incision, orientation and leg deformation. Another geometry creates a new target. Measure force, grip motion, crack front and local opening simultaneously. Advance the crack in **0,5 mm** segments and observe at least five segments after initiation. All area is specified in reference configuration with **one** geometric fracture area, A_crack=t·Δa for a planar through-crack; two new surfaces do not automatically give a factor of two in G.

For each segment record:

`G_total = (ΔW_external − ΔΨ_recoverable − ΔD_legs/bulk − ΔW_contact)/ΔA_crack`.

`G_total = Gamma_severance + Σ_j gamma_j·ΔA_newly_detached,j/ΔA_crack + G_other`.

Measure reversible leg work through simultaneous unload/reload sections before propagation and a separate intact-leg control. Measure grip slip and machine compliance. Bulk cycle work and local interface work must not be subtracted/added twice. No recovery at a single time point does not prove all work is terminal dissipation; the hold/relaxation curve is needed. Measure delamination/bridge work as *alternative resolutions of the same ledger*, not two extra contributions.

Proposed instrument level: force accuracy ≤1 mN for cutting and small bridges; choose a separate load cell handling tearing force and ≤0,5 % of signal. Log 1 kHz at 2/20 mm/min and 2 kHz at 200 mm/min. Crack-front position must have verified absolute error ≤5 µm per 0,5 mm step; thickness error ≤1 % (15 µm in the example). Cameras at least 100 frames/s at 2/20 and 1000 frames/s at 200 mm/min; actual front movement between frames ≤5 µm is the acceptance requirement. Calibrate pixel scale, time synchronization and motion blur under load.

At t=1,5 mm and Δa=0,5 mm, ΔA=7,5·10⁻⁷ m². An error of 5 J/m² corresponds to **3,75 µJ** in segment work or **7,5 mN** in steady cutting. The latter does not approve a 7,5 mN sensor: contact, released energy and area must also share the budget. 2 % of the Mode III reference corresponds to 0,309 mJ per segment, a proposed budget for global work accounting, not identifying the much smaller cutting intercept.

## 2. Native DEJ and dermal interfaces to full separation

DEJ: 10 mm wide T-peel strip, small mechanically initiated separation, **no** enzyme, salt, heat or chemical pretreatment weakening the interface. If native DEJ cannot be initiated without other damage, record it and classify the test as assay failure; do not replace the value with tissue-cultured DEJ. Post-test histology/TEM must locate the fracture path. Measure F, width b, peel angles, crack length a, both leg lengths and grip separation Δ. Calculate `gamma_eff=(∫F dΔ − ΔΨ_arms − D_arms)/ΔA_int`, with observed ΔA_int. Ideal T-peel `2F/b` is a control, not automatic conversion of all published N/m values. Proposed peel-force accuracy ≤0,1 mN, width error ≤0,05 mm and measurable leg compliance. This can detect low J/m² levels without the instrument floor dominating.

Dermis: locate actual fiber/bundle contacts with microscopy. If no repetitive lamellar planes exist, stop calling the specimen interlaminar. Use microtension of a marked fiber bridge or a local interbundle split specimen. Measure the **full** force–local-opening curve until force returns to baseline, reloading/relaxation and initial contact/embedding area. Archive separation versus fiber rupture as distinct fracture paths. Microdissection can create damage; document initial damage and run adjacent preparation controls.

For an initial fiber-bridge specimen, 0,3 mm ligament width and 10×30 mm specimen can be used per Tang2025, with undamaged adjacent compliance control. Their normalization concerns ligament cross-section, not interface area. Require local opening accuracy **≤0,8 µm** to distinguish 5 J/m² at assumed traction of 6 MPa. That traction is an earlier synthetic calculation example; replace with actual measured traction and recalculate the requirement. A camera with 10 µm resolution can map millimeter bridging but cannot alone distinguish a skin intercept of 150 J/m² with this budget. Integrate the measured traction curve; assume no triangular curve or friction constant.

## 3. Mikrostrukturens aktiva areacensus

Measure `M_j=ΔA_newly_detached,j/ΔA_crack` **independently of G_tear** in the same load/history process and both modes. For each contact store initial reference area, orientation, contact identity, opening/slip and first irreversible separation. A surface is paid once on separation; remaining work/history follows the contact. Slip without new area needs its own measured work port. Macroscopic DEJ, interbundle and interfibrillar areas have different denominators and remain separate.

For micrometer fibers, SHG/TPEF with calibrated **≤0,5 µm lateral and ≤1 µm axial resolution** is proposed, registered steps before/after loading and from both sides. Place the volume in reference coordinates and test at least two independent segmentations. This resolution target is an instrument specification, not verified full-skin imaging. **SHG does not automatically see an entire 1,5 mm thick dermis**: measure actual imaging depth and report the unseen core. The tearing-tip region's field-of-view must cover at least 5×5 mm in-plane, full accessible thickness and control regions; expand until it does not truncate newly detached area. Patch width 5 mm is a starting size, not a process-zone bound.

Fibrils around 82 nm require electron tomography/serial blockface with target **≤20 nm lateral and ≤40 nm axial sampling** in small registered volumes, e.g. 20×20×20 µm, stratified across depth and crack distance. Optical layer images do not count as resolved nanocontact area. Preparation/fixation may alter spacing and create false cracks; undamaged controls and correlative registration are needed. Destructive specimens at different arrest points are separate specimens giving an ensemble, not observed continuous history in one specimen.

Require total **≤10 % deterministic area error** from segmentation, registration, unobserved core and sampling extrapolation. Sample size is determined from the pilot's measured heterogeneity; six animals do not prove 10 %. If unseen area lacks a usable upper bound, M becomes UNKNOWN. Diameter 2,23 µm, thickness/d or number of visible intersected fibers must not replace this census.

## 4. Separate blade contact and transition

On adjacent matched specimens: edge radii 1, 10, 50, 100 and 500 µm, same documented wedge angle and prestrain. Measure radius before/after in calibrated profile/SEM; for the 1 µm edge target ≤0,1 µm radius error. Choose 20 mm/min as primary comparison and 2/200 mm/min as speed control. Record local blade arrival, opening, contact pressure/friction or independent contact work, and reinsertion work curve. Reinsertion alone does not identify all previous contact/damage work.

From INCISION R2: `Gamma0=F_cut/t+G_release−G_contact−W_before`. All terms must share area convention and specimen history. For Γ0≈150 J/m² and 20 % total budget, all contributions together must stay ≤30 J/m². Earlier synthetic acquisition cutoff±10, contact coefficient±5 and bridge work±5 already gave ±15,125 J/m²; remaining uncertainty allowance was14,875. This is a design warning from the exact algebraic model, not sensor validation. Define cutting/tearing transition radius in advance, e.g. radius at 50 % of separately measured extra bridge work. Freeze definition and test radii before fitting; hold out at least one radius and one animal.

## 5. Thermal intervention and cell port

On matched specimens: 55/60/65 °C with locally measured temperature history and 0,1/1/10/100 s duration where the instrument can actually deliver these histories; reject unmeasured nominal pulses. Record shrinkage, water loss, SHG/DSC and full mechanical curves before/after on separate matched specimens. Proposed temperature precision ≤0,15 °C and time synchronization ≤1 % of exposure time. Earlier chicken SHG kinetics at60 °C gives d ln k/dT≈0,52/K; precision is a design reference, not skin calibration. Require joint A/Ea-covariance and residual mechanics; native fraction alone is insufficient per R2.

Cell survival requires separate typed viability measurement in epidermis/dermis at registered local load/temperature, with undamaged controls and specified time after load. SHG loss or ECM delamination is not cell membrane rupture. Numerical cell thresholds remain null until such joint acquisition exists.

## Acceptance, control and raw data contract

Primary energy/area gate: measured upper `Σ gamma_j M_j + Gamma_cut` below measured lower G_tear falsifies *the specified interface-only ledger*. No larger M may be rescued by fitting. If fitting requires different M_I/M_III, the difference must have been observed before prediction. Missing native gamma, active areas or separate bulk/fiber work gives UNKNOWN.

Freeze cross-prediction target at ≤25 % error for both Mode I and III on the same cohort with heldout specimens. For dominant product work, 10 % error each in gamma and M gives up to21 % product error. Together with20 % error in the small cut term and3 % other budget, total becomes23,98 % of both earlier mean references: by design within25 %, **outside10 %**. These error bounds are requirements, not achieved uncertainties. Inadequate model closure, species transfer and confounding are not included as a free remainder.

Strongest control: an ordinary mode/speed/history-dependent cohesive law gets the same work curves, microscopy and heldout blocking. Same area ledger is algebraic TIE. Count laboratory acquisition, prep, imaging, storage, setup, fit, refinement and fallback; no cost advantage demonstrated.

Deliver raw force/displacement/time/temperature/crack-front in CSV/HDF5, calibrations, image volumes with pixel scale, reference registration and unique contact masks; specimen ID, animal block, layer, orientation, fluid and speed in JSON. Then deliver joint parameter/posterior or deterministic intervals with identities/covariance and area/energy convention. Data goes into a new local results directory and independent review; nothing is ordered or admitted by this specification.
