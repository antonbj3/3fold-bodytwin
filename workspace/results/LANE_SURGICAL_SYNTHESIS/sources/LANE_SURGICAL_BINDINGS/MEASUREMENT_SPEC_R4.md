# Measurement specification after round 4

Status: PROPOSED_MEASUREMENT / PENDING_INDEPENDENT_REVIEW. No measurement, order or data collection has been started. Consumer: LANE_SURGICAL_INCISION; parents GRAPH_INVERSE, MULTIPHYSICS and BT-CTX-SURG-COLLAGEN. SOURCES_R4.json distinguishes measured analogues from missing skin parameters.

**The decision that the measurement should enable:** separate interface work from fiber breakage, bulk-/benarbete and blade contact, and test if the same measured condition gives both Mode I/III-rivarbete and an unadapted cutting prediction. First test if delamination has sufficient *newly added area*. A fiber count photograph does not identify this.

## Sample and blocking

For a first pilot: fresh, unfixed juvenile porcine whole skin from **six independent animals**, same anatomical region, documented age, slaughter time and longitudinal/transverse orientation. Six is ​​a chosen pilot number, not a statistical power calculation. Transport on ice; no freeze cycles; trial period after slaughter is recorded. PBS and 22±1 °C are suggested default choices. Document water content, thickness, bias and time in the bath; path before/after. Measure actual thickness in each working segment, including epidermis and dermis separately. 1,5 mm is only a calculation example, no forced sample thickness.

Take matched, adjacent samples for (a) full-skin tear Mode I and III, (b) isolated dermis, (c) DEJ-peel, (d) local fiberbrygga/interface, (e) sheet series, and (f) heat treatment. A spoiled sample cannot be reused as a virgin control. Same animal is not same sample: keep position and common animal blocking and measure differences between adjacent samples. No art-/donormarginaler may be treated as a common independent posterior.

Primary mechanical speeds: 2, 20 and 200 mm/min, as a decided series to expose speed-dependent, not transferred skin parameters from Tang2025. Record actual crack rate separately from grip-/bladhastighet. Primary scenario is unheated tissue; heat is a separate intervention.

## 1. Matched work ledger and crack area

Use published Mode I/III-geometri of Pissarenko2020 when reproducing their goals; document grip dimensions, cut, orientation and deformation of the legs. If a different geometry is used, it becomes a new target. Measure force, grip movement, crack front and local opening simultaneously. Allow the crack to advance in segments about **0,5 mm** and observe at least five segments after initiation. All area is given in reference configuration with **en** geometric fracture area, A_crack=t·Δa for plane through crack; two newly formed surfaces do not automatically give a factor of two in G.

For each segment, post:

`G_total = (ΔW_external − ΔΨ_recoverable − ΔD_legs/bulk − ΔW_contact)/ΔA_crack`.

`G_total = Gamma_severance + Σ_j gamma_j·ΔA_newly_detached,j/ΔA_crack + G_other`.

Measure reversible leg work by simultaneous av-/pload parts before propagation and separate intact leg control. Measure grip slip and machine compliance. Bulk cycle work and local interface work must not be drawn av/adderas twice. Absence of recovery at a single point in time does not prove that all work is terminal dissipation; the hold/relaxation curve is needed. Measure delaminerings-/bryggarbete as *alternate resolutions of the same ledger*, not as two additional increments.

Suggested instrument level: force accuracy ≤1 mN for cutting and small bridges; select separate load cell that can withstand the tearing force and ≤0,5 % of the signal. Log 1 kHz at 2/20 mm/min and 2 kHz at 200 mm/min. The crack front position must have verified absolute error ≤5 µm per 0,5 mm step; thickness error ≤1 % (15 µm in the example). Cameras at least 100 bilder/s at 2/20 and 1000 bilder/s at 200 mm/min; actual front movement between images ≤5 µm is the acceptance requirement. Calibrate pixel scale, time sync and motion blur under load.

At t=1,5 mm and Δa=0,5 mm, ΔA=7,5·10⁻⁷ m². An error of 5 J/m² corresponds to **3,75 µJ** in the segment work or **7,5 mN** in steady cutting. The latter is not an approval of a 7,5 mN sensor: also contact, released energy and area must share the budget. 2 % of Mode III-referensen corresponds to 0,309 mJ per segment and is a suggested budget for global work accounting, not for identifying the much smaller cutting intercept.

## 2. Native DEJ and dermal interfaces to full separation

DEJ: 10 mm wide T-peel strip, small mechanically initiated separation, **no** enzyme, salt, heat or chemical pretreatment that weakens the interface. If native DEJ cannot be initialized without other damage, record it and classify the attempt as assay failure; do not replace the value with tissue cultured DEJ. Histology/TEM after the attempt must locate the breaking path. Measure F, width b, peel angles, crack length a, both leg lengths and grip separation Δ. Calculate `gamma_eff=(∫F dΔ − ΔΨ_arms − D_arms)/ΔA_int`, with observed ΔA_int. Ideal T-peel `2F/b` is a check, not an automatic conversion of all published N/m values. Suggested peel force accuracy ≤0,1 mN, width error ≤0,05 mm and measurable bone compliance. So even low J/m² levels can be detected without the instrument floor dominating.

Dermis: locating true fiber-/buntkontakter by microscopy. If no repetitive lamellar planes are present, stop calling the sample interlaminar. Use micropulling of marked fiber bridge or a local interbundle gap test. Measure **full** force–local opening curve until force returns to baseline, reload/relaxation and initial kontakt-/inbsaturation area. File separation versus fiber breakage as different breakage paths. A microdissection can create damage; document the starting damage and run adjacent preparation checks.

For an initial fiber bridge test, 0,3 mm ligament width and 10×30 mm test can be used according to Tang2025, with undisturbed adjacent compliance control. Their normalization applies to the cross-section of the ligament, not the area of ​​the interface. Require local opening accuracy **≤0,8 µm** if 5 J/m² is to be separated at an assumed traction of 6 MPa. That traction is an earlier synthetic calculus example; replace with actual measured traction and recalculate the requirement. A camera with 10 µm resolution can map millimeter bridging but alone cannot distinguish a skin intercept of 150 J/m² with this budget. Integrate measured traction curve; assume no triangular curve or friction constant.

## 3. Mikrostrukturens aktiva areacensus

Measure `M_j=ΔA_newly_detached,j/ΔA_crack` **independence of G_tear** in the same last-/historiefrun and in both modes. For each contact initial reference area, orientation, contact identity, opening/sliding and first irreversible separation are stored. An area is paid once upon separation; remaining arbete/historia follows the contact. Sliding without new area needs its own measured working port. Macroscopic DEJ, interbundle and interfibrillar area have different denominators and are kept separate.

For micrometer fibers, SHG/TPEF with calibrated **≤0,5 µm lateral and ≤1 µm axial resolution ** is proposed, recorded steps before/after load and from both sides. Put the volume in reference coordinates and try at least two independent segmentations. This resolution target is an instrument specification, not verified full skin imaging. **SHG does not automatically see an entire 1,5 mm thick dermis**: measure actual imaging depth and account for unseen core. The field of view of the tear tip region shall include at least 5×5 mm in the plane, the entire accessible thickness and control areas; expand the field of view until it does not cut off the newly detached area. Patch width 5 mm is a starting size, not a process zone limit.

Fibrils around 82 nm require elektrontomografi/seriell blockface with targets **≤20 nm lateral and ≤40 nm axial sampling** in small recorded volumes, for example 20×20×20 µm, stratified over depth and distance to the crack. No optical layer images count as resolved nanocontact area. Preparering/fixering can change distance and give false cracks; non-destructive controls and correlative registration are needed. Destructive samples at different arrest points are separate samples and provide an ensemble, not observed continuous history in the same sample.

Require total **≤10 % deterministic area error** from segmentation, registration, unobserved kernel, and sample extrapolation. Sample size is determined from the pilot's measured heterogeneity; six animals is not proof of 10 %. If unseen area lacks a useful upper bound, M becomes UNKNOWN. Diameter 2,23 µm, tjocklek/d or number of visible cut fibers must not replace this census.

## 4. Separate blade contact and transition

On adjacent matched samples: edge radii 1, 10, 50, 100 and 500 µm, same documented wedge angle and bias. Measure radius before/after in calibrated profil/SEM; for 1 µm edge target ≤0,1 µm radius error. Select 20 mm/min as primary comparison and 2/200 mm/min as speed control. Record local blade arrival, opening, kontakttryck/friktion or independent contact work, and the reinsertion work curve. Reinsertion does not alone identify all previous kontakt-/skadearbete.

From INCISION R2: `Gamma0=F_cut/t+G_release−G_contact−W_before`. All terms must have the same area convention and test history. For Γ0≈150 J/m² and 20 % total budget, all contributions together must hold ≤30 J/m². Previous synthetic purchase cutoff±10, contact coefficient±5 and bridging work±5 already gave ±15,125 J/m²; left for other uncertainty var14,875. It is a design caveat from the exact algebraic model, no sensor validation. The radius of the cutting/tearing transition is defined in advance, for example the radius at 50 % of the separately measured extra bridging work. Freeze definition and test radii before fit; endure at least one radius and one animal.

## 5. Thermal intervention and cell gate

On matched samples: 55/60/65 °C with locally measured temperature history and 0,1/1/10/100 s duration where the instrument can actually deliver these courses; reject nominal pulses not measured. Record shrinkage, water loss, SHG/DSC and full mechanical curve before/after on separate matched samples. Suggested temperature accuracy ≤0,15 °C and time synchronization ≤1 % of exposure time. The former chicken SHG-kinetiken at 60 °C yields d ln k/dT≈0,52/K; the precision is a design reference, not a skin calibration. Require joint A/Ea-covariance and residual mechanics; native fraction alone is insufficient according to R2.

Cell survival requires separately typed viability measurement in epidermis/dermis at registered local belastning/temperatur, with undamaged controls and specified time after load. SHG loss or ECM-delaminering is not cell membrane disruption. Numeric cell thresholds are left null until such a common purchase exists.

## Acceptance, control and raw data contracts

Primary energi-/areagrind: measured upper `Σ gamma_j M_j + Gamma_cut` below measured lower G_tear traps *the specified interface-only-ledger*. No major M must be saved by fit. If the fit requires different M_I/M_III, the difference must have been observed before prediction. Missing native gamma, active areas or separate bulk-/fiberarbeten gives UNKNOWN.

Freeze the cross-prediction target to ≤25 % error for both Mode I and III on the same cohort with holdout samples. For dominant product work, 10 % error each in gamma and M gives up to 21% product error. Together with 20% error in the small cut entry and 3% remaining budget, the total is 23,98% of both previous average references: design wise within 25%, **outside 10%**. These error limits are requirements, not achieved uncertainties. Lack of model closure, species transfer and confounding are not included as free residuals.

Strongest control: a regular mode-/hastighets-/historieberoende cohesive team gets the same work curves, microscopy and holdout blocking. The same area ledger is algebraic TIE. Count laboratory purchases, prep, imaging, storage, setup, fit, refinement and fallback; no cost benefit has been demonstrated.

Deliver raw force/displacement/time/temperature/crack-front in CSV/HDF5, calibrations, pixel-scaled image volumes, reference registration, and unique contact masks; sample ID, animal blocks, layers, orientation, fluid and velocity in JSON. Then deliver joint parameter/posterior or deterministic intervals with identiteter/covariance and area-/energikonvention. Data to go to a new local results folder and independent review; nothing is ordered or admitted by this specification.
