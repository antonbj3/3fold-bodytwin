# After round 4: measurement handoff, lane closed in this form

Status: LANE_SCOPE_CLOSED_DATA_ACQUISITION_REQUIRED / PENDING_INDEPENDENT_REVIEW. Parent GRAPH_INVERSE/MULTIPHYSICS, BT-CTX-SURG-COLLAGEN. Consumer LANE_SURGICAL_INCISION. Native joint skin prediction UNKNOWN; innovation gate FAIL; equally informed energy ledger TIE. Read the latest steering before this material. No automatic fifth synthetic design is proposed.

Capability → the same independently measured microstructure/history should provide Mode I/III, cutting work and thermal change.
Conflict → large tearing energy requires irreversible full-fracture work and activated area, but the interface values found belong to other tissues/assays.
Obstacle → native dermal gamma, newly activated interface area and separate blade-transfer-work are missing; physical interface count does not give area-ratio.
Changed operation → implement MEASUREMENT_SPEC_R4.md or acquire a publicly available matched raw-data package filling exactly the same fields.

## Start directly

1. Read CHECKPOINT_R4_FINAL.json, SOURCES_R4.json, INTERFACE_BUDGET_RESULTS_R4.json, PORTS_R4_FINAL.json and MEASUREMENT_SPEC_R4.md. Preserve all r1–r4, especially original optimizer/serialization/heat failures and access failures. New acquisition/execution must get a new suffix, not overwrite r4.
2. Prioritize matched unfixed juvenile porcine skin/dermis: same animal, anatomical position, orientation, hydration, speed and measured full-separation-work; cutting/tearing and active interfaces must share lineage. Adjacent specimens are not the same specimen. The six-animal proposal is a pilot, not verified power or precision.
3. A new acquisition must fill `native_gamma_interface`, `DeltaA_newly_detached/DeltaA_crack`, `full_workcurve`, `recoverable_leg/bulk_work`, `fiber_failure_work` and `blade_contact/opening_before_arrival`, with uncertainty and shared specimen metadata. G and fiber counts are insufficient. Consider searching for the identified Tobaruela2016 table/UPM-thesis44197 as a separate analogue; access yielded three bytes, no numbers were acquired.
4. Sources to reuse: Wang2021 native porcine aortic media106,42±25,82/183,78±58,08 J/m² at0,2mm/s with elastic leg correction; Wu2005/2006 SC1–18 J/m² across different regimes; Sierra-Sánchez2023 tissue-cultured DEJ F/b1,47±0,47/2,84±2,48 N/m; Kumar2009 amnion F/b1,7–26,2 N/m. Force/width is not intrinsic work. T-peel2F/b requires its specific leg/area geometry. Gupta2013 antler0,18–0,27 J/m² is a mineralized NCP interface and must not be copied to wet dermis.
5. Tang2025 offers the concrete full-fracture measurement: direct ligament pullout with local compliance correction plus `G(delta_tail)=G_A+integral sigma ddelta` under small-scale bridging. Fresh bovine pericardium gave95,9 kJ/m² extra between reported onset/global failure; it is composite bridge-work, not gamma per microsurface. Use the protocol as a starting point, not the numbers as skin parameters.10µm camera resolution is insufficient for a150J/m² intercept at assumed6MPa traction.
6. Perform an actual active-area-census. Dermis has3D bundle structure; no published active lamellar count was acquired. For a through crack `A_crack=t*da` and parallel delamination `A_int=w*da`: `M=sum(w)/t`, not the number of anatomical interfaces. Follow reference identities and pay for newly detached area once. A full crack-plane cohesive tractioncurve and microinterface work from the same pullout must not be added.
7. Acceptance before fitting: a measured upper `sum(gamma*M)+Gamma_cut` below the measured lower tearing work rejects **the specified interface-only ledger**. Without native bounds: UNKNOWN. The same gamma/cut requires observed mode-specific area or other work; at the old target means, area-ratioI/III1,478–1,484 is required. This is diagnostic and must not be fed back as parameter fitting.
8. The strongest control receives exactly the same full curves/imaging/history and mode/speed freedom. An algebraic energy sum is TIE. Freeze a ≤25% held-out cross-prediction target for both modes and at least one blade radius/one animal. Charge physical acquisition, setup, certification, rebuild and fallback. No larger cost difference has been shown.

## Final ports and negative bounds

Empirical Γ_cut per layer, transition radius and thermal change are null. Dermis150–380 is **synthetic** severance/wedge-work, not a published cutting measurement or CI. R2 thermal cut75/75 and tear5042,80/744,48 are scenarios with unmeasured joint marks. R2-bending57,5µm and INCISION-bridge732µm describe different closures.

A single equally sized analogue interface misses all tearing targets. SASS nominal M3600/5322 and aortic media Long111/164 are requirement diagnostics, not dermal counts. Generous cylsurface with phi0,3, d82nm, borrowed gamma≤0,30 and assumedV/A1,5mm gives G≤6965J/m²; d2,23µm gives≤622J/m². This rejects the frozen sensitivity, not all delamination: skin thickness/DIC-ahead length does not automatically boundV/A. R3 partial cycle0,46MJ/m³ does not measure terminal fracture work.

Best mechanism picture: recruitment/reorientation and long fiber bridging/pullout with terminal interface/fiber rupture; a quantitative skin explanation is missing. Cell viability, residual mechanics after heating, native cut and shared uncertainty are still UNKNOWN. No further mesh/LP/cache variant resolves the missing acquisition.

## Reproduktion

COMMANDS_R4.md describes the frozen source/PREREG ledger, analytical budget, figure and verification. Scripts write r4 results exclusively; for new reproduction create a new lane-local subfolder, copy scripts and their frozen inputs there and run there. acquire_r4.py retrieves only public primary sources but has several preserved failures; already acquired material is reused. NCBI BioC files underliterature_r4 are the successful full texts for Wang/Kumar/Fang. The filename's ARORA was acquisition naming; the actual source is Kumar2009.

Only lane/shared-data writable; no subagents/cloud/queue/graph/product writes. GRAPH_FEEDBACK_R4.json is local material for the coordinator's binding on definition/review, no admission. The research goal is open but this lane round is closed with a measurable prerequisite.
