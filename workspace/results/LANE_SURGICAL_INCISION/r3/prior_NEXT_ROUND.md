# Next supporting operation after round 2

PENDING_CONSTRUCTION / PENDING_INDEPENDENT_REVIEW. Parents GRAPH_INVERSE, MULTIPHYSICS, BT-CTX-SURG-INCISION and BINDINGS R1; no graphs/sources have changed. Empirical gate FAIL; best same-information control TIE; skin Γ0 and joint cross-prediction UNKNOWN. New steering always takes precedence over this handoff.

Capability → cutting/tearing force and layer gap from a shared measured bond/fiber state.
Conflict → a radius series is easier to measure than separate irreversible work terms, but the intercept contains both.
Obstacle → Γ0/B0 has an exact null direction. 150J/m² is the synthetic upper ceiling of short-fibril pullout; nominal12, not chemical cleavage. 15kJ/m² misses both skin tear modes.
Changed operation → same-sample independent work on fiber bridges until blade arrival and contact/prestress energy; then an unfitted held-out blade prediction.

## Start directly

1. Read CHECKPOINT_R2_FINAL.json, SOURCE_INTERPRETATION_CORRECTIONS_R2.md, DECOMPOSITION_R2.json and r2/acquisition_results.json. Preserve all R1/R2 files. Do not continue Green/residual certification as the main track.
2. Primary acquisition: published **skin**, same sample/cohort and orientation, cut/tear work or complete force curves plus thickness, radius/angle, local opening and traction or another independent process-work observation. Measure subsequent reinsertion for contact but do not assume it subtracts constant pre-cut damage or restores prestress. Data from muscle/gel, scissor/needle/scalpel and other species are kept separate. The Smith2013 PDF had access errors; no numerical data were acquired. The Barnett file with .pdf suffix is an HTML error response, not evidence.
3. Freeze new world ordinates before model/fit. Report the exact missing observation if matched open data are not found. Sources are in SOURCE_TARGETS_R2.json; Pissarenko Table1 is already verified and is not a new held-out test.
4. Execute the work port's measurability test: `Gamma0=F_cut/h + G_release - G_contact - W_before`. R2 synthetic precision: cutoff±10J/m², contact coefficient±5, bridge-work±5 gaveΓ134.875..165.125 atR1µm. With bridge-work±20 it gave119.875..180.125 and FAIL. ForΓ150 with20% budget,14.875J/m² remains for ALL release/closure/geometry errors. Nominal T6MPa requires opening uncertainty≤0.83µm for5J/m² if traction were exact; real precision UNKNOWN. Use deterministic intervals or measured joint covariance, no independent marginals as an invented posterior.
5. Mechanical consumer: `layered_skin_r2.py` has Layers(nx,ny,angle,coupling), energy(U,hessian=True), solve([a_epi,a_dermis,a_subcutis],stretch=1.21). It is full sparse Newton in synthetic compressible2D layers with±theta fibers and reversible slip. Prescribed incision and additiveΓ mean post-cut gap alone does not identifyΓ. No warm capability free of global state. Existing final results r2/layers_v2_n96 and r2/layer_extension_results.json may be reused; importing the module writes nothing.
6. Prioritize actual blade arrival `delta_b(R,wedge,friction,prestretch)` and recruitment/damage history over a larger mesh. R2's explicit `delta_b=2R` gives half the pullout work atR732µm,73 fiber diameters: fiber diameter alone does not derive the transition. Replace this closure with an independently measured/local contact law and test held-out sharpness without a new Γ fit. If it does not exist, parameterized intervals with UNKNOWN world gate.
7. The matched control gets the same microscopy/work/contact data and a mode-separated cohesive law with the full work curve. Measure prediction error and the number of **independent** empirical acquisitions, not just training fit. Algebraic identity is TIE. Real bond capability requires also surpassing this strong control contract without extra free fitting.

## Numerical raw data and limitations

Derma16→20mm over intact subcutis releases0.21639N in all layers; Γ150*1.5mm=0.225N leaves0.00861N before unknown process/contact. The 64→96 release difference1.040J/m² is a sensitivity, no rigorous continuous or empirical error bound. Epi-only96→128 gap change0.963%,dermis8mm0.962%; stop routine mesh increases. Final verification is in VERIFICATION_R2.json. Initial raw files, oblique-family/length convention and mis-signed energy report are preserved; corrected_energy_check is bookkeeping, no biological validation.

Reproduce new data in a NEW output directory: `OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2 NUMEXPR_NUM_THREADS=2 python results/LANE_SURGICAL_INCISION/layered_skin_r2.py --nx 64 --out results/LANE_SURGICAL_INCISION/r3_layers_probe`. Do not run acquisition/finalize over R2; they write exclusive R2 files. For a new acquisition round, copy the module/PREREG to a new suffix and change the destination before running. COMMANDS_R2.md has the original commands.

Biological cell, vessel, bleeding and hemostasis ports are UNKNOWN. No clinical admission. Only this lane/shared-data is writable; no subagents/cloud/queues/source-graph changes. Graph feedback GRAPH_FEEDBACK_R2.json is local material for the coordinator's regular review binding.

## Final addition

Full-depth8mm0°/45°64→96 missed the2% gate; preregistered128 cases gave at most1.421% change. Final verification/cost is in VERIFICATION_R2_FINAL.json/COST_ACCOUNTING_R2_FINAL.json. BINDINGS R2 has been read: its bending-rupture radius57.5µm and our accumulated bridge-work radius732µm are different mechanisms and both conditional. Read PEER_CONTEXT_R2.md and BINDINGS CHECKPOINT_R2_FINAL/PORTS_R2_FINAL before shared coupling; RESPONSE still had no NEXT_ROUND or numerical night-round at the final check.
