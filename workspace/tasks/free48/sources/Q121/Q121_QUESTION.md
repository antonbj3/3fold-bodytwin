### Q121 — Colon transit and gas

**Question:** Can a model distinguish gas production, absorption and retention as causes of measured pressure and intestinal volume?

**Resolution target:** Colon folds and wall activity, heterogeneous contents, gas bubbles, dissolution and flow obstructions.

**Physical basis:** Compressible or regime-valid incompressible flow, heat/mass transport and momentum/energy exchange with the wall.

**Mechanism:** Couple colon motility to separate amounts of gas phase and dissolved gases. Production, uptake into blood, transfer between phases and distal gas release need their own fluxes; proximal return transport requires an explicit path.

**Mathematical start:** `dM_c/dt=J_ileal−M_c/τ_c; dN_g/dt=J_prod+J_upstream−J_rect−J_abs−J_diss. M_c [g], J_ileal [g/h], N_g [mmol], gas-J [mmol/h], τ_c [h]. J_upstream is measured gas transfer from earlier GI segments, not swallowed air directly.`

**Detailed reference to develop:** 3D flow in measured moving geometry, with tissue, transport and acoustics for the question at hand; space/time refined separately.

**Laws/parameters to determine:** Tissue/mucus rheology, turbulence description, fluid interface, perfusion and activation.

**Output:** colon transit and fecal flow, h and g/h respectively; gas volume, mL/min; gas partial pressure, amount and rectal gas flow, mmHg and mL/min respectively

**Reuse:** `src/bodytwin/cells/gastrointestinal/gi_motility_slow_waves.py`; `src/bodytwin/cells/gastrointestinal/gut_microbiome.py`

**Next work/data:** The anchors give movement and substrate outputs but no coherent colonic transit time, gas dissolution or rectal gas opening.

**Comparison:** Constant colon transit and gas formation proportional to substrate.

**Test that can reject the proposal:** The candidate does not distinguish retention, gas formation and dissolution on held-out transit changes or misses measured gas/mass outflows.

**Cross-connections to test:** Q049, Q057, Q059, Q120. Mechanistic development proposal without a specific external source in this review.

