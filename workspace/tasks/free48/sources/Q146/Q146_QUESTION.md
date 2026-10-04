### Q146 — Nasal route with swallowed and systemic share

**Question:** Can local nasal geometry improve the prediction of both local and systemic exposure when the swallowed fraction is followed separately?

**Resolution target:** Regional nasal microgeometry resulting from both mucosal uptake, swallowing and organ/cell distribution.

**Physical basis:** Common control volumes and identities for mass, element, charge, force/work and energy across each interface.

**Mechanism:** Deposition distributes substance on mucosa; dissolution, epithelial transport and mucus clearance compete. Swallowed fraction goes to the GI model and its first pass.

**Mathematical start:** `dN_nose/dt=I_dep−J_epi−F_swallow−F_external−R_met. The same F_swallow [mol/s] enters the GI and J_epi enters regional tissue/blood after any local metabolism.`

**Detailed reference to develop:** Co-run of the resolved sub-models with common time, interface and substance identity; each part has its own measuring anchor.

**Laws/parameters to determine:** The material and reaction laws of the submodels, time delays and uncertain boundary/initial conditions.

**Output:** Regional residual quantity [mol], swallowed fraction [1], free local concentration [mol/m³], plasma AUC [mol·s/m³].

**Reuse:** `src/bodytwin/cells/respiratory/mucociliary_clearance.py`; `src/bodytwin/cells/gastrointestinal/hepatic_clearance.py`

**Next work/data:** Formulation-specific dissolution and deposition, regional clearance and same substance PK. Direct pathway to brain is a separate mechanism test.

**Comparison:** Route split compartment model with common absorbed fraction but same measurement and calibration budget.

**Test that can fail the proposal:** Geometry gain disappears on held-out anatomy/formulation or plasmafit hides wrong local uptake.

**Cross connections to try:** Q009, Q013, Q015, Q132, Q134. Source inputs: S03, S05.

