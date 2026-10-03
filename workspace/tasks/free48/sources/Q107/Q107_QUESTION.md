### Q107 — Stomach mechanics and emptying

**Question:** When does fundus, antrum and pylorus geometry improve the prediction of gastric emptying for liquid and solid food?

**Resolution target:** Stomach wall layers and electrical waves, food phases and the pylorus opening/closing process.

**Physical basis:** Mass, momentum and energy; contact work and thermodynamically consistent damage/viscoelasticity.

**Mechanism:** The wall's electrical activation and contractile waves shape pressure, mixing, particle retention and pylorus opening. A reduced hydraulic model gives a control case for coupled wall–fluid–particle mechanics.

**Mathematical starting point:** `dV_g/dt=−Q_pyl+Q_in; Q_pyl=C_D·A_p·√(2·max(ΔP_pyl,0)/ρ) for an orifice-like inertial regime. V_g [m³], Q [m³/s], A_p [m²], ΔP [Pa], ρ [kg/m³], C_D [-]. For viscous/particle-bearing emptying, a separately calibrated hydraulic law is used.`

**Detailed reference to develop:** Convergence-tested 3D contact/tissue solution with measured microstructure and coupled fluid where needed.

**Laws/parameters to determine:** Anisotropy, fracture/friction laws, active stress, fluid permeability and the individual's microstructure.

**Output:** gastric volume, mL; pylorus flow and pressure, mL/s and mmHg respectively; solute and residence-time profiles, min

**Reuse:** `src/bodytwin/cells/gastrointestinal/gi_motility_slow_waves.py`; `src/bodytwin/cells/endocrine/glucose_meal_dallaman2007.py`

**Next work/data:** The anchors document slow waves and k_empt respectively but not the proposed moving wall or solid–liquid separation.

**Comparison:** Calibrated emptying model with separate solid/liquid fractions and the same meal, volume and time data.

**Test that can reject the proposal:** Rejected if the model cannot match held-out volume and retention curves better than the same measurements with a scalar k_empt.

**Cross-couplings to test:** Q012, Q049, Q057, Q059. Source entry points: S01.

