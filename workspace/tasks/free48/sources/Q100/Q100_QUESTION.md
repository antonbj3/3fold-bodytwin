### Q100 — Vocal fold mechanics and self-oscillation

**Question:** When do layered tissue, contact and surface fluid improve the prediction of phonation compared with an established reduced vocal fold model?

**Resolution target:** Vocal fold epithelium and lamina propria, fiber orientation, contact surface and glottal airflow at every oscillation.

**Physical basis:** Compressible or regime-valid incompressible flow, heat/mass transport and momentum/energy exchange with the wall.

**Mechanism:** The tissue's elastic and viscous response, air loading and contact can create self-oscillation. A thin surface layer can change contact and damping. Test layer detail against established reduced dynamics in the same pressure and image sequence.

**Mathematical start:** `m ü+c u̇+k u=F_air(t,u,u̇)+F_contact+F_meniscus. m [kg], u [m], c [N·s/m], k [N/m], all F [N]. F_air is calculated from the coupled airflow and F_contact from tissue contact; meniscus force is included only where such fluid geometry exists. This is a reduced modal building block; the reference solves moving airflow and layered tissue together and checks energy exchange.`

**Detailed reference to develop:** 3D flow in measured moving geometry, with tissue, transport and acoustics for the question at hand; space/time refined separately.

**Laws/parameters to determine:** Tissue/mucus rheology, turbulence description, fluid interface, perfusion and activation.

**Output:** prephonatory pressure [Pa]; mucosal wave speed [m/s]; contact fraction [%]; fundamental frequency [Hz]

**Reuse:** No direct code anchor mapped in this scoped review.

**Next work/data:** Layer geometry and material laws, high-speed video, air pressure/flow and sound, as well as contact and fluid-film data where needed.

**Comparison:** Calibrated two-mass vocal fold model with the same subglottal pressure and measurement budget.

**Test that can refute the proposal:** Layer and contact detail does not improve held oscillation frequency, glottal area and threshold pressure compared with the two-mass model.

**Cross-couplings to test:** Q049, Q050, Q055, Q101. Source entry points: S02.

