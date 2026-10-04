### Q156 — Perivascular, CSF and interstitial exchange

**Question:** Which measured transport contributions between perivascular spaces, CSF and brain interstitium are required for a given solute’s time course?

**Resolution target:** Micrometres in gaps to organ region; seconds to days.

**Physical basis:** Conservation of substance in porous spaces with explicit boundary flow to CSF and blood.

**Mechanism:** Compare diffusion, possible advective flow and binding as alternative mechanisms; a general bulk-flow mechanism is not assumed. The existing Peclet fragment constrains a small substance in parenchyma and provides no general proof of CSF/glymphatic transport.

**Mathematical start:** `∂(ε c_s)/∂t=∇·(ε D_eff,s∇c_s−q c_s)+R_s. ε [-] is fluid fraction, c_s [mol/m³ fluid], D_eff [m²/s], q [m/s] is surface-averaged volume flow and R_s [mol/(m³ total volume·s)]. Domains and boundary flows are specified separately. Pure diffusion q=0 is a comparison, with the same ε and boundary conditions.`

**Detailed reference to develop:** 3D moving perivascular/CSF/ISF domains with time-resolved tracer data and explicit uncertainty in the flow hypothesis.

**Laws/parameters to determine:** Effective diffusion/tortuosity, porosity, possible directed flow, sleep/pressure dependence and binding.

**Outputs:** Regional concentration [mol/m³]; Exchange flow [mol/s]; Effective spreading time [s]

**Reuse:** `src/bodytwin/cells/nervous/glymphatic_peclet_bound.py`

**Next work/data:** Define measurable boundary conditions and alternative transport hypotheses before choosing a mechanism.

**Comparison:** Diffusion in a fixed homogeneous interstitial volume.

**Test that can refute the proposal:** An added advective path is not identified by held-out time and spatial data or violates volume balance.

**Cross-connections to test:** Q153, Q155, Q168. Mechanistic development proposal without a specific external source in this review.

