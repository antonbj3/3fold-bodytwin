### Q090 — Tooth, ligament and jaw joint under chewing load

**Question:** How does the compliance of the tooth, periodontal ligament and jaw joint change the measured load distribution?

**Resolution target:** Root, ligament fibres, bone microstructure and jaw joint layers; local fluid–fibre load is coupled to the movement of the whole jaw.

**Physical basis:** Mass, momentum and energy; contact work and thermodynamically consistent damage/viscoelasticity.

**Mechanism:** Local compliances can be used as a reduced load path; the actual force paths of the jaw joint and tooth roots are a network. Series springs are a limited control case, not the topology of the whole jaw.

**Mathematical starting point:** `C_eff=C_en+C_PDL+C_TMJ; δ=F·C_eff. C [m/N], δ [m], F [N]. Applies only to a defined one-dimensional series load path with small deformations; general jaw joint mechanics require a coupled stiffness matrix.`

**Detailed reference to develop:** Convergence-tested 3D contact/tissue solution with measured microstructure and coupled fluid where needed.

**Laws/parameters to determine:** Anisotropy, fracture/friction laws, active stress, fluid permeability and the individual's microstructure.

**Outputs:** contact load F [N]; tooth displacement δ [µm]; PDL pressure [kPa]; jaw joint angle [°]

**Reuse:** No direct code anchor mapped in this limited review.

**Next work/data:** Individual layer thicknesses, ligament fibre direction, bone and muscle coupling and simultaneous force and motion measurement are missing from the evidence.

**Comparison:** Rigid tooth with a single PDL spring and imposed jaw joint motion.

**Test that can reject the proposal:** If systematic compliance changes do not shift measured load and motion outcomes within the measurement uncertainty, the separate series model is not needed.

**Cross-links to test:** Q050, Q053, Q054. Mechanistic development proposal without a specific external source in this review.

