### Q160 — Biofilm EPS mechanics and transport

**Question:** How does the extracellular polymer matrix change local diffusion, flow resistance and detachment of microbial communities?

**Resolution target:** Nanometre pores to micrometre colonies; seconds to days.

**Physical basis:** Water, substance and cell balance and force balance in a porous deformable matrix.

**Mechanism:** Follow cells, EPS, water and dissolved substances as separate volume fractions; mechanical deformation and diffusion share moving geometry.

**Mathematical starting point:** `q=−(κ/μ)∇p, with surface-averaged fluid flux q [m/s], κ [m²], μ [Pa·s], ∇p [Pa/m]. J_s=−ε D_eff,s∇c_s+q c_s [mol/(m²·s)], ε [-], c_s [mol/m³ fluid]. Storage is ε c_s per total volume. This Darcy control case assumes a stationary matrix and flow with negligible inertia; deformation requires relative flow and matrix balance.`

**Detailed reference to develop:** 3D poroelastic cell/EPS biofilm with measured flow, deformation and tracer penetration.

**Laws/parameters to determine:** EPS elasticity, permeability, adhesion, shear-dependent detachment and binding.

**Outputs:** Biofilm thickness [m]; Effective diffusion [m²/s]; Detachment flux [celler/s]

**Reuse:** No direct code anchor mapped in this limited review.

**Next work/data:** Define EPS amount and distinguish biofilm from free luminal bacteria in Q159.

**Comparison:** Fixed biofilm volume with constant diffusivity.

**Test that can reject the proposal:** Coupled mechanics do not improve held-out penetration and detachment curves.

**Cross-links to test:** Q149, Q159, Q161. Mechanistic development proposal without a specific external source in this review.

