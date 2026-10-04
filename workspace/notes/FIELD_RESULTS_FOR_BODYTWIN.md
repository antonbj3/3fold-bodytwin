# Field engine results as a basis for BodyTwin: transfer and hypotheses

Agent O3, 2026-09-22. Read-only in `romi_collab` and source projects. No lanes, queues or computations have been run. I read the numbers below in the reports. I have not reproduced them.

**Status labels.** *Confirmed* means that an independent audit reproduced the result within the stated condition. *Failed* means that the result was rejected by an audit or a locked gate. *Hypothesis* is my own proposal and has not been tested. Where the audit and producer differ, the audit’s scope applies.

**Sources.** Research root `R=../3fold-motion-engine/_private/romi_collab`, files `R/build/<LANE>/RESULTS.md` and `remaining_obligations.json` for A185, A194, A201 and A204. Also `~/HANDOVER_FIELD_2026-09-22.md`, `~/research/field_handover_20260922/LANE_INVENTORY.md`, `R/lanes/{INNOVATION_OBJECTIVES_2026-09-22.json,SEED_DECOMPOSITION.md}`, `~/research/inference_training_20260921/SHARED_GEOMETRY_MATH.md`, `~/research/COLLABORATOR_BODYTWIN_CONTEXT_2026-09-22.md` and `STARTUP_MESSAGE.md` (the section "Our geometric relationships"). For A196 (U280) and A198 (U282), I have only read the summary. They are included as side notes.

---

## 1. Five cross-cutting lessons

Every BodyTwin experiment based on the field results needs to take these into account.

1. **First check that the coarse model actually fails.** In U262/A178, the coarse model already met the error requirement of 0.015 mm. Then there is no gain to obtain however correct the refinement is. Every resolution or adaptivity hypothesis therefore needs to first show that the coarse model misses the task’s requirement.
2. **Simple local rules often win.** In U264/A180, the adjoint-based selector chose the same cell as the rule "refine the output cell" in 5/5 cases, but cost 30–50 ms more. In U290/A206, an anchor at a declared location reached the same outcome as the exact task planner without planning cost, whereas the planner cost 113 982 evaluations. In U265/A181, warm-start PCG beat the Woodbury adapter. Such simple rules are the strong baselines.
3. **Sum checks do not see partition errors.** In A194, `V1+V2=|cell|` still held when the distribution between materials was wrong by +3200 mm³. In A204, a hidden tunnel did not appear in the component graph, although β1 changed from 0 to 1. The check metric must start from the recipient: volume per material, orientation and handles.
4. **Best approximation and the Galerkin solution are different things.** This was shown in U285/A201 and A185. A space’s ability to represent the field does not mean that the energy-minimizing solution meets the requirement. Finite-order quadrature does not give a lower bound either: Q2 was underestimated by 24 % in A185.
5. **Self-consistency is not physical validation.** In U270, U286 and U290, the "truth" was the same CSG model as the engine, with a deviation of 1e-18 to 1e-19 m. The accuracy numbers thus say that the model is consistent with itself. In A202, the cost gain also depended on a currency that was not charged. It is no measured gain.

---

## 2. Results by track

For each result, what was confirmed, what failed, which BodyTwin recipient it may have and what a transfer requires are stated. The requirement is given as a mathematical mapping and the conditions that must hold.

### 2.1 Local elasticity and adaptive DOF

**U262/A178: changed degrees of freedom with an exact global Schur update**
- *Confirmed:* Q1 plus zero-trace bubbles in interior elements. The DOF count changes 192 → 195 → 216 → 195. Global Schur agrees with the direct solution within ≤6e-16. 27.7 % of `A⁻¹B` lies outside B’s support, so the nonlocal coupling must be retained.
- *Failed:* No total cost gain at the locked requirement of 0.015 mm, because the coarse model already met it (error 0.0085 mm). The result applies to a static manufactured test case, not general h-refinement.
- *Recipient:* A locally refined region in a frozen bone geometry, for example an implant seat or a ligament attachment, where the rest of the bone is solved coarsely.
- *Transfer:* The block system `[A B; Bᵀ C][u;a]=[f;g]` gives `S=C−BᵀA⁻¹B`. It requires frozen geometry and frozen material in A, that the enrichment has zero trace toward neighboring elements and that `A⁻¹B` is retained globally. Before this is meaningful, one must show that the coarse model misses the requirement for the biomechanical quantity.

**U264/A180: DOF enrichment selector driven by the output’s adjoint**
- *Confirmed:* The adjoint score predicts the change from a single bubble exactly (<1e-14 mm). The fallback via a residual guard works.
- *Failed:* The selector chose the same cell as the rule "output cell" in 5/5 cases and met 0.005 mm in 1/5 cases. It cost 32–57 ms more. Uniform 10³ passed 5/5.
- *Recipient:* Choosing where to refine for a point quantity, for example strain at an attachment.
- *Transfer:* The score `|(l−BᵀA⁻¹q)ᵀS⁻¹(g−Bᵀu₀)|` measures the change from a single enrichment, not the remaining error. It cannot be used as a stopping criterion.

**U265/A181: Woodbury reuse for a local material change**
- *Confirmed:* The algebra for `K'=K+W diag(λ)Wᵀ` holds (≤8e-11). Against a new factorization, a local change gave 1.2–7.5 times lower cost.
- *Failed:* In topopt, 97.4 % of the elements change per iteration, so the cache is invalidated at every step (ratio 1.0043). Warm-start PCG beats the adapter in all bounded rows, for example patch8 at 0.32 s versus 0.60 s. The rank cap 256 became binding.
- *Recipient:* Repeated small material changes on a frozen base. One example is a few elements in a lesion or a local remodeling zone.
- *Transfer:* Requires `rank(ΔK_free)` to be small. An entire tissue region, such as PDL or trabecular bone, is not small in rank. A parameter family `K(s)=K0+s·K1` with large support for K1 gives no low-rank gain. Warm-start PCG is a mandatory comparison.

**U269/A185: interior bubbles compared with shared edge and face modes**
- *Confirmed:* Interior bubbles do not meet the volume requirement 0.06 (0/5), and the space is too small (projection error 0.039–0.106). A185 rejected U269’s stronger conclusion that h-refinement is necessary: shared edge modes lower the error 0.0975 → 0.021, and full Q2 passes 5/5 at the same DOF as uniform 10³.
- *Failed:* The claim that h-refinement is uniquely necessary. Uniform 10³ also fails at corners (0.0672) and at edge-edge (0.0620).
- *Recipient:* Choosing a basis space when local accuracy in stress or strain is wanted.
- *Transfer:* The enrichment must be compatible higher order, i.e. shared edge and face modes. Bubbles are insufficient. The cost is not measured algorithmically because the element assembly is written in Python.

**U285/A201: output-centered radius-2 patch with Q2 trace modes**
- *Confirmed:* The patch with 25 cells and 252 modes meets both requirements on the locked cases (5/5 + 1/1) with exact conformity. At the same accuracy, the patch is cheaper than the uniform refinement required, for example 0.225 s versus 1.267 s at n=16.
- *Failed:* The frozen rule with a ring residual expanded in 1/6 cases. Generalization was rejected: the corner case `adv_corner18` only passed at radius 3, because the patch was truncated against the boundary. `fresh_corneredge12` fails for all trace modes on 5³ and is a real limit requiring h-refinement. No speed gain arose among the cases that met the requirements (0/13).
- *Recipient:* Repeated questions about stress at a point on frozen geometry, for example for many load cases during a gait cycle.
- *Transfer:* The patch is centered on the output point, with mirrored handling near the boundary (open hypothesis in A201). Three conditions apply. Features narrower than coarse h require h-refinement. The solution must be smooth within the patch, which is interrupted by material boundaries between cortical and trabecular bone (not tested). More right-hand sides with uniform factorization is a baseline that has not been tested.

### 2.2 Thin materials, topology and conservation

**U266/A182: local topology guard when editing an occupancy grid**
- *Confirmed:* No false accepts among 573 409 edits. This is empirical and not a proof. The refusing guard gives a factor of 0.198–0.212 against a global check.
- *Failed:* With full coverage and fallback, the method loses (1.17–1.70 times global). The state diverges in 2/17 cases. 33 false rejections are conservative.
- *Recipient:* Segmentation editing or morphing that must preserve topology. Examples are the marrow cavity (cavity), foramina (tunnel) and the joint space, where two cartilage surfaces must not grow together.
- *Transfer:* The condition is that β0, β1 and β2 are unchanged in a window of binary occupancy with connectivity (26,6) and an explicitly declared outer boundary. The result only applies to a discrete grid, not a continuous surface, and is no certificate. Speed requires accepting refusal.

**U272/A188: typed handoff contract between the field and BodyTwin**
- *Confirmed:* Round-trip, orientation (error 0 of 6984), bit-exact conversion mm → m and 15/15 header mutations are rejected.
- *Failed:* The provenance flag is not connected. The registry is not enforced (M=9 passes). UNKNOWN is never output: queries against a lost subregion gave region 1 in 968/968 cases. Region 3, a thin wall, has 0 cells, the source mesh is non-manifold and `transfer_ok=true` despite `valid=false`.
- *Recipient:* The contract for a *geometry manager* with stable region, body and material identities, units and coordinate frames. This is closest to the collaborator's expressed needs.
- *Transfer:* Four things must exist: explicit UNKNOWN when a subregion is lost, an enforced registry, provenance from the right key and a conservation check per region. Transfer is not the same as conservation.

**U278/A194: conservative subcell moments `(V, S)` per material**
- *Confirmed:* The arithmetic is exact (≤1e-15) and conservative. Affine covariance holds: `V'=V`, `S'=RS+tV`. A 1 mm thin layer with pitch 3.125 mm is preserved exactly, while majority labeling gives −2.7 % to +144 %. Queries against stored values are about 190 times faster than new clipping.
- *Failed:* The signed summation is unsound when slabs cross (ADV_CROSS): +3200 mm³ in material 2 and 482 cells throwing errors. There is no memory gain. The method is known as moment-of-fluid (MoF).
- *Recipient:* The segment’s mass and center of mass in a multimaterial voxel. Material fraction per cell during FEM assignment.
- *Transfer:* The partition must be disjoint (see U288). `(V,S)` gives no topology, layer order or stiffness: two configurations with the same `(V,S)` had 2 and 1 components respectively, and the bending stiffness differed by a factor of 1.95. Second moments are not implemented.

**U288/A204: disjoint moments plus an interface graph**
- *Confirmed:* The disjoint partition repairs ADV_CROSS and gives 252800/259200. The moments distinguish components. The consumer `lastfalt_v1_fem` gives the same E as the series formula, but this is an algebraic identity at ν=0 and no validation.
- *Failed:* The graph nodes are not derived from the interface. They coincide with a separate union-find over the tetrahedra. The hidden tunnel is not detected. Non-manifold contact at edges and corners is accepted. The contraction gate is tautological. UNKNOWN and the registry are missing. `lastfalt_v1_fem` has no refusal path.
- *Recipient:* Thin barriers where a hole must not be treated as sealed. Examples are cortical shell, cartilage, PDL and septa.
- *Transfer:* A β1/Euler witness quantity for the material-restricted complex (A204’s next hypothesis) and a refusal path in the FEM consumer are required.

### 2.3 Geometrirepresentation

**U267/A183: adaptive interval enclosure with UNKNOWN**
- *Confirmed:* The definitive labels are sound in real arithmetic, with 0 misclassifications on 4000 unseen points. The uniform baseline errs on 187–240 points. The encoding fix works.
- *Failed:* The floating-point enclosure is up to 4.24e-15 mm too narrow and needs an outward margin. The exact expression is 12–336 times smaller than the adaptive source, and the adaptive source plus exact fallback is strictly dominated by the expression alone. Coverage without fallback is 72–93 %.
- *Recipient:* Queries about minimum thickness or distance where abstention is allowed, for example cortical thickness ≥ t_min or the gap between an implant and a bone canal.
- *Transfer:* For analytical CAD, such as implants and dental manufacturing, use the exact expression. Anatomical geometry is a mesh or SDF and has no expression. Rigorous interval bounds for a mesh SDF are not demonstrated, so they are a gap.

### 2.4 OED and geometry measurement

**U270/A186: rank-first plan with a coverage requirement per feature**
- *Confirmed:* Negative result. Rank 6/6 is reached at 9 truths, but the corrected full cost 34296 is higher than adaptive OED (29958).
- *Failed:* The gate against false accepts is tautological (E2). The median is calculated incorrectly (E3). The planning cost was overcharged by 28704 (E1).
- *Recipient:* Choosing where to measure on the surface (landmarks or points) to estimate shape parameters.
- *Transfer:* The condition is a real contract for observation, access and noise. No physical validation has been done.

**U286/A202: reusable analytical Jacobian dictionary for OED**
- *Confirmed, partly:* Rank 6 within the entire trust box and an error of 2.15 µm were reproduced on this test case. This is self-consistency.
- *Failed:* The currency of 0 evaluations: with counted differences, the planning cost becomes 19008–32832 and the gate fails. The validity radius is no certificate, because 56/1056 probes change branch in the box’s corners. Fallback and guard were not connected.
- *Recipient:* Same as U270. The result also shows that boundary changes (topology or contact) invalidate local Jacobians.
- *Transfer:* The Jacobian’s rows must be paid for or derived from the model’s own expression tree. Validity must be checked in every corner of the parameter box.

### 2.5 Geometric alternatives and shared evidence

**U290/A206: OED that distinguishes discrete geometric worlds**
- *Confirmed:* Negative result. An anchor at a declared location suffices. Local Fisher within a branch is blind to the branch, with det = 0.
- *Failed:* The contract with 15 probes was implemented with 13. w0 and δ are only identified as a sum, rank 1 (E2). Expected loss was unnormalized (E4). Undeclared tracks gave confidently wrong decisions in all arms.
- *Recipient:* Discrete anatomical variants, for example presence or absence of an osteophyte, bifid structure or accessory bone. Also calibration parameters that only appear as a sum.
- *Transfer:* A state combining branch and nuisance. Refusal upon deviation for undeclared features is an open hypothesis.

**U291/A207: fusion accounting for factors and provenance**
- *Confirmed:* A fusion retaining the joint distribution over (m, θ) reproduces exact Bayes, with conditional error 0.026 versus 0.467 for a naive product. The stronger baseline in A207, where θ is handled correctly, gives 0.19.
- *Failed:* The prior in the code (7/10) differed from the contract (3/4). The provenance guard was a stub. Byte-level dedupe is defeated by re-encoding. An underlying layer crossing a cell boundary breaks the whole-cell output (conditional error 0.478).
- *Recipient:* Fusion of CT, MR, surface scanning and atlas prior. A common atlas registration or segmentation model is the shared factor.
- *Transfer:* The shared calibration must be included in the state. The provenance key should be the identity of the raw acquisition plus algorithm version, not a hash of bytes.

**U292/A208: compression of worlds with exact equivalence**
- *Confirmed:* The grouping is exact and commutes.
- *Failed:* It gives no gain in evaluations or storage compared with a memoized baseline. Reused factor IDs lose evidence. The same mean and covariance can give different futures (F3), so there is no moment closure. Underlying layers require 925 repetitions in a binary interface, while a continuous residual guard refuses (F4).
- *Recipient:* Representation of uncertainty in shape, for example a Gaussian SSM versus an ensemble.
- *Transfer:* Branches in shape uncertainty require retaining the members. A Gaussian SSM cannot carry branches.

### 2.6 Side notes (only the summary read)

- **A196 (U280):** Material sensitivity (μ, e) through a contact event agrees within 4.7e-9, in a planar rigid model with few events. The frozen event time gives errors up to 3.9. A guard against a kink in the final time is missing. Recipients are muscle paths that start or stop wrapping around a surface and joint contact that is engaged. The scope in A191 applies.
- **A198 (U282):** Quadratic boundary closure wins on wall widths that do not fall on the grid (6/6), while halfway wins on the locked cases. Recipients may be flow in vessels and airways. The track is outside what is prioritized here.

### 2.7 Seed: `SHARED_GEOMETRY_MATH.md`

- *What transfers is mathematics:* the local Fisher tensor `G(z)=Σ_q w_q J_qᵀ diag(1/P_q) J_q` and its kernel. The kernel consists of directions that do not change any measured prediction at first order. This is local observability, not global (§5). Pair distances in prediction space are zero for worlds that the query family does not distinguish (§3).
- *What does not transfer:* A prediction distance, Hellinger or Fisher–Rao, is not an anatomical distance in mm. An anatomical application requires a typed observation mapping `θ_form → y` with units and noise, known correspondences (landmarks and region identity) and comparison with a baseline in mm, for example Procrustes/PCA. Pair optimization without anchors allows labels to be permuted (§4). This must be closed with typed outcomes.

---

## 3. Negative results that every transfer must respect

| # | Negative result | Consequence for BodyTwin |
|---|---|---|
| N1 | U262/A178: the coarse model already met 0.015 mm | First show that the coarse geometry or model misses the task’s requirement. |
| N2 | U264/A180: adjoint scanning = the "output cell" rule | Local refinement at the output point is a mandatory baseline. |
| N3 | U265/A181: warm-start PCG beats Woodbury, and broad changes invalidate the cache | Speed through reuse applies only to low-rank changes on a frozen base. |
| N4 | U269: bubbles are insufficient. A201: the patch fails at corners and for narrow features | Compatible higher order, boundary-aware patch, h-refinement below coarse h. |
| N5 | A194: the sum check does not see a +3200 mm³ partition error | Check volume per material and per region, never only total volume. |
| N6 | A204: tunnels are invisible, non-manifold is accepted, the contraction gate is tautological | The topology witness must include β1. Manifoldness should be checked via edge and corner incidence. |
| N7 | A188: UNKNOWN is never output, and a thin subregion disappears silently (0 cells) | A lost subregion should give UNKNOWN or refusal, never a neighbor label. |
| N8 | A183: the floating-point enclosure is 4.24e-15 mm too narrow, and the exact expression dominates | Outward margin. Use the exact expression when one exists (CAD). |
| N9 | A186/A202: tautological gates and uncharged cost currency | Every cost is counted in paid evaluations. Shape gates must not be satisfied by definition. |
| N10 | A206: the anchor suffices, local Fisher is blind to branch, rank 1 for a sum parameter | Discrete variants require a branch state. Sum parameters should be reported as unidentifiable. |
| N11 | A207: a naive product gives confidently wrong decisions, and byte dedupe is defeated | Shared calibration or registration must be in the joint state. |
| N12 | A208: the same mean and covariance give different futures | Gaussian shape uncertainty is insufficient at branches. |
| N13 | A191/A196: frozen events and kinks in the final time give errors | Sensitivity through events, such as wrapping and contact onset, requires an event term and a guard. |
| N14 | U270/U286/U290: the "truth" is the engine’s own model | Neither µm accuracy nor planning says anything about anatomical validity. |

---

## 4. Hypoteser

All hypotheses below are *hypotheses*. None has been tested. Each hypothesis is stated with idea, mechanism, equation, operation, representation, assumptions, expected outcome, strong baseline, counter-test, cost and data. The mapping to the tracks is stated as **J** (the collaborator: shape parameter → biomechanical quantity) and **B** (broad BodyTwin: material, tissue and uncertainty).

### H1. Task-driven geometric error bound: what resolution does this particular computation require? (J+B, priority 1)

- **Idea:** The geometric accuracy required is determined by the shape derivative of the biomechanical quantity, not by a global surface metric or volume error. Different quantities on the same bone should require tolerances differing by orders of magnitude.
- **Mechanism:** Hadamard’s shape derivative. A normal displacement ε(x) of the surface gives the first-order change `δQ = ∫_∂Ω g_Q(x) ε(x) dA`.
- **Equation:**
  - Mass: `g_m = ρ`.
  - Center of mass: `g_c = ρ(x−c)/m`.
  - Inertia about the center of mass: `g_I = ρ(|x−c|² 𝟙 − (x−c)(x−c)ᵀ)`.
  - Moment arm `r(θ)=∂ℓ/∂θ` for a straight path or via-point path: g_r is concentrated at the locations of attachments and via-points. The surface’s position affects it only through how the attachment is anchored to the surface, or through the surface normal during wrapping.
  - Budget: with the requirement `|δQ| ≤ τ_Q`, the allowed local tolerance is a weighted bound. The worst case is `sup|ε| ≤ τ_Q / ∫|g_Q| dA`, and the local distribution is given by `ε(x) ∝ 1/|g_Q(x)|` under a cost budget.
- **Operation:**
  1. Take a real bone mesh.
  2. Calculate Q exactly from the mesh: polyhedral mass and inertia via the divergence theorem, moment arm from a declared path.
  3. Voxelize or convert to SDF at pitch h ∈ {4, 2, 1, 0.5} mm.
  4. Measure δQ(h) and compare with the prediction `∫ g_Q ε_h dA`, where ε_h is the measured distance between the surfaces.
  5. Then build a refinement weighted by g_Q and compare it against uniform refinement at the same cell count.
- **Representation:** Surface mesh (truth) → field SDF, the existing `compose_mesh_to_field.py` → voxel or octree. Landmarks are carried with stable identity.
- **Assumptions:** Homogeneous density per region. The moment arm’s path is declared (straight line or via-point). First order applies for small ε, but not near wrapping transitions (see H5).
- **Expected outcome:**
  - Mass and center of mass already pass at coarse h, because the errors in ε cancel (mean near 0).
  - Inertia requires finer h at surfaces far from the axis.
  - The moment arm is barely affected by voxel resolution if the attachment is a landmark point, but strongly if the attachment is defined as the nearest point on the surface.
  - If this holds, the conclusion is that the resolution decision per quantity can be made before running.
- **Strong baseline:**
  1. Uniform refinement with Richardson extrapolation.
  2. The rule "refine around the attachment or output", as in U264/A180.
  3. The volume error metric used today, 0.0975 % for the kidney mesh.
- **Counter-tests that reject the hypothesis:**
  1. If the coarse pitch 2 mm already gives all Q within τ (N1), the hypothesis fails for these Q.
  2. If `∫ g_Q ε_h dA` does not predict δQ(h) within a factor of 2 across all h, the mechanism fails.
  3. If refinement weighted by g_Q does not beat the local rule at the same cell count, the mechanism is unnecessary (N2).
  4. Mandatory: attachments near bone edges and corners, given A201’s corner results.
- **Cost:** CPU in minutes. No new solver is needed, only the existing mesh → SDF chain and polyhedral integrals.
- **Data:** A public bone mesh with a license. D2/O2 maps what is available locally; for me it is UNKNOWN. τ_Q is taken from published variation between individuals, which must be sourced separately, otherwise it is declared in advance.

### H2. Identifiability: which shape and material parameters can be distinguished, and which measurement decides most? (J+B, priority 2)

- **Idea:** Many combinations of shape, material and calibration appear only as sums or products in available observations. Examples are scale versus density in mass measurement and registration shift versus attachment location. Their kernel can be calculated and should be reported. The next measurement is then chosen to reduce variance in the *task*, not a global logdet.
- **Mechanism:** Local Gauss–Newton or Fisher information (seed §5, A206 E2). The task variance is `Var(Q) ≈ g_Qᵀ G⁺ g_Q` only if `g_Q ⟂ ker G`. Otherwise Q is unidentifiable and must be reported as UNKNOWN. The pseudoinverse must not be used as a guess.
- **Equation:** θ = (b₁..b_k SSM/scaling parameters, ρ_kortikal, ρ_trabecular, E_kortikal, E_trabecular, R, t registration, σ-calibration). The observation model is `y = h(θ) + e`, `e~N(0,Σ)`, `G = JᵀΣ⁻¹J` with `J = ∂h/∂θ`. Candidate measurement j gives `G_j = G + j jᵀ/σ_j²`. The choice is made with the criterion `argmin_j g_Qᵀ G_j⁺ g_Q`, with charged cost.
- **Operation:**
  1. Build a synthetic but anatomically parameterized model: shape modes, two material regions and a registration.
  2. Calculate J with counted finite differences (N9: no uncharged currency).
  3. Eigendecompose G for each set of observations: {surface points}, {surface points + landmarks}, {+ total mass}, {+ CT-HU profile}, {+ a load–displacement measurement}, {+ a natural frequency}.
  4. Report the kernel and task variance for Q = moment arm, inertia and bending stiffness.
  5. Check globally with multistart and mirror branches, because local rank is not global identifiability (seed §5).
- **Representation:** A parametric shape model with correspondence. Material regions via disjoint moments (U288), because a partition error otherwise contaminates J (N5).
- **Assumptions:** Gaussian noise with declared σ; a missing σ must not be replaced with a default value. Small perturbations. Known correspondences.
- **Expected outcome:**
  - With surface observations alone, (E_kortikal, E_trabecular) lie entirely in the kernel.
  - With mass, ρ and a total scale factor are partly confounded until a density measurement (CT-HU) is added.
  - Registration shift and attachment location along the surface are confounded without a landmark.
  - Task-based selection gives at least as low Var(Q) as D-optimal at the same cost.
- **Strong baseline:**
  1. D-optimal (logdet).
  2. An anchor at a declared location: measure directly at the attachment or output (the A206 winner).
  3. Random selection.
  4. "Measure the quantity directly", where possible.
- **Counter-tests:**
  1. If the anchor gives the same Var(Q) as task-based selection, the mechanism fails (N10).
  2. If the predicted kernel directions are not confirmed in a global multistart (are there different optima with the same residual?), the local analysis is insufficient.
  3. If rank changes within the parameter box (branch change, compare A202’s 56/1056), the local G is invalid there.
- **Cost:** CPU. The cost is dominated by J with k+8 parameters times the number of observations. Small.
- **Data:** Synthetic model. Later a real SSM with a license (D2). Noise and σ for each modality must be sourced. If missing, it becomes UNKNOWN and no default values are used.

### H3. Thin regions and material partition through conversion and deformation, checked via the recipient’s stiffness (B, also J via inertia; priority 3)

- **Idea:** For a layer thinner than the pitch, such as cortical shell, articular cartilage or PDL, neither material fraction nor total volume suffices. The recipient is FEM stiffness, and it depends on fraction, the layer’s normal and continuity. A cell representation carrying disjoint fraction and first moment per material, plus a β1 witness, should reproduce the fine reference’s stiffness. Majority labeling and isotropic mixing should not.
- **Mechanism:** Laminate theory for a mixed cell with layers parallel to the normal n:
  - Along the layers (Voigt): `E_∥ = Σ f_m E_m`.
  - Across the layers (Reuss): `E_⊥ = (Σ f_m/E_m)⁻¹`.
  - The normal is estimated from the difference between the materials’ centers of mass within the cell: `n ∝ S_2/V_2 − S_1/V_1` (U278/U288). This gives a transversely isotropic cell stiffness.
  - A shell’s bending stiffness `D ∝ E t³` amplifies thickness errors: `δD/D ≈ 3 δt/t`.
- **Equation, deformation:** Under a morph with local Jacobian F, volume scales by `det F` and the layer’s area by `det F·|F⁻ᵀn|`. The new thickness therefore becomes `t' = t/|F⁻ᵀn| = t·(nᵀC⁻¹n)^{-1/2}`, with `C=FᵀF`. The new normal is `n' ∝ F⁻ᵀn`. For affine segments, this is exact. Moments transform affinely (confirmed in U278, `S'=RS+tV` in the rigid case). The affine case is `V'=|det F|V` and `S'=|det F|(F S + t V)`.
- **Operation:**
  1. Use a real or synthetic shell geometry, for example a cylindrical cortical tube of thickness t ∈ {0.3, 1, 2} mm with a trabecular core, or a tooth root with a PDL layer.
  2. Run fine conforming tetrahedral-mesh FEM as truth.
  3. Convert to a coarse grid, h = 1–3 mm, with four variants: (a) majority label, (b) fraction + isotropic Voigt or Reuss, (c) fraction + normal + laminate, (d) same as (c) after an RBF morph.
  4. Measure bending and torsional stiffness and axial stiffness.
  5. Measure volume *per material and region* (N5) and β0/β1 for the thin material (N6).
- **Representation:** `u288.mixedstate` (disjoint V and S) plus a β1 witness (A204’s next hypothesis) plus UNKNOWN when a layer disappears (N7).
- **Assumptions:**
  - The layer is locally planar across a cell, i.e. the radius of curvature is much larger than h.
  - Linear elasticity with E values declared per tissue. The tissue values are a model closure, not validated.
  - Perfect bonding, because sliding between the interfaces is missing from the contract (U272).
- **Expected outcome:**
  - (a) misses the shell or makes it too thick, with stiffness errors of tens to hundreds of percent (compare the layer test in U278: +144 %).
  - (b) overestimates bending stiffness across the layer.
  - (c) lies within the declared tolerance when h is less than about half the radius of curvature.
  - (d) preserves volume per material within the FEM tolerance during morphing.
  - The β1 witness distinguishes SOLID from TUNNEL.
- **Strong baseline:**
  1. Conforming tetrahedral-mesh FEM at the coarse level, with a layer at least one element thick.
  2. Isotropic Voigt and Reuss bounds.
  3. A refined grid where h < t.
- **Counter-tests:**
  1. If isotropic Reuss or Voigt already gives stiffness within tolerance, the normal is unnecessary.
  2. If (c) fails where the radius of curvature is about equal to h, for example at the femoral neck’s cortex, the hypothesis is restricted to planar regions.
  3. If the normal from S is unstable under a shift of half a pitch, the representation fails.
  4. If β1 is identical for SOLID and TUNNEL, the topology part fails (A204).
- **Cost:** CPU. The fine reference FEM dominates and takes minutes to hours depending on the mesh.
- **Data:** Synthetic first. Then real cortical geometry from a bone mesh with an inner surface (UNKNOWN whether one exists locally), and for PDL a tooth root from the dental track. Coordinate with dental and do not move their data.

### H4. Moment arm sensitivity through wrapping events: shape parameter → moment arm (J, priority 4)

- **Idea:** First-order sensitivity `dr/db` is cheap and correct far from wrapping transitions. At a transition, where the muscle path starts or stops touching a wrapping surface, the curve has a kink. A derivative with a frozen contact location gives errors, analogous to A196’s error of 3.9 with frozen event time.
- **Mechanism:** The gap function `g(θ,b)` = the distance from the path to the wrapping surface. The event angle θ* is defined by `g=0`, and `dθ*/db = −(∂g/∂b)/(∂g/∂θ)`. Moment arm sensitivity is split into two cases: continuous within each regime and with a jump or kink at θ*.
- **Equation:** `r(θ;b) = ∂ℓ/∂θ`. Within a regime, `dr/db = ∂²ℓ/∂θ∂b` holds. The event’s contribution follows from ℓ being continuous but ∂ℓ/∂θ changing slope. For a fixed angle, the difference between regimes is `Δr(θ*)`, and it is shifted by `dθ*/db`.
- **Operation:**
  1. Use a two-dimensional or three-dimensional model of a joint with a cylindrical or spherical wrapping surface whose radius or position is controlled by b.
  2. Compare at 50 angles and 5 b values. (i) FD over b serves as truth. (ii) The derivative with a frozen regime. (iii) Event-aware derivative with a guard, which refuses when |θ−θ*| is less than a threshold (A196’s missing guard).
- **Representation:** The path as via-points and wrapping surfaces with stable identities carried through the morph (U272 contract).
- **Assumptions:** Geodesic path on the wrapping surface. No muscle forces, only kinematics.
- **Expected outcome:** The event-aware derivative agrees with FD outside the guard band. The frozen derivative has large errors near θ*. The guard refuses in a narrow band.
- **Strong baseline:** FD over b. It is cheap when b is low-dimensional, so a cost gain is required only if b is large.
- **Motprov:**
  1. If no wrapping transitions occur in the physiological angle interval, the hypothesis is irrelevant for this muscle.
  2. Om den frysta derivatan har samma fel som den eventmedvetna, faller mekanismen.
  3. If the guard band must be so wide that coverage falls below 80 %, the method is unusable.
- **Kostnad:** CPU i sekunder till minuter.

### H5. Exact inertia for multimaterial segments through second moments and affine deformation (J+B, priority 5)

- **Idea:** U278 and U288 carry `(V, S)` exactly, which suffices for mass and center of mass. Inertia requires the second moment `T = ∫ x xᵀ 1_m dV`. If T is stored disjointly per material, the segment’s inertia is preserved exactly during resolution changes and piecewise affine morphing. Resampling a density field (CT-HU → ρ) to FEM resolution does not.
- **Equation:** Under `x' = Fx + t`:
  - `V' = |det F| V`
  - `S' = |det F|(F S + t V)`
  - `T' = |det F|(F T Fᵀ + F S tᵀ + t Sᵀ Fᵀ + t tᵀ V)`
  - Inertia about the center of mass: `I = ρ(tr T_c · 𝟙 − T_c)`, where `T_c = T − S Sᵀ/V`.
  - The sum over child cells is exact when merging. Check per material, not total (N5).
- **Operation:**
  1. Take a synthetic density field with sharp boundaries between cortical and trabecular bone, and marrow if needed.
  2. Calculate I exactly at source resolution.
  3. Compare at coarse h: (a) stored `(V,S,T)`, (b) trilinear resampling of ρ, (c) majority label with ρ per region.
  4. Apply a piecewise affine morph and an RBF morph, and compare the error in I.
- **Assumptions:** Density is constant per material within the cell. If it is not (HU gradient), ρ-weighted moments are required, and exactness then applies only on the source grid. RBF is not affine and gives an error of order `‖∇F − F̄‖·h`.
- **Expected outcome:** (a) gives machine-precision exactness during resampling and affine morphing. (b) and (c) give a percentage error in I that increases with h. Under RBF, (a) follows the error bound.
- **Strong baseline:** Polyhedral mass and inertia per region mesh (Mirtich) are exact and cheap. The hypothesis only adds something when no region mesh exists, i.e. when the material comes from CT voxels or a field operation.
- **Counter-tests:**
  1. If (b) already gives I within a tolerance smaller than the variation between individuals (sourced), the gain is irrelevant.
  2. If a region mesh always exists in the collaborator's or BodyTwin’s chain, the baseline suffices and the hypothesis is unnecessary.
  3. The signed summation must not be used (A194). This is a prerequisite, not a counter-test.
- **Cost:** CPU in minutes.
- **Data:** Synthetic. Then a CT volume with HU calibration, if a license exists (UNKNOWN).

### H6. Fusion of multiple modalities that refuses for undeclared anatomical variants (B, priority 6)

- **Idea:** A boundary location, for example cartilage thickness or cortical thickness, is estimated from CT, MR and atlas prior. The sources share factors: the same atlas registration and the same version of the segmentation model. Naive fusion counts them twice and becomes confident but wrong (A207: conditional error 0.47 versus 0.03). Declared model families make confidently wrong decisions for undeclared features (A206). In A208, a continuous residual guard refused an underlying layer.
- **Mechanism:** Joint posterior over (thickness, shared calibration θ_reg, θ_seg) with each distinct factor counted once. The provenance key is (raw acquisition ID, algorithm version, registration ID). A residual guard refuses when `r = ‖y − ŷ(θ̂)‖²_Σ⁻¹` exceeds a χ² quantile calibrated *only on development cases*.
- **Operation:**
  1. Use synthetic observations from real meshes with implanted variants: osteophyte, focal cartilage defect and cortical perforation.
  2. Compare four arms: naive product, corrected for shared information, full joint distribution, and full joint distribution plus guard.
  3. Measure the number of confidently wrong decisions, coverage and Brier score.
- **Expected outcome:** The joint distribution plus guard gives few confidently wrong decisions on the variants and retains at least 90 % coverage on declared cases.
- **Strong baseline:**
  1. Anchor at a declared location (A206).
  2. The corrected naive fusion with correct θ that A207 added (conditional error 0.19).
  3. Only a confidence bound on the posterior.
- **Counter-tests:**
  1. If the guard refuses in more than 10 % of normal development cases, or misses more than half the variants, it fails.
  2. If the anchor reaches the same number of confidently wrong decisions, the mechanism is unnecessary.
  3. If re-encoded data with the same raw source are treated as independent, the provenance key is wrong (A207 `ADV_PROVENANCE_REENCODE`).
- **Cost:** CPU in minutes. The model is finite and can be computed exactly with fractions as in A207.
- **Data:** Synthetic on real meshes. Noise models per modality must be sourced.

### H7. A task metric in shape space may differ from anatomical variance ordering (J+B, priority 7)

- **Idea:** The seed’s prediction geometry gives a metric on shape space: the pullback `G_Q(b) = J_Qᵀ Σ_Q⁻¹ J_Q`, where `J_Q = ∂Q/∂b` and Q are the biomechanical outputs. It measures how much a shape change affects the task. It is not an anatomical distance in mm. The PCA ordering in an SSM is an mm metric. If the two orderings differ, an SSM truncated by task relevance may reach the Q tolerance with fewer modes.
- **Operation:**
  1. Calculate `J_Q` via FD for Q = {m, I, r(θ) for some muscles}.
  2. Rank the modes by `σ_k²·‖J_Q e_k‖²_Σ⁻¹` and compare with σ_k² (PCA).
  3. Compare the error in Q after truncation to K modes for both orderings, on individuals held outside the estimation.
- **Expected outcome:** For inertia and mass, the ordering follows PCA, because the large modes are also scale. For moment arms, small local modes at attachments may dominate.
- **Strong baseline:** PCA truncation at the same K. Also the rule "use all modes", if K is small.
- **Counter-test:** If the Spearman correlation between PCA ordering and task ordering is ≥ 0.9, PCA suffices and the hypothesis fails. Compare U290, where predictive variance ranked as task value with Spearman 0.994. If truncation does not reduce the omitted error in Q at the same K, the hypothesis also fails.
- **Condition:** A real SSM with correspondences is required. Without it there is no mapping, and the hypothesis cannot be tested. Label the result as a mathematical identity plus an empirical question, not as a new finding.
- **Cost:** CPU in minutes once an SSM exists.
- **Data:** SSM with a license (D2). Status UNKNOWN.

### H8. Boundary-aware Q2 trace patch for repeated local strain queries in bone with a material boundary (B, priority 8, most engine-oriented)

- **Idea:** Extends U285/A201 to a tissue recipient. The output is strain at a point, for example cortical strain in the femoral neck or an implant seat, for K load cases. The geometry is frozen, and the material boundary between cortical and trabecular bone crosses the patch.
- **Mechanism:** `V = Q1 ⊕ span{edge and face modes within radius 2}`, mirrored at the boundary. The coarse factor is reused via Schur over 252 modes, with fully charged cost.
- **Counter-tests:**
  1. If uniform refinement plus one factorization for all K right-hand sides is cheaper at the same accuracy, the hypothesis fails. This is the strongest baseline and it has not been tested in A201.
  2. If the patch fails when the material boundary cuts it (smooth solution assumed), the hypothesis is restricted.
  3. Corner cases according to `adv_corner18`.
  4. Features below coarse h fail as predicted in A201.
- **Expected outcome:** Gain only at K ≥ K*. It may be absent entirely, because A201 had 0/13 speed gains among cases that met the requirements.
- **Cost:** Requires vectorized assembly before any cost conclusion (A185). Therefore lowest priority.

---

## 5. Prioritering

| Prio | Hypothesis | Track | Why first | Depends on |
|---|---|---|---|---|
| 1 | H1 task-driven error bound | J+B | Directly answers "what resolution is required" and is cheap. Requires no new solver. N1 and N2 are built in from the start. | A real bone mesh (O2/D2) |
| 2 | H2 identifiability and task-based OED | J+B | Answers "what can be distinguished and which measurement decides". The lessons from A206, A202 and A186 apply directly. | Synthetic model; σ sourced |
| 3 | H3 thin layers via stiffness | B (+dental) | Reformulates the lessons from A194 and A204 into a recipient-based check. Shared with the dental track via PDL. | U288 code (read-only), fine FEM |
| 4 | H4 event-aware moment arm | J | The most collaborator-specific and cheap. Transfers A196 with an explicit condition. | Wrapping model |
| 5 | H5 second moments for inertia | J+B | A strong baseline already exists (polyhedral). Relevant only without a region mesh. | — |
| 6 | H6 fusion with refusal | B | Important for uncertainty but depends on noise models. | Noise per modality |
| 7 | H7 task metric versus PCA | J+B | Requires an SSM. | D2 |
| 8 | H8 Q2 patch | B | Engine-oriented and cost-limited by Python. | Vectorized assembly |

## 6. Open gaps affecting the hypotheses

- Which bone meshes, SSMs and CT volumes with licenses exist locally: UNKNOWN (O2/D2).
- Sourced tolerances τ_Q for biomechanical quantities (variation between individuals, measurement noise): missing. They must not be replaced with default values.
- Rigorous interval bounds for a mesh SDF, unlike an analytical expression: not demonstrated.
- Contact with sliding between tissues, fiber directions and an uncertainty channel in the handoff contract: missing (U272’s capability map).
- Refusal path in `lastfalt_v1_fem`: missing (A204).
