# Motion Engine for BodyTwin: what's usable and what's missing (M1, 2026-09-22)

The engine reviewed is `3fold-motion-engine` on aae29ed (public checkout, read but not changed). Two components are tested on real human gait (gait2392 subject01, force plate). Numbers and deviations are in `results/M1/README.md`, the criteria in `results/M1/PREREG.md`.

## Svaret i korthet

- **Inverse dynamics works.** The motor's URDF path (`pin.buildModelFromUrdf` followed by `pin.rnea`) gives OpenSim's joint moment of 3·10⁻⁹ N·m. The engine's own batchable RNEA (`rnea_warp`) kernel gives the same result of 10⁻⁵ N·m in float32. The requirements: the model is converted with coupled joints, gravity is applied along y, and GRF is applied as an external wrench.
- **The identifiability gate works on humans.** It shows that a gait cycle with force plate only identifies:
  - total mass (72,50 ± 0,18 kg, model says 72,6);
  - The first part of the thighs and torso;
  - torso inertia and parts of thigh inertia.

  No individual segment mass and no lower leg or foot parameters can be identified.
- **Contact, events, mesh collision and IK are not tested.** The map below shows what each part requires.

## Map: engine component → link in the chain → K1 leaf

The Status column indicates:
- **tested M1**: run on real time, with number;
- **executable**: the code exists and fits, but has not been run here;
- **requires building**: the code is there, but the human core is missing;
- **not relevant**.

| Engine Component | File:Function | Biomechanical Link | K1 Leaf | What Needs to Be Adapted | Status |
|---|---|---|---|---|---|
| URDF loading (one code path) | `motion_stack.py:MotionStack.__init__` → `pin.buildModelFromUrdf` | segment model: joints, frames, inertia | `jnt.anatomical_frame` , `inr.*` (consumer) | osim→URDF generator (`results/M1/m1_opensim_ref.py:write_urdf`):<br>• ball joint = 3 revolute with massless intermediate links<br>• knee coupled translation = 2 prismatic joints, prescribed via the chain rule and projected on θ<br>• gravity (0, −9,81, 0) instead of pinocchio's z<br>• m and kg·m² (H5 and D1 works in mm and kg·mm²)<br>The `MotionStack` class itself is unusable: it requires collision meshes and effort boundaries, prunes "always colliding" pairs, and derives with `np.gradient` | **tested M1** (A1) |
| Inverse dynamics, CPU | `motion_stack.py:MotionStack._effort` (`pin.rnea`) | segment parameters → joint moments | `id.newton_euler`, `id.cut_chain` (`data.f` = net joint reaction without muscles) | external forces (fext on calcn); floating base like 3P+3R; projection for coupled coordinates; derivation outside the engine | **tested M1**: 12 moments and 6 residuals, RMS ≤ 3·10⁻⁹ N·m |
| Inverse dynamics, GPU core | `rnea_warp.py:_build_kernel`, `WarpRNEA` | batched parameter sweeps (H5 rules, joint center offset K1) | `id.newton_euler`, `jnt.center_shift` | The core only supports serial chain, revolute around z, fixed base, no external forces and fixed `cuda:0`. In M1 solved with axis alignment (jR), pelvis acceleration as gravity per frame, −JᵀF and CPU-device. For trees, prismatic joints and fect, a new kernel is required | **tested M1** (A2): max 1,6·10⁻⁵ N·m, right leg only |
| Inertia regressor | is missing from the engine; pinocchio's `computeJointTorqueRegressor` was used | identification of mass and inertia from motion + GRF | `inr.*` (control), H2 | of the engine "identification on Panda" is friction on top of URDF inertia (`scripts/moment_gate.py`, `friction_aware_effort_backend.py`) | pinocchio, **tested M1** (B) |
| Identifiability gate | `identifiability_governor.py:governor_gate` | which segment parameters are distinguishable (H2) | H2, `inr.affine_rule` | Column scaling before the gate; the gate says UNVERIFIABLE as soon as a single direction is zero, so the response per parameter (projection on the null space) must be calculated separately | **tested M1**: null space 51 by 120; 23 individual parameters identifiable |
| System identification by contact | `ncp/sysid_gn.py` (adjoint-GN, Talos + box) | mass and μ out of contact impulses, i.e. foot–ground instead of hand–box | `id.grf_predict` | The scene is hardcoded for Talos and box. The foot must be modeled as a contact body. The measurement is joint moment, which man lacks; GRF would have to be the measurement | requires construction |
| Contact solver | `contact_engine.py:SplitImpulseEngine`, `ncp/ncp_ref.py:solve_ncp_*`, GPU-variants | foot–ground contact, predicted GRF | `id.grf_predict` (predicted branch) | foot geometry (heel/toe spheres or SDF), stiffness, shoe–floor friction; validation against measured GRF (found here) | requires build |
| Contact sensitivity and events (P2) | `ncp/ncp_ref.py:sensitivity`, `sensitivity_fd` | heeling and toe off as events, sensitivity through the contact change | is missing in K1 (events do not exist as leaves) | In M1 the events were taken with a threshold of measured GRF (Fy > 20 N), not with the engine. P2 tool requires a contact problem (G, b, μ) from the foot | requires build |
| Support and load regions (P4) | `reproducibility/p4/lp_reference.py:solve_lp`, `friction_polygon.py` | static support surface and load limit (balance, the collaborator's box lift with feet in floor contact) | `jr.reaction` (statics) | foot contact points +  friction cone; joint torque limits as caps | drivable |
| Exact mesh collision along motion | `motion_stack.py:MotionStack.verify_path_exact` (coal); `scripts/collision_free_trajectory_cert_swept_ccd_plus_margin_binding_segment_witness.py` | bone–bone penetration after morphing (D1/P1), joint boundaries | `reg.surface_nonrigid` → `jnt.*` (control) | bone meshes per segment in URDF (m, segment frame); knee and hip are **always** in contact, so cropping in `MotionStack` removes exactly the pairs to be checked; a distance or penetration depth measure is needed, not just a binary collision | requires construction |
| SDF composition with the field engine | `examples/compose/field_sdf_to_contacts_to_graph.py` | soft tissue/bone as field → contact → graph node | `inr.soft_tissue`, `jr.contact_fe` | SDF out of BodyTwin geometry (mm → m) | drivable |
| Multistart-IK (0,1 mm) | `scripts/movement_ik.py:ik_multistart`, `ik_solve` | markers → leads | `kin.ik_lsq` | Solves an end-effector pose, not weighted least-squares over ~40 markers. Requires marker residual, weights and ball joints. The multistart is useful against local minima | requires construction |
| Derivation and filtering | missing (`_effort` uses `np.gradient`) | q → q̇, q̈ | `kin.diff` (partial) | In M1 used OpenSim's chain: pad, IIR 6 Hz, GCV-spline degree 5 | missing |
| Torque gate | `dynamic_feasibility.py:feasibility`, `scripts/moment_gate.py` | need versus capacity (joint torque or muscle strength) | `mus.*`, `rec.*` | angle- and speed-dependent capacity instead of constant effort limit | operable with adaptation |
| Friction models (8) | `scripts/friction_tournament_harness.py` etc. | joint friction (negligible in man); shoe–floor | – | – | not relevant |
| Humanoid gait cells | (README "Humanoid") | gait as a reference | – | the gates do not pass | not useful as a reference |
| Exact particle count (P3) | reproducibility/p3 | – | – | – | not relevant |

## Measured in real time

| Test | Resultat |
|---|---|
| A1: engine against OpenSim-ID (same model, kinematics and GRF) | 12 joint moments and 6 residuals: RMS 2,7–3,1·10⁻⁹ N·m |
| A1k: knee as usual URDF-revolute | knee 1,9 N·m RMS (3,3–4,1 % of top), hip flexion up to 7 N·m, pelvic residuals up to 9,7 % |
| A2: `rnea_warp` core against pinocchio | ≤ 1,6·10⁻⁵ N·m |
| A3: H5 scale rules (thigh +10 % length) vs exact affine transport | peak Δ hip moment: 2,15 N·m (isotropic, mass conservation), 5,14 (isotropic s = ls), 4,85 (no update). Prediction was 2–14. Measured swing-α 18,8 rad/s² vs H5's assumed 50–200 |
| B: identifiability from movement + GRF | see below |

Details for B:
- Governor gave UNVERIFIABLE, with rank 69 of 120.
- The corridor reaches the entire structural rank. Four of the directions come from the mm translation of the knee (65 with fixed knee).
- No segment mass is individually identifiable.
- Total mass ± 0,18 kg.
- Least-squares in the row space lowers the basin residuals by 19–64 %.

## What it takes for the engine to become a BodyTwin consumer

1. **ModelAdapter.** An osim or the reference model-to-engine adapter that carries coupled coordinates, such as a `coupled` entry next to URDF with spline and derivative, plus gravity direction and SI. M1's generator is a working first version for OpenSim CustomJoint.
2. **Inverse dynamics input.** A ID input that takes (q, q̇, q̈), external wrenches per body and fluid base, and gives generalized forces in the source model coordinates. That's about 60 rows (`m1_engine.py:run_id`, `gen_forces`). The derivation must be outside the engine or be a separate, tested part.
3. **General GPU-RNEA.** `WarpRNEA` needs to be generalized to tree topology, arbitrary axes, prismatic joints, fext and a device parameter. Then thousands of parameter variations can be run on GPU in one call: H5 rules, K1's joint center shifts, and Monte Carlo over segment parameters. Today, one call is required per picture and leg.
4. **Inertial regressor as a motor component.** It is combined with `governor_gate` and the response per parameter, and then becomes H2's tool: which movements or measurements make the segment masses separable. B shows that more gait cycles do not help.
5. **For contact and events.** Foot model (geometry, stiffness, μ), GRF-prediction validated against measured GRF, and P2 sensitivity by heel placement. It is the largest unbuilt step.
6. **For bone-to-bone collision.** Bone meshes in segment frames (m), no pruned joint pairs and a penetration depth instead of a binary collision.
