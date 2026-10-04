# CX-CONNECTIVE — the graph node MSK-CONNECTIVE: tendon/ligament tension and cartilage load-bearing in the knee, validated against public cadaver data

Background: `python3 ~/projects/bodytwin/scripts/anchor_graph_tools.py --graph data/MECHANISM_ANCHOR_GRAPH.json priority` ranks MSK-CONNECTIVE highest in biomechanics (reach 236, status ASSUMED, evidence 0). Claim: "Tendon-tension/cartilage load-bearing via mechanobiology/tensile/ultrasound legs vs cadaveric knee-FE ligament-force model." Anton (24/9): the knee track should go here. This is ONLY musculoskeletal biomechanics. The graph's private soft-tissue cluster (TIS-*) must NOT be touched.

## Build on (read first)
- The graph node plus its neighbours, and `~/projects/bodytwin/docs/MECHANISM_TENDON_*`, `MECHANISM_CONTACT_MUSCLE_DECOMP.md`, `MECHANISM_JOINT_FORCE_SCORECARD.md`.
- `~/projects/bodytwin/scripts/msk` (grep ligament, cartilage, tendon, openknee, contact).
- Today's results:
  - A312/A324/A358/A362–A368: the knee chain falls; the patella and operating range are identified as causes;
  - CX-ROLLBACK (running);
  - CX-LUNGEID;
  - Field's F-8: contact from measured kinematics is ill-posed, so use force-controlled contact.
- Field's cartilage contact U384 `cartilage_contact(R,t,E,ν,δ)` (64/64 within 2 % of FE, not audited): `../3fold-motion-engine/_private/romi_collab/build/U384/code/` (contact_api.py, regime.py, poroelastic.py). Use it via an adapter with a source reference, only in its valid regime.

## Facit (public)
- Look for cadaver knee data with ligament forces and cartilage contact, starting with Open Knee(s) (SimTK). The lane runner sandbox has network access. Also check local datasets (`tasks/index/DATASETS.json`, grep openknee/ligament/cadaver). Download to `external_media` and record the source + licence.
- Literature values with DOI for tendon/ligament stiffness and cartilage thickness/modulus.

## Tasks (PREREG.md + sha256 first)
1. **Build** a knee module: ligaments (ACL/PCL/MCL/LCL as nonlinear springs with slack), the patellar tendon, and cartilage contact via U384 (medial/lateral), with force-controlled contact (Field's lesson).
2. **Test** against the cadaver facit: ligament forces and contact force/area under the published load cases (e.g. anterior drawer, varus/valgus, flexion under compression). Criterion: errors within stated tolerances. Null model: rigid contact without ligaments.
3. **Connect** it to the knee chain: does the module change L1's early-stance overestimate (A312: +0.7–1.0 BW) or the medial/lateral split (compare with the F-8 facit only through Field's sealed procedure — do not open eTibia yourself)?
4. Update the graph node's claim/evidence as a proposal in `results/CX-CONNECTIVE/GRAPH_UPDATE.md`. Do not write in the graph (read-only).

Resources: lane runner has full permissions in the workspace. Internal TLEM/the collaborator data stays local or on OVH. Write in `results/CX-CONNECTIVE/`. `RESULTS.md` starting with `# CX-CONNECTIVE`, plus results.json and code with pytest.
