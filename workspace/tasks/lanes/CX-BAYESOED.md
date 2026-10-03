# CX-BAYESOED — Bayesian, decision-relative OED: which measurements make the individual's knee contact force (N_med/N_lat) DETERMINED enough, per unit cost?

## Background (read these)
- The register A1384 (BODYGRAPH: gait2392 rank 17, gauge 75), A1385 (U401: the collaborator nullity 1126), A1389 (U405: stiffness+EMG do not close it), A1394 (U413: patellar tendon + separate medial/lateral hamstring tensiometry closes it, rank 26, width 449/402 N vs 590/848 N), A1392 (CX-SETVALUED: the hip spread is 37.6 BW), and A1387/A1390 (EMG adds no robust information in the GC data).
- Field's code: `the public staging tree/3fold-motion-engine/_private/romi_collab/build/cloud_results/CLOUD-F-BODYGRAPH/`, `build/U405/`, `build/U413/` (gait2392, sensor models, posterior). Reuse them with a source reference; Field owns the rank analyses, and you own the OED/protocol.
- The graph session's input (the principle for this lane): `results/CX-BAYESOED/GRAPH_INPUT.md` (copied below).

## The graph session's principles (mandatory)
1. **c-optimal, not D-optimal:** the contact force is ONE output c=gᵀθ in the gauge space. Minimise Var(gᵀθ|S)=gᵀ(Λ0+Σ_{o∈S}F_o)⁻¹g per unit cost. Close only the part of the gauge that g projects onto.
2. **An exact check first:**
   - `GE/observability_sigma_min_cert.py`: the output is identifiable from S ⇔ g ∈ rowspace(J_S). Right singular vectors name the remaining combinations.
   - Operator gauge test (`C24/lanes/S3_identifiable_patches/branch_b.py`): a sensor class that cannot reach g gets excluded before any Bayes step.
   - Output-only closure: find the SMALLEST sensor set with g ∈ rowspace(J_S) (an exact combinatorial question).
3. **Machinery:**
   - `GE/identifiability_oed.py` (greedy on ΣF_o, with the objective swapped for gᵀΣ_post g);
   - `C24/lanes/XD3_lowrank_det_oed/lowrank_oed.py` (rank-k updates, 10–32× faster; for a single output, Sherman–Morrison on gᵀΣg);
   - `GE/plan_value.py` (cost, reliability, a sequential stopping rule);
   - `GE/decision_cert.py` (p_flip);
   - `GE/margin_net.py` (shares/N_eff);
   - `GE/alarm.py` (a sequential residual test catching a misspecified model).
   Here GE=~/projects/graph_workspace/pub/src/graph_engine and C24=~/research/inference_training_20260921/claude_24h.
4. **Lessons:**
   - Correlated sources are not independent: gait cycles, steps and channels from the same person, EMG sharing electrodes/normalisation, tensiometry sharing calibration. Declare shares, report N_eff, and make the certificate unit = PERSON.
   - FC: ICC + marginal error does not give the risk floor, so specify the full joint error model.
   - Rank at ONE linearisation point is local: run a multistart collision test (distinct θ with identical sensor responses) before calling a set identifying.
   - At one-sided constraints (F ≥ 0, ≤ Fmax): NO Laplace/Fisher-only. Use constrained sampling (truncated Gaussian or polytope HMC with reflection) plus coverage checks on simulated truth (cf. Field U403/OED-CALIB).
   - **The prior's share:** along g, report Var_post/Var_prior, the width with each prior component removed (physiology, population, EMG), and that data-only width = ∞ when g ∉ rowspace. If the gauge is not closed, write plainly "the width is 100 % prior".
5. **Angles:**
   - **Decision-relative OED:** value = the expected reduction in decision loss (e.g. "medial overload N_med > threshold", or a surgical choice), via plan_value + decision_cert, not reduced width.
   - **Sequential/adaptive:** cheap sensors first, expensive ones (tensiometry/ultrasound/implant) only if the decision is still open → a personal minimal protocol.
   - **Two-tier population:** expensive sensors on a few persons to learn a hierarchical prior over the gauge, cheap ones on everyone. OED over PERSONS × sensors.

## Tasks (PREREG.md + sha256 first)
1. Models:
   - gait2392 (Field's BODYGRAPH/U413, public);
   - the collaborator's system (john_WR: `/mnt/shared_data/bodytwin_work/share_field/john_WR/`, INTERNAL, local only);
   - the Grand Challenge persons, where the implant force exists as facit for the decision's truth (NOT jw_lungef1).
2. Candidate sensors with cost and noise, from literature with DOI: surface EMG per muscle, fine-wire EMG, shear-wave tensiometry per tendon (patellar/Achilles/semimembranosus/biceps femoris), joint stiffness perturbation, IMU, ultrasound fascicle length, force plate, the implant (as reference).
3. Deliverables:
   - (a) the output-only minimal sensor set, exactly;
   - (b) a Pareto front of cost against gᵀΣ_post g, with the prior's share per point;
   - (c) a decision-relative ranking for "medial overload" with p_flip;
   - (d) a sequential protocol per person with a stopping rule;
   - (e) a two-tier population design;
   - (f) a coverage check on simulated truth, person-level.
4. Everything is a node that expands (FRAMING_NOTE): report where the protocol cannot close g, and which measurement node would be needed.

Resources: lane runner has full permissions in the workspace; the graph and Field code is read-only (copy with a source reference). Heavy sampling goes to OVH (`tasks/cloud_run.sh`, BodyTwin ≤ 12 vCPU, finish by 07:30) or UpCloud; internal the collaborator data stays local/OVH only. Write in `results/CX-BAYESOED/`. `RESULTS.md` starting with `# CX-BAYESOED`, plus results.json and pytest.
