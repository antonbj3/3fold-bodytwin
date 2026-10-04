# CX-D1PARITY — D1 (plan §11/§10): Field's exact shape derivatives (U380) through BodyTwin's solver → dF/du with parity, and the effect against uncertainty

Read first:
- `results/CX-WHATIF2/RESULTS.md` (A369: what-if curves, implicit KKT derivative 0.00035 % at the peak step, thigh inertia frozen).
- `results/CX-WHATIF2/U380_FROM_FIELD.md` (format, paths, caveats).
- `results/CX-FIELDSHARE/femur_edit.py`.
- `external_research_path` §10.1 and §11 D1.

## Tasks (PREREG.md + sha256 before the first run)
1. Adapter: U380 `shape_massprop` (bone V, c, I, dV, dc, dI about the centre of mass, density 1 → multiply by cortical/trabecular ρ, state the choice) → the thigh segment's mass properties in N40b/CX-SOLVER2's inverse dynamics. Bone is only part of the segment. Add the bone's change to the segment's soft-tissue share (keep the soft tissue fixed; state it).
2. Parity: dF_hip/du (u = twist, varus, lengthening) through OUR solver, analytically (U380 + our implicit KKT derivative on a frozen active set) against central FD with a full rebuild. Four femurs (C01LFE, C01RFE, C02LFE, z001) or those that have a matching N40b individual; state the mapping. Criterion: ≤ 1e-4 relative where the active set is stable (D1's own fall criterion).
3. Size: the inertia effect's share of the total dF/du, compared with the attachment path from A369. Hypothesis: < 1 % of the force effect (bone ≪ soft tissue, A181). If so, say so plainly — that is a result, not a failure.
4. First optimisation step (D1 bold): which u minimises the peak hip force in the collaborator's box lift within ±10° / ±10 mm? Use the gradient + projection. Verify the optimum with a full rebuild (not Taylor, Field's caveat). Report whether the optimum lies inside the model's uncertainty band (N2b A363).
5. Counter-tests: permuted u axes; another person's femur.

Resources: lane runner has full permissions. Internal the collaborator/restricted model data stays local or on OVH (`tasks/cloud_run.sh`, BodyTwin ≤ 12 vCPU, finish by 23:59), not Modal. Read-only in `3fold_staging` (copy U380's code with a source reference). Write in `results/CX-D1PARITY/`. `RESULTS.md` starting with `# CX-D1PARITY`, plus results.json, code, and pytest.
