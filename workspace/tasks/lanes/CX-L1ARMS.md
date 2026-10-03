# CX-L1ARMS — is L1's MUSCLE GEOMETRY (knee moment arms A, contact projection cj) the reason the measured force lies below Fmin?

## Background (all other geometric explanations have been tested)
- A1438 CX-KNEECENTER: the functional knee centre/tibial axis gives ΔC0 ≈ 0 N, and the frames below Fmin INCREASE (20 %→30 %).
- A1435 CX-EPSBAND: C0 = −F_knee·y_tibia, with no independent geometric correction in the data.
- A1425 U447 (Field): the frames below Fmin are repaired by C0 −312 N OR by a 3.9 % change in contact arm; F0 widening repairs only 13/106.
- A1418 CX-MINCONTACT: at the LP minimum, rectus femoris carries 54–66 % of the knee moment, and gastroc 0 %.
- L1 code: `results/L1/code/` (how A, cj, Mkx, Lmt are built from the model: which OpenSim model, which muscle paths/wrapping), plus `results/L1/prep/*.npz`.

## Tasks (PREREG.md + sha256 BEFORE comparing against the implant force)
1. **Inventory:** which musculoskeletal model and which knee muscle paths L1 uses. Tabulate L1's knee flexion moment arm per muscle group as a function of knee angle:
   - rectus femoris, vasti, hamstrings (semimem, semiten, BF long/short head), gastroc, gracilis, sartorius, TFL.
   Compare against published references: e.g. Buford 1997, Arnold 2010, Rajagopal 2016, Herzog & Read 1993, Spoor & van Leeuwen 1992. State the source and the table.
   Do the same for cj, each muscle's contribution to tibiofemoral compression per unit force: the projection of the muscle force onto the tibial axis at the knee. Is it reasonable?
2. **Swap test with no implant fitting:** replace L1's knee-crossing moment arms and projections with a public model's (Rajagopal 2016 or Arnold 2010, via OpenSim/opensim-python if available locally; otherwise the literature's angle curves, frozen before scoring). Rebuild A/cj on L1's 113 trials, keeping everything else as it is. Compute the LP lo/hi and count the frames below Fmin.
   Frozen criterion: the share of frames below Fmin falls ≥50 % without the share above hi growing, AND the person-median of (meas−lo) moves towards a constant.
3. Sensitivity: which single muscle group's arm/projection, scaled ±10 %, moves lo the most? That identifies the node.
4. Counter-test: a random ±10 % perturbation of the arms must not give the same improvement.

Deliver RESULTS.md starting with `# CX-L1ARMS`, results.json, the scripts, and pytest. Run the LP under bigmem.lock with 2 threads, or on OVH via `tasks/cloud_run.sh` (finish by 07:30). Internal data stays local. Exclude jw_lungef1. lane runner has full permissions in the workspace; `~/projects/bodytwin` is read-only. Every outcome is a node that expands.
