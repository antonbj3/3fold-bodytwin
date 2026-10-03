# CX-MUSCLE-CT2 — T2 on all 30 VSD persons (run on OVH, which has internet and 47 GB of free disk)

Read first: `results/CX-MUSCLE-CT/RESULTS.md` (A360, pilot on 5 persons: shape did not beat the BMI mass law, 27.1 % vs 24.4 %). Reuse analyze.py, the PREREG criteria and the HU thresholds unchanged. This is the same question at the preregistered size.

## Tasks
1. New PREREG.md + PREREG.sha256 (same criteria as CX-MUSCLE-CT; state that this is the full sample and that the pilot has been seen).
2. Run on OVH: `ssh -i ~/.ssh/hunt_20260923 -o IdentitiesOnly=yes -o UserKnownHostsFile=~/research/sol6_recovery_20260923/CLOUD_HUNT_20260923/infra/known_hosts rocky@54.37.4.45`. Use sudo systemd-run units via `tasks/cloud_run.sh` (BodyTwin ≤ 12 vCPU in total, and the CX-WHATIF shards already use some; check with `status` and keep ≤ 4 vCPU for this). Download the VSD archives from Zenodo 8302449 directly on OVH into /opt/bt/data/vsd/, which is public data (CC BY-NC-SA). Process there and bring home only the per-thigh tables (not CT volumes).
3. Also add M-f: M-e + sex (from VSD metadata). A359 showed that the mass laws are biased by sex (men underestimated, women overestimated).
4. LOSO on 30 persons (60 thighs). Report the median absolute error per model, sex residuals, and the permutation counter-test.
5. Finish by 22:30 (OVH lease 23:00). Clean up /opt/bt/data/vsd afterwards.

Write in `results/CX-MUSCLE-CT2/`. `RESULTS.md` starting with `# CX-MUSCLE-CT2`, plus results.json and code with pytest.
