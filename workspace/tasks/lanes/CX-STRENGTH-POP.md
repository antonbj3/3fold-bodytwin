# CX-STRENGTH-POP — T1: which scaling law predicts an individual's measured strength? (NHANES, ~5,000 people)

Facit: NHANES 1999–2002, muscle strength (MSX: isokinetic knee extension 60°/s, peak torque) + DXA (total and leg lean mass, fat %) + height/weight/age/sex. Download from CDC.
Models (predict peak knee extension torque per person; CV 10-fold + held-out survey cycle):
- L0 population median;
- L1 the ScalingLengthMass baseline (strength ∝ mass);
- L3 PCSA ∝ M^(2/3);
- L4 ∝ height² (as in A121);
- L5 fat-free mass (DXA);
- L6 DXA leg lean mass (= "the individual's own composition", the M4 idea at population level);
- L7 L6 + thigh length (moment-arm scaling).
Report RMSE and R² per model, the fitted exponent p in strength ∝ mass^p, residuals by sex/age/BMI, and where the reference model's laws go most wrong (e.g. high BMI). Counter-test: permuted DXA. Criterion: L6 beats L1 and L2 by ≥ 10 % RMSE in CV and in the held-out cycle.

## Common
- PREREG.md + PREREG.sha256 BEFORE the first run: facit, criteria with numbers, counter-tests (permuted inputs), and what counts as a fall. Symmetric skepticism.
- The lane runner sandbox has network access (public data may be downloaded; save it on external_media and record the URL + sha256). Locally: nice, 2 threads, ≤ 60 s per test; heavier work via Modal (`tasks/modal_run.py`) for public data, or OVH (`tasks/cloud_run.sh`, BodyTwin ≤ 8 vCPU) for internal data.
- Write only in `results/CX-STRENGTH-POP/`. `~/projects/bodytwin` is read-only. `RESULTS.md` starting with `# CX-STRENGTH-POP`, plus results.json and code with pytest.
