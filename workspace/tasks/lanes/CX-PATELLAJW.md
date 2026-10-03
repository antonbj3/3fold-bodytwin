# CX-PATELLAJW — does L1's patella solver (ρ clipped at 1.3) give a wrong quadriceps arm, and what happens when it is corrected?

Background: A1450 CX-QUADARM. For JW at ~50° knee flexion, L1's effective quadriceps arm is 4.40 mm against Rajagopal's 36.68 mm. The median ρ is clipped at 1.3, and the derived patellar tendon arm is only 4.20 mm. JW has 48.7 % of frames under lo (the worst person).
Read: `results/L1/code/n12_model.py` (the patella solver, rho, bounds 0.4–1.3), `results/CX-QUADARM/` (mechanism.json, geometry.py), `results/CX-L1ARMS2/`, and the lane briefs `tasks/lanes/CX-L1ARMS-3D.md` and `CX-L1ARMS-RECTUS.md` (proposed by an earlier lane).

## Tasks (PREREG.md + sha256 first)
1. **Diagnose the solver per person/angle:**
   - how often ρ hits 0.4/1.3;
   - where the patellar tendon arm is < 15 mm;
   - whether it is the solver branch, the knee axis, or the JW/DM transplant geometry (the tibial tuberosity/patellar tendon attachment).
   Compare against the literature: patellar tendon arm 35–50 mm at 30–60° (Herzog–Read, via results/CX-RQANGLE); ρ from Grood/van Eijden.
2. **Correct WITHOUT fitting to the implant:**
   - (a) remove the clipping and solve the patella equilibrium properly;
   - OR (b) replace ρ(θ) with a frozen literature curve (van Eijden 1986/Grood);
   - (c) both.
   Rebuild A[3] and cj for the quadriceps elements consistently (ρ affects BOTH the arm and the patellofemoral/tibiofemoral projection).
3. **Frozen criteria:**
   - JW gait: under lo ≤ 25 % (from 48.7 %);
   - all 4 persons: under lo ≤ 13 %;
   - JW4 non-gait: a reduction ≥ 40 %.
   Also report the ε=10 % coverage, and lo+c/stress-2 against N1g.
4. Deliver `knee_geometry(model='l1_patella_fixed')` + pytest.

Deliver RESULTS.md starting with `# CX-PATELLAJW` and results.json. Run under bigmem.lock with 2 threads. Internal data stays local; no jw_lungef1. lane runner has full permissions in the workspace; `~/projects/bodytwin` is read-only. Every outcome is a node that expands.
