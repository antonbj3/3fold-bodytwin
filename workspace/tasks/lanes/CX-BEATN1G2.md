# CX-BEATN1G2 — round 2: CERTIFIED correction of N1g + the leads from round 1

Read first:
- `results/CX-BEATN1G/RESULTS.md` (A324: H1 wins in DM/SC/PS, but JW falls 0.35→0.83).
- `results/CX-KNEENET/RESULTS.md` (A325: activity×GRF interaction 0.367 < N1g 0.379; knee angle has the highest measurement value; MAP wins in DM/SC/PS).
- `results/CX-EARLYSTANCE/RESULTS.md` (A326: the L1 ID moment differs 49 % from OpenSim ID in DM6; contact-weighted recruitment −48 %).
- `results/CX-BEATN1G/HINT_FROM_GRAPH.md`: the graph engine's margin_net/decision_cert/alarm, and HMM marginals instead of hard phases.
- The code in `~/projects/graph_workspace/pub/src/graph_engine/`.
- Reuse the code in `results/CX-BEATN1G`, `results/CX-KNEENET` and `results/L1`. Do not redo what is there.

## Hypotheses (new PREREG.md + sha256 before the first run; same LOPO, 113 trials and windows as A324)
- C1 certified clipping: H1's feasible set as margin_net requirements with a declared shared model error (`shares`) between moment, capacity and geometry. Clip N1g only where decision_cert gives p_flip ≤ α (α = 0.05 and 0.2, both preregistered). Report the fraction of frames corrected, the error on exactly those frames, and N_eff for the physics terms. Hypothesis: JW's collapse disappears because JW's uncertain geometry (DM transplant) gives wide margins.
- C2 correct ID moment: H1 and C1 with OpenSim ID moments (`/mnt/shared_data/bodytwin_work/X1b/*_id.sto`, where they exist) instead of L1's b[:,3]. Does the feasible set become right?
- C3 preregistered activity×GRF interaction (K6: 0.367) plus knee angle as an additional term (K5), as a new null model candidate N2. Report it as a stronger null model if it holds, not as a physics win.
- C4 the best of C1–C3 combined, only if each part holds separately.
- Placebo for each hypothesis (shifted moments/limits/angle). Criterion as in A324: ≥ 5 % below N1g in ≥ 3/4 LOPO rotations and not worse in early stance. Symmetric skepticism.

## Resources
- Locally: nice, 2 threads, ≤ 60 s per test.
- Write only in `results/CX-BEATN1G2/`.
- `RESULTS.md` starting with `# CX-BEATN1G2`, plus `results.json` and code with pytest.
