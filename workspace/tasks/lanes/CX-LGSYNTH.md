# CX-LGSYNTH — compile the test of the parameter-free knee law on 50 trials (BT-LG) + build the next The swarm wave

## 1. Compile
Compile `results/BT-LG-*/` (current run, together with `_run1/` where it exists; the current run has explicit implant side/BW/force column, and _run1 guessed the plate).
- Per trial: person, activity, chosen plate/side, number of frames, k_law vs k_meas, peak error (BW, %), RMSE, correlation, and N1g for comparison.
- Check each result against its own data. Recompute k and the peak from the inputs yourself for at least 10 trials, spread across persons.
- Mark the trials where the agent chose the wrong column, the wrong side or the wrong plate. For those, compute them yourself from the inputs with the same method, and state which ones you replaced.
- Table per person × activity group (gait / bouncy/crouch / squat/chair rise/calf raise/1-leg stand).
- Where does the law hold (peak within ±15 %, k within ±20 %) and where not? Which pattern does the deviation follow: knee flexion angle, two-leg vs one-leg support, dynamics (bouncy), person?
- Compare against A1391/A1398 (frozen: the gait peak median −2.1 %, k_stairs/k_gait 1.337).

## 2. The next The swarm wave
Write `QUEUE_PROPOSAL.txt`. The coordinator reads it before anything is queued. Create 30–60 packets in `results/BT-LG2-*/` that follow from the compilation's deviations. Examples:
- r_q as a function of knee angle (a literature curve) instead of a constant;
- co-contraction in two-leg activities;
- separating the dynamics term for bouncy;
- the ipsilateral GRF share in 2-leg activities.

Each packet must test against measured implant force, with a frozen criterion in the BRIEF.

Rules:
- Only real files in inputs/ (NO symlinks), < 20 MB per packet.
- The BRIEF contains explicit metadata: implant side (JW R, DM R, SC L, PS L), BW (JW4 66.7, DM6 70.0, SC3 78.4, PS5 75.0 kg), and the force column (eTibia norm(Fx,Fy,Fz)·4.4482216152605; JW eKnee = (PM+AM+AL+PL)·4.4482216152605).
- Exclude jw_lungef1.
- No reanalysis of our own results; every packet computes something new against the measured force.
- Copy tasks/NIGHT_PREAMBLE.md into inputs/.

Write `results/CX-LGSYNTH/RESULTS.md`, starting with `# CX-LGSYNTH`, plus results.json and the aggregation script. lane_runner has full permissions in the workspace; `~/projects/bodytwin` is read-only. No verdict words; every deviation is a node that expands.
