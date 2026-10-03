# CX-RQANGLE — does the knee law hold with a LITERATURE curve for the quadriceps arm r_q(θ)? + why does the law overestimate JW's bouncy?

## Background
- A1428, BT-LG2 (The swarm, unverified): the crude r_q(θ)=0.045·[1+0.25cosθ] gives ≥10 % RMSE gain in crouch (DM 22–24 %, SC 14 %) and DM mtpgait2 (12.9 %). The inertia correction is ≈0 everywhere, the ipsilateral GRF share ≈0, co-contraction mixed.
- A1419, CX-LGSYNTH: the law is within ±15 % peak / ±20 % k in only 9/50 non-gait trials. It underestimates squat/chair rise/one-leg stand in k by 27–38 % and overestimates JW bouncy by +31 to +110 % in peak.
- Read: `results/CX-LGSYNTH/aggregate.py` (the shared method, reuse it), `results.json`, `results/BT-LG2-*/` (their scripts), `results/CX-LAW2`, `results/CX-ANKLETERM`.

## Tasks (PREREG.md + sha256 BEFORE scoring; choose the curves from the literature BEFORE looking at the implant force)
1. **Freeze 2–3 published r_q(θ) curves** for the effective quadriceps/patellar-tendon arm about the tibiofemoral contact. Candidates: Kellis & Baltzopoulos, Herzog & Read, Buford, or the patellofemoral ratio curve (van Eijden). Give the source, the table/equation, and the conversion to "effective quadriceps force arm", including the patellar tendon/quadriceps ratio. Also take a literature curve for the hamstring arm if the law's flexion moment is used.
2. **Apply the curves with NO fitting** to ALL 50 BT-LG trials + L1's 108 gait trials (the aggregate.py method). Score against the measured force: RMSE, peak, and k per person × activity group, against the constant law and N1g. Freeze in advance the criterion ≥10 % person-median RMSE gain in non-gait without worse gait.
3. **The JW bouncy overestimate:** neither dynamics nor the plate share explains it (LG2). Test the following, and report which one explains the most with no implant fitting:
   - (a) the sign of the knee moment: a flexion moment gives hamstring force, not quadriceps, but the law uses |M|;
   - (b) that JW's eKnee sum (PM+AM+AL+PL) is a different quantity from eTibia (a regression/calibration factor, per CX-CONTACTDEF);
   - (c) the marker knee centre in bouncy;
   - (d) the ankle term 0.29|M_ankle|/0.045 at high plantarflexion moments.
4. Every deviation is a node that expands. Report which activity or angle range still falls outside after the curves are applied.

Deliver RESULTS.md starting with `# CX-RQANGLE`, results.json, the script, and pytest. It is cheap computation; run locally with 2 threads. Exclude jw_lungef1. lane runner has full permissions in the workspace; `~/projects/bodytwin` is read-only.
