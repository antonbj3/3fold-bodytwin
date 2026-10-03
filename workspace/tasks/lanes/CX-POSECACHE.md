# CX-POSECACHE — speed for real: a tolerance-based pose cache + warm start for the contact band (GAUGE2D's next node)

Background: A1446 CX-GAUGE2D. The gauge is exactly 2-D, but the exact null-space formulation is 5–24× SLOWER (dense N); A changes every frame. Field U478 attacks the problem with a projection/Fourier–Motzkin (read `the public staging tree/3fold-motion-engine/_private/romi_collab/build/U478/` if it exists; read-only).
Read `bodytwin_core/contact_band_fast.py`, `results/CX-EPSBAND/eps_band.py`, `results/CX-GAUGE2D/`.

## Tasks (PREREG.md + sha256 first)
1. **Warm start:** the LP/QP for frame t starts from frame t−1's basis/solution (HiGHS basis reuse; Clarabel warm start where supported). Measure frames/s against a cold start.
2. **Pose cache:** cluster frames by pose (knee/hip/ankle angle, GRF direction). Compute the band exactly for the cluster representatives. For the other frames, use a first-order correction (sensitivity of the LP optimum: the dual values give d(lo)/d(A,b)). Always report the error against the exact result per frame, and a certified fallback when the error bound exceeds 0.01 BW.
3. **Batching:** run many frames in parallel via vectorised numpy/jax where possible.
4. **Target:** ≥ 20× on the whole L1 mask (11,736 frames) with max band error ≤ 0.01 BW, verified against the frozen bands in CX-EPSBAND/INVERSEOC.
Deliver `contact_band_batch(trials, eps, tol)` in bodytwin_core + pytest, RESULTS.md starting with `# CX-POSECACHE`, and results.json with the timing. Run under bigmem.lock with 2 threads. lane runner has full permissions in the workspace; `~/projects/bodytwin` is read-only.
