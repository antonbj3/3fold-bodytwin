# CX-GAUGE2D — speed: the contact band per frame as a 2-D problem (particular solution + the contact-relevant gauge directions) instead of a 166-dim LP/QP

Background:
- Field U442/U447/U459 (A1417, A1425, A1440): only ~2 of 75 gauge directions touch the knee contact; hamstrings carry Δf.
- U467 (A1441): ℓᵀz is determined ⇔ ℓ ⊥ ker[J A].
- CX-EPSBAND (A1435): eps_band.py solves a QP/SOCP per frame (Clarabel).
- U456 (A1433): exact rational certification, ~0.8 s/frame.

## Tasks (PREREG.md + sha256 first)
1. Per frame from L1 (A, b, F0, cj, C0, med/lat):
   - a particular solution f_p (min-E);
   - a basis N of the null space of A;
   - the contact-relevant subspace = the directions in N that change cj·f (and med/lat). Report the dimension per frame, and whether it is 1–2 as Field says, or larger in L1 (166 elements).
2. **A reduced solver:**
   - the LP band = min/max of cj·f over {f_p + N z : 0 ≤ f ≤ F0, med/lat ≥ 0}. Project onto the contact-relevant coordinates: the non-contact coordinates only determine feasibility, so handle them exactly (Fourier–Motzkin/projection, or an LP in reduced dimension);
   - the ε-band = the same, plus E(f_p + N z) ≤ (1+ε)E_min, where E is quadratic in z.
   Verify EXACTLY against CX-EPSBAND/INVERSEOC on all 11,736 frames (max difference in BW; target ≤ 1e-6 BW).
3. **Timing:** frames/second for the reduced solver against Clarabel/HiGHS, single-threaded, batched numpy, and optionally GPU (warp/jax if available; not required). Target ≥ 100× faster with identical bands.
   - Option: the offline part (N, projection) is precomputed per trial pose if A varies slowly; report reuse across frames.
4. **Deliver:** `contact_band_fast(A, b, F0, cj, C0, med, lat, eps)` in bodytwin_core-compatible form + pytest (identical to the reference on a synthetic system and on 3 real trials).

Deliver RESULTS.md starting with `# CX-GAUGE2D`, results.json, and the benchmark script. Run locally with 2 threads under bigmem.lock for large runs. Internal data stays local. lane runner has full permissions in the workspace; `~/projects/bodytwin` is read-only.
