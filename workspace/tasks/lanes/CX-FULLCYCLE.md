# CX-FULLCYCLE — IK/ID and L1 operators through the WHOLE gait cycle, swing included (a prerequisite for FDK, A1782)

Background: A1782 CX-FDK1. L1's operators/IK cover only 60–65 % of Marra's window (heel strike → heel strike) on PS5. Late stance and swing are missing, and that is exactly where Marra et al.'s FDK has lower reported error (A1781).
Read: `results/L1/code/` (l1_prep.py, n12_io.py, n12_model.py: why the frames are cut off, e.g. the GRF threshold, plate coverage or marker gaps), `results/CX-FDK1/results.json` (the missing intervals), `results/CX-ANYBODYWIN` (the window).

## Tasks (PREREG.md + sha256 first)
1. Diagnose why L1 excludes the frames: a GRF threshold, a missing plate, marker gaps, or the accept mask. For each cause, give the share of frames it removes.
2. Build IK + ID + the L1 operator (A, b, F0, cj, C0, Mkx, mx, d) for the WHOLE cycle on PS5_ngait_og_ss1 and PS5_rightturn6. In swing, GRF = 0 and gravity+inertia are what remain. Marker gaps: interpolate ≤ 50 ms, otherwise mark the frame. Use the quadriceps-arm geometry as an option (CX-QUADARM geometry.py).
3. Cross-checks:
   - on the frames L1 already had, the new operator must reproduce L1's to within 1e-9 (the same A/b);
   - LP feasibility over the whole cycle;
   - lo/hi and N1g in swing against the implant force (report it; no criterion).
4. Extend to all 113 L1 gait trials if time allows (batch). Save to external_media (at most 1 GB; the disk is tight, so check `df` first).

Deliver RESULTS.md starting with `# CX-FULLCYCLE`, results.json, and pytest. Run under bigmem.lock with 2 threads. Internal data stays local. lane_runner has full permissions in the workspace; `~/projects/bodytwin` is read-only.
