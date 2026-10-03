# CX-JWPATGEOM — JW's patellar tendon geometry (attachment points, neutral length, knee axis at 40–55°) against INDEPENDENT image/implant geometry

Background:
- A1450/A1452: JW's L1 quadriceps arm is 3–6 mm at 40–60° (the literature has 35–50 mm). The ρ ratio does not explain it; the tendon-length constraint is off by 4.6–8.4 mm.
- A1702 (pooled): swapping only the quadriceps arm lowers JW gait from 48 % to 21 % under lo.
- The morning summary's next node 1: measure JW's patellar attachments, neutral tendon length and knee axis against independent geometry.
Read:
- `results/CX-JWGEOM/` (bones.py, bounds_attachments_JW.npz, bounds_axis_JW.npz, archive_inventory.json; the GC's JW implant/CT/fluoro documentation, see its RESULTS/docs);
- `results/L1/code/n12_model.py` (how L1 places JW's patellar tendon/tibial tuberosity; transplanted from DM?);
- `results/CX-QUADARM/mechanism.json`;
- `results/CX-PATELLAJW/` (branch_probe).

## Tasks (PREREG.md + sha256 first; no implant force for choices)
1. Locate the tibial tuberosity/patellar tendon insertion, the patella's thickness/length, and the knee axis from JW's own image/implant geometry. Use the implant CAD/fluoro/CT in the GC archive if it exists; otherwise the bounds from CX-JWGEOM. Compare against L1's placement, in mm, together with the uncertainty.
2. Rebuild JW's patella/quadriceps geometry with the independent values. Compute the effective quadriceps arm over 0–90° against the literature, and the LP share under lo on JW's gait mask + JW4 non-gait.
   Frozen criterion: JW gait under lo ≤ 25 % with NO curve from another model (only JW's own geometry).
3. If the archive lacks independent geometry: report exactly what is missing, and which measurement (a lateral X-ray, fluoro frame, CT) would close it.

Deliver RESULTS.md starting with `# CX-JWPATGEOM`, results.json and pytest. Run under bigmem.lock with 2 threads. Internal data stays local; no jw_lungef1. lane runner has full permissions in the workspace; `~/projects/bodytwin` is read-only.
