# CX-IMPLREG — register the implant mesh ↔ femur/tibia marker segments for PS5 (and JW4, DM6, SC3): a prerequisite for FDK and the JW repair

Background:
- A1782 CX-FDK1 and A1703 CX-JWPATGEOM: there is no verified rigid transform from the implant STL (femoral component, tibial insert) to the marker segments in the gait trials.
- CX-CT2MARKER found 0 accepted PS5 transforms (a 28.2 mm residual).
- The GC archive contains fluoroscopy-derived tibiofemoral kinematics, CT bones and static trials.
Read: `results/CX-CT2MARKER`, `results/CX-FLUOROLINK`, `results/CX-JWGEOM` (archive_inventory.json, bones.py), and the GC documentation (/media/anton/sdc1-tmp/bodytwin/CX-JWGEOM/docs/).

## Tasks (PREREG.md + sha256 first; the implant force is NEVER used for the registration)
1. Inventory per person what exists:
   - CT bone models with the implant;
   - fluoro poses (which trials/frames, and in which frame);
   - static marker trials with anatomical landmarks (epicondyles, malleoli);
   - the implant's CAD/STL frame.
2. Choose a registration chain that has independent support. Examples:
   - (a) static trial landmarks ↔ the CT bone's landmarks (a rigid fit plus a residual);
   - (b) fluoro pose ↔ markers in a shared trial/frame;
   - (c) a combination of these.
   Report the residual in mm and the uncertainty (bootstrap over landmarks). Freeze an acceptance threshold (≤ 5 mm residual, ≤ 3° rotation).
3. Validate the registration against something it was NOT fitted on:
   - the fluoro kinematics in another frame;
   - OR the tibial insert's position against the knee centre from SCoRE (CX-KNEECENTER);
   - OR the implant's joint line against the marker-based flexion axis.
4. Deliver `implant_to_segment(person)` with a transform + uncertainty, and pytest.
If it is impossible with the archive, say EXACTLY which measurement is missing and how it would be acquired.

Deliver RESULTS.md starting with `# CX-IMPLREG` and results.json. 2 threads, at most 300 MB of intermediate files. Internal data stays local. lane runner has full permissions in the workspace; `~/projects/bodytwin` is read-only.
