# CX-LUNGESYNC — BodyTwin's delivery to Field's in vivo validation F-8 (U388): knee angle, time sync, and the N1g null model for JW's lunge jw_lungef1

Consumer: the Field lane (anton-12), U388. It computes contact force between the femoral component and the insert from fluoroscopy poses and compares against the measured eTibia force. Joint priority F-8 in `external_research_path` §14.

## Data
- The Grand Challenge 4th competition: `external_media` → "Synchronized Motion Data.zip" → jw_lungef1:
  - fluoro kinematics segments 1–5 (84 frames, 30 Hz);
  - markers 120 Hz;
  - GRF/force plates;
  - eTibia.
- The frame table and times: `results/CX-FLUOROLINK/results.json` → JW4.
- Extract only what is needed to `external_media`.

## Deliverables in `results/CX-LUNGESYNC/`
1. `lunge_frames.csv`, one row per fluoro frame (84):
   - time (s, the synced time base);
   - knee flexion angle from the fluoro pose (state the convention);
   - the same from the markers, as a check (difference in °);
   - GRF vector + |GRF| (N, and in BW with JW's body weight — state the source);
   - measured eTibia total force + medial/lateral, if they are in the archive (otherwise UNKNOWN; do not open the facit for Field before they freeze their PREREG — deliver eTibia in a separate file `facit_etibia.csv` that they read AFTER PREREG).
2. The N1g null model per frame: k(activity)·|GRF| with k from the existing LOPO fits (results/L1/scores.json, BT-B122/BT-B128, N39). The lunge is probably not among N1g's activity groups: state which k is used (the nearest group, or the overall N1 k) and give both N1 and N1g.
3. `SYNC.md` (≤ 30 lines): time base, the offset between fluoro/markers/GRF/eTibia (measured, e.g. cross-correlation of the load), known gaps (152 missing frame numbers between segments), units, and frames.
4. `RESULTS.md` starting with `# CX-LUNGESYNC`: what was delivered, checks (fluoro vs marker angle), and N1/N1g error against eTibia **after** Field has frozen its PREREG. You may compute N1g's error against eTibia immediately, since it is our own null model; say so.

Resources: local, nice, 2 threads, ≤ 60 s per step. lane runner has full permissions in the workspace; the zips are read-only. When done, send a short note to Field: write `results/CX-LUNGESYNC/READY_FOR_FIELD.md` with the paths.
