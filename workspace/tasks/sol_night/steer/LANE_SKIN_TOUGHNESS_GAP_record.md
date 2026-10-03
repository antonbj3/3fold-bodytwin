# Steering LANE_SKIN_TOUGHNESS_GAP — round 2 (the coordinator, 1/10 03:00)

Review of round 1: good narrowing, with numbers (verified MECHANISM_TABLE_R1_FINAL.json, r1/axis_corrected_results.json):
- Poroelastic dissipation ≤ 0,11 kJ/m² also generous → ruled out as main mechanism.
- Single ripple correction ≤ 0,14 kJ/m² → excluded.
- **Contiguous mm bridges** (L ≈ 5 mm, strain 0,5) gives 19–25 kJ/m² — right order of magnitude, only candidate left. But the orientation ratio trans/long becomes 0,78 to measured 1,30 (40 % error) under affine fiber orientation, and mode III is unknown. Good that the first orientation meeting was canceled after the shoulder straightening.

## Obstacles and changed operation

The brew needs two independent anchors, no more parameters:
1. **Bridge zone length and opening at bridge failure from published images** — Yang 2015 (Nat Commun 6:6649) and similar in situ studies (SEM/SAXS/video of tearing) show fibers spanning the crack. Measure or quote bridge zone length and crack width at failure and compare to required (triangular opening 2,6–4,4 mm for mode I). Freeze before comparison.
2. **Nonaffine reorientation**: Yang shows that fibers rotate against the loading direction before bearing. Replace fixed affine orientation with a reorientation law (fibers rotate against the tensile direction with a measured angular change) and predict trans/long-kvoten **without** adapting to 1,30. Then predict mode III out of the same zone (shear loading gives different bridge opening) and compare with 20,6 kJ/m².

Falsifiers: if the bryggzon/opening from the images is far from what is required, or if the reorientation cannot give the correct sign and magnitude of the anisotropy, the bridging explanation fails — then the gap is truly unexplained and it should be clear.

Provide updated mechanism table and PORTS (bridging zone length, G_bridge per orientation and mode) for RESPONSE (wound strength also carried by bridging fibers).
