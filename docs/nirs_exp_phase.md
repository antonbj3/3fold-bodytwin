# Unique-input exponential measurement

Existing reduced-chain profile: flat exponential63.48/64.85ms, accumulation49.88/50.11ms. Temporary Modal observer measures256 unique path-dot/exponential/collection operations; full-chain timing is not inferred.

| Fixture | Batch | In-place | First ms | Second ms |
| --- | ---: | --- | ---: | ---: |
| photon_curved_detector_paths_l4.npz | 1 | False | 14.712 | 13.406 |
| photon_curved_detector_paths_l4.npz | 1 | True | 13.365 | 13.332 |
| photon_curved_detector_paths_l4.npz | 4 | True | 14.419 | 14.712 |
| photon_curved_detector_paths_l4.npz | 16 | True | 14.139 | 14.133 |
| photon_curved_detector_paths_l4.npz | 64 | True | 13.945 | 14.008 |
| photon_detector_paths_l4.npz | 1 | False | 52.562 | 15.225 |
| photon_detector_paths_l4.npz | 1 | True | 17.006 | 15.257 |
| photon_detector_paths_l4.npz | 4 | True | 17.152 | 17.010 |
| photon_detector_paths_l4.npz | 16 | True | 16.636 | 16.218 |
| photon_detector_paths_l4.npz | 64 | True | 21.311 | 22.862 |

All20 measured output collections have zero bit mismatches. No candidate: on the flat fixture in-place15.257/17.006ms fails the strict paired speed criterion against15.225/52.562ms. Batch16 needs1054464bytes per stacked working array, batch64 needs4217856; output retention and stacking costs are included in these measurements, not peak memory. Original control calls remain independent. Evidence: reports/nirs_exp_phase.json.
