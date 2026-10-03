BT-HX-Q049

## Result
First executable model is a 1D layered poroelastic consolidation model: interfacial film, superficial zone, bulk and calcified/barrier. `model.py` includes Darcy flow, pressure storage, series compliance, boundary conditions and unit control; `test_model.py` passes five tests.

**Verified anchors:** Park et al. 2003, DOI `10.1016/S0021-9290(03)00231-8`, Table 1: bovine 94 ± 4% at the surface and 71 ± 8% in deep zone. Krishnan et al. 2003, DOI `10.1115/1.1610018`, Tables 1–2: 0.93–0.98 `Wp/W` upon contact. These are references, not calibration data for the model.

## Key figures (`results.json`)
- Full: `psi_s(1 s)=0.05322`, `psi_s(10 s)=0.001017`; surface tension `1942.6 Pa` and `37.10 Pa` respectively.
- Without movie: `psi_s(1 s)=0.000522`; film-delta `0.05270`.
- With the surface merged in bulk: `psi_s(1 s)=0.04930`; bone flow at 10 s differs relatively with `0.4345`.
- Without barrier: leg flow `4.5071e-9 m/s` against full `2.0826e-7 m/s`; barrier effect `0.9784`.
- Analytical error `9.97e-16`; mass error `2.37e-14`; 16/32 resolution diff `1.74e-4` for `psi_s` and `3.66%` for cumulative leg flux diff.

## Assessment
Film and barrier are transport relevant in this model scenario; superficial zone affects bone flow. The robust game is signature stable for all three ±50% sweeps. The core test is **FAIL** only for the plausibility anchor: the full model's 1 s-value is below the frozen 0.90–1.00 range. No measurement data has been adjusted; film and barrier parameters are assumptions. The model therefore does not alone answer the question of biological necessity.

## Next step
Measure film conditional conductance and bone boundary transport in matched contact, supplement with lateral contact field and cyclic compression, and rerun preregistered ablations.
