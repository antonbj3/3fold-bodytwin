# LANE_SURGICAL_SYNTHESIS

The surgical hunt of the night (1/10) has yielded four lanes of ports, findings, and measurement specs. This lane makes them **a coherent, runnable and tested BodyTwin capability** — what Anton asked for: high resolution where the scalpel opens the body. Results folder `results/LANE_SURGICAL_SYNTHESIS/`.

## Background (read everything, build on it, don't redo)

- `results/LANE_SURGICAL_INCISION/` — PORTS_R4_FINAL.json, MEASUREMENT_SPEC_R4.md: needle toughness tool specific (same tool over speeds ≤ 2,8 % error, new tools 22–55 %), held needle geometry 11 %, prospective Γ0 ≲ 29 J/m², stored bias model.
- `results/LANE_SURGICAL_BINDINGS/` — PORTS_R4_FINAL.json, MEASUREMENT_SPEC_R4.md: Γ_cut scenario 150–380 J/m², interfacial energies, micro budget ≤ 7 kJ/m².
- `results/LANE_SKIN_TOUGHNESS_GAP/` — MECHANISM_TABLE_R2_FINAL.json, PORTS_R2_FINAL_V2.json, MEASUREMENT_SPEC_R2.md: poroelasticity ≤ 0,11 and ripple ≤ 0,14 kJ/m² excluded; mm bridging 19–25 kJ/m² right size; mode split unexplained.
- `results/LANE_SURGICAL_RESPONSE/` — RESPONSE_PORTS_R3.json and later: U/I/M collagen inventories give Levenson strength RMSE 8,31 pp without strength fit; HP42 −52 %; round 4 is in progress (read its results when available).
- Swarm's SURG result (BT-FW48-AUTO-* with SURG_* in source_keys under external_mount): e.g. bleeding is borne by radius distribution (number-matched control underestimates 25×), wall shear in platelet plug.
- Private code (read-only) `~/projects/bodytwin/scripts/msk/` and HX Q033/Q036/Q049.

## Delivery

1. **An executable chain** `surgical_chain/` (Python package in results folder): tools + cut path + tissue stack → damage zone and gap → bleeding/hemostasis → oxygen edge → healing inventories → strength over time. Each step takes a documented port (device, uncertainty, source, status: MEASURED / SYNTETISKT / UNKNOWN) and no number hides that it is unknown.
2. **Regression tests** that lock the verified numbers of the night (those in the NIGHT_LOG and the laners' VERIFICATION files) so that future changes don't silently break them; analytical limit cases per step.
3. **A unified measurement specification** that merges the four lanes' MEASUREMENT_SPEC into a prioritized list: which single measurement unlocks the most of the chain (value per cost).
4. **Graph binding**: one `./graph dispatch` + `./graph feedback` per target (BT-CTX-SURG-INCISION, -COLLAGEN, -HEMOSTASIS, -HEALING), kind review/define, PENDING_INDEPENDENT_REVIEW, with negative results preserved.
5. RESULTS.md for a reader: ability → what holds against the world (with file) → what is unknown → next measurement. No process narration.

No internal data. No change in canonical BodyTwin code or source graph; all in the results folder. Strongest control for the chain as a whole: separately calibrated curves per observable — honestly report where the chain predicts something they can't and where it doesn't.
