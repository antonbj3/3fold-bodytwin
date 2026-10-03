# LANE_SURGICAL_SYNTHESIS

Tonight's surgical search (1/10) has produced four lanes with ports, findings and measurement specifications. This lane turns them into **one coherent, executable and tested BodyTwin capability** — what Anton asked for: high resolution where the scalpel opens the body. Results directory `results/LANE_SURGICAL_SYNTHESIS/`.

## Supporting material (read everything, build on it, do not repeat it)

- `results/LANE_SURGICAL_INCISION/` — PORTS_R4_FINAL.json, MEASUREMENT_SPEC_R4.md: tool-specific needle toughness (same tool across speeds ≤ 2,8 % error, new tools 22–55 %), held-out needle geometry 11 %, prospective Γ0 ≲ 29 J/m², stored prestrain model.
- `results/LANE_SURGICAL_BINDINGS/` — PORTS_R4_FINAL.json, MEASUREMENT_SPEC_R4.md: Γ_cut scenario 150–380 J/m², interface energies, microbudget ≤ 7 kJ/m².
- `results/LANE_SKIN_TOUGHNESS_GAP/` — MECHANISM_TABLE_R2_FINAL.json, PORTS_R2_FINAL_V2.json, MEASUREMENT_SPEC_R2.md: poroelasticity ≤ 0,11 and crimp ≤ 0,14 kJ/m² excluded; mm bridging 19–25 kJ/m² has the right magnitude; mode split unexplained.
- `results/LANE_SURGICAL_RESPONSE/` — RESPONSE_PORTS_R3.json and later: U/I/M collagen inventories give Levenson strength RMSE 8,31 pp without strength fitting; HP42 −52 %; round 4 underway (read its results when available).
- Warmness SURG-resultat (BT-FW48-AUTO-* med SURG_* i source_keys under /mnt/games-240/research/bunny48_20260926/bodytwin/): e.g. bleeding is carried by radius distribution (antal-matchad kontroll underskattar 25×), wall shear in platelet plug.
- Private code (read-only) `~/projects/bodytwin/scripts/msk/` and HX Q033/Q036/Q049.

## Leverans

1. **An executable chain** `surgical_chain/` (Python package in the results directory): tool + incision path + tissue stack → damage zone and gap → bleeding/hemostasis → oxygen edge → healing inventories → strength over time. Each step takes a documented port (unit, uncertainty, source, status: MEASURED / SYNTHETIC / UNKNOWN) and no number hides its unknown status.
2. **Regression tests** locking tonight's verified numbers (those in NIGHT_LOG and the lanes' VERIFICATION files) so future changes do not silently break them; analytical limiting cases per step.
3. **A consolidated measurement specification** merging the four lanes' MEASUREMENT_SPEC into a prioritized list: which single measurement unlocks the most of the chain (value per cost).
4. **Graph binding**: one `./graph dispatch` + `./graph feedback` per goal (BT-CTX-SURG-INCISION, -COLLAGEN, -HEMOSTASIS, -HEALING), kind review/define, PENDING_INDEPENDENT_REVIEW, with negative results preserved.
5. RESULTS.md for a reader: capability → what holds against the world (with file) → what is unknown → next measurement. No process narration.

No internal data. No changes to canonical BodyTwin code or the source graph; everything in the results directory. Strongest control for the chain as a whole: separately calibrated curves per observable — report honestly where the chain predicts something they cannot and where it does not.
