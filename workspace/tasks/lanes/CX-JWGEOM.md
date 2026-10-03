# CX-JWGEOM — make each Grand Challenge person's knee geometry individual → does the physics bound become correct for JW?

Lead from CX-BEATN1G2 (`results/CX-BEATN1G2/RESULTS.md`): the certified physics bound corrects N1g correctly in SC (the error on the corrected frames falls 1.91 → 0.17 BW), but in JW it gets it badly wrong (0.49 → 3.32 BW). The model's lowest possible contact force is too high for JW. JW's COMAK/model is a DM.osim transplant (A326), so JW has the wrong geometry. Hypothesis: with individual knee geometry (patellar tendon lever arm, muscle attachments/paths, knee centre, contact location), JW's F_min/F_max becomes correct, and the certified correction of N1g holds in ≥ 3/4 persons.

## Build on (read first)
- `results/CX-BEATN1G2` (C1 code, bounds_*.npz, margin/decision_cert wiring) and `results/CX-BEATN1G`, `results/CX-KNEENET`, `results/CX-EARLYSTANCE`, `results/L1` (code/, prep/, tlem_arch.json).
- The geometry manager `results/N7c/geomgr` (registration op_register_surface, attachment transfer, certificate) and `results/CX-FIELDSHARE/femur_edit.py`.
- `~/projects/bodytwin/scripts/msk`: grep `knee_jw_transplant`, grand_challenge, comak, tibia, patella. The graph's KNEE-CELL node in `~/projects/bodytwin/data/MECHANISM_ANCHOR_GRAPH.json`.
- The Grand Challenge competition zips `/media/anton/sdc1-tmp/*Competition-latest.zip`: find CT/MR, implant CAD, and registered bone geometry per person (JW, DM, SC, PS). Extract only what is needed to `/media/anton/sdc1-tmp/bodytwin/CX-JWGEOM/`.

## Tasks (PREREG.md + PREREG.sha256 before the first run)
1. Inventory which individual geometry exists per person: CT/implant/landmarks. Say plainly what is missing.
2. Build individual knee geometry for JW (and for the other three, same method) with geomgr: register the TLEM femur/tibia/patella to the person's bone/implant geometry, transfer attachments and the knee centre, and compute the patellar tendon lever arm per knee angle.
3. Rebuild L1's A and b (moment arms, contact projection) with the individual geometry. Rerun H1/C1 (the feasible set + decision_cert) exactly as in CX-BEATN1G2.
4. Criteria:
   - JW's error on corrected frames < N1g's on the same frames;
   - C1 median ≥ 5 % below N1g in ≥ 3/4 LOPO rotations and not worse in early stance.
   Counter-test: another person's geometry (swap JW↔DM↔SC) must lose the gain.
5. Report how much of the JW error the geometry explains (the lever arm alone, attachments alone, the knee centre alone).

## Resources
- Locally: nice, 2 threads, ≤ 60 s per test.
- Write only in `results/CX-JWGEOM/`.
- `RESULTS.md` starting with `# CX-JWGEOM`, plus results.json and code with pytest.
