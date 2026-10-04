# ORTHO_rate: does PDL pressure from the reference FE order orthodontic movement rate?

Cell `DENT-CELL-ORTHO-RATE-01`, node `DENT-IF-TPB-REMODELING`. PREREG was frozen before production runs and before comparison with the anchors
(`PREREG.md`, sha256 `37a60742…de27`; code and anchors in `PREREG.code_sha256`). ATOMS: `atoms_runner --execute evidence.json` →
**ALLOW/GREEN 28/28** (`atoms_run.txt`). Node rows: `notes/expansion/ortho_remodeling.jsonl`.

## Answer
**No.** The FE-pressure model plus the Schwarz/Reitan law passes 5 of 8 held contrasts. Constant rate and uniform stress F/A pass
5 of 8 as well, and pure force scaling 1 of 8. H1 fails according to the preregistered falsifier (5 ≤ 5), and the FE distribution adds
nothing beyond F/A (H3 fails).

The reason is that even the lowest clinical force, 4 kPa (0,28 N on the canine), gives a median PDL pressure on the compression side of
3,7 kPa. Then 69 % of the compression side is above capillary pressure, 2,2 kPa. Across all clinical forces, the model therefore becomes practically
constant, with ratios between 0,96 and 1,00. Data, however, rise weakly: ×1,43 from 4 to 13 kPa (Nickel), ×1,5 from 50 to 200 cN
(Owman-Moll) and significantly from 50 to 300 g (Yee).

## Setup
- **FE.** A227's surface-refined PDL S2L4, E 0,735 MPa, ν 0,45.
  - Teeth: OFJ P12/P14, 4 canines (distal force) and 4 first premolars (buccal force).
  - Two unit loads per tooth: 1 N at the bracket and 1 N·mm force couple. Tipping = force only. Translation = force plus M/F = r_b, where r_b gives zero tipping.
  - Pressure p = −n·σ·n per surface column.
  - Solver CalculiX PARDISO, locally via `tasks/heavy_run.sh`, because the Modal account has reached its spending limit.
  - ITERATIVE CHOLESKY was rejected: it gave 0,9 % error in pressure under force load and 5 % under force couple (`solver_check/`).
- **Law M.** φ = p/p_cap up to p_cap = 2,2 kPa (Schwarz 20–25 g/cm², via Asiry 2018). Above p_cap, φ = 1 − t_lag/T,
  with t_lag = 21 d. The tooth's rate is obtained by fitting a planar rigid-body motion to φ on the compression side and read at the bracket.
- **Baselines.** B0 constant, B1 ∝ F, B2 = same law on uniform stress F/A_proj (A/2 at tipping).

## Against held anchors (median over 4 teeth; ratio v_high/v_low)
| | source | band (95 %) | obs | M | B0 | B1 | B2 |
|---|---|---|---|---|---|---|---|
| N1 4→13 kPa, translation | Nickel 2014 | 1,02–2,01 | 1,43 | 1,00 F | 1,00 F | 3,25 F | 1,00 F |
| N2 13→26 | ” | 0,97–1,61 | 1,25 | 1,00 P | P | 2,00 F | P |
| N3 26→52 | ” | 0,87–1,34 | 1,08 | 1,00 P | P | 2,00 F | P |
| N4 52→78 | ” | 0,93–1,37 | 1,13 | 1,00 P | P | 1,50 F | P |
| O1 50→100 cN, premolar tipping | Owman-Moll 1996a | 0,77–1,30 | ns | 0,97 P | P | 2,00 F | P |
| O2 50→200 cN | Owman-Moll 1996b | 1,07–1,93 | 1,50 | 0,96 F | F | 4,00 F | F |
| Y1 50→300 g, rabbit sliding mechanics | Yee 2009 | > 1,10 (sign only) | sign. | 0,98 F | F | 6,00 P | F |
| E1 50→150 g, month 1 | El-Salam 2026 (table 3, read after freezing) | 0,85–1,87 | 1,26 | 1,00 P | P | 3,00 F | P |
| **Total** | | | | **5** | **5** | **1** | **5** |

- **Sensitivity.** M gets 5/8 at p_cap 1,96 and 2,45, t_lag 20, measurement point in the crown's center or tip, 6-DOF fitting and A = 44,7 mm².
  At t_lag 30, M gets 4/8. M never ends up above the baselines.
- **Material.** E 0,694 and 0,782 change the ratios ≤ 0,002 in ln, because pressure under force control does not depend on E.
  ν 0,49 changes ≤ 0,052 in ln (one tooth).
- **H2, direction (secondary).** Tipping/translation at the bracket becomes 1,63–1,69 (M) and 1,72 (FE-lin), against Nakano 2014 (rat) ≥ 2.
  H2 fails.
- **Secondary Iwasaki series (same cohort).** M and B0 give log-RMSE 0,24–0,76. B1 gives 0,80–2,2.

## Stress point: center of resistance
- For the canine under distal force, CR in the FE is only **2,8–4,2 mm apical to the bracket**, thus near the alveolar crest
  (1,8 mm at ν 0,49).
- For the premolar under buccal force, CR is 6,8–8,8 mm apical to the bracket.
- Clinical canine translation requires M/F 9–13 mm (Iwasaki 2000). The almost incompressible linear PDL therefore places CR
  too cervically, or the clinical M/F also includes the appliance's play.
- M/F for translation is a published measure that can identify ν, the same parameter F1, F2 and VAL F could not determine
  (node `DENT-PHYS-PDL-CR-TRANSLATION-MF`).

## Post hoc (not preregistered)
- The fit v ∝ F^b gives b = 0,26 on Nickel 2014 only.
- It passes all 4 independent held contrasts: O1 1,20, O2 1,43, Y1 1,59 and E1 1,33.
- It lowers log-RMSE in Iwasaki 2017 from 0,66 to 0,10 (`posthoc_powerlaw.json`).
- Data therefore suggest a weak force dependence without threshold in the range 0,3–5 N, not a capillary pressure threshold.
- The next step is to freeze b = 0,26 and test against new anchors outside the Iwasaki cohort (node `DENT-BIO-ORTHO-RATE-FORCE-POWERLAW`, OPEN).

## Rerun and gaps
- **Rerun.** P12:5 S2L4 was rerun in a new process and gave a bit-identical npz. The analysis was rerun in a new process and gave the same
  sha256 for `ortho_summary.json`.
- **Resolution gap.** S3L4 (1:9) could not be run:
  - systemd-oomd killed two attempts under memory pressure from other sessions,
  - the Modal account has reached its spending limit,
  - OVH had only 8 GB free.

  The replacement S1L4 → S2L4 changes ratios ≤ 0,029 in ln and CR by 0,02 mm. N2 changes verdict on one individual tooth at the band edge.
  The median over 4 teeth does not change direction. M remains flat and the H1 verdict rests on that.
- **Other caveats.**
  - Iwasaki stress is mapped via our own A_proj.
  - Yee has only sign; the floor 1,10 is my choice.
  - O1's and O2's bands are derived from ranges (SD ≈ range/4).
  - The direction anchors come from rat.
- **Bookkeeping.** `DENT-BIO-ORTHO-RATE-PDL-PRESSURE` is set to REFUTED (for the Schwarz law with FE in the clinical force range).
  `DENT-IF-TPB-REMODELING` remains OPEN and now has falsifiers, tolerances and anchors.
- **Files.** FE output is in `fe_out/` (npz, 13 MB). OFJ data is in `external_media` (78 MB).
  Transient solver files were in `/dev/shm` and are deleted.
