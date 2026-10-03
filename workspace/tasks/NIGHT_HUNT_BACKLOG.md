# Breakthrough hunt night 24/9 — ready tasks to start when resources/instructions arrive

Each mission: read `tasks/NIGHT_PREAMBLE.md` first. Shoulder number = "Based on:" = start point that MUST be read. Model selection is set when the lane runner-lanner instructions for reserve_worker/space swarm arrived. Already running: N1 (axel 1), N2 (3), N3 (6), AU6, LV3.

| Id | Axis | Question (falsifiable) | Builds on | Ground truth / counter-test |
|---|---|---|---|---|
| N4 | 2 | How large are the moment arm bands from TLEM attachments under atlas uncertainty, and which muscles fall outside the literature? | `scripts/msk/moment_arm_lit_compare.py`, X1b, H2/H2b/AT1 | literature moment arms; null model gait2392/Rajagopal (lane file `tasks/lanes/BT-N4.md`) |
| N5 | 5 | Does an expanded certificate reject gonial +60 and adversarial deformations without false rejections? | geometry_validity_v1 (clone b95a8dc), D1/X1/X1b, C1b | false approvals/rejections before/after (lane file `BT-N5.md`) |
| N9 | 3+10 | Error-bound-driven resolution (LV3 + U359 stopping rule): how much computation is saved with unchanged decision in the entire chain geometry → moment arm → force? | LV3 results, field U357/U359, FIELD_RESULTS H1 (Hadamard shape derivative) | full-resolution reference; decision flip = error |
| N10 | 1+7 | the collaborator's own idea: conditional likelihood over population — which minimum measurement set (X-ray, landmarks, star movement, IMU) reaches 5 mm HJC/attachment? Information value per measurement | IM1/IM3, P1, K5 VOI pattern, FIELD_RESULTS H2 (identifiability) | LOSO; cost per measurement declared in PREREG |
| N11 | 4 | The whole collaborator's NATO chain with error budget: activity → tissue load → fatigue → risk, which link dominates and which measurement is bought first | F1, C2/C3/C4, K5, `scripts/fatigue_life.py`, `docs/MECHANISM_FEMORAL_NECK_STRESS.md`, `bone_stress.py` | null model; variance per link sums to total |
| N12 | 6 | Knee cell with the day’s geometry as a new independent leg: more Grand Challenge subjects, COMAK on DM/JW/SC/PS | KNEE-CELL, `scripts/msk/knee_jw_transplant/`, `comak_*`, X1b T4 (−346 N bias), MAP2-A | measured eTibia force; null model 104 %BW |
| N13 | 5 | Differentiable chain: gradient of joint load with respect to shape parameters (shape derivative) versus finite difference — enables optimization/calibration through shape derivatives | FIELD_RESULTS H1, rnea_warp, N2’s batch solver | FD reference; sensitivity ranking against N3’s error decomposition |
| N15 | 2 | Attachment regions instead of points (region fields on bone surface) — does this reduce moment arm error against point attachments at the same data? | H2/AT1/JS1, LHDL, field region representation (U288 mixedstate) | literature attachment surfaces; counter-test: point = region centroid |
| A-N* | all | Independent review of each N result before "holds" | respective results/<id> | recalculation of key numbers |

Order when resources exist: N7 (broad coverage of all geometry manager interpretations, possibly split among multiple agents per interpretation) and N6 first, then N14 (parity), N8, N10, N12, N9, N13, N11, N15; N4/N5 as lanes. Ongoing review.
