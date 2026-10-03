# REMODEL_crestal: marginal bone loss around implants over time

Node `DENT-BIO-CRESTAL-BONE-REMODELING` with four more nodes in `notes/expansion/bone_remodeling.jsonl`.
- PREREG was frozen 2026-09-24 14:08 (sha256 `d4f8780d…`). Code and predictions were hashed in `PREREG.code_sha256`. The anchors' abstracts were read only after freezing.
- Additions A1–A4 in `AMENDMENTS.md` were made after reading: multiple systems in A1, human or animal in A5, grouping rule and anchors without CI.
- ATOMS: 33/33 GREEN (`atoms_run.txt`).
- Rerun in a new process: bit-identical on two primary cases (`rerun_compare.json`), with a check that the vectorised law = bodytwin's `mechanostat_rate`.
- Code: `cells/physics/crestal_remodeling.py`, `cells/design/design_implant_v2_mbl.py` (DESIGN_implant_v2 is untouched).

## Model
- **L2 bio.** Three conditions determine the crest: biological width, microgap and surface boundary. The bone drops to the deepest of three levels:
  - h_BW − t_muc (h_BW 3,3 mm, Berglundh 1996)
  - s + max(0, d_mg − δ) (d_mg 1,5 mm, Hermann 2000/2001; δ = platform switching along the surface path)
  - s + h_pol (polished neck)
- **L1.** L2 plus a mechanostat on harmonic axisymmetric FE:
  - The FE model: implant, abutment, crestal cortex from BONE_PV (1,65 mm, density 0,70) and cancellous bone 500 MPa.
  - Load: 350 N at 30° (LOAD_spectrum M1, exceeded 36 times per day).
  - Stimulus = ε_implantat / ε_tand, where ε_tand is the strain when a tooth with PDL (F1) carries the same load (site-specific homeostasis, Huiskes).
  - The law: bodytwin `mechanostat_rate` with a lazy zone ±0,75. Overload resorption (Li 2007) above 3,5× was added to it.
- **Verification.** The cantilever against beam theory: the end effect scales as h² laterally (18,7 against 18,4) and is constant axially (1,000). Halving mesh size changes MBL by −0,06 mm.
- All parameters have a source (DOI/PMID in `sources/`) or stand as UNKNOWN and were swept: 37 scenarios.

## Before anchors (pilot, in PREREG)
- **The Frost bands do not work absolutely.** With 50/1500/3000 µε, the FE gives at chewing load p50 6 600 µε in peri-implant cancellous bone. The bone around a natural tooth is at 5 800 µε. The result is collapse to the apex.
  - With the generic bands, even the tooth is therefore not in homeostasis.
  - A site-specific reference is therefore used. The alveolar setpoint level is UNKNOWN (DENT-IF-BI-MECHANOSTAT).
- **The tooth-referenced mechanics is near a cliff.** Collapse already at 1× load in 8 of 37 scenarios: κ 1,5, s 0,35, mesial/distal angle, wall 4 mm, R_t 3, additive bio and the two absolute variants.

## Against held-out anchors
| | anchor (test − control, mm) | L1 full | L2 bio | baseline |
|---|---|---|---|---|
| level, 5 years from loading | Laurell 2011: 0,49 (systems 0,24–0,75) | 0,56 PASS | 0,29 PASS | Albrektsson 2,3 FAIL |
| (a) thin − thick mucosa | Suárez-López 2016: +0,80 [0,42; 1,18] | −0,01 FAIL | +0,27 FAIL (factor 2,9) | 0 |
| (b) platform switching − matched | Atieh 2010: −0,37 [−0,55; −0,20]; Chrcanovic −0,29; Strietzel −0,52 | +0,01 FAIL | **−0,27 PASS (3/3)** | 0 |
| (c) microthreads − without | Niu 2017: −0,09 [−0,18; −0,01] | 0 FAIL | 0 FAIL | 0 |
| (d) overload | Heitz-Mayfield 2004 (dog): no effect; Isidor 1996 (monkey): 5/8 lose osseointegration | collapse: 1 of 2 | 0: 1 of 2 | – |
| (d) bruxism | PMID 36896080: ns (p 0,60) | 0 PASS | 0 PASS | – |
| (e) subcrestal 1 mm − equicrestal, PS | Valles 2018 human −0,18 [−1,31; 0,95] ns; animal −0,45 [−0,66; −0,24] | +0,30 PASS (human) / FAIL (animal) | +0,87 PASS / FAIL | – |

**Verdict according to PREREG: OPEN for both branches.** L1 passes 1 of 5 contrasts and L2 2 of 5. Both pass the level.
The mechanics term is **REFUTED as a driver of the crestal level** (`DENT-BIO-MECH-TERM-CRESTAL`): it resorbs the cortical remainder in the saddle down to around t_c and thereby erases the effects of platform switching and mucosa shown by three meta-analyses.

## Steelman (`steelman.json`)
- PS success in L2 depends partly on choosing the surface-path rule. Prior knowledge of the order of magnitude was reported as A0. With a spherical zone, the effect becomes −0,03.
- L2 predicts a constant PS effect over time, but Chrcanovic reports that it grows.
- The mucosa anchor is matched with h_BW 3,9 from the same source (0,79, post hoc). The anchor thus identifies h_BW and the composition rule.
- Subcrestal PASS rests only on a broad human CI. Animal data and CBL relative to the shoulder (−1,09) go the other way. The model allows no bone above the platform.
- Level from placement: the model gives 1,55–1,93 mm. The patent lane's BQ3 extraction, read after freezing and only descriptive, gives a median of 0,50 mm at 1 year and 0,75 mm at 5 years. Hermann's d_mg 1,5 from dogs is probably too high for modern connections.

## Design coupling (`design_mbl.json`)
- **The bio optimum lies in x_bio and not in the diameter:** maximum platform switching, equicrestal placement and rough neck.
  - Subcrestal placement always increases loss as long as platform switching is less than d_mg.
  - Thin mucosa (< h_BW) cannot be compensated by any implant parameter in the model.
- **Diameter matters only through possible platform switching.** With a minimum abutment radius 1,5 mm (UNKNOWN), δ_max = D/2 − 1,5.
- **DESIGN_implant_v2 with MBL as the fourth ε target:** the R choice changes at 0 % of sites with a matched platform and at 2,7 % with max PS and t_muc 3. The rule chooses robustness (P̂) first.
- **The cost of R's narrow choice** (mean-D 3,6 mm) with max PS: 1,21 mm predicted MBL against 0,94 for rule of thumb B1. That is 0,27 mm more.

## Next stress point
- **The bio branch:** a preregistered dose–response for mismatch (δ 0,2/0,4/0,6) and a time course against held-out RCTs. A human histometric value for d_mg at a conical connection. A rule for bone above the shoulder with PS + subcrestal.
- **The mechanics:** in vivo strain in alveolar bone during chewing, giving the absolute setpoint. Also a formation response that precedes damage resorption, and 3D ridge geometry instead of axisymmetry.

Raw data: `/media/anton/usb-stage/dental_scratch/REMODEL/runs/mech.jsonl` (274 FE-runs). Reproduktion:
`run_remodel.py pred|design` → `analyze_remodel.py` → `compare_anchors.py` → `rerun_compare.py` → `make_evidence.py` → atoms_runner.
