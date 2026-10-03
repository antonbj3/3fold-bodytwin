# PREREG — BT-HX-Q146 (frozen before first run)

Nasal route with swallowed and systemic fraction. First executable mechanistic model of
the core mechanism: deposition → (mucosal uptake | mucus clearance to GI | external loss)
→ local metabolism → plasma, with the same F_sw [mol/s] fed into the GI model.

Freeze date: written BEFORE the first run of `model.py`.
Run: `python3 model.py` (writes `results.json`), `python3 test_model.py`.

---

## 1. Hypotes (H1)

A regional, geometrically resolved nasal model where **absorption rate and
mucus clearance compete locally** (dimensionless number Pe = k_epi·a / v_clear) and where
the swallowed fraction enters a GI model with published first-pass, predicts
both (a) systemic exposure (AUC/Cmax/F) and (b) locally remaining amount
from **one** parameter set, without site-specific fitting.

## 2. Reference answers (outputs according to the question)

| Storhet | Enhet |
|---|---|
| N_nose, regional remaining amount | mol |
| F_swallow, swallowed fraction of the dose | 1 |
| C_free, maximum free local concentration in a mucosal strip | mol/m³ |
| AUC_plasma | mol·s/m³ |
| F_abs, absolut bioavailableitet (iv-refererad) | 1 |

## 3. Reference values (published, looked up with web search before the run)

**R1 (PRIMARY).** Johansson C-J, Olsson P, Bende M, Carlsson T, Gunnarsson PO.
"Absolute bioavailability of nicotine applied to different nasal regions."
*Eur J Clin Pharmacol* 1991;40(6):581–588. DOI **10.1007/BF00314989**.
PMID 1815971. Verified: publication page (link.springer.com + doi.org) read 2026-09-25.
- 8 healthy volunteers; 1 mg nicotine intranasally (dropper on septum or conchae,
  and spray pump 100 µL/activation) against 0,7 mg nicotine iv as a 2 min infusion.
- **F_abs = 60–75 % for the nasal applications** (means: septum 76 %, conchae 64 %,
  spray 58 %; n = 8; individual values e.g. spray 17–85 %).
- CL_iv = 820 ml/min = 49,2 L/h; V = 2,8 l/kg; t½ ≈ 2 h; t_max ≈ 11–13 min
  (spray 19,2 min); C_max_nasal/C_max_iv = 54–60 %.
- The authors’ own conclusion: *"there was no significant difference in the rate or
  extent of absorption between the different nasal treatments."*

**R2 (GI-ben, F_oral).** "Nicotine Population Pharmacokinetics in Healthy Smokers
After Intravenous, Oral, Buccal and Transdermal Administration." *Clin Pharmacokinet*
2020/2021. DOI **10.1007/s40262-020-00960-5** (PMC8016787). Verified 2026-09-25.
- F_oral = 40 % (published interval 17–69 %); swallowed fraction in buccal
  formulations 55–69 %; K_gi (oral) = 1,55 h⁻¹; transit K_tr 3,5–5,5 h⁻¹;
  CL = 67 L/h; V = 4,3 l/kg. Route-separated model: oromucosal absorption F = 100 %,
  swallowed fraction via GI with F = 40 % — **this is exactly our null model**.

**R3 (hepatic extraction).** *Hukkanen/Benowitz*, "Nicotine Chemistry, Metabolism,
Kinetics and Biomarkers" (PMC2953858): total CL ≈ 1200 ml/min; non-renal CL ≈
70 % of liver flow ⇒ approximately 70 % extraction per liver passage. Verified 2026-09-25.
The model here uses well-stirred: E_H = 1 − CL/Q_H with Q_H = 90 L/h ⇒ E_H ≈ 0,23–0,26.

**R4 (slemclearance).** Caponnetto P et al. "Saccharin test for the study of
muciliary clearance: reference values for a Spanish population." PMID **19006634**:
n = 249, nasal MCT median 16 min, medel 17,17 ± 8,43 min, 2,5–97,5 percentile 6–36 min.
Hastighet: 5–8 mm/min (Puchelle F et al., *Acta Oto-Laryngol* 1981;91:297–303,
DOI 10.3109/00016488109138511; jmf. popline-abstract 1990: 5,0–8,0 mm/min).

**R5 (GEOMETRI).** Ličen A, Grmaš J, Gubič Š, Trontelj J, Gomboc T, Hriberšek M et al.
"Development, Fabrication and Application of a Sectioned 3D-Printed Human Nasal Cavity
Model for In Vitro Nasal Spray Deposition Studies." *Biomedicines* 2026;14(2):329.
DOI **10.3390/biomedicines14020329**. Verified 2026-09-25.
- Nasal volume 15–20 ml, **total area ≈ 150 cm²**, length 12–14 cm.
- **Vestibule: "unsuitable for systemic absorption"**; respiratory region
  **≈ 130 cm²** (turbinates, "the main site for systemic uptake"); olfactory region ≈ 15 cm².
- Five regions: vestibule, olfactory, middle/upper turbinate, lower turbinate, nasopharynx.
- **Measurement (3D-printed nasal model, two spray devices):** deposition in olfactory region and
  nasopharynx **< 1 %**; "majority of the dose accumulated in the vestibule and lower
  turbinate regions"; wider plume + smaller droplets ⇒ more vestibular deposition.
- Support: Larimi M et al., *Pharmaceutics* 2023;15(3):406, DOI 10.3390/ph16030406:
  **> 90 % of 16–30 µm droplets deposit in the vestibule/nasal valve**; ~20 % of < 5 µm.

**R6 (local metabolism, qualitative).** CYP2A is active in human olfactory epithelium
(Chen Y et al., *Biochem Pharmacol* 2003;66:1245–1251) and cotinine is measured in
nasal lavage (Özdenes MH et al., *Int J Environ Res Public Health* 2009;5(1):11,
DOI 10.1186/1617-9625-5-11). ⇒ R_met is physically necessary but its magnitude is
**not** published for humans → assumption, not source.

## 4. Frozen acceptance criteria

- **C1 (primary, absolute level).** |F_abs_pred(nasal spray) − 0,58| within a factor of 2
  (0,29–1,16). Met = yes/no.
- **C2 (all sites, one parameter set).** F_pred for septum, conchae and spray
  all within a factor of 2 of the respective 0,76 / 0,64 / 0,58.
- **C3 (structural ceiling — the trap).** The model’s **ceiling** F_max = max over free
  parameters of F_pred must be ≥ measured F for each site. If F_pred fits within
  a factor of 2 but F_max < F_measured, the level is "hidden by plasma fit" ⇒ C3 = NOT MET.
- **C4 (spread).** Predicted difference F(septum) − F(spray) within a factor of 2 of
  measured 0,18, with positive sign. Note: the measured values are **not significantly different**
  according to R1, so C4 tests the model’s **content**, not the measurement.
- **C5 (identity/mass balance).** Σ(fractions) = 1 to < 1e-9 relative; all
  flows have the correct dimension; the analytical limiting case is reproduced to < 1e-8 relative.
- **C6 (null model/placebo).** A route-separated pooled model with ONE common absorbed
  fraction (the comparison baseline in the question) is reported on the same three data points.
  Geometry content counts only if the model with the same number of free parameters gives
  **smaller maximum error than the null model** on the three sites.

## 5. What counts as an error

- Fitting k_epi, p_vest or χ per site to hit R1.
- Reporting F within a factor of 2 without also reporting C3.
- Using n = 8 individual values with the spread 17–85 % as "support" for geometry.
- Inventing measured data, or calling an assumption a source.
- Falling back on a null model if the geometry model fails.

## 6. Builds on / reuse

- Builds on: `inputs/QUESTION.md` (Q146), `inputs/NIGHT_PREAMBLE.md` §2 (null models
  required), §3 (precision, source vs derivation vs hypothesis, UNKNOWN).
- Sibling `results/BT-HX-Q132` (empty, brief only) and `results/BT-HX-Q134` exist but
  lie outside the working directory and were **inaccessible** in this run — nothing
  is reused from there.
- The question points to `src/bodytwin/cells/respiratory/mucociliary_clearance.py` and
  `.../gastrointestinal/hepatic_clearance.py`. **These files do not exist** in this
  workspace (checked with glob) ⇒ the pointer is dead; the physics is rewritten from first
  principles here. No `import *`.

## 7. Resurser

1 thread, < 1 GB, no batch/sweep, no cloud, no GPU. The computation is an
ODE integration at a couple of hundred time points.
