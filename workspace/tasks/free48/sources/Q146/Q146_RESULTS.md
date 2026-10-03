BT-HX-Q146

# Q146 — Nasal pathway with swallowed and systemic share: first workable mechanistic model

**Conclusion:** deposition in 5 nasal regions → competition between epithelial
uptake and mucus clearance → swallowed portion to GI (first passage) → plasma.
Total control volume per region; the same F_sw [mol/s] enters the GI model.
Running: `model.py` (→ `results.json` ), `test_model.py` (12/12). PREREG
frozen before execution: `PREREG.md` + `PREREG.sha256` (fed66c26…).

## Prediction against reference

Referens **R1** = Johansson C-J m.fl., *Eur J Clin Pharmacol* 1991;40:581–588,
DOI 10.1007/BF00314989 — 1 mg nikotin intranasalt (septum/conchae/spray) mot
0,7 mg iv, n = 8. GI-ben **R2** = *Clin Pharmacokinet* 2020/21,
DOI 10.1007/s40262-020-00960-5 (F_oral = 0,40). Geometri **R5** = Ličen m.fl.,
*Biomedicines* 2026;14(2):329, DOI 10.3390/biomedicines14020329.

| site | F_pred | F_measured (R1) | quota | ceiling F_max | F_sw | C_free,max [mol/m³] | AUC_plasma [mol·s/m³] | t_max [min] |
|---|---|---|---|---|---|---|---|---|
| septum  | 0,584 | 0,76 | 1,30 | 0,873 | 0,903 | 44,7 | 0,184 | 77 |
| conchae | 0,576 | 0,64 | 1,11 | 0,808 | 0,914 | 92,5 | 0,182 | 78 |
| spray   | 0,440 | 0,58 | 1,32 | 0,643 | 0,688 | 184,9 | 0,139 | 77 |

The model underpredicts all three by 11–32 %, **without site-specific adaptation**.

## Kriterier (frysta i PREREG §4)

| | kriterium | utfall |
|---|---|---|
| C1 | F(spray) within factor 2 of 0,58 | **YES** (1,32) |
| C2 | all three sites within factor 2 | **YES** (1,30 / 1,11 / 1,32) |
| C3 | model ceiling F_max ≥ F_measured (structure test) | **YES** (0,873/0,808/0,643 ≥ 0,76/0,64/0,58) — but the spray ceiling is only 11 % above the measurement |
| C4 | dispersion F(septum)−F(spray): sign and factor 2 against 0,18 | **YES** (+0,144, ratio 1,25) |
| C5 | identities: Σ shares = 1; mass balance; 12/12 tests; dimensional analysis | **YES** (1,1e-16; 1,6e-11 mol; 15 expressions) |
| C6 | the geometry model should beat pooled null model (1 free parameter) on the same 3 data points | **NOT** |
| — | t_max (not in PREREG, still reported) | **NOT**: 77 min vs. 11–13 min (R1) |

## What can actually be shown — identifiability (independent of k_epi)

The measurement F_in = 0,58–0,76 (R1) is **above F_oral = 0,40** (R2). The swallowed
path can therefore give at most 0,40. Conclusion without model parameters:

- F_nasal_rutt ≥ F_in − F_oral = **+18 pp (spray), +24 pp (conchae), +36 pp (septum)**
  must come from the nasal mucosa. Model's own contributions: 0,165 / 0,211 / 0,222 — compatible.
- ΔF = 18 pp between septum and spray occurs with **same molecule,
  same dose, same GI leg**; only the place of deposit differs.
  Thus, geometry is *proven to exist* as a cause in principle.
- The model's requirement for k_core: Pe = k_epi·a/v must be **3,2 (septum), 1,5
  (conchae), 4,5 (spray)**, i.e. k_epi ≈ 0,16–0,49 mm/s against my assumption 0,01
  mm/s (passive permeability). The requirement is **consistent across sites** — it
  is a measurable requirement for carrier-mediated uptake, not a free speech.
- Structurally fresh: F_max(spray) = 0,643. If more than **≈5 %** of the spray dose
  lands in the non-absorbing vestibule, the model shape is falsified by R1.

## What failed

- **C6: the geometry is not separable from a single pooled number.** The null model (F
  = 0,66 for all three) gives a maximum error of 0,10; the geometry model gives a maximum
  error of 0,18. With 3 data points, 1 dex uncertainty in k_epi and n = 8 this is the expected
  outcome — and R1 himself says that the differences between the sites are **not significant**
  (spray individual values ​17–85 %). F146's geometry finding is thus **not proven
  by existing measurement**; only the existence condition above applies.
- **t_max misses by a factor of 4–6.** The model makes the nose step
  ~0,08 s; the measurement shows ~10–20 min. Either uptake is carrier mediated
  (slower than membrane), or deposition/delivery is time-prolonged.
  Unsolved, and it is the single most important resolved issue.
- Sensitivity ±50 % (spray, order): **p_vestibulum ±17,5 %** > v_clear +7,6/−2,8 %
> k_epi −4,3/+3,9 %. The outcome is controlled by **where** the dose lands, not how fast
  it is taken up.
- External loss 11,3% of dose (vestibular runoff)
  is an assumption (χ = 0,15), not a source.

## Next resolution step (data, geometry, measurement)

1. **Measure F_sw per site** (gamma machine scintigraphy or ^14C butter sample, Ličen method on
   3D printed nose model) for exactly the same pump and droplet volume as R1. This gives
   p_vestibulum and p_olfactory directly and determine the falsification limit of 5 %.
2. **Measure k_epi for nicotine on human nasal mucosa** (Ussing chamber, cut from
   turbinatectomies, pH 7,4). Required ~0,2–0,5 mm/s for R1 to apply; about that
   turns out 10^-5 mm/s is R1's 58–76 % inexplicable with F_oral = 0,40 and
   the question must be rephrased.
3. **Plasma with 30–60 s sampling 0–30 min** to distinguish the nasal step
   time constant from GI passage's; The t_max separation is cheap and decisive.
4. **Paired subject×site-data** (n ≥ 12, same person, all three places) — without
   pairing, the geometry cannot be separated from person-to-person variation.
5. `k_epi`, `f_met` and χ should be quoted or measured, not assumed.

##files
`PREREG.md` , `PREREG.sha256` , `model.py` , `test_model.py` , `results.json` (all numbers).
Builds on: `inputs/QUESTION.md` , `inputs/NIGHT_PREAMBLE.md` §2–3. The reuse pointers
in the question (`src/bodytwin/cells/respiratory/mucociliary_clearance.py`, `.../gastrointestinal/hepatic_clearance.py`
) **do not exist** in this workspace — the physics is rewritten from
first principles. Sibling `results/BT-HX-Q134` exists but was inaccessible
outside the working directory; nothing is reused from there.
