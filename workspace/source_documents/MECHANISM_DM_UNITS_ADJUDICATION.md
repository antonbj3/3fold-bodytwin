# DM eTibia units adjudication — Newtons or pounds-force?

**Confidence tier: in-vivo-anchored** (published instrumented-knee force literature + an on-disk,
unit-unambiguous Newton cross-check, independent of that literature).
**Evidence JSON:** `docs/MECHANISM_DM_UNITS_ADJUDICATION_evidence.json`
**Independence note:** derived fresh from raw data + live literature fetches this session; a coordinator's
prior view was deliberately not consulted.

## Verdict (lead, so it isn't buried)

The raw `Fx,Fy,Fz` (and `Tx,Ty,Tz`) columns in `data/external/simtk_kneeloads_DM_measured_force/DM_*_knee_forces.csv`
are **pounds-force (lbf), not Newtons**. DM's true peak tibiofemoral contact force in level overground gait is
in the **~250–330 %BW bucket** (mean ≈ 287–305 %BW depending on body-mass assumption) — close to the
**~296 %BW** figure named in the task prompt — **not** the 58–68 %BW bucket. Both bucket names in the prompt
were tested; the sub-100 %BW bucket is falsified by three independent checks below, one of which needs no
literature and no body-mass assumption at all.

---

## 1. Raw peak resultant force, both unit hypotheses

Computed `sqrt(Fx²+Fy²+Fz²)` per row, per trial, from raw CSVs (script:
`/tmp/coordinator-1000/-home-anton/6f480b09-3647-4c07-8360-2f41db78e3c2/scratchpad/plateau_and_final.py`).

| trial set | n trials | peak resultant (raw units) |
|---|---|---|
| `DM_ngait_og{1,2,3,4,5,6,7,9}` (overground, the literal glob given) | 8 | 386.8 – 482.2 (global max in `og4`) |
| broader `DM_ngait*` (+ treadmill/transition) | 15 | up to 572.5 (in `ngait_tm_ss1`) |

Body mass: `DM.osim` (23-body scaled full-body model, bilateral-symmetric leg masses) sums to **69.98 kg**
— about 3–4 kg above the task's 66–68 kg prior. Both are reported; the verdict is robust across the whole
66–70 kg range (shown below).

**%BW at mass = 68 kg (BW = 666.9 N), 8-trial og-set:**

| interpretation | mean | range |
|---|---|---|
| raw = Newtons | **66.6 %BW** | 58.0 – 72.3 %BW |
| raw = lbf (×4.4482216) | **296.1 %BW** | 258.0 – 321.6 %BW |

At mass = 66 kg: N→68.6 %BW (59.8–74.5), lbf→305.1 %BW (265.8–331.4).
At mass = 69.98 kg (on-disk model mass): N→64.7 %BW (56.4–70.3), lbf→287.7 %BW (250.7–312.5).
The qualitative split (well under 100 %BW vs. ~250–330 %BW) does not change anywhere in this mass range.

## 2. Decorrelated external literature anchor (live-verified)

**Primary, same-lineage anchor** — D'Lima DD, Fregly BJ, Patil S, Steklov N, Colwell CW Jr. "Knee joint
forces: prediction, measurement, and significance." *Proc Inst Mech Eng H.* 2012;226(2):95–102.
DOI `10.1177/0954411911433372`, PMCID **PMC3324308** (fetched live this session). This is the Scripps/D'Lima
e-knee team's own summary of their instrumented-knee measurements — the same lab/device lineage the
Grand-Challenge "DM" (Sixth Competition) dataset comes from. Direct quotes:

> "forces transmitted across the knee joint during normal walking range between **2 and 3 times body weight**"

> Table 1 — Walking (Laboratory Floor): **2.5–2.8×BW**; Treadmill walking 1–3 mph: 2.1±0.2×BW; Power walking
> 4 mph: 2.8±0.4×BW.

No DM-specific (Sixth-Competition-only) published %BW figure was found live this session (checked the paper's
own text and searched for a COMAK/JAM validation paper naming "DM" — none found); the anchor is the cohort
figure from the same measurement lineage.

**Decorrelated cross-lab anchor** — Kutzner I, Heinlein B, Graichen F, Bender A, Rohlmann A, Halder A, Beier
A, Bergmann G. "Loading of the knee joint during activities of daily living measured in vivo in five
subjects." *J Biomech.* 2010;43(11):2164–2173. DOI `10.1016/j.jbiomech.2010.07.036`. This is a fully
independent lab (Berlin/Charité "OrthoLoad"), independent instrumented-implant device, independent subject
cohort. Its primary text was **not** directly fetchable this session (Elsevier redirect required JS; a wrong
PMC-ID guess returned an unrelated spine paper) — verified only **indirectly**, via a live-fetched tertiary
review: Wu XD et al., "Sensor-integrated hip and knee prostheses," *Front Bioeng Biotechnol.* 2025;13:1721499,
DOI `10.3389/fbioe.2025.1721499`, PMCID **PMC12907417** (fetched live), which quotes:

> "peak contact forces of **3.3 times the body weight** during level walking" (Kutzner et al., 2010)

**Literature consensus: ~2.1–3.3×BW (210–330 %BW) peak tibiofemoral force in level walking — decisively not
~0.6×BW.**

## 3. Physical-consistency check (statics) — can peak JCF be <1×BW in single-limb stance?

**No.** During single-limb stance nearly the whole GRF is transmitted axially up the limb to the knee; active
muscle forces crossing the joint (needed to balance the moment the GRF creates about the knee center, since
muscles have much shorter moment arms than the GRF's lever arm) can only **add** compression, never subtract.
So peak JCF must be **≥** peak GRF during that same stance phase. Textbook peak vertical GRF in walking is
~1.0–1.2×BW — already >1×BW before any muscle contribution.

This was **not** taken on faith — it was checked directly against DM's own force-plate data:

- `experimental_data/motion_analysis/DM_ngait_og{1,3,4,5,6,7,9}_grf.sto` (OpenSim-JAM mirror, standard OpenSim
  `.sto` Newton convention, unambiguous) — peak vertical GRF **736–762 N**, mean 751.9 N →
  **109.6–117.8 %BW** at mass 66–70 kg. DM's own measured GRF confirms the textbook ~1.0–1.2×BW figure
  directly, not just generically.

**The decisive, mass-independent, literature-independent check:** for the 7 trials that have both an
independently-sourced `_knee_forces.csv` (SimTK archive) and `_grf.sto` (OpenSim-JAM GitHub mirror), computed
the ratio peak-resultant-JCF(raw) / peak-vertical-GRF(N). This ratio needs **no** body-mass assumption and
**no** literature citation — only internal consistency between two independently-provenanced files for the
same subject.

| trial | peak GRF (N) | peak JCF (raw) | ratio if raw=N | ratio if raw=lbf |
|---|---:|---:|---:|---:|
| og1 | 758.7 | 437.9 | 0.577 | 2.567 |
| og3 | 757.5 | 386.8 | 0.511 | 2.271 |
| og4 | 759.1 | 482.2 | 0.635 | 2.825 |
| og5 | 745.3 | 456.9 | 0.613 | 2.727 |
| og6 | 762.2 | 450.2 | 0.591 | 2.627 |
| og7 | 736.1 | 415.6 | 0.565 | 2.512 |
| og9 | 744.1 | 442.4 | 0.595 | 2.645 |
| **mean** | | | **0.584** | **2.596** |

Mechanics requires this ratio >1 (JCF exceeds GRF). **Under raw=N, all 7 trials give a ratio <1 — mechanically
impossible** (would mean the knee transmits less force than the ground reaction passing through the same
limb). **Under raw=lbf, all 7 trials land in [2.27, 2.83]**, tightly matching the independently-sourced
~2–3× literature expectation, across every trial in the diverse instance-space checked.

*Caveat forced honestly:* the CSV (SimTK telemetry clock) and `.sto` (OpenSim-JAM mocap/force-plate clock) are
two independently-sourced files for the "same" named trial; a spot-check found GRF=0 at the CSV's own peak-force
timestamp, i.e. **the two clocks are not aligned**. The ratio argument above uses trial-level peaks (not
same-instant samples) specifically so it survives this: GRF genuinely reaches ~750 N at some point in every
trial, which mechanically necessitates a comparably large JCF at that same real-world instant, and under raw=N
the CSV's own reported trial-peak never reaches that magnitude (falls short by ~1.6–2×) — the contradiction
does not depend on exact timestamp pairing.

## 4. Secondary check (reported honestly — mixed, not decisive)

`DM_1legstand{1,2}_knee_forces.csv` (quasi-static single-leg balance): isolating the genuine low-variance
plateau (Fz > 70 % of trial max, knee angle constant at ~6° for 6+ seconds) gives 546.6–615.9 raw. Under
raw=N at 68 kg: ~82–92 %BW (plausible for sustained single-limb support, slightly under 1×BW, explainable by
handrail support/sensor offset). Under raw=lbf: ~364–410 %BW sustained for 6+ seconds at a passively-stable,
near-fully-extended knee angle — not corroborated by any activity-matched literature value (D'Lima's review
has no "one-leg stance" row to check against). This is an **honest, unresolved tension**, not swept under
the rug — but it carries less evidential weight than checks 1–3 (no literature anchor exists for this specific
activity, and the exact clinical protocol/dynamics of "1legstand" — quiet balance vs. a deliberately
challenging test — is not independently confirmed), so it does not overturn the primary verdict.

## Adjudication

- **Units: pounds-force (lbf).** DM's true peak resultant tibiofemoral contact force in level overground gait
  ≈ **287–305 %BW** (mean across mass assumptions), individual-trial range **250–330 %BW** — matching the
  ~296 %BW bucket, not 58–68 %BW.
- **Falsifier:** this verdict is overturned if a primary source (original SimTK `kneeloads` codebook/README, or
  a Grand-Challenge/D'Lima methods paper) explicitly labels the eTibia telemetry columns as Newtons. No such
  explicit label was found this session — not in `PROVENANCE.md`, not in `DM.osim`, and not on the live SimTK
  project page (checked directly: "the page does not specify units for force data files"). This is an
  acknowledged gap; absent a direct label, the verdict rests on three independent, mutually-reinforcing
  indirect checks (§1 magnitude-vs-literature, §3 statics-impossibility, §3 mass/literature-independent
  JCF/GRF ratio) that all agree, plus one honestly-flagged unresolved tension (§4).
- **Residual uncertainties:** exact DM mass (66–68 kg prior vs. 69.98 kg on-disk model — doesn't flip the
  verdict); CSV↔`.sto` clock alignment unestablished (argument built to survive it); no DM-specific published
  %BW found; Kutzner primary text verified only indirectly via a tertiary review; §4 remains an open tension.

## Paths

- Data: `source_repository/data/external/simtk_kneeloads_DM_measured_force/DM_ngait_og*_knee_forces.csv`,
  `source_repository/data/external/simtk_kneeloads_DM_measured_force/DM_1legstand{1,2}_knee_forces.csv`
- Cross-check GRF: `source_repository/data/external/opensim_jam_grand_challenge_DM/experimental_data/motion_analysis/DM_ngait_og{1,3,4,5,6,7,9}_grf.sto`
- Body mass source: `source_repository/data/external/opensim_jam_grand_challenge_DM/DM.osim`
- Evidence JSON: `source_documents/MECHANISM_DM_UNITS_ADJUDICATION_evidence.json`
- Scratch scripts (not part of the deliverable, kept for repro): `/tmp/coordinator-1000/-home-anton/6f480b09-3647-4c07-8360-2f41db78e3c2/scratchpad/{compute_peak_force.py,plateau_and_final.py,grf_vs_kneeforce.py}`
