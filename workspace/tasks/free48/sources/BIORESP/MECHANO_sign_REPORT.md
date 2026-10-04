# MECHANO_sign: does the tooth referenced mechanostat give the right sign against load proxies?

Node `DENT-IF-BI-MECHANOSTAT` (direction test). Node rows: `notes/expansion/mechanostat_sign.jsonl`.
- PREREG frozen before the test arms and before reading the anchors (sha256 `4edbeed3…`, code in `PREREG.code_sha256`). AMENDMENTS: A1 encoding of "ns without CI", A2 fallback according to the rule.
- ATOMS 20/20 GREEN (`atoms_run.txt`). Own overtaking in new process bit-identical in three cases. The load additions checked: axial resultant × 1,70, moment exactly F_ax·7 mm.

## Setup
REMODEL L1 (THETA0 unchanged, G0 matched platform, t_muc 3) run for check and test per proxy, 6 scenarios, 39 FE runs (11,5 min, 230 MB). The baseline is
constant loss, so 0 for each contrast. Rule: |ΔMBL 5 years| ≥ 0,10 mm counts as power, and the sign of the armature counts only when the CI excludes 0.

## Results (primary)
| proxy | model ΔMBL 5 year | Δ stimulus q_p50 | held-out anchor (test − control) | L1 | baseline |
|---|---|---|---|---|---|
| C/I 0,7 → 1,5 | +7,19 (kollaps) | +0,81 | ns, p 0,155 (PMID 33072508, 4350 implantat) | FAIL | PASS |
| smal 3,3 − regular 4,1 | +0,25 | +0,23 | −0,15 [−0,32; 0,01] (PMID 37085829) | FAIL | PASS |
| cantilever (liver × 1,7) | 0,00 | +0,24 | ns, trend + (PMID 39570465) | PASS | PASS |
| bruxism × 1,5 | +9,5 (collapse) | – | UNKNOWN: none of three selected MA had MBL-analys | – | – |

**Verdict according to PREREG: STRUKEN also as a modifier.** L1 gets 1/3 PASS and the baseline 3/3. No held-out anchor has a significant sign, so the model's sign cannot be confirmed. Same verdict in all 6 scenarios (s_lazy 0,35, κ 3, ℓ 0,1, θ-mean, platform switching 0,35).

## Interpretation
- **The direction of the stimulus is physically plausible:** ε_imp/ε_tand rises with all three proxies. The outcome in the leg is still controlled by a rock. The answer is either 0 (cantilever) or collapse (C/I, bruxism, cantilever as moment +9,2). There is no graded response to which the clinic can be compared.
- **The diameter anchor goes the other way** (narrow implants tend towards less loss, p 0,07). The model says more loss.
- **Bruxism causes collapse at 1,5 × load.** The non-held-out PMID 36896080 (ns) and REMODEL's overload tankers do not point there.
- Together with REMODEL (the mechanics REFUTED of the level) there is thus nothing left for the mechanics term in MBL-kedjan. Bio leg L2 becomes sole predictor. The load proxies do not affect L2, which is consistent with three ns anchors.

## Disclaimer
- The contrast C/I 1,5 is right on the split boundary of the anchor (> 1,5 vs < 1,5).
- The cantilever anchor blends front and rear areas.
- AAbstract and a full text (PMC11582258) were read. Forest plot numbers were not read (they are only available as an image).
- A0: approximate prior knowledge of C/I- and the cantilever literature is reported in PREREG. The selected articles were unread.
- An ns anchor cannot distinguish "no effect" from "too small an effect". The test imposes the term on the rule (no gain against the baseline), not on the loss being proven to be independent of load.

## Next point of stress
A held-out series with significant load dependence MBL (eg measured occlusal force or EMG-bruxism vs. MBL). Then the sign can be tested for real. A graded mechanostat without a rock is only tested afterwards.

Reproduction: `mechano_sign.py run` → `compare.py` → `rerun_compare.py` → `make_evidence.py` → atoms_runner. Raw data `external_media`.
