# MECHANISM Fmax vs PCSA VALIDATION — is the twin over-strong? (2026-07-21)

Directly executes `docs/MECHANISM_CROSS_SUBJECT.md` Sec.6 gap #4's own named next step ("would need a
controlled Fmax-rescaling experiment across subjects to test directly") and checks whether the
twin's ~1.4-2.0x over-prediction of in-vivo knee/hip contact force (`docs/
MECHANISM_JOINT_FORCE_SCORECARD.md`, `docs/MECHANISM_CROSS_SUBJECT.md`) is plausibly explained by
muscles that are simply too strong. Anchor is Handsfield et al. 2014 MRI muscle-volume regressions —
**independent of the OrthoLoad contact-force anchor** every other joint-force cert in this repo is
validated against.

**Revised headline (the union-aggregate first suggested "matches"; a decisive causal test caught that
this was too coarse — see Sec.4):** the twin's **knee-crossing prime movers are systematically,
one-directionally over-strong** relative to Handsfield-PCSA-implied Fmax (12/12 muscles above 1.0,
mean ratio 1.146 at specific_tension=60N/cm², implied specific tension 60.5-79.3 N/cm²) — and
**causally correcting them shrinks subject2's own knee over-prediction by 17.1%** (1.515x → 1.255x
vs OrthoLoad), closing roughly half the gap to a ratio of 1.0. The **hip-crossing prime movers show
a genuinely mixed pattern** (ratios 0.649-1.446, roughly as many muscles under 1.0 as over) that nets
to **no material effect when corrected** (hip ratio 1.412x → 1.417x, effectively unchanged). So:
**strength calibration IS a real, substantial, causally-confirmed partial contributor to the KNEE
over-prediction specifically, and is RULED OUT for the hip** — a joint-differentiated answer, not a
single up/down verdict for "the twin."

## 0. Method, anchor, and a citation this task's own prompt got wrong

**PCSA identity** (unchanged from `scripts/msk/subject_specific_scaling.py`, which used the same
formula with a cruder bone-geometry volume proxy): `PCSA = muscle_volume / optimal_fiber_length`;
`Fmax_pcsa = PCSA x specific_tension`. Volume here comes from **real in-vivo MRI data on actual
muscle tissue** (Handsfield et al. 2014) rather than either of that prior doc's two bone/mass proxies
(explicitly disclosed there as "neither proxy ever measures muscle tissue at all").

**Citation correction (forced, not silently propagated):** the task specified **PMID 24581800** for
Handsfield et al. 2014. Verified LIVE via raw `curl` of PubMed's own `citation_*` meta tags (not a
summarizer model — which, on a first attempt, *also* independently hallucinated a wrong paper) that
this PMID belongs to an unrelated sodium-channel-subunit paper (Liu et al., *C R Biol* 2014). The
correct citation, independently re-derived and cross-confirmed via NCBI eutils + Semantic Scholar's
API, is: **Handsfield GG, Meyer CH, Hart JM, Abel MF, Blemker SS. "Relationships of 35 lower limb
muscles to height and body mass quantified using MRI." *J Biomech.* 2014 Feb 7;47(3):631-8.
doi:10.1016/j.jbiomech.2013.12.002. PMID 24368144.** Full text is paywalled (Unpaywall:
`is_oa:false`; Semantic Scholar: `openAccessPdf.status:"CLOSED"`; not in PMC). The per-muscle
regression table used here is a disclosed **third-party transcription** of the paper's own
Supplementary Table 3 (`github.com/johnjdavisiv/scale-muscle-strength`), internal-consistency-checked
before use (Sec.1). Full provenance: `data/external/handsfield2014/PROVENANCE.md`.

**Specific-tension band:** pre-registered at the task's own requested **[45, 60, 100] N/cm^2**
(low/mid/high). This is deliberately **charitable to the null** (biased against finding
over-strength): live-verified from the actual primary source of this exact model lineage's Fmax
values (Arnold, Ward, Lieber & Delp 2010, *Ann Biomed Eng*, PMC2903973, fetched + grepped directly),
the "true" single-fiber/mammalian literature values are **lower** — "higher than the range of values
(11-47 N/cm^2) reported previously, and larger than the experimentally measured value for mammalian
muscle of 22.5 N/cm^2" (quoted directly from that paper). Using the lower, more-standard values
instead would shrink the PCSA-implied-Fmax denominator and make any over-strong finding **larger**,
not smaller.

**Subject anthropometrics** — machine-verified live this session from each subject's own
`sessionMetadata.yaml` (not recalled from prior docs' prose): subject2 (this task's primary target)
1.96m/78.2kg/M; subject3 1.69m/63.5kg/F; subject4 1.68m/62.6kg/F (both bonus cross-subject checks,
reusing `docs/MECHANISM_CROSS_SUBJECT.md`'s own cohort); Ward/Arnold cadaver-mean 1.684m/82.7kg
(control point, re-typed directly from Arnold 2010's own Methods text: "average height of
168.4+-9.3cm and weight of 82.7+-15.2kg").

## 1. Verification gates (machine-checked, all PASS)

Run: `scripts/msk/fmax_pcsa_validation.py` (real run, exit 0).

| gate | result |
|---|---|
| Handsfield CSV per-muscle-group `volume_proportion` sums to 1.0 (29 groups, right side) | **PASS**, all == 1.0 |
| XML-parse vs OpenSim-API-parse agreement on model Fmax/Lopt (80 muscles x 2 fields = 160 comparisons) | **PASS**, 0 mismatches |
| Model muscle-name-set (n=80) vs Handsfield-table rows (n=80) | **EXACT MATCH** |
| Spot-check vs `docs/MECHANISM_MUSCLE_AUDIT.md`'s already-published table (gasmed_r, soleus_r, bflh_r, glmax1_r) | **PASS** — same model, same values, fresh script |
| Bilateral (L/R) Fmax exact-match (40 pairs) | **PASS**, all exact (re-confirms the structural "Fmax never rescaled" fact live, this script) |
| Predicted muscle volume positive, all subjects x 40 muscles (160 checks) | **PASS**, 0 non-positive (psoas_r/tfl_r have negative regression intercepts but stay positive at every H×M tested) |

## 2. Per-muscle table — subject2 (primary target), sorted by model Fmax

`PCSA` = Handsfield-predicted volume at subject2's own H=1.96m/M=78.2kg, divided by the model's own
optimal fiber length. `ST_impl` = the specific tension (N/cm^2) that would make model Fmax exactly
equal the PCSA-implied value — the single cleanest per-muscle summary number.

| muscle | Fmax_model (N) | PCSA (cm^2) | Fmax@45 | Fmax@60 | Fmax@100 | ratio@60 | ST_impl (N/cm^2) |
|---|---:|---:|---:|---:|---:|---:|---:|
| soleus_r | 6194.8 | 93.16 | 4192.2 | 5589.7 | 9316.1 | 1.108 | 66.5 |
| vaslat_r | 5148.8 | 71.05 | 3197.1 | 4262.8 | 7104.6 | 1.208 | 72.5 |
| gasmed_r | 3115.5 | 41.67 | 1875.3 | 2500.5 | 4167.4 | 1.246 | 74.8 |
| vasmed_r | 2747.8 | 39.23 | 1765.3 | 2353.8 | 3923.0 | 1.167 | 70.0 |
| semimem_r | 2201.0 | 27.75 | 1248.8 | 1665.0 | 2775.1 | 1.322 | 79.3 |
| recfem_r | 2191.7 | 35.22 | 1585.0 | 2113.4 | 3522.3 | 1.037 | 62.2 |
| tibpost_r | 1730.2 | 25.09 | 1129.0 | 1505.4 | 2509.0 | 1.149 | 69.0 |
| vasint_r | 1697.4 | 21.82 | 981.9 | 1309.2 | 2182.0 | 1.296 | 77.8 |
| gaslat_r | 1575.1 | 20.95 | 942.6 | 1256.8 | 2094.7 | 1.253 | 75.2 |
| psoas_r | 1426.8 | 27.16 | 1222.3 | 1629.7 | 2716.2 | 0.875 | 52.5 |
| glmax2_r | 1406.0 | 24.35 | 1095.6 | 1460.8 | 2434.6 | 0.963 | 57.8 |
| bflh_r | 1313.2 | 20.39 | 917.8 | 1223.7 | 2039.5 | 1.073 | 64.4 |
| tibant_r | 1227.5 | 18.71 | 842.1 | 1122.8 | 1871.4 | 1.093 | 65.6 |
| perlong_r | 1115.4 | 17.33 | 780.0 | 1040.0 | 1733.4 | 1.072 | 64.3 |
| glmed1_r | 1093.5 | 18.96 | 853.2 | 1137.5 | 1895.9 | 0.961 | 57.7 |
| piri_r | 1029.8 | 15.35 | 690.9 | 921.2 | 1535.3 | 1.118 | 67.1 |
| iliacus_r | 1021.1 | 17.99 | 809.6 | 1079.4 | 1799.0 | 0.946 | 56.8 |
| glmax1_r | 983.8 | 18.63 | 838.3 | 1117.8 | 1863.0 | 0.880 | 52.8 |
| glmax3_r | 947.8 | 14.94 | 672.2 | 896.3 | 1493.9 | 1.057 | 63.4 |
| addlong_r | 916.8 | 15.07 | 678.2 | 904.3 | 1507.2 | 1.014 | 60.8 |
| fhl_r | 907.8 | 15.25 | 686.2 | 915.0 | 1525.0 | 0.992 | 59.5 |
| glmed3_r | 871.2 | 14.42 | 648.8 | 865.1 | 1441.9 | 1.007 | 60.4 |
| glmed2_r | 765.1 | 12.03 | 541.5 | 722.0 | 1203.3 | 1.060 | 63.6 |
| addbrev_r | 625.8 | 9.72 | 437.4 | 583.2 | 972.0 | 1.073 | 64.4 |
| edl_r | 603.5 | 9.61 | 432.4 | 576.6 | 960.9 | 1.047 | 62.8 |
| addmagDist_r | 597.3 | 7.62 | 343.1 | 457.5 | 762.5 | 1.306 | 78.3 |
| addmagIsch_r | 597.3 | 8.75 | 393.5 | 524.7 | 874.5 | 1.138 | 68.3 |
| addmagMid_r | 597.3 | 9.67 | 435.2 | 580.3 | 967.2 | 1.029 | 61.8 |
| addmagProx_r | 597.3 | 12.75 | 573.6 | 764.8 | 1274.7 | 0.781 | 46.9 |
| semiten_r | 591.3 | 9.45 | 425.3 | 567.1 | 945.1 | 1.043 | 62.6 |
| bfsh_r | 557.1 | 9.14 | 411.4 | 548.5 | 914.2 | 1.016 | 60.9 |
| perbrev_r | 521.2 | 8.99 | 404.6 | 539.5 | 899.1 | 0.966 | 58.0 |
| glmin3_r | 446.8 | 11.48 | 516.7 | 688.9 | 1148.1 | 0.649 | 38.9 |
| fdl_r | 423.2 | 6.25 | 281.2 | 374.9 | 624.8 | 1.129 | 67.7 |
| tfl_r | 411.2 | 7.59 | 341.6 | 455.5 | 759.2 | 0.903 | 54.2 |
| glmin2_r | 394.8 | 5.39 | 242.6 | 323.5 | 539.2 | 1.220 | 73.2 |
| glmin1_r | 374.0 | 4.31 | 194.0 | 258.7 | 431.2 | 1.446 | 86.7 |
| ehl_r | 285.9 | 4.27 | 192.1 | 256.1 | 426.8 | 1.116 | 67.0 |
| grac_r | 281.3 | 4.33 | 194.7 | 259.6 | 432.6 | 1.084 | 65.0 |
| sart_r | 249.4 | 4.12 | 185.6 | 247.5 | 412.5 | 1.008 | 60.5 |

Plausibility spot-check (well-established, non-controversial anatomy, used only as a forensic sanity
check, not as ground truth): predicted PCSA magnitudes are physiologically sane — soleus (93.2 cm^2)
and vastus lateralis (71.1 cm^2) as the largest, sartorius (4.1 cm^2) and gracilis (4.3 cm^2) among
the smallest despite sartorius being the body's longest muscle (near-zero pennation, textbook-known
low-PCSA paradox) — no garbage values.

## 3. Aggregate — does the model systematically exceed PCSA-implied strength?

| subject | muscle set | n | mean ratio@60 | median ratio@60 | min–max | mean ST_impl (N/cm^2) |
|---|---|---:|---:|---:|---:|---:|
| **subject2 (primary)** | all 40 | 40 | 1.071 | 1.066 | 0.649–1.446 | 64.3 |
| **subject2 (primary)** | **knee-crossing (12)** | 12 | **1.146** | 1.126 | **1.008–1.322 (12/12 > 1.0)** | 68.8 |
| **subject2 (primary)** | **hip-crossing (25)** | 25 | 1.040 | 1.037 | **0.649–1.446 (genuinely mixed)** | 62.4 |
| subject2 | union prime-movers (31) | 31 | 1.070 | 1.057 | 0.649–1.446 | 64.2 |
| Ward/Arnold cadaver-mean (control) | union (31) | 31 | 1.159 | 1.150 | 0.695–1.549 | 69.5 |
| subject3 (cross-subject bonus) | union (31) | 31 | 1.431 | 1.402 | 0.829–1.849 | 85.8 |
| subject4 (cross-subject bonus) | union (31) | 31 | 1.454 | 1.422 | 0.840–1.873 | 87.2 |

**Answer to "by what factor" — joint-differentiated, not single-number:** across ALL 40 lower-limb
muscles, subject2's model exceeds Handsfield-PCSA-implied Fmax by a modest ~7% on average (mean
1.071). That modest overall number **hides a real asymmetry**: the **12 knee-crossing muscles are
uniformly above 1.0** (every single one, mean +14.6%, implied specific tension a tight 60.5-79.3
N/cm^2 band) while the **25 hip-crossing muscles are a genuine mix** (0.649-1.446, roughly as many
below 1.0 as above) that averages out close to 1.0. A single "union" or "all-40" aggregate, taken
alone, would have reported "matches" and missed the knee-specific signal — caught only by
disaggregating per anatomical crossing-set (Sec.4 shows this is not just an artifact of averaging:
the disaggregation predicts, and a causal test then confirms, materially different downstream
effects at the two joints).

**Geometric governor identified:** because Fmax is a frozen, never-rescaled-per-subject quantity
(re-confirmed live, Sec.1) while Handsfield-implied PCSA scales with each subject's own
height x mass product, the ratio is a monotonically DEcreasing function of a subject's H×M. Subject2
(H×M=153.3, the tallest/heaviest of the 4 points tested — even slightly above the Ward/Arnold
cadaver-template's own 139.3) sits closest to a ratio of 1.0 not because the model was tuned for
subject2, but because subject2 happens to be close in body size to whatever implicit body the frozen
Fmax values best describe. Subject3/4 (H×M=107.3/105.2, much smaller) sit further out on the same
curve (ratio 1.40-1.45) — this is the SAME governing relationship, not a separate phenomenon, and
explains why subject3/4 show a larger mismatch without needing a different mechanism per subject.

## 4. The decisive causal test (forced, not left as a dangling correlation)

A striking parallel emerged from Sec.3 that this analysis is obligated to force, not admire from a
distance (OODA discipline): subject3/4 show BOTH a larger Fmax/PCSA mismatch (this doc) AND worse
joint-contact-force over-prediction (`docs/MECHANISM_CROSS_SUBJECT.md`: knee 1.68-1.99x vs subject2's
1.515x). Two Handsfield-PCSA-corrected models were built (specific_tension=60 N/cm^2, each subject's
own height/mass) and run through the IDENTICAL, unedited Static-Optimization + self-computed +
JointReaction pipeline (`scripts/msk/static_opt_knee.py`, `scripts/msk/validate_hip_force.py`, reused
via `scripts/msk/subject_specific_scaling.py`'s own `run_joint_forces_for_model`/
`compute_joint_forces_from_so_output`, and `scripts/msk/cross_subject_validation.py`'s
`configure_for_subject` for subject-path monkeypatching — **zero new force/geometry code**).

**Pre-registered threshold:** >=10% relative drop in the ratio-vs-OrthoLoad counts as a real,
quantified contributor.

| subject | joint | baseline ratio | corrected ratio | relative change | verdict |
|---|---|---:|---:|---:|---|
| **subject2 (primary)** | **knee_r** | 1.515x | **1.255x** | **−17.1%** | **REAL, substantial contributor** |
| subject2 (primary) | hip_r | 1.412x | 1.417x | +0.3% (self) | no material effect |
| subject3 (worst cross-subject case) | knee_r | 1.990x | 1.786x | −10.2% | REAL, quantified contributor (clears pre-registered 10% bar) |
| subject3 | hip_r | 1.664x (self) / 1.745x (JR) | 1.629x / 1.738x | −2.1% / −0.4% | does not explain the spread |

Per-muscle correction direction, subject2 (why the knee moves and the hip doesn't): every major
knee-crossing muscle's Fmax was corrected DOWNWARD — `gasmed_r` 3115.5→2500.5N (x0.803, the single
muscle already established in `docs/MECHANISM_CMC_SECOND_SOLVE.md`/`docs/MECHANISM_SUBJECT_SPECIFIC.md`
as supplying 45.7-47.1% of total knee-crossing force), `vaslat_r` 5148.8→4262.8N (x0.828), `vasmed_r`
x0.857, `vasint_r` x0.771, `gaslat_r` x0.798, `soleus_r` x0.902, `semimem_r` x0.756 — a consistent,
compounding, one-directional reduction. Hip muscles moved in BOTH directions: `psoas_r` 1426.8→1629.7N
(x1.142, UP), `glmax1_r` x1.136 (UP), `iliacus_r` x1.057 (UP), `glmed1_r` x1.040 (UP, the muscle
already established as the single largest hip-crossing contributor, ~23-24% of total) — netting to
approximately zero change in total hip contact force.

**Quality gates on both corrected runs** (not just trusting the output number): SO `convergence_pass`
subject-agnostic checks PASS both; `muscles_over_fmax_ratio_1p5={}` and `joint_reserve_leaning_hard={}`
for both (no degenerate/runaway solve despite reduced knee-muscle capacity); anatomical crossing-set
anchors 0 missing at both joints, both subjects; hip negative control (`bfsh_r`) correctly excluded
both; self-computed vs JointReaction agree to <0.01% at the knee both subjects (unaffected by the
known hip edge-truncation confound already diagnosed in `docs/MECHANISM_CROSS_SUBJECT.md` Sec.2.3,
which reproduces identically here, not a new artifact).

**Interpretation, precisely stated:** correcting knee-crossing Fmax to Handsfield-implied values
closes **~50% of subject2's own gap** between its baseline knee ratio (1.515x) and a ratio of 1.0
((1.515−1.255)/(1.515−1.0) = 0.505). It does **not** eliminate the over-prediction — a ratio of 1.255
is still a 25.5% over-prediction — so strength calibration is a real, substantial, but **partial**
contributor to the knee specifically, not the sole explanation. For the hip, the same causal test
rules strength calibration out entirely (<0.5% effect, direction not even consistent between
subjects).

## 5. Symmetric QC — the task's own named adversary, forced

**The adversary:** "LaiArnold itself was partly built from PCSA data, so some agreement is expected."

**Forced, not assumed:** verified live (Arnold, Ward, Lieber & Delp 2010, PMC2903973, raw HTML grepped
directly) that this model lineage's Fmax values ARE literally `PCSA x 61 N/cm^2`, where PCSA came from
**Ward et al. 2009's 21-cadaver dissection study** (mean height 168.4±9.3cm, mass 82.7±15.2kg;
`Ward SR, Eng CM, Smallwood LH, Lieber RL. Clin Orthop Relat Res. 2009;467:1074-1082`, independently
confirmed via its own reference-list entry in Arnold 2010). The repo author of the Handsfield CSV used
here independently notes specific tension "same as Rajagopal (60), Arnold (61)" — cross-confirming
the 60-61 N/cm^2 convention from a second angle. So this comparison is **PCSA-based vs PCSA-based**
(some methodological-family agreement genuinely IS expected, and is not itself informative) but
**decorrelated on**: (a) data source (cadaver dissection vs in-vivo MRI), (b) subject pool (21 older
cadaver donors vs 24 young healthy MRI volunteers), (c) subject2's own actual body vs either source
population's average.

**The control-point device that isolates the genuine part:** the Ward/Arnold cadaver-mean control
point (Sec.3, row 5) shows ratio 1.159/ST_impl 69.5 — i.e. even evaluated at the population the
model's Fmax literally came from, Handsfield's independently-collected MRI volumes don't reproduce it
exactly (a 16% gap even for the "expected-to-agree" case) — a reminder that "same methodological
family" does not mean "identical," and that the ~15% knee-specific gap actually measured for
subject2 is not simply an artifact of the comparison methodology, since the same methodology applied
at the "home" population still shows daylight. The subject2-vs-cadaver-control delta (ratio 1.057 vs
1.150 for the union set, i.e. subject2 sits CLOSER to a match than the cadaver control itself) is
the part attributable to subject2 happening to be a similar body size to the source population, not
evidence of a specially-good fit.

**One caveat on the control point, disclosed, not hidden:** Handsfield's regression was fit on YOUNG
HEALTHY volunteers; applying it at the cadaver population's height/mass implicitly assumes
age-invariance of the height×mass-to-volume relationship, which sarcopenia likely violates (Ward's
donors were presumably older). This makes the cadaver-control comparison an imperfect, secondary
sanity check, not a fully rigorous test — flagged, not treated as decisive. Subject2's own age is not
recorded in `sessionMetadata.yaml` (an honest gap, not fabricated), but the LabValidation cohort this
subject was drawn from is a typical mocap-validation population, presumptively a better age-match to
Handsfield's own young-healthy cohort than 2009-era cadaver donors — this is inference, not confirmed.

## 6. Pre-registered falsifiers — resolved explicitly

- **"Model Fmax for knee-crossing prime movers systematically exceeds PCSA-implied values, beyond
  specific-tension uncertainty, and a corrected Fmax measurably shrinks the over-prediction ratio"**:
  **PARTIALLY CONFIRMED, precisely scoped to the knee.** 12/12 knee-crossing muscles exceed PCSA-implied
  Fmax (mean +14.6%); a decisive causal correction test (not just a threshold comparison) shows a real
  17.1% (subject2) / 10.2% (subject3) drop in the knee ratio — clears this doc's own pre-registered
  >=10% "real effect" bar in both subjects tested. Whether this counts as "beyond specific-tension
  uncertainty" in the strict sense is genuinely borderline (implied ST 60.5-79.3 N/cm^2 sits inside the
  wide defensible literature band) — but since this model's own builders used ONE uniform specific
  tension for all muscles (Arnold 2010, quoted above), not a joint-varying value, there is no
  principled reason internal to the model's own construction to explain away a knee-specific-only
  bias via specific tension. The parsimonious reading is a genuine knee-muscle PCSA/volume
  representation gap, not a specific-tension artifact.
- **"Model Fmax matches PCSA-implied — strength calibration is NOT the source"**: **CONFIRMED for the
  hip** (ratio 1.040 mean, genuinely mixed sign, <0.5% causal effect when corrected in both subjects
  tested) — a clean, valuable negative for that joint. **NOT confirmed for the knee** (see above) —
  an earlier pass through this analysis, using only the coarser "union" aggregate (median 1.057,
  inside the pre-registered [0.7,1.3] "matches" band), would have wrongly reported "matches" for
  BOTH joints; the decisive per-subject correction test caught that this aggregate was too coarse
  before it was reported as final. Recorded here as a disclosed near-miss in this analysis's own
  process, not swept away.
- **Does NOT explain the twin's baseline over-prediction floor.** Even fully corrected, subject2's
  knee ratio (1.255x) remains a 25.5% over-prediction, and subject3's (1.786x) remains 78.6% — so
  something else (SO's minimum-effort recruitment strategy vs real co-contraction, no RRA, moment-arm/
  geometry, or another mechanism already flagged elsewhere in this repo, e.g. `docs/
  MECHANISM_EMG_DRIVEN.md`'s real-EMG co-contraction finding) remains the dominant, unresolved driver
  of the over-prediction floor itself. This analysis narrows, and partially quantifies, one
  contributor — it does not close the whole question.

## 7. Confidence tier and honest gaps

**Confidence tier: cadaveric/MRI published-plausibility, with an added third-party-transcription
ceiling.** Handsfield 2014's volumes are real in-vivo MRI data (a genuine, independent-of-OrthoLoad
external anchor, as the task required) — but:
1. **Single-source risk on the regression table itself.** The primary paper is paywalled (verified
   live, Sec.0); the per-muscle coefficients used here are a third-party transcription
   (`data/external/handsfield2014/PROVENANCE.md`), not independently re-derived from Supplementary
   Table 3. Mitigated (not eliminated) by: exact 80/80 muscle-name match to this model lineage,
   internally-consistent volume-proportion sums, independently-corroborated specific-tension detail,
   and physiologically sane predicted magnitudes — but no second independent transcription was found
   to fully decorrelate this risk.
2. **The multi-compartment volume-split proportions** (e.g. gluteus maximus → glmax1/2/3 at
   0.294759/0.421277/0.283964) are the CSV author's own convention for dividing Handsfield's
   combined-muscle volumes across this model's finer compartments — their derivation (e.g. whether
   drawn from Arnold 2010's own per-compartment cadaver PCSA ratios) was not verified this session,
   an honest, disclosed unknown.
3. **n=1 formal primary-subject test (subject2) + n=3 descriptive cross-subject bonus** (subject2/3/4)
   — sufficient to demonstrate the mechanism and its geometric governor (Sec.3), insufficient to
   claim a population-level magnitude.
4. **Specific tension remains a single scalar choice**, even though this doc used the model-family's
   own documented convention (60 N/cm^2) rather than an arbitrary pick — a different, still-defensible
   choice within [45,100] would shift every ratio uniformly (not the knee-vs-hip asymmetry, which is
   robust to any single global specific-tension choice since it is a JOINT-CONTRAST, not an absolute
   level).
5. **Only the specific-tension-corrected knee/hip re-solve was tested** for subject2 and subject3
   (not subject4, not left side, not other trials) — same scope caveat inherited from every cert in
   this family.
6. **Cadaver-control comparison (Sec.5) has an age-mismatch caveat** (Handsfield=young healthy,
   Ward/Arnold cadavers presumptively older) not resolved here.

## 8. Files

- `scripts/msk/fmax_pcsa_validation.py` — Part 1-3 pipeline (load/verify, per-muscle PCSA-implied
  Fmax at 4 anthropometric points, aggregate). Real run, exit 0.
- `scripts/msk/fmax_pcsa_correction_test.py` — the decisive subject3 causal correction test (builds
  corrected model, reruns SO+JR, compares vs cached baseline). Real run, exit 0.
- `scripts/msk/fmax_pcsa_correction_test_subject2.py` — the same test for subject2 (this task's
  primary subject), reusing `build_corrected_model` from the script above and
  `subject_specific_scaling.py`'s `compute_baseline_fresh()`/`run_joint_forces_for_model`. Real run,
  exit 0.
- `data/external/handsfield2014/muscle_scaling_coefficients.csv` — the Handsfield 2014 regression
  table (third-party transcription, see provenance file).
- `data/external/handsfield2014/PROVENANCE.md` — full citation correction + source-verification
  chain for the CSV above.
- `data/msk_smoketest/subject2_walking1/fmax_pcsa_validation/fmax_pcsa_validation_results.json` —
  every number in Sec.1-3, full 80-muscle x 4-subject detail.
- `data/msk_smoketest/subject2_walking1/fmax_pcsa_correction_test/fmax_pcsa_correction_test_subject2_results.json`
  — subject2 correction-test full detail (Sec.4).
- `data/msk_smoketest/subject3_walking1/fmax_pcsa_correction_test/fmax_pcsa_correction_test_results.json`
  — subject3 correction-test full detail (Sec.4).
- `data/msk_models/subject2_scaled_handsfield_fmax_corrected.osim`,
  `data/msk_models/subject3_scaled_handsfield_fmax_corrected.osim` — the two corrected models built
  and tested (new files, Fmax overwritten only, all other parameters untouched).
- Reused, unedited: `scripts/msk/validate_joint_force.py`, `scripts/msk/static_opt_knee.py`,
  `scripts/msk/validate_hip_force.py`, `scripts/msk/subject_specific_scaling.py`,
  `scripts/msk/cross_subject_validation.py`.
- Prior docs this builds on: `docs/MECHANISM_CROSS_SUBJECT.md`, `docs/MECHANISM_MUSCLE_AUDIT.md`,
  `docs/MECHANISM_SUBJECT_SPECIFIC.md`, `docs/MECHANISM_CMC_SECOND_SOLVE.md`.

Isolation respected throughout: bodytwin only; LabValidation model/session data read in place on the
read-only `mechanism_data` NTFS drive; external web fetches were read-only (PubMed, NCBI eutils,
Semantic Scholar, Unpaywall, GitHub raw content) with no credentials and no write access. No git
commit, no git push, no git add performed.
