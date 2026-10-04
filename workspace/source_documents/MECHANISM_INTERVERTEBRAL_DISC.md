# MECHANISM INTERVERTEBRAL DISC — nucleus pulposus pressure, disc stiffness, and creep: a first falsifiable compressive-mechanics model (2026-07-22)

**Question.** `docs/MECHANISM_SPINE_GAIT_VBR.md` and `docs/MECHANISM_SPINE_FORCE.md` compute lumbar
compressive **force** (a lumped bone/facet/disc contact residual) but treat the disc as a rigid
contact — no pressure, no viscoelasticity, no height loss. This build adds the first layer of
disc-specific physics: **P_NP = k·(F/A)** (Nachemson's nucleus-pressure/compressive-stress relation),
a disc compressive-stiffness → instantaneous height-loss model, and a first creep (sustained-load
height-loss) model — each anchored to live-fetched, PMID-cited literature, and used to run a
genuine, pre-registered falsifier against the existing gait cert's own already-machine-measured
walking force.

## Headline result

| | value | source |
|---|---:|---|
| **Falsifier: predicted P_NP at twin's walking Tier-2b peak ∈ [0.53, 0.65] MPa (Wilke 1999 in-vivo walking)?** | **FAIL** (1.18–1.36 MPa, k=1.3–1.5) | this build, machine-measured |
| Ratio over Wilke's upper bound (k=1.5 central) | **2.09×** | |
| **Forced adversary** (generous k, BW-scaled area, even an extreme literature-unsupported combo) | still FAILS (0.65 MPa vs 0.65 MPa boundary — the *only* combination that doesn't fail requires k<1, contradicted by every cited source) | this build |
| **Refined, fraction-free over-prediction estimate** (via Wilke's pressure, sidesteps the VBR-fraction confound) | **1.8–2.6×** — tighter than the existing doc's raw VBR-force ratio | this build, new |
| Existing doc's raw VBR-force ratio (`MECHANISM_SPINE_GAIT_VBR.md`) | 3.94–5.39× | prior cert |
| Refined implied VBR load-share fraction (own-median / Damm2017) | **0.46–0.65 / 0.34–0.48** | this build, new — more biomechanically plausible than the existing doc's naive 0.19–0.25 |
| Disc compressive stiffness (human, ex-vivo) | 1734 ± 446 N/mm (1.73±0.45 MN/m) | Beckstein 2008, PMID 18344845 |
| Twin's own instantaneous elastic height loss at walking peak | 0.94 mm (0.75–1.26 mm) | this build — matches literature "step displacement" 0.90±0.12 mm |
| Empirical Nachemson factor k (back-derived from an independent MSK-modeling paper) | **1.524** (1.6% above the task's stated 1.3–1.5 prior) | this build, cross-check |

**Confidence tier: in-vivo-anchored (Wilke intradiscal pressure) for the pressure-load sub-model.**
The stiffness and creep sub-models are **in-vitro/ex-vivo-anchored** (human cadaveric disc testing)
— a genuinely different, lower tier, stated explicitly, not blended into the headline tier.

The falsifier **FAILS as stated** — but this is not a new defect: it is a **decorrelated
confirmation**, via a completely different in-vivo measurement modality (an implanted pressure
transducer, not a VBR-implant force telemetry), of the **same** SO-driven over-prediction mechanism
`MECHANISM_SPINE_GAIT_VBR.md` already documented (§6 there). Because Wilke's pressure anchor is
**fraction-free** (a transducer inside the nucleus reads the disc's true total load directly, with no
VBR-device load-sharing ambiguity), it **tightens** the existing raw 3.9–5.4× ratio down to a
disclosed 1.8–2.6×, and **refines** the existing doc's own open VBR-fraction question toward a more
biomechanically plausible ~46–65% (vs. the previous, fraction-confounded 19–25%) — §6 below.

---

## 1. Pre-registration (stated before running the falsifier)

- **Falsifier**: predicted nucleus-pulposus pressure at the twin's own walking Tier-2b
  (muscle+ligament-subtracted, PRIMARY) gait-cycle peak force (1628.10 N,
  `data/msk_smoketest/subject2_walking1/spine_gait_force/spine_gait_force_results.json`), converted
  via `P = k·(F/A)`, must fall inside Wilke et al. 1999's measured in-vivo walking range
  **[0.53, 0.65] MPa** for **both** ends of the task's stated k range (1.3 and 1.5) to PASS.
- **Adversary I am tempted to skip** (this claim leans negative — predicted pressure over-shoots):
  the strongest fair rescue is "the conversion PARAMETERS (k, A) are wrong in the PASS-favoring
  direction," not "the twin's force estimate is wrong." Forced via a sweep in §5: even a
  literature-unsupported, maximally generous parameter combination is required to approach a PASS,
  and even that combination lands almost exactly on the boundary, not comfortably inside it.
- **Symmetric-QC target** (task's own instruction): verify Nachemson's ~1.3–1.5× pressure/stress
  factor **direction** before using it — a factor error in the wrong direction would flip every
  downstream conclusion. Resolved in §3 via an independent empirical back-derivation from an
  unrelated paper's own reported numbers (not from memory/recall).
- **Anchor discipline**: every constant below is cited with a PMID, fetched live this session via
  NCBI E-utilities (`eutils.ncbi.nlm.nih.gov`) and the Europe PMC REST API — never a tautology gate,
  never recalled without a citation.

## 2. Method — the geometry (not rote algebra)

The nucleus pulposus behaves as a confined, pressurized, quasi-hydrostatic gel: axial compressive
load F applied over the disc's cross-sectional area A does not distribute as a uniform stress F/A —
curved vertebral endplates concentrate load centrally (where the nucleus sits), and the annulus
carries a disproportionate share via hoop tension rather than direct axial stress. Nachemson's
classic in vitro calibration (cadaveric discs under known load, nucleus pressure measured directly)
found the nucleus pressure runs **higher** than the naive average stress by a factor k ≈ 1.3–1.5:

```
P_NP = k · (F / A)                      k ≈ 1.3–1.5  (Nachemson; verified §3, not assumed)
```

This is the SAME geometric relationship used, in the reverse direction, to estimate compressive
FORCE from a measured in-vivo pressure (`F = P·A/k`) — the method `docs/MECHANISM_SPINE_GAIT_VBR.md`'s
own VBR comparison could not use (no device-specific load-share fraction exists, §6 there) but that
Wilke's direct nucleus-pressure measurement sidesteps entirely (§6 below).

Two further sub-models, each geometrically distinct:
- **Instantaneous (elastic) compression**: `Δh = F / k_disc`, disc treated as a linear spring at the
  scale of a single loading event (a first-order simplification — real IVD stiffness is load-
  dependent/nonlinear, disclosed in §8).
- **Creep (time-dependent, under SUSTAINED load)**: a stretched-exponential (KWW) relaxation,
  `x(t) = d_∞·(1 − exp(−(t/τ)^β))` — the standard viscoelastic form for IVD creep, reflecting a
  superposition of at least two relaxation processes at different rates (fast: nucleus/endplate
  fluid flow; slow: annulus fibrosus), not a single arbitrary curve fit.

## 3. Literature anchors (all PMIDs live-verified this session, verbatim quotes preserved)

| # | citation | PMID | what it gives |
|---|---|---|---|
| 1 | Wilke HJ, Neef P, Caimi M, Hoogland T, Claes LE. "New in vivo measurements of pressures in the intervertebral disc in daily life." *Spine*. 1999;24(8):755-62. | **10222525** | N=1, 45yo, 70kg, non-degenerated L4-L5, telemetric transducer, ~24h. Verbatim: relaxed standing 0.5 MPa; standing flexed forward 1.1; sitting unsupported 0.46; sitting max flexion 0.83; lifting 20kg round-back 2.3, flexed-knees 1.7, close-to-body 1.1; lying prone 0.1; night 0.1→0.24. |
| 2 | Wilke HJ, Neef P, Hinz B, Seidel H, Claes L. "Intradiscal pressure together with anthropometric data — a data set for the validation of models." *Clin Biomech*. 2001;16(Suppl 1):S111-26. | **11275349** | SAME subject/disc as #1 (reproduces relaxed-standing = 0.5 MPa exactly — internal consistency confirmed). Anthropometric companion dataset; disc cross-sectional area **A = 1800 mm²** cited from this dataset by an independent downstream paper (#8). |
| 3 | Sato K, Kikuchi S, Yonezawa T. "In vivo intradiscal pressure measurement in healthy individuals and in patients with ongoing back problems." *Spine*. 1999;24(23):2468-74. | **10626309** | Independent (different lab/subjects/instrument), n=8. Computed spinal LOAD directly: prone 144N, lateral 240N, **standing 800N**, sitting 996N (p<0.0001). |
| 4 | Nachemson A, Morris JM. "In Vivo Measurements of Intradiscal Pressure." *J Bone Joint Surg Am*. 1964;46:1077-92. (via secondary review PMC9000064) | 14193834 | n=8. NP pressure: prone/supine 0.091–0.539 MPa; seated 0.46–1.33; **standing 0.5–0.87 MPa**. |
| 5 | Nachemson A, Elfcurrent G. "Intravital dynamic pressure measurements in lumbar discs..." *Scand J Rehabil Med Suppl*. 1970;1:1-40. | 4257209 | Confirmed real citation (monograph, pre-index, no abstract available). The original "walking/movements" dynamic-pressure study. |
| 6 | Dreischarf M, Rohlmann A, Zhu R, Schmidt H, Zander T. "Is it possible to estimate the compressive force in the lumbar spine from intradiscal pressure measurements? A finite element evaluation." *Med Eng Phys*. 2013;35(9):1385-90. | **23570899** | Tests the `P=k·(F/A)` method itself against FE ground truth: model-specific k gives <4% error (except extension, >27%); **generic/non-model-specific k gives up to 44% error**; in-vivo standing force from IDP estimated at 430–600N depending on k choice. Directly informs this build's uncertainty budget. |
| 7 | Roman-Liu D, Kamińska J, Tokarski T. "Differences in lumbar spine IDP between standing and sitting postures: a comprehensive literature review." *PeerJ*. 2023;11:e16176. | 37872945 (PMC10590571, open access, full text fetched) | Compiles Nachemson/Sato/Wilke *relative* pressure data; used for cross-orientation, did not itself contain the walking absolute value or the k factor. |
| 8 | "The effect of walking speed on spinal loads and trunk muscle forces using subject-specific musculoskeletal modeling..." | (PMC12584331, open access, full text fetched) | **Independent** MSK-model estimate: L4-L5 peak compression ≈**638N (92.9%BW)** during normal-speed walking, 70.4kg avg subject; "assuming a disc cross-sectional area of 1800 mm² [Wilke 2001]" → derived IDP ≈**0.54 MPa**; explicitly states this "fell within the measured intradiscal pressure range (**0.53–0.65 MPa**) reported by Wilke et al. 1999 ... during walking." → the walking-range figure the task asked for, cross-corroborated by a second, independent paper (below) and used to back-derive k empirically (§3.1). |
| 9 | "Viewing low back pain through the lens of spinal evolution..." | (PMC12810782, open access, full text fetched) | Independently states, in a paragraph discussing Wilke et al.: "Walking yields **0.53–0.65 MPa**." Second, decorrelated confirmation of the same range. |
| 10 | Adams MA, Hutton WC. "The mechanical function of the lumbar apophyseal joints." *Spine*. 1983;8(3):327-30. | **6623200** | Verbatim: facet joints "share in resisting the intervertebral compressive force, **but only in lordotic postures**." In flexion (this gait trial's own measured posture throughout, `MECHANISM_SPINE_GAIT_VBR.md` §1), facets do not meaningfully share axial compression → supports treating ≈100% of the Tier-2b residual as disc-transmitted here. |
| 11 | Beckstein JC, Sen S, Schaer TP, Vresilovic EJ, Elliott DM. "Comparison of animal discs used in disc research to human lumbar disc: axial compression mechanics and glycosaminoglycan content." *Spine*. 2008;33(6):E166-73. | **18344845** | HUMAN lumbar disc, ex-vivo, standardized protocol: compressive stiffness **1734±446 N/mm**; step (instantaneous) displacement 0.90±0.12mm; creep displacement 0.55±0.03mm; ROM 1.21±0.18mm. (Secondary review PMC9000064 in-text calls this "Jesse et al." — the matching reference-list entry is Beckstein et al.; cited by the verified reference-list entry, the likely-safer of the two, discrepancy disclosed not silently resolved.) |
| 12 | O'Connell GD, Jacobs NT, Sen S, Vresilovic EJ, Elliott DM. "Axial creep loading and unloaded recovery of the human intervertebral disc and the effect of degeneration." *J Mech Behav Biomed Mater*. 2011;4(7):933-42. | **21783103** | HUMAN lumbar discs, ex-vivo: **1000N creep load for 4h**, unloaded recovery 24h. Two-component viscoelastic split of equilibrium displacement: fast (NP+endplate) = 10–15%; slow (annulus fibrosus) = 40–70%. Recovery is **3–4× slower** than loading. τ≈14h (via secondary review citing this paper); >80% of equilibrium within 4h; full equilibrium ~12h. |
| 13 | Tyrrell AR, Reilly T, Troup JD. "Circadian variation in stature and the effects of spinal loading." *Spine*. 1985;10(2):161-4. | **4002039** | Whole-spine (NOT per-disc) context: mean circadian stature loss 19.3mm (1.1% of stature); **54% occurs in the first hour** after rising; static shoulder loads (2.5–40kg) → nonlinear shrinkage-rate increase. |
| 14 | Dowzer CN, Reilly T, Cable NT. "Effects of deep and shallow water running on spinal shrinkage." *Br J Sports Med*. 1998;32(1):44-8. | **9562163** | Whole-spine context: 30-min treadmill running → 4.59±1.48mm stature loss. |
| 15 | Rohlmann A, Dreischarf M, Zander T, Graichen F, Bergmann G. "Loads on a vertebral body replacement during locomotion measured in vivo." *Gait Posture*. 2014;39(2):750-5. | 24211089 | Already used in `MECHANISM_SPINE_GAIT_VBR.md`; reused here unchanged: VBR walking force = 171% of standing. |

All PMIDs above were resolved live this session via `eutils.ncbi.nlm.nih.gov/entrez/eutils/e{search,summary,fetch}.fcgi` and `www.ebi.ac.uk/europepmc/webservices/rest/` (full-text XML fetched for open-access papers #7–9), not recalled from training memory — WebSearch quota was exhausted for this session/account (a pre-existing, shared-account constraint also logged in `MECHANISM_SPINE_FORCE.md`'s own honest gaps), so literature discovery this session used direct PubMed/Europe PMC REST calls instead, which is at least as strong a provenance chain (machine-parseable JSON/XML, not a search-snippet paraphrase).

### 3.1 Verifying the Nachemson factor's direction and magnitude (task's explicit ask)

The task asked to verify "nucleus pressure ~1.3–1.5× the applied axial stress" rather than assume
it — a real risk, since a reciprocal misstatement (pressure = stress/1.5) would flip every
downstream conclusion. Source #8 (an independent MSK-modeling paper, unconnected to this build)
reports its OWN three numbers together: **F=638N, A=1800mm², derived P=0.54 MPa**. Back-solving:

```
k_empirical = P / (F/A) = 0.54 / (638/1800) = 1.524
```

This **confirms both the direction** (pressure exceeds nominal stress, k>1 — not the reciprocal) **and
the magnitude** (1.524, 1.6% above the task's stated 1.3–1.5 upper bound, i.e. consistent with the
prior) — from an independently-computed, real paper's own numbers, not from memory. Machine-checked:
`k_empirical_within_task_prior_range_strict: false` (1.524 > 1.5 strictly) but
`k_empirical_pct_above_task_prior_upper_bound_1.5: 0.0157` — a 1.6% overshoot, reported plainly, not
padded away.

A second, independent cross-check (Sato #3 vs. Nachemson #4, two different primary in-vivo studies,
neither aware of the other's exact numbers): Nachemson's standing pressure range (0.5–0.87 MPa,
midpoint 0.685) converted via the same k=1.5, A=1800mm² gives an implied load of 822N — Sato's own,
independently measured/computed standing load is 800N. **2.75% difference** (well inside a 15%
tolerance, the same convention this repo's cert family already uses elsewhere) — a genuine
decorrelated consistency check between two 1960s–1990s independent measurement programs.

## 4. The falsifier, machine-computed

Reusing (not retyping) this repo's own already-machine-measured gait-cert force
(`spine_gait_force_results.json`, Tier-2b PRIMARY, 158-frame gait cycle, t=0.09s peak):

| quantity | value |
|---|---:|
| F_walk (Tier-2b, PRIMARY) | 1628.10 N |
| A_disc | 1800 mm² |
| Nominal stress F/A | 0.9045 MPa |
| **Predicted P_NP, k=1.5** | **1.357 MPa** |
| Predicted P_NP, k=1.3 | 1.176 MPa |
| Predicted P_NP, k=1.524 (empirical) | 1.378 MPa |
| Wilke in-vivo walking range | **0.53–0.65 MPa** |
| Ratio over upper bound (k=1.5 / k=1.3) | **2.09× / 1.81×** |
| **Falsifier verdict** | **FAIL** (both k=1.3 and k=1.5 land above the range) |

## 5. Forced adversary — does the FAIL survive generous, fair parameter choices?

Per the discipline, a leaning-negative claim must force the adversary to its strongest fair form
before being accepted. The adversary here: "the FAIL is a parameter artifact (wrong k or A), not a
genuine finding about the twin's force estimate."

| adversary move | result | still fails? |
|---|---:|---|
| k needed to hit Wilke's upper bound (0.65 MPa), holding A=1800 | **k=0.719** — below 1, contradicted by every one of 5 independent sources found (all say k>1) | adversary FALLS |
| A needed to hit 0.65 MPa, holding k=1.5 | **A=3757 mm²** — ~2× any cited lumbar disc area (typically 1200–2000mm²) | adversary FALLS |
| A scaled up for subject2's greater body mass vs. Wilke's subject (78.98/70=1.117×) → A=2031mm², k=1.5 | 1.203 MPa | still FAILS (1.85×) |
| Same BW-scaled area, k at the task's own low end (1.3) | 1.042 MPa | still FAILS (1.60×) |
| **Extreme, literature-unsupported combo** (A=2500mm², well above any citation; k=1.0, contradicted by every source) | **0.6512 MPa** | **still FAILS — but only by 0.2%, essentially at the boundary** |

Even the most generous **fair** combination (BW-scaled area + task's own low-end k) still over-predicts
by 1.6×. Only an **extreme, literature-unsupported** combination (k below every cited value, area
roughly double every citation) approaches the boundary, and even that combination does not cross it.
**Verdict: the FAIL is robust to parameter uncertainty — it is not a conversion-constant artifact.**

## 6. Orienting the FAIL (OODA, not a one-shot "honest negative") — a real, useful refinement

A bare "FAIL, honest negative" would be premature surrender. Orienting on *why*: `MECHANISM_SPINE_GAIT_VBR.md`
already diagnosed and disclosed (its own §6, §10 gap #1) that this **exact** Tier-2b force
(1628.10N) over-predicts OrthoLoad's VBR-implant walking telemetry by a **raw** 3.94–5.39× — and
flagged, as an open, unresolved question, that this raw ratio conflates two unknowns: (a)
architecture-driven over-prediction (same mechanism as knee/hip, 1.4–2.0×) and (b) the VBR device's
own unknown anterior/posterior load-share fraction. That doc could not decompose the two (a
dedicated literature search found no VBR-device-specific fraction).

**This build's pressure anchor is fraction-free** — Wilke's transducer sits inside the nucleus and
reads the disc's real total compressive load directly, with no VBR-device load-sharing ambiguity to
disentangle. That lets this build compute an independent estimate of the *true* total in-vivo
compressive force during walking, `F_true = P_Wilke·A/k`, and compare it directly to the twin's
own Tier-2b estimate:

| | value |
|---|---:|
| F_true, implied from Wilke's walking pressure range (k=1.3–1.5) | **636–900 N** |
| F_twin (Tier-2b) | 1628.1 N |
| **Pressure-derived over-prediction ratio** | **1.81–2.56×** |
| (vs. existing doc's raw, fraction-confounded VBR ratio) | 3.94–5.39× |

This is a **materially tighter** estimate of the pure over-prediction factor than the existing raw
VBR ratio — because it removes the unknown VBR-fraction confound entirely. Using it to REFINE
(not replace) the existing doc's own "implied VBR fraction" question:

| | naive (existing doc, assumes Tier-2b = true force) | refined (this build, fraction-free F_true) |
|---|---:|---:|
| Implied fraction, own-median OrthoLoad (413.4N) | 0.254 | **0.459–0.650** |
| Implied fraction, Damm 2017 (302.1N) | 0.186 | **0.336–0.475** |

A VBR load-share fraction of ~46–65% is considerably more biomechanically plausible for an
**anterior-column** load-bearing device (which mechanically replaces the vertebral body/disc space
itself) than the naive 19–25% — clinical VBR+posterior-instrumentation constructs are generally
understood to leave the anterior column carrying a substantial share, not a small minority, of axial
load. **This does not fully resolve `MECHANISM_SPINE_GAIT_VBR.md`'s open question (no independent
VBR-device-specific fraction was found this session either — same disclosed gap), but it narrows it
with a genuinely new, decorrelated, fraction-free data point.**

**Bracket sanity** (Tier-1, pure kinematic, SO-independent): F=433.47N → predicted P=0.313–0.361 MPa
(k=1.3–1.5), which **under**-shoots Wilke's range (expected: Tier-1 omits SO-modeled muscle
co-contraction, which physically ADDS compression). The true in-vivo value (implied 636–900N) sits
between Tier-1 (433N) and Tier-2b (1628N), closer to Tier-1 — consistent with "SO substantially
over-adds co-contraction compression, pure kinematics substantially under-adds it," the same
mechanism already named (not newly discovered) in the existing gait cert's §6.

## 7. Secondary, decorrelated cross-check: stoop-lift (a third activity)

`MECHANISM_SPINE_FORCE.md`'s own static 10kg stoop-lift Tier-2 force (2984.50N with_box, 2274.46N
no_box) converted the same way (k=1.5): **2.49 MPa / 1.90 MPa**. Wilke's 20kg-lift range is
1.1–2.3 MPa — the twin's 10kg prediction sits **at/above** the range for a **20kg** lift (half the
external load). Caveat stated plainly: 10kg vs 20kg is not a clean match, so this is directional
context, not a second pass/fail gate — but it shows the **same** qualitative over-prediction
direction on a third, independent activity (lifting, not gait), reinforcing §6's diagnosis rather
than contradicting it.

## 8. Disc compressive stiffness → instantaneous height loss

Using Beckstein 2008's human ex-vivo stiffness (1734±446 N/mm) as a **first-order linear spring**
(disclosed simplification — real IVD stiffness is load-dependent/nonlinear; not modeled here):

| scenario | F (N) | Δh, central (mm) | Δh, range (mm) |
|---|---:|---:|---:|
| Walking, Tier-2b peak | 1628.1 | **0.94** | 0.75–1.26 |
| Stoop-lift, with_box | 2984.5 | 1.72 | 1.37–2.32 |
| Stoop-lift, no_box | 2274.5 | 1.31 | 1.04–1.77 |

The twin's own computed elastic height loss at its walking peak (0.94mm central) lands within the
independently-measured human "step displacement" literature figure (0.90±0.12mm, same paper) —
**a coincidental-but-supportive order-of-magnitude match**, explicitly NOT claimed as an independent
proof (the literature figure's own applied load was not confirmed identical to 1628N in the
extracted text).

## 9. Creep model — a first falsifiable form, explicitly lower-confidence than §4–6

`x(t) = d_∞·(1 − exp(−(t/τ)^β))` (KWW stretched exponential, standard for IVD creep). Primary
anchor: O'Connell 2011 (PMID 21783103), 1000N/4h human ex-vivo protocol, τ≈14h (same paper, via
secondary review), fast(10–15%)/slow(40–70%) equilibrium-displacement split, recovery 3–4× slower
than loading.

**A real inconsistency, forced and disclosed, not hidden**: the literature also states ">80% of
equilibrium reached within 4h" (same secondary review, same underlying paper). Solving for the
stretch parameter β that would make BOTH `τ=14h` and `f(4h)=0.80` hold under the KWW form gives
**β≈−0.85 — unphysical** (β must be in (0,1]). Machine-checked, not assumed:

```
fraction_at_4h_with_tau14_beta1  = 0.249   (NOT 0.80)
fraction_at_12h_with_tau14_beta1 = 0.576
fraction_at_24h_with_tau14_beta1 = 0.820
```

**Conclusion: these two literature facts come from different fitted protocols and must NOT be
force-combined into one curve** — reported as a genuine, disclosed lower-confidence gap in this
sub-model, not silently reconciled. An illustrative (simple-exponential, τ=14h, β=1 — likely
**understating** early-time creep, since the true response is faster-initial per the literature's
own >80%@4h claim) curve, applying the Beckstein 2008 creep-displacement magnitude (0.55mm) on top
of the instantaneous step (0.90mm):

| t | fraction (τ=14h, β=1) | creep-only Δh (mm) | total Δh incl. step (mm) |
|---:|---:|---:|---:|
| 1h | 0.069 | 0.038 | 0.938 |
| 4h | 0.249 | 0.137 | 1.037 |
| 8h | 0.435 | 0.239 | 1.139 |
| 12h | 0.576 | 0.317 | 1.217 |
| 24h | 0.820 | 0.451 | 1.351 |

This table is **illustrative first-order context**, not a validated per-disc in-vivo creep curve —
flagged explicitly as such. Whole-spine (not per-disc) in-vivo creep context (Tyrrell 1985: 19.3mm
diurnal stature loss, 54% in the first hour; Dowzer 1998: 4.59mm after 30min treadmill running)
confirms the same qualitative fast-initial/slow-tail shape at a different anatomical scale, but was
**not apportioned to a single lumbar disc** (no cited apportionment fraction found this session) —
disclosed as a real limitation, not silently extrapolated.

**Confidence tier for this sub-model: in-vitro/ex-vivo-anchored (human cadaveric disc testing) —
explicitly LOWER than the in-vivo-anchored pressure-load sub-model (§4–6).** A genuine future
falsifier for this sub-model specifically: if a per-disc in-vivo (not cadaveric) creep measurement
becomes available (e.g. MRI-based disc-height tracking over a work shift), check whether the
modeled curve's magnitude at t=4h/8h/24h matches — not done this session, forward-looking and
explicitly flagged as such.

## 10. Machine-checked gates

| gate | threshold | measured | verdict |
|---|---|---:|---|
| Headline falsifier | P_NP(walk, Tier-2b) ∈ [0.53,0.65] MPa, both k=1.3 and k=1.5 | 1.176–1.357 MPa | **FAIL** |
| Forced adversary (extreme combo) still fails† | P > 0.65 MPa even at A=2500mm², k=1.0 | 0.6512 MPa | **PASS**† |
| Tier1/Tier2b bracket sane | Tier1 < Wilke_lo < Tier2b(k=1.3) | 0.313 < 0.53 < 1.176 | **PASS** |
| Sato-vs-Nachemson standing cross-check | \|diff\| < 15% | 2.75% (800N vs 822N) | **PASS** |
| Empirical k vs. task prior | report exact overshoot, no padding | k=1.524, +1.57% over 1.5 | **reported (marginal, disclosed)** |
| Script reproducibility | bit-for-bit across 2 independent runs | identical | **PASS** |
| **Overall** | falsifier is the headline claim; robustness/consistency gates support the diagnosis | | **FALSIFIER FAILS (robustly); diagnosis PASSES; net contribution: a tighter, fraction-free over-prediction estimate** |

†The adversary-row "PASS" means the *robustness check* passed, i.e. it correctly confirms the FAIL
is not a parameter artifact (the extreme adversary still could not rescue a PASS on the headline
falsifier) — it is NOT a claim that the headline falsifier itself passed. Flagged explicitly to
avoid the two senses of "PASS" being conflated.

## 11. Honest gaps (full list)

1. **The headline falsifier FAILS** — the twin's Tier-2b walking force, converted to pressure,
   over-predicts Wilke's in-vivo walking range by ~1.8–2.6× (fraction-free estimate) to ~2.1× (raw,
   k=1.5). This is the primary, load-bearing finding of this build, not softened.
2. **Disc area (1800mm²) is Wilke's own subject's value** (70kg individual), not subject2's own
   measured disc area (no CT/MRI on subject2 exists in this repo) — a BW-scaling sensitivity check
   was run (§5) but a subject-specific disc area was not independently obtained.
3. **The Tier-2b "bone/facet contact" residual is assumed ≈100% disc-transmitted** in this trial,
   supported by Adams & Hutton 1983's qualitative finding (facets share compression only in
   lordosis, this trial stays in flexion) but not by a precise, trial-specific facet-engagement
   measurement (this model has no separately-modeled facet joints at all — inherited limitation,
   same as `MECHANISM_SPINE_GAIT_VBR.md` gap #5).
4. **The k=1.3–1.5 Nachemson factor itself carries real uncertainty**: Dreischarf 2013 (PMID
   23570899) shows generic (non-model-specific) correction factors can carry up to 44% error, and
   model-specific factors, while accurate (<4%) for most postures, degrade badly (>27% error) in
   extension specifically — this gait trial stays in flexion (not extension), so the better-behaved
   regime, but this is a directional argument, not an independently re-verified error bound for
   THIS specific model/posture.
5. **Creep sub-model parameters (τ, β, d_∞) are NOT jointly fit to a single consistent curve** —
   disclosed and forced in §9, not hidden. The illustrative curve provided is first-order context.
6. **Creep sub-model is entirely ex-vivo/cadaveric-anchored**, no in-vivo per-disc human creep
   measurement was found/used this session (in-vivo diurnal-stature data exists but is whole-spine,
   not per-disc — §9).
7. **WebSearch quota was exhausted for this session/account** (pre-existing, shared-account
   constraint, also logged in `MECHANISM_SPINE_FORCE.md`'s own gaps) — literature discovery used
   direct NCBI E-utilities / Europe PMC REST calls instead (arguably a stronger provenance chain,
   machine-parseable, but a different discovery method than a semantic web search would have been).
8. **Beckstein 2008's stiffness/displacement numbers are not tied to one specific applied load** in
   the text extracted this session (only the paper's scope/abstract were independently verified;
   the specific numbers came via a secondary review) — the "coincidental match" noted in §8 is
   explicitly flagged as not independently load-matched.
9. **Single lumped lumbar joint, single trial, single subject** — same scope inherited from every
   spine build in this repo (`MECHANISM_SPINE_GAIT_VBR.md` gaps #6, #9).
10. **The "Jesse et al." vs "Beckstein et al." attribution discrepancy** in the secondary review's
    own in-text citation vs. its reference list (§3, source #11) — resolved by citing the verified
    reference-list entry and disclosing the discrepancy, not silently picking one without comment.
11. **subject2's mass appears as both "78.2 kg" (prose in the two existing spine docs) and
    "78.976 kg" (this session's live JSON re-read of `total_mass_kg`)** — a small (~1%), likely
    rounding/provenance difference; this build used the JSON's own machine value (internally
    consistent with its own `body_weight_N`), not the prose figure — disclosed, not silently
    reconciled.

## 12. Files

- `scripts/msk/model_disc_pressure.py` — the full, self-contained model (stdlib only, no OpenSim/
  `.venv-msk` dependency; re-runnable with plain `python3`). Every literature constant is a named,
  commented module-level variable with its PMID; reads (does not retype) this repo's own existing
  gait-cert and stoop-lift-cert JSON outputs.
- `data/msk_smoketest/disc_pressure_model/disc_pressure_model_results.json` — every number in this
  document, machine-written, bit-for-bit reproducible across 2 independent runs (verified this
  session).
- Reused, unedited: `data/msk_smoketest/subject2_walking1/spine_gait_force/spine_gait_force_results.json`,
  `data/msk_smoketest/subject2_spine_stoop_lift/spine_force_validation_results.json`.
- Compared against: `docs/MECHANISM_SPINE_GAIT_VBR.md` (the walking gait force cert this build's
  falsifier targets directly, and whose own §6 open question this build refines), `docs/
  MECHANISM_SPINE_FORCE.md` (the static stoop-lift cert used for the secondary cross-check).
- No git commit, no git push performed (isolation respected, per `COORDINATOR.md` §1 and the task).
