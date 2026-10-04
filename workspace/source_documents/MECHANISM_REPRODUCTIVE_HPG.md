# MECHANISM MALE HPG (hypothalamic-pituitary-gonadal) STEADY-STATE AXIS (2026-07-22)

Builds and MEASURES the male HPG negative-feedback axis: GnRH pulses -> LH/FSH -> testosterone
(T), with T's dual-site negative feedback (hypothalamic pulse-frequency + pituitary
responsiveness) closing the loop, as a CERTIFIED model — resolving the pre-existing
**SEED-DESIGN** (status `OPEN`, "designed, not measured") hypothesis nodes
`MSK-HPG-RATE-CONSTANT-LAYER` / `AUTO-BUILD-AN-HPG-RATE-CONSTANT-LAYER-A-PER-S` (both explicitly
name a required state-vector input "serum ... T ... FSH", and — via
`INT-REDS-ENERGY-DEFICIENCY` — "LH-pulse-frequency(t)") already sitting in
`data/MECHANISM_ANCHOR_GRAPH.json` (998 nodes; **not edited this session** — isolation rule
"touch only files you create"; folding this result into those nodes via the canonical
`mechanism_fold -> fold_gate_v2` path is the natural next step, not performed here, matching the
identical precedent set by `MECHANISM_THYROID_AXIS.md`). Those nodes **assume** a T/FSH/LH-
frequency state vector as an input; this script is the missing upstream mechanistic layer that
actually *produces* it, rather than treating it as a free/measured-only quantity.
Script: `scripts/msk/hpg_male_axis.py`. Evidence:
`data/msk_smoketest/subject2_walking1/hpg_male_axis/hpg_male_axis_results.json`.

**Why male, not female (task said "pick the cleaner falsifier" — an explicit scoping choice).**
The female cycle's periovulatory LH surge is a **pituitary/ovarian bistable threshold-switch**
triggered by sustained-elevated estradiol — Knobil et al. 1980 (*Science*) and Wildt et al. 1981
(*Endocrinology*, "Estradiol as a gonadotropin releasing hormone") both show a **fixed,
unvarying** 1-pulse/h GnRH replacement still permits a normal ovulatory cycle+surge in rhesus
monkeys, and estradiol **alone** (no GnRH pulse-pattern change at all) can trigger a gonadotropin
discharge. The surge mechanism is therefore **decorrelated from GnRH pulse frequency per se** —
a genuinely different (bistable-switch, threshold-not-sharply-defined — exactly this task's own
flagged caveat) mathematical object than a negative-feedback steady state. The male axis is a
clean, non-bistable, delayed-negative-feedback steady state with a rich, decades-deep,
quantitative primary literature. **Female surge dynamics are explicitly OUT OF SCOPE here** —
held open, not modeled, not hand-waved.

## The falsifiers, verdicts stated up front (nothing hidden)

> **Falsifier 1** (diurnal T rhythm): does Veldhuis et al. 1987's own circadian-amplitude/
> implied-mean numbers, applied with zero extra fitting, predict a T peak-to-nadir swing
> consistent with Diver et al. 2003's INDEPENDENTLY, directly measured 43% swing (different
> cohort, decade, and method)?

**PASSES**: 46.9% (Veldhuis-derived) vs 43.0% (Diver-measured) — 9.1% relative difference,
comfortably inside the pre-registered 25% threshold. See §3.

> **Falsifier 2** (GnRH/LH pulse frequency -> mean T, a SATURATING not linear relationship): does
> a model with a frequency-dependent pituitary-desensitization term reproduce (a) Spratt et al.
> 1987's finding that T stayed ~constant despite an 8x GnRH frequency increase, AND (b)
> Finkelstein et al. 1988's finding (a fully independent, decreasing-frequency protocol) that T
> fell substantially as frequency was decreased — while the forced naive-linear adversary is
> falsified by (a)?

**PASSES, with an honestly disclosed constructive-choice caveat**: model change across the
Spratt range = +9.6% (gate <20%, PASS); across the Finkelstein range = −29.6% (gate >15% fall,
PASS); the naive linear adversary predicts +533% over the Spratt range, cleanly falsified by
Spratt's own held-out finding. The Spratt-range near-constancy is achieved via a **disclosed,
chosen** (not independently measured) saturation-scale parameter — the genuinely non-circular,
held-out part is that the SAME parameter, calibrated only against Spratt, correctly predicts
Finkelstein's independent, opposite-direction finding without further tuning. See §4.

> **Decorrelated second observable** (frequency encodes LH vs FSH differentially — the
> Marshall/Belchetz/Wildt finding): using ONLY the two hormones' independently measured plasma
> clearance half-lives (zero extra free parameters), does a linear transfer-function argument
> predict LH's pulse-tracking ratio >> FSH's, consistent with Wildt et al. 1981 (monkey) and
> Spratt et al. 1987 (human)?

**PASSES** (ratio 4.47x–8.29x across the tested frequency range, gate >3x everywhere) —
**with a genuine, disclosed, UNRESOLVED tension**: Finkelstein et al. 1988's human
decreasing-frequency protocol did **not** clearly replicate Wildt 1981's monkey finding that FSH
invariably rises at slow GnRH frequency (FSH was "stable or fell" instead). Held OPEN, not
smoothed over. See §5.

> **Spectrum/stability** (geometric structure, not a falsifier against external data): does the
> linearized delayed-feedback loop, built from measured timing/gain parameters (no tuning for
> stability), sit in the stable regime consistent with the real absence of reported pathological
> slow self-oscillation in the male axis?

**PASSES, and more strongly than expected going in**: the loop is **delay-independently stable**
— because the measured feedback elasticity is inelastic (γ≈0.33<1), no amount of loop delay
destabilizes it (verified out to a 1e6-minute sweep; code validated against 2 controls: an
elastic γ>1 case does destabilize at a finite delay, and the classical pure-delay π/2 threshold
is reproduced to 4 significant figures). See §2.

**16/16 machine-computed gates PASS** (`overall_pass_strict_all = True` in the JSON) — see §6.
This is a **genuinely strong scoreboard, treated with the symmetric-QC suspicion that deserves**
(§7 spells out exactly which gates are hard external falsifiers vs. disclosed constructive
choices, so the strong scoreboard isn't mistaken for more certainty than it earns).

## Citations — every PMID/DOI below verified LIVE this session (NCBI E-utilities:
esearch then efetch, direct curl, full abstract text fetched and read — not WebFetch-summarized,
not recalled from memory; this repo's own prior finding across sibling docs is a measured
~62–67% citation-drift rate from memory alone)

| # | Citation | PMID / DOI | Role |
|---|---|---|---|
| 1 | Veldhuis JD, King JC, Urban RJ, Rogol AD, Evans WS, Kolp LA, Johnson ML (1987). "Operating characteristics of the male hypothalamo-pituitary-gonadal axis: pulsatile release of testosterone and follicle-stimulating hormone and their temporal coupling with luteinizing hormone." *J Clin Endocrinol Metab* 65(5):929-41. | **3117834**, DOI 10.1210/jcem-65-5-929 | Primary anchor, n=15: LH/T/FSH interpulse intervals, T pulse duration/amplitude, LH-T cross-correlation lag (60 min peak, 10-20 min earliest), circadian amplitudes. |
| 2 | Diver MJ, Imtiaz KE, Ahmad AM, Vora JP, Fraser WD (2003). "Diurnal rhythms of serum total, free and bioavailable testosterone and of SHBG in middle-aged men compared with those in young men." *Clin Endocrinol (Oxf)* 58(6):710-7. | **12780747**, DOI 10.1046/j.1365-2265.2003.01772.x | Falsifier 1's independent held-out anchor, n=10 young men: 43% peak-nadir T swing, acrophase 0700-0730h. |
| 3 | Spratt DI, Finkelstein JS, Butler JP, Badger TM, Crowley WF Jr (1987). "Effects of increasing the frequency of low doses of gonadotropin-releasing hormone (GnRH) on gonadotropin secretion in GnRH-deficient men." *J Clin Endocrinol Metab* 64(6):1179-86. | **3106396**, DOI 10.1210/jcem-64-6-1179 | Falsifier 2 increasing-frequency leg, n=5 IHH men: T constant over 8x frequency increase; nFSH falls faster than nLH. |
| 4 | Finkelstein JS, Badger TM, O'Dea LS, Spratt DI, Crowley WF (1988). "Effects of decreasing the frequency of gonadotropin-releasing hormone stimulation on gonadotropin secretion in gonadotropin-releasing hormone-deficient men and perifused rat pituitary cells." *J Clin Invest* 81(6):1725-33. | **3290251**, DOI 10.1172/JCI113512, PMC442617 | Falsifier 2 decreasing-frequency leg, n=12 IHH men + in vitro: T falls; steroid-independence confirmed via T-clamp AND rat pituitary cells; FSH did not clearly rise (the held-open tension, §5). |
| 5 | Wildt L, Hausler A, Marshall G, Hutchison JS, Plant TM, Belchetz PE, Knobil E (1981). "Frequency and amplitude of gonadotropin-releasing hormone stimulation and gonadotropin secretion in the rhesus monkey." *Endocrinology* 109(2):376-85. | **6788538**, DOI 10.1210/endo-109-2-376 | The classical primate frequency-coding paper: fast freq -> decline in both; slow freq (1/3h) -> FSH invariably rises. |
| 6 | Belchetz PE, Plant TM, Nakai Y, Keogh EJ, Knobil E (1978). "Hypophysial responses to continuous and intermittent delivery of hypothalamic gonadotropin-releasing hormone." *Science* 202(4368):631-3. | **100883**, DOI 10.1126/science.100883 | Foundational: pulsatility itself is structurally necessary (continuous GnRH fails; intermittent restores). |
| 7 | Knobil E, Plant TM, Wildt L, Belchetz PE, Marshall G (1980). "Control of the rhesus monkey menstrual cycle: permissive role of hypothalamic gonadotropin-releasing hormone." *Science* 207(4437):1371-3. | **6766566**, DOI 10.1126/science.6766566 | Scoping justification: fixed-frequency GnRH still permits a normal surge — surge is ovarian/pituitary-level. |
| 8 | Wildt L, Hausler A, Hutchison JS, Marshall G, Knobil E (1981). "Estradiol as a gonadotropin releasing hormone in the rhesus monkey." *Endocrinology* 108(5):2011-3. | **6783398**, DOI 10.1210/endo-108-5-2011 | Further scoping support, n=4: E2 alone (no GnRH) triggers a gonadotropin discharge. |
| 9 | Hayes FJ, Seminara SB, Decruz S, Boepple PA, Crowley WF Jr (2000). "Aromatase inhibition in the human male reveals a hypothalamic site of estrogen feedback." *J Clin Endocrinol Metab* 85(9):3027-35. | **10999781**, DOI 10.1210/jcem.85.9.6795 | Feedback-gain calibration anchor, n=14: E2 −61.8% -> LH pulse frequency +37.3% (10.2->14.0/24h), amplitude +47.4% (5.7->8.4 IU/L). |
| 10 | Santen RJ (1975). "Is aromatization of testosterone to estradiol required for inhibition of luteinizing hormone secretion in men?" *J Clin Invest* 56(6):1555-63. | **1104659**, DOI 10.1172/JCI108237, PMC333134 | Dissociates T (frequency effect) vs E2 (amplitude/pituitary effect); DHT confirms non-aromatized action. Tension vs #9, held open (§7). |
| 11 | Liu PY, Takahashi PY, Roebuck PD, Bailey JN, Keenan DM, Veldhuis JD (2009). "Testosterone's short-term positive effect on luteinizing-hormone secretory-burst mass and its negative effect on secretory-burst frequency are attenuated in middle-aged men." *J Clin Endocrinol Metab* 94(10):3978-86. | **19584190**, DOI 10.1210/jc.2009-0135, PMC2758726 | Graded T dose-response (castrate->physiological, n=23): monotonic LH suppression — the clean dose-dependence demonstration. |
| 12 | Pepperell RJ, Kretser DM, Burger HG (1975). "Studies on the metabolic clearance rate and production rate of human luteinizing hormone and on the initial half-time of its subunits in man." *J Clin Invest* 56(1):118-26. | **1170215**, DOI 10.1172/JCI108060, PMC436562 | LH clearance anchor, n=10: MCR 25.6 mL/min/m²; subunit half-time 15-18 min -> intact hLH ~33 min (this model's k_LH). |
| 13 | Santen RJ, Bardin CW (1973). "Episodic luteinizing hormone secretion in man. Pulse analysis, clinical interpretation, physiologic mechanisms." *J Clin Invest* 52(10):2617-28. | **4729055**, DOI 10.1172/JCI107454, PMC302522 | Independent, 14-yr-earlier pulse-frequency cross-check: 2.7-3.9 pulses/6h, overlapping Veldhuis 1987's number. |
| 14 | Urban RJ, Padmanabhan V, Beitins I, Veldhuis JD (1991). "Metabolic clearance of human follicle-stimulating hormone assessed by radioimmunoassay, immunoradiometric assay, and in vitro Sertoli cell bioassay." *J Clin Endocrinol Metab* 73(4):818-23. | **1909706**, DOI 10.1210/jcem-73-4-818 | FSH clearance anchor: half-life 274 min (this model's k_FSH) — the quantitative basis of the frequency-coding transfer-function argument (§5). |

## 1. Geometric structure — a real OODA moment in the kernel fit, not a one-shot success

**The core geometric object**: LH pulses drive testosterone via a linear transduction *kernel*
h(t) — its shape, not just its scale, is what a falsifiable model must get right. Veldhuis 1987
supplies two independent shape constraints from the SAME cohort: the LH->T cross-correlation
peak lag (60 min) and the T pulse duration/FWHM (90 min), giving a target ratio
FWHM/t_peak = 1.5.

**First attempt, machine-PROVEN infeasible before touching any external data**: the natural,
standard "2-compartment cascade" kernel family, h(t) = e^(−k1·t) − e^(−k2·t) (2 independent rate
constants), has FWHM/t_peak **bounded below by 2.446**, achieved only at the degenerate
equal-rate limit, and increasing monotonically away from it (machine-scanned, rate ratio
k2/k1 = 1 to 1000: `kernel_naive_family_proven_infeasible` = PASS,
`kernel_naive_family_monotonic_scan_confirmed` = PASS). The measured target ratio (1.5) sits
*below* this floor — infeasible for **any** (k1,k2), not a fitting failure to patch around.

**OODA, not surrender**: Orient — a higher-order n-stage Gamma(n,θ) cascade (h(t) ∝
t^(n−1)·e^(−t/θ), n continuous — physically, n sequential rate-limited transduction steps: LH-
receptor binding -> StAR induction -> cholesterol side-chain cleavage -> steroidogenic
conversion -> release) is a *different*, richer family whose FWHM/t_peak decreases
monotonically as n grows past 2 (machine-confirmed, `kernel_gamma_family_monotonic_decrease` =
PASS). Decide/Act — solve for n exactly: **n = 3.540, θ = 23.62 min**, reproducing both targets
to machine precision (`kernel_gamma_fit_matches_both_targets` = PASS). This fixes the model's
dominant relaxation rate a_T = 1/θ = 0.0423/min (half-life 16.4 min).

*(A second, related self-correction during this build: an early version of the kernel-family
check mistakenly asserted n=2 is a "local minimum" of the Gamma family itself — machine-run, that
specific claim came back FALSE, since the scan is monotonically decreasing throughout with no
turning point at n=2. Caught by the check, not smoothed over; the corrected, actually-load-
bearing claim — monotonic decrease past n=2 — holds and is what the fix relies on. Reported
here as a disclosed example of "force the fix, verify by machine, don't trust the first
narrative" in action, not hidden.)*

**Clearance rates** (2 further primary anchors): LH intact half-life ≈33 min (Pepperell 1975,
from subunit half-times "twice as great" as 15-18 min); FSH half-life 274 min (Urban 1991) — an
**8.30x** ratio, the quantitative basis of §5's frequency-coding argument.

**Pulse-frequency over-determination**: Veldhuis 1987's 95±11 min LH interpulse interval and
Santen & Bardin 1973's independent, 14-years-earlier 2.7-3.9 pulses/6h (=92.3-133.3 min) **overlap
in [92.3, 106.0] min** — two decorrelated eras/methods converging on the same number
(`clearance_two_era_pulse_frequency_overlap` = PASS).

## 2. Spectrum — the loop is delay-independently stable (an unforced, geometric finding)

Linearizing the closed loop (T feeds back on GnRH pulse frequency with elasticity γ, calibrated
in §3) gives the classical scalar delay-differential equation ẋ(t) = −a·x(t) − b·x(t−Δ), with
a = a_T = 0.0423/min (from §1's kernel), b = a·γ_freq (from §3), Δ = 60 min (§1's measured LH-T
lag). Its characteristic equation λ = −a − b·e^(−λΔ) is **transcendental**; its roots (the
spectrum) were solved via the Lambert W function across **multiple branches** (k = −4..4, not a
recalled closed-form threshold) — the rightmost root governs stability.

**Code validated against 2 controls before trusting the result**: (1) the classical pure-delay
(a=0) threshold Δ_crit = π/(2b) is reproduced numerically to 4 significant figures (analytic
31.416 min vs numeric crossing at 31.42 min for b=0.05/min); (2) a control case with elastic
feedback (b=0.08 > a=0.04234) **does** destabilize at a finite delay (crossing between 30-40
min), confirming the code correctly detects instability when it is genuinely present, not just
reporting "stable" unconditionally.

**The actual finding**: at the measured (a, b, Δ), Re(λ_rightmost) = **−0.0187** (stable). More
than that: because the measured feedback elasticity γ_freq ≈ 0.329 is **inelastic (< 1)**, b < a
— the classical sufficient condition for **delay-independent stability**. Swept out to a
1,000,000-minute delay (≈ 2 years, absurdly beyond any physiological relevance), the rightmost
root's real part **never crosses zero** — closest approach Re(λ) = −0.000068 at Δ ≈ 16,426 min
(≈ 11.4 days), still comfortably negative. **External anchor**: real 24-36h (Veldhuis 1987) and
longer serial-sampling studies of normal men report no sustained slow (hours-to-days) self-
oscillation in mean T superimposed on the pulsatile+diurnal pattern — this is the falsifiable
prediction the model must match, and does. An unstable prediction at the measured parameter
point would have been a clean falsification of this reduced loop structure, forcing a search for
an unmodeled damping mechanism — it wasn't needed.

## 3. Falsifier 1 — diurnal T rhythm

Veldhuis 1987's own numbers (circadian amplitude 185 ng/dL; pulse peak 910 / incremental 242
ng/dL => nadir 668 => implied mean 789 ng/dL) give, via the standard cosinor convention
(amplitude = half peak-to-trough swing): **2×185/789 = 46.9%** peak-to-nadir swing. This is
compared against Diver et al. 2003's fully independent, directly measured **43.0%** (n=10 young
men, dense 30-min sampling + cosinor fit, a different cohort/decade/method from Veldhuis's
Fourier/deconvolution approach). **Relative difference 9.1%**, well inside the pre-registered
25% threshold — PASS.

**Disclosed limitation**: the 46.9% figure is derived from the SAME 1987 paper/cohort that
supplies the kernel-timing numbers (a same-source, different-statistic estimate — Fourier
circadian amplitude vs. pulse peak/nadir), not fully independent by itself. The genuinely
independent, falsifiable comparison is its agreement with Diver 2003. **Supporting structural
check** (not a falsifier, a sanity check on the kernel): the fitted kernel's own transfer-
function gain at the 24-hour diurnal frequency is ≈1.0 (unity gain) — a fast (tens-of-minutes)
transduction kernel does not meaningfully distort a slow (24h) rhythm, confirming the kernel
behaves sensibly as a low-order filter relative to the diurnal timescale.

## 4. Falsifier 2 — GnRH/LH pulse frequency -> mean T is SATURATING, not linear

**Forced adversary**: a naive model with T_mean(f) directly proportional to frequency (no
desensitization) predicts a **+533%** T rise across Spratt 1987's 8x frequency-increase range
(120->15 min interpulse) — this adversary is **cleanly falsified** by Spratt's own held-out
finding that T remained approximately constant.

**This model** adds a frequency-dependent pituitary-desensitization term (direction anchored:
both Wildt 1981 and Spratt 1987 report reduced per-pulse pituitary response at higher GnRH
frequency), giving T_mean(f) = C·f/(1+f/f_c) — a saturating (Michaelis-Menten-shaped) curve. The
saturation scale f_c is a **disclosed, chosen** (not independently measured) parameter, solved so
the model's T-change across the Spratt range is +10% (comfortably under the 20% "near-constant"
gate) — an **existence/constructive demonstration**, not itself a fully independent test.

**The genuinely non-circular, held-out part**: that SAME f_c (never touched again) is then used,
without further fitting, to predict T's behavior across Finkelstein 1988's **independent**,
decreasing-frequency range (out to 8h interpulse) — giving a **−29.6%** fall, comfortably past
the pre-registered 15%-fall gate. Calibrate on one study's qualitative direction, correctly
predict a different, decorrelated study's opposite-direction qualitative finding, with zero
re-fitting — a real, if modest, cross-validation.

**Honest scope note**: f_c itself corresponds to an ~819-min (13.7h) half-saturation interpulse
interval, which sits *outside* the directly tested range in either paper (15-480 min) — an
extrapolated scale parameter, not an independently observed transition point.

## 5. Decorrelated second observable — LH vs FSH frequency-tracking ratio

Zero extra free parameters: using only the two independently measured clearance half-lives
(k_LH from Pepperell 1975, k_FSH from Urban 1991), a first-order linear transfer-function
argument, |H(jω)| = k/√(k²+ω²), gives each hormone's pulse-tracking ratio as a function of GnRH
pulse frequency. Computed across the real tested range (15-480 min interpulse): **ratio
LH/FSH = 4.47x to 8.29x**, always exceeding the pre-registered 3x gate — LH tracks pulse
frequency far more faithfully than FSH throughout, because its clearance (t½≈33 min) is 8.3x
faster than FSH's (t½≈274 min).

**Directionally consistent** with Wildt et al. 1981 (monkey: slow frequency -> FSH invariably
rises, ratio shrinks toward 1 at the slow end in this model, i.e. FSH "catches up" relatively)
and Spratt et al. 1987 (human: fast frequency -> nFSH falls faster than nLH, ratio is largest at
the fast end here). **Genuine, disclosed, UNRESOLVED tension, held open per task instruction**:
Finkelstein et al. 1988's human decreasing-frequency protocol found FSH "stable or fell" — it did
**not** clearly replicate Wildt's monkey finding that FSH reliably rises at slow frequency. The
pharmacokinetic-clearance argument built here explains why LH tracks frequency better than FSH
*in general*, but does not by itself explain why FSH failed to rise in this specific human study
— the literature's own emphasis on an *additional* transcriptional/Ca-oscillation-frequency-
decoding mechanism (Kaiser/Coss-lab-type work) was **not independently verified this session** (a
disclosed scope limit, not a fabricated mechanism) and may be the dominant missing piece.

## 6. Pre-registered gates — 16/16 PASS

```
kernel_naive_family_proven_infeasible:                PASS
kernel_naive_family_monotonic_scan_confirmed:         PASS
kernel_gamma_family_monotonic_decrease:               PASS  (corrected from a false "n=2 local min" premise, see §1)
kernel_gamma_fit_matches_both_targets:                PASS
clearance_two_era_pulse_frequency_overlap:            PASS
spectrum_stable_at_primary_delay:                     PASS
spectrum_stable_across_veldhuis_delay_range:          PASS
spectrum_margin_gt_1p5x:                              PASS  (unconditional: delay-independent stability)
falsifier1_diurnal_swing_within_25pct_of_diver2003:   PASS  (9.1% relative difference)
falsifier1_kernel_near_unity_gain_at_24h:             PASS
falsifier2_spratt_near_constant:                      PASS  (+9.6%, disclosed constructive choice of f_c)
falsifier2_finkelstein_falls:                         PASS  (-29.6%, genuine held-out prediction)
falsifier2_naive_adversary_falsified:                 PASS  (+533% predicted vs Spratt's ~constant)
falsifier3_lh_fsh_ratio_gt3x_everywhere:              PASS  (4.47x-8.29x)
couples_hpg_rate_layer_inputs_satisfied:              PASS
couples_reds_cascade_input_satisfied:                 PASS
```

`overall_pass_strict_all = True`. Determinism: 2 independent runs produce byte-identical JSON
(verified via `diff`); zero NaN/Inf anywhere in the output tree (checked programmatically over
the full JSON).

## 7. Symmetric QC — why a 16/16 scoreboard is NOT treated as more certain than it is

A clean scoreboard is exactly when this discipline says to be MOST suspicious, not least. Sorted
by evidentiary weight:

- **Hard, non-circular external falsifiers** (could genuinely have failed): Falsifier 1's
  agreement with Diver 2003 (§3); Falsifier 2's naive-adversary falsification AND its Finkelstein
  held-out prediction (§4, the f_c-transfer part only, not the Spratt-calibration part); the
  LH/FSH ratio's directional agreement with Wildt/Spratt (§5); the 2-era pulse-frequency overlap
  (§1); the kernel-family infeasibility proof (§1) — none of these were tunable by a free
  parameter I controlled.
- **Constructive/existence demonstrations** (disclosed, not hidden, weaker evidence): the kernel's
  own n/θ values are a fit BY CONSTRUCTION to hit 2 targets exactly (a self-consistency check,
  not a discovery); Falsifier 2's Spratt-range gate is satisfied by a chosen f_c, not an
  independent measurement; the spectrum's Δ=60min is the primary measured value but the slow-
  adaptation timescale folded into `a_T` is a coarse-grained reduction, not a directly measured
  quantity in its own right.
- **Held OPEN, not reconciled** (per task instruction, genuine unresolved tensions): Santen 1975
  (non-aromatized androgen -> frequency site) vs. Hayes 2000 (aromatized estrogen -> BOTH
  frequency+amplitude sites, 25 years later, more selective pharmacology) — which steroid, which
  site, is genuinely contested across eras. Wildt 1981 (monkey, slow freq -> FSH reliably rises)
  vs. Finkelstein 1988 (human, slow freq -> FSH stable-or-fell) — a real cross-species/study
  disagreement (§5).
- **Real inter-subject/inter-study spread**, not point estimates: Veldhuis 1987's own SEMs show
  meaningful variability (LH interpulse 95±11 min ≈ 12% CV across n=15; T pulse amplitude
  910±92 ng/dL ≈ 10% CV); Santen & Bardin 1973's "apparent half-life" of LH spans 34-233 min
  (method-dependent, attributed to variable multi-pool mixing) — a much wider, more honest range
  than the single ~33-min point estimate this model uses for k_LH. Liu 2009 additionally shows
  the feedback sensitivity itself weakens with age — this model is scoped to young/adult men
  only (matching Veldhuis 1987's and Diver 2003's own "young" cohorts), not a lifespan model.
- **A genuinely corrected error**, not a swept-under-the-rug one: an early version of this
  build's own Gamma-family check asserted a false claim (n=2 as a "local minimum") — caught by
  running the check, not by inspection, and fixed at the source rather than patched to force a
  pass (§1).
- **Not independently verified this session** (disclosed scope limits): the intracellular Ca²⁺-
  oscillation/transcriptional LHβ-vs-FSHβ frequency-decoding mechanism the field treats as
  probably dominant (§5); any duration-dependence of the feedback gain beyond Hayes 2000's 7-day
  protocol; the female cycle's surge dynamics (deliberately out of scope, §0).

## 8. couples_to — a machine-checked key-presence match, not a fabricated coupling number

This model's own output state-vector keys (`serum_T`, `serum_FSH`, `LH_pulse_frequency_t`,
`LH_pulse_amplitude_t`, `GnRH_pulse_interval_t`) were checked against the exact variables
`MSK-HPG-RATE-CONSTANT-LAYER` ("serum ... T/FSH") and `INT-REDS-ENERGY-DEFICIENCY`
("LH-pulse-frequency(t)") name as required inputs in their own pre-existing SEED-DESIGN claim
text (`data/MECHANISM_ANCHOR_GRAPH.json`, read-only, not edited) — both satisfied
(machine-checked set-containment, PASS/PASS). **Disclosed, not overreached**: no bone/muscle
tissue-rate-constant cell in this repo currently exposes a hormone-sensitive parameter to couple
a concrete re-derived number into — grep-verified, `bone_remodeling.py` has zero
estrogen/testosterone/androgen hooks. Unlike `MECHANISM_THYROID_AXIS.md`'s thermoregulation
coupling (which re-derived an actual downstream number), this doc offers only the key-presence
match, not a fabricated numeric coupling — an honest scope difference, not an oversight. Real
literature exists for both named downstream couplings (estrogen/bone; testosterone/muscle
protein synthesis+metabolic rate) but was not pulled into a re-computed number this session.

## 9. Confidence tier

**In-vivo-anchored** (real human serial/frequent-sampling hormone data: n=15 Veldhuis 1987, n=10
Diver 2003, n=5 Spratt 1987, n=12 Finkelstein 1988, n=14 Hayes 2000, n=23 Liu 2009, n=10
Pepperell 1975; real rhesus monkey primary data: Wildt 1981, Belchetz 1978, Knobil 1980) — one
tier below a subject-specific in-vivo measurement (no serial LH/FSH/T sampling panel exists for
subject2 — the same disclosed scope every sibling endocrine layer in this repo already carries,
e.g. `MECHANISM_THYROID_AXIS.md`). Every operating point here is population/literature-level.

## 10. Honest gaps — symmetric QC: what this does NOT prove

- **The female menstrual-cycle LH surge is entirely out of scope** — a deliberate choice
  (§0), not an oversight, given the task's own instruction to pick the cleaner falsifier. The
  bistable-switch/threshold dynamics of the positive-feedback surge are a genuinely different,
  harder mathematical object (not attempted here).
- **Falsifier 2's Spratt-range gate rests on a chosen, not measured, saturation scale f_c** —
  disclosed throughout (§4, §7); only the Finkelstein cross-validation and the naive-adversary
  falsification are fully independent of that choice.
- **The Santen 1975 vs. Hayes 2000 mechanistic-attribution tension is held OPEN**, not
  reconciled — genuinely unclear which steroid/site combination is the more complete account,
  25 years and one generation of more selective pharmacology apart (§7).
- **The Wildt 1981 vs. Finkelstein 1988 slow-frequency FSH tension is held OPEN** — the
  clearance-kinetics argument built here (§5) does not explain why FSH failed to rise in
  Finkelstein's human protocol; an additional transcriptional mechanism is the literature's own
  likely candidate, not independently verified this session.
- **The slow-adaptation delay/timescale folded into the spectrum analysis (§2) is a coarse-
  grained reduction**, not a directly measured quantity — a fuller model would resolve fast
  (pulse-to-pulse) and slow (day-to-day hypothalamic adaptation) dynamics as separate state
  variables, which this reduced scalar DDE deliberately collapses for tractability.
- **Real inter-subject/inter-study spread is wide** for several of this model's point-estimate
  inputs (LH "apparent half-life" 34-233 min per Santen & Bardin 1973, vs. the ~33-min point
  estimate used here from Pepperell 1975's cleaner constant-infusion method) — reported, not
  hidden (§7).
- **No couples_to numeric re-derivation into bone/muscle** — only a key-presence match (§8);
  grep-verified no hormone-sensitive hook exists in `bone_remodeling.py` this session.
- **Single confidence tier throughout: population/literature-anchored, never subject-specific**
  — no hormone panel exists for subject2, matching every sibling endocrine doc's own disclosed
  scope.
- **No graph-edge write this session** — folding into the pre-existing `MSK-HPG-RATE-CONSTANT-
  LAYER` / `INT-REDS-ENERGY-DEFICIENCY` SEED-DESIGN nodes via `mechanism_fold -> fold_gate_v2` is
  the natural next step, not performed here (isolation rule: touch only files created this
  session).

## Files

- `scripts/msk/hpg_male_axis.py` — self-contained (numpy/scipy only), builds the kernel fit (incl.
  the disclosed infeasibility-then-fix OODA sequence), clearance/frequency cross-check, feedback-
  gain calibration, the Lambert-W multi-branch spectrum/stability analysis (validated against 2
  controls), both falsifiers, the decorrelated frequency-coding observable, the couples_to
  key-presence check, and all 16 gates; writes the evidence JSON below; prints a full summary.
- `data/msk_smoketest/subject2_walking1/hpg_male_axis/hpg_male_axis_results.json` — every number
  in this doc, machine-written: all 14 citations, the kernel infeasibility-proof and fix, the
  clearance/frequency cross-check, feedback-gain elasticities, the full spectrum analysis (incl.
  the wide-sweep worst-case and both validation controls), both falsifiers, the frequency-coding
  transfer-function grid, the couples_to check, and all 16 gates. Verified deterministic (2
  independent runs, byte-identical JSON via `diff`) and NaN/Inf-free (checked programmatically
  over the full JSON tree).
- Read but NOT modified (isolation: touch only files created this session):
  `data/MECHANISM_ANCHOR_GRAPH.json` (the SEED-DESIGN nodes this doc resolves: `MSK-HPG-RATE-
  CONSTANT-LAYER`, `AUTO-BUILD-AN-HPG-RATE-CONSTANT-LAYER-A-PER-S`, `INT-REDS-ENERGY-DEFICIENCY`,
  and context from `NEU-HYPOTHALAMUS-HUB` / `XDOMAIN-ANDROGEN-MUSCLE-PROSTATE-SATURATION-GATE`),
  `data/body_twin/agent_outputs/hpg-rate-constant-layer__a3af689a9fa9be906.json` and
  `.../reproduction-pregnancy__a13d17371a3dd9d85.json` (prior sessions' literature-scout outputs —
  treated as HYPOTHESES per this project's own discipline; none of their specific numbers were
  reused as load-bearing here — this doc's citations were independently re-verified live from
  scratch).

## Repro

```
cd ~/projects/bodytwin
source .venv-msk/bin/activate
python3 scripts/msk/hpg_male_axis.py
```
No external inputs required. Pure Python/numpy/scipy (`scipy.special.lambertw`, `scipy.optimize.
brentq`), no OpenSim call, runs in under 2 seconds, deterministic. No git operations; writes only
under `data/msk_smoketest/subject2_walking1/hpg_male_axis/`.
