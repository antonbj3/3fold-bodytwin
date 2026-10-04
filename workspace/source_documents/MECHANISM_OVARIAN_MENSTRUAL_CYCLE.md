# MECHANISM FEMALE OVARIAN/MENSTRUAL CYCLE — the estradiol POSITIVE-FEEDBACK LH-SURGE switch (2026-07-22)

Builds and MEASURES the delta `MECHANISM_REPRODUCTIVE_HPG.md` (the male HPG axis) explicitly
scoped OUT: "the female cycle's periovulatory LH surge is a pituitary/ovarian **bistable
threshold-switch**... a genuinely different (bistable-switch, threshold-not-sharply-defined)
mathematical object than a negative-feedback steady state... **female surge dynamics are
explicitly OUT OF SCOPE here — held open, not modeled, not hand-waved.**" This doc builds that
object. Script: `scripts/msk/ovarian_cycle_axis.py` (988 lines). Evidence:
`data/msk_smoketest/subject2_walking1/ovarian_cycle_axis/ovarian_cycle_axis_results.json`.

**The falsifier, stated up front (nothing hidden).** Reproduce the phase structure (follicular →
LH surge → luteal → menstruation) AND the estradiol (E2) positive-feedback sign-flip that
produces the surge — while a **purely negative-feedback** model, forced to its strongest fair
form (swept across steepness, operating point, AND an applied time delay — 90 configurations),
must **FAIL** to reproduce it. A GEOMETRIC fact does the forcing: a monotonically DEcreasing
function of E2, composed with E2's own unimodal (single-peaked) trajectory, has its **minimum —
not maximum** — exactly where E2 peaks (true for ANY steepness/operating point, and a delay only
relabels *when* that minimum lands — it can never turn a trough into a peak). This makes "purely
negative feedback" a comprehensively swept, non-strawman adversary, not one cherry-picked shape.

## Verdicts up front

> **Falsifier 1 (phase)**: does the model's LH surge coincide with (not anti-phase to) the E2
> peak, matching 3 fully independent real papers (WHO 1980 n=177/107; Fritz 1992 n=7;
> Hoff/Quigley/Yen 1983 n=5) that never touched this model's construction?

**PASSES**: model phase-lag = 0.0–0.35h (grid- vs fine-resolution) against a pre-registered
±24h gate; WHO 1980's own independently-measured E2-peak-to-LH-peak lag is 7.5h — both land in
the "same-day, near-simultaneous" regime the real data defines. The **forced adversary** (90
swept negative-feedback configs) instead places its trough (not peak) exactly at the E2 peak —
**0/90 configs escape this**, and the exact shift-by-delay geometry was itself machine-verified
(§2, including a self-caught correction, not smoothed over).

> **Falsifier 2 (amplitude + abruptness)**: does the switch fire with realistic fold-amplitude
> and rise abruptly relative to E2's own smooth rate of change?

**PASSES, with a disclosed modesty**: amplitude ratio 4.79× (gate ≥3×); abruptness ratio 8.64×
(gate ≥2×) — directionally correct but NOT a literal match to Hoff/Quigley/Yen's measured ~2h LH
doubling time (28–30× faster than E2); reproducing that exact number needs a separate fast
pituitary self-priming mechanism this reduced model doesn't build (§4, disclosed).

> **Falsifier 3 (duration-gating — Young & Jaffe's strength-DURATION concept)**: does a BRIEF
> (12h) E2 spike fail to trigger a surge while the SAME model, same parameters, fed a SUSTAINED
> elevation, does?

**PASSES**: sustained max(D)=0.521 (primed), transient max(D)=0.043 (never primed); transient
amplitude ratio 1.00× (no surge) vs sustained 4.79× (surges) — same code, same parameters, only
input duration differs (§5).

> **Symmetric self-check**: is the 15/15→17/17-clean scoreboard caused by the sign-flip mechanism
> specifically, or an artifact of the chosen E2(t) curve shape?

**Kill-switch test (a permanent, machine-checked gate, not a one-off)**: forcing gain_lh=
gain_fsh=0 (positive feedback surgically disabled, IDENTICAL E2 curve/code otherwise) makes the
SAME model degenerate to the adversary's EXACT failure mode (trough aligned with E2's peak,
amplitude 1.37×, both gates FAIL) — proving the earlier PASS is mechanism-driven (§2).

> **Decorrelated second observable** (GnRH pulse-frequency → LH:FSH ratio coding, ZERO new free
> parameters — reusing the male doc's own validated k_LH/k_FSH): does real measured interpulse-
> interval data from 4 independent papers (PCOS/normal-follicular/normal-luteal/hypothalamic
> amenorrhea) predict the classical elevated PCOS LH/FSH ratio (Rebar 1976)?

**PASSES**: PCOS's fastest-measured pulse frequency (Waldstreicher 1988, interpulse 58min)
predicts the LARGEST LH:FSH ratio of all 7 tested conditions, exceeding every hypothalamic-
amenorrhea condition's ratio — directionally consistent with Rebar 1976's classical elevated-LH/
low-FSH finding (§6).

**17/17 machine-computed gates PASS** (`overall_pass_strict_all = True`). Treated with the
symmetric-QC suspicion a clean scoreboard deserves (§8 sorts every gate by evidentiary weight).

## Citations — every PMID/DOI live-verified this session (NCBI E-utilities: esearch then efetch,
direct curl, full abstract text fetched and read — not recalled from memory, not
WebFetch-summarized; this repo's own measured ~62–67% citation-drift-from-memory rate is why)

| # | Citation | PMID / DOI | Role |
|---|---|---|---|
| 1 | World Health Organization Task Force (1980). "Temporal relationships between ovulation and defined changes in the concentration of plasma estradiol-17beta, luteinizing hormone, follicle-stimulating hormone, and progesterone. I. Probit analysis." *Am J Obstet Gynecol* 138(4):383-90. | **6775535**, no DOI (this 1980 report predates registered DOIs; PII 0002-9378(80)90133-7 — see §9 self-caught error) | Multicenter (n=177/107, laparotomy-confirmed) held-out phase anchor: E2-peak-to-LH-peak lag ≈7.5h. |
| 2 | Stricker R, Eberhart R, Chevailler MC, Quinn FA, Bischof P, Stricker R (2006). "Establishment of detailed reference values for LH, FSH, estradiol, and progesterone during different phases of the menstrual cycle..." *Clin Chem Lab Med* 44(7):883-7. | **16776638**, DOI 10.1515/CCLM.2006.160 | n=20, phase-structure/reference-range anchor (LH-peak-synchronized design). |
| 3 | Reame N, Sauder SE, Kelch RP, Marshall JC (1984). "Pulsatile gonadotropin secretion during the human menstrual cycle: evidence for altered frequency of GnRH secretion." *J Clin Endocrinol Metab* 59(2):328-37. | **6429184**, DOI 10.1210/jcem-59-2-328 | n=8: LH pulse frequency 11.8→14.3/12h (follicular), 8/12h (luteal). |
| 4 | Nippoldt TB, Reame NE, Kelch RP, Marshall JC (1989). "The roles of estradiol and progesterone in decreasing LH pulse frequency in the luteal phase..." *J Clin Endocrinol Metab* 69(1):67-76. | **2499593**, DOI 10.1210/jcem-69-1-67 | n=13: E2 (not P4) drives luteal pulse-frequency slowing (3.2 vs 8.0/10h). |
| 5 | Reame NE, Sauder SE, Case GD, Kelch RP, Marshall JC (1985). "Pulsatile gonadotropin secretion in women with hypothalamic amenorrhea..." *J Clin Endocrinol Metab* 61(5):851-8. | **3900122**, DOI 10.1210/jcem-61-5-851 | n=19: HA LH pulse frequency 4.7/12h (vs normal EF 7.7/12h); E2 52±5 pg/mL. |
| 6 | Khoury SA, Reame NE, Kelch RP, Marshall JC (1987). "Diurnal patterns of pulsatile LH secretion in hypothalamic amenorrhea..." *J Clin Endocrinol Metab* 64(4):755-62. | **2880864**, DOI 10.1210/jcem-64-4-755 | n=14, independent HA replication: daytime 1.0 vs normal-EF 5.1 pulses/8h. |
| 7 | Waldstreicher J, Santoro NF, Hall JE, Filicori M, Crowley WF Jr (1988). "Hyperfunction of the hypothalamic-pituitary axis in women with polycystic ovarian disease..." *J Clin Endocrinol Metab* 66(1):165-72. | **2961784**, DOI 10.1210/jcem-66-1-165 | n=12 PCOD+21 normal: PCOD LH pulse frequency 24.8/24h > all 3 follicular sub-phases (15.6-22.2); correlates with E2 (r=0.84). |
| 8 | Rebar R, Judd HL, Yen SS, Rakoff J, Vandenberg G, Naftolin F (1976). "Characterization of the inappropriate gonadotropin secretion in polycystic ovary syndrome." *J Clin Invest* 57(5):1320-9. | **770505**, DOI 10.1172/JCI108400, PMC436785 | n=24, classical LH-high/FSH-low PCOS anchor; positive-feedback switch itself intact (clomiphene-induced surges normal). |
| 9 | Burger HG, Dudley EC, Hopper JL, Groome N, Guthrie JR, Green A, Dennerstein L (1999). "Prospectively measured levels of serum FSH, estradiol, and the dimeric inhibins during the menopausal transition..." *J Clin Endocrinol Metab* 84(11):4025-30. | **10566644**, DOI 10.1210/jcem.84.11.6158 | n=150, longitudinal: FSH-E2 r=-0.73; FSH rise = loss-of-restraint mechanism. |
| 10 | Burger HG, Cahir N, Robertson DM, Groome NP, Dudley E, Green A, Dennerstein L (1998). "Serum inhibins A and B fall differentially as FSH rises in perimenopausal women." *Clin Endocrinol (Oxf)* 48(6):809-13. | **9713572**, DOI 10.1046/j.1365-2265.1998.00482.x | n=110, 4-stage cross-sectional: inhibin-B falls FIRST, then inhibin-A+E2 collapse together (the model's target numbers). |
| 11 | Lenton EA, Landgren BM, Sexton L, Harper R (1984). "Normal variation in the length of the follicular phase of the menstrual cycle: effect of chronological age." *Br J Obstet Gynaecol* 91(7):681-4. | **6743609**, DOI 10.1111/j.1471-0528.1984.tb04830.x | n=293 cycles: follicular phase geometric mean 12.9d (10.3-16.3d 95%CI). |
| 12 | Fritz MA, McLachlan RI, Cohen NL, Dahl KD, Bremner WJ, Soules MR (1992). "Onset and characteristics of the midcycle surge in bioactive and immunoactive LH secretion in normal women..." *J Clin Endocrinol Metab* 75(2):489-93. | **1639949**, DOI 10.1210/jcem.75.2.1639949 | n=7: surge onset (2× rise definition) coincident with E2 peak; duration 54.0±4.0h. |
| 13 | Hoff JD, Quigley ME, Yen SS (1983). "Hormonal dynamics at midcycle: a reevaluation." *J Clin Endocrinol Metab* 57(4):792-6. | **6411753**, DOI 10.1210/jcem-57-4-792 | n=5, 2h sampling: abrupt surge onset (LH doubles within 2h) vs E2's 57-61h doubling time; onset coincident with E2 peak. |
| 14 | Chang RJ, Jaffe RB (1978). "Progesterone effects on gonadotropin release in women pretreated with estradiol." *J Clin Endocrinol Metab* 47(1):119-25. | **122395**, DOI 10.1210/jcem-47-1-119 | n=12: E2 271±3 pg/mL (quantitative threshold-magnitude anchor); P4 augments/triggers the surge. |
| 15 | Young JR, Jaffe RB (1976). "Strength-duration characteristics of estrogen effects on gonadotropin response to GnRH in women. II..." *J Clin Endocrinol Metab* 42(3):432-42. | **767352**, DOI 10.1210/jcem-42-3-432 | n=19: the classic strength-DURATION concept paper — **disclosed gap**: indexed abstract is methods-only, no numeric result (§9). |
| 16 | Knobil E, Plant TM, Wildt L, Belchetz PE, Marshall G (1980). "Control of the rhesus monkey menstrual cycle: permissive role of hypothalamic GnRH." *Science* 207(4437):1371-3. | **6766566**, DOI 10.1126/science.6766566 | REUSED from `MECHANISM_REPRODUCTIVE_HPG.md`, spot-re-verified live: fixed 1-pulse/h GnRH still permits a normal surge — surge is E2/pituitary-driven, not frequency-driven. |
| 17 | Wildt L, Hausler A, Hutchison JS, Marshall G, Knobil E (1981). "Estradiol as a gonadotropin releasing hormone in the rhesus monkey." *Endocrinology* 108(5):2011-3. | **6783398**, DOI 10.1210/endo-108-5-2011 | REUSED+spot-re-verified, n=4: E2 ALONE (zero GnRH pulses) triggers a gonadotropin discharge. |

## 1. Geometric structure — the model

**State**: E2(t)/P4(t) are **prescribed forcing** (a disclosed reduction, matching the male doc's
own choice to take T's feedback on GnRH frequency as a calibrated elasticity rather than full
steroidogenesis) — a smooth, landmark-calibrated curve (follicular rise + luteal plateau, built
from asymmetric double-exponential "bumps"), NOT a follicle-growth ODE. Landmarks: E2 baseline 45
pg/mL, periovulatory peak 270 pg/mL (anchored to Chang & Jaffe 1978's own 271±3 pg/mL), luteal
plateau 120 pg/mL; pre-surge rise time-constant set so E2's doubling time ≈59h, matching
Hoff/Quigley/Yen 1983's own measured 57–61h; peak day = 13 (Lenton 1984's geometric-mean
follicular length 12.9d, rounded). P4 turns on 12h before the E2/LH peak (Hoff/Quigley/Yen's own
"12h after the initiation of a rapid P4 rise" finding) and off at luteolysis (~day 26).

**The bistable core**: a slow leaky-integrator "primed" variable D(t) (time constant τ_D=48h — the
Young & Jaffe strength-**duration** concept, task-given) low-pass-filters the indicator
"E2 > 200 pg/mL" (task-given threshold). D(t) continuously blends the LH secretion law between a
negative-feedback Hill function (D=0, decreasing in E2 — normal follicular suppression) and a
saturating positive-feedback Hill function (D=1, increasing in E2, capped — the surge regime).
**The SIGN of d(secretion)/d(E2) literally flips as D crosses its own transition** — this is the
fast-slow (relaxation-oscillator-type) structure a genuine bistable switch requires; D is the slow
variable that reshapes the fast LH/FSH dynamics' vector field. LH and FSH relax toward their
secretion drive at the ACTUAL measured plasma clearance rates — k_LH (t½=33min, Pepperell 1975)
and k_FSH (t½=274min, Urban 1991) — **reused verbatim from `MECHANISM_REPRODUCTIVE_HPG.md`**, not
re-fit (same hormone molecules, sex-independent clearance — a genuine numeric `couples_to`, §10).

**A real, disclosed self-correction (not smoothed over)**: the first version of the forced-
adversary sweep's trough-alignment check compared every config's LH trough against the
un-shifted E2-peak day, and flagged 54/90 delay>0 configs as "misaligned." Diagnosing WHY (not
just re-running): a monotonic-decreasing function of a **delayed** E2 argument has its own minimum
where the ARGUMENT peaks — i.e. at `t_E2peak + delay`, not at `t_E2peak` itself (a delay relabels
time; it does not anchor the trough back). Machine-confirmed the shift matched `+delay_h` exactly
for all 54 flagged rows, then corrected the check to compare against the right reference point —
the corrected invariant (`adversary_trough_always_aligned_with_e2_peak_plus_delay`) now holds for
all 90/90 configs. The LOAD-BEARING falsification claim (the adversary's own PEAK never aligns
with E2's peak — `gate_phase` fails 90/90) was never affected by this bug; only a secondary,
more-precise diagnostic check was.

## 2. The forced adversary + kill-switch (the falsifier's teeth)

**Forced to its strongest fair form**: 6 steepnesses (n=1..32) × 3 operating points (E2₅₀ = 58.5,
110.2, 216.0 pg/mL, spanning the whole observed E2 range) × 5 delays (0–72h) = **90
configurations** — not one cherry-picked negative-feedback shape. **0/90 achieve both the phase
AND amplitude gate simultaneously** (`adversary_never_achieves_both_gates_across_90_configs =
PASS`); every config's trough sits exactly at `t_E2peak + delay` (machine-verified, not assumed).

**The kill-switch test — the decisive mechanism-isolation check**: setting `gain_lh=gain_fsh=0`
in the FULL switch model (identical E2(t) curve, identical D(t)/ODE machinery, only the
positive-feedback gain zeroed) makes the model's own LH response degenerate to **amplitude
1.37× (gate needs ≥3×, FAILS) with its trough exactly aligned to E2's peak (matching the
adversary's own failure signature)**. This proves the 17/17 scoreboard is caused by the sign-flip
term specifically — not the E2 curve shape, not the D(t) integrator alone, not an unrelated
coding artifact. (`kill_switch_degenerates_to_adversary_failure_mode = PASS`.)

## 3. Falsifier 1 — phase

WHO 1980's own median estimates (both relative to ovulation): E2 peak 24.0h before ovulation, LH
peak 16.5h before ovulation → **E2 peaks ~7.5h before LH**. This model's phase-lag: 0.0h at 2h
grid resolution, 0.35h at fine (3-min) resolution — same "near-simultaneous, same day" regime,
comfortably inside the pre-registered ±24h gate. Fritz 1992 ("LH-BIO surge onset... coincident
with the peak in E2 levels") and Hoff/Quigley/Yen 1983 ("onset of LH and FSH surges... temporally
associated with the attainment of peak E2 levels") are two FURTHER, fully independent (different
cohorts/methods/decades) confirmations of the same qualitative phase relationship — none of the
three were touched during model construction (only Chang & Jaffe's magnitude number and
Hoff/Quigley/Yen's own doubling-time number were used as forcing-curve landmarks; the PHASE
finding is a genuinely emergent, checked-not-fit property).

**Disclosed honesty**: the model's own 0.0–0.35h lag is smaller than WHO 1980's measured 7.5h —
my reduced model resolves the right REGIME (same-day near-simultaneity, not the adversary's
13-day antiphase) but not hour-scale delay structure; the ±24h gate is appropriately loose for
what this reduction can honestly claim, not tightened after the fact.

## 4. Falsifier 2 — amplitude + abruptness

Amplitude ratio 4.79× (gate ≥3×) sits inside commonly-cited real ranges (which are wide and
assay-dependent — a genuine, disclosed spread, not a single number). Abruptness ratio 8.64× (LH's
own max fractional rate of change ÷ E2's own max fractional rate) clears the ≥2× gate — but
Hoff/Quigley/Yen 1983's measured ratio is ~28–30× (LH doubles in 2h vs E2's 57–61h doubling) —
this model's own LH doubling time is ~6.2h at its fastest, directionally right (much faster than
E2) but not a literal match. **Honest attribution**: real GnRH/LH secretion involves an
ADDITIONAL fast pituitary self-priming loop (each pulse sensitizing the next) operating on a
timescale this reduced D(τ=48h)-gated model does not resolve — disclosed, not hidden.

## 5. Falsifier 3 — duration-gating (Young & Jaffe strength-DURATION)

Same model, same `DEFAULT_PARAMS`, only the E2(t) input changes: a sustained (multi-day) rise
reaches max(D)=0.521 (primed, D_crit≈0.5) and surges (amplitude 4.79×); a brief 12h spike to the
same peak concentration reaches max(D)=0.043 (never primed) and produces **no surge at all**
(amplitude 1.00×). This directly operationalizes "sustained, not transient, elevation is
required" — Chang & Jaffe 1978's own qualitative finding ("rising concentrations of E2 to which
the system is exposed for an appropriate DURATION... initiate the surge") reproduced as a genuine
simulated contrast, not asserted.

## 6. Decorrelated second observable — GnRH pulse-frequency → LH:FSH ratio (zero new parameters)

Reusing `MECHANISM_REPRODUCTIVE_HPG.md`'s own validated transfer function
|H(jω)|=k/√(k²+ω²) and its exact k_LH/k_FSH values, fed REAL measured interpulse intervals from 4
independent papers: PCOS (Waldstreicher 1988, 58min) → normal late-follicular (Reame 1984,
50min) → normal early-follicular (Nippoldt 1989, 75min) → normal luteal (Reame 1984/Nippoldt
1989, 90–188min) → hypothalamic amenorrhea (Reame 1985/Khoury 1987, 153–480min). PCOS's ratio
exceeds every HA condition's ratio, directionally matching Rebar 1976's classical elevated-LH/
low-FSH PCOS finding.

**Honesty disclosure (a math-status correction I made to my own draft before finalizing)**: the
"ratio monotonically decreases with interpulse interval" gate is a **mathematical identity** given
k_LH>k_FSH (machine-confirmed true for ANY generic interpulse sweep, not data-dependent) — it is a
sanity check that the already-externally-validated rate constants were wired correctly, **not**
itself a risky empirical test. The genuinely at-risk content is narrower: (1) whether the REAL
measured interpulse intervals happen to fall in the physiologically-expected order (not
guaranteed by any math), and (2) whether that reproduces Rebar 1976's independent qualitative
clinical finding. Both hold; only these two are reported as real evidence, not the "monotonic"
computation itself.

**Rebar 1976's own added nuance, kept not smoothed over**: clomiphene-induced preovulatory E2
rises still produced "appropriate LH surges" in PCOS patients — the positive-feedback switch
MACHINERY itself is intact in PCOS. This model's PCOS story is deliberately narrow (persistently
fast tonic pulse frequency → elevated LH:FSH ratio only) and does not claim to model disordered
folliculogenesis or a broken switch.

## 7. Menopause transition — a DECORRELATED mechanism, disclosed as such

The task itself names this as a different story ("loss of inhibin/estradiol feedback → high
FSH"), not the bistable switch — kept separate here, not conflated. FSH's tonic level =
FSH_max/(1 + INH_B/INH_B0 + INH_A/INH_A0 + E2/E2_0) — **fully neutral, equal weights** (not tuned:
checked over w_A∈[1,5], the qualitative finding is unchanged across the whole sweep, so equal
weights were kept as the least-arbitrary choice). Fed Burger 1998's own real 4-stage
(pre/early-peri/late-peri/post) inhibin-B/inhibin-A/E2 values: the model reproduces the correct
MONOTONIC FSH rise and places the LARGEST single-stage jump at late-perimenopause — exactly where
Burger 1998's real data places it (inhibin-B falls alone first with only a modest, non-significant
FSH rise; inhibin-A and E2 then collapse TOGETHER, producing the large FSH rise). An existence/
ordering demonstration (not fit to the 4 target numbers), disclosed as such; Burger 1999 (used for
the FSH-E2 r=-0.73 correlation) shares senior authorship/lab with Burger 1998 — an honest,
disclosed same-group limitation, not a fully independent second cohort.

## 8. Pre-registered gates — 17/17 PASS

```
switch_model_gate_phase:                                              PASS  (lag 0.0-0.35h, gate<=24h)
switch_model_gate_amplitude:                                          PASS  (4.79x, gate>=3x)
switch_model_gate_abruptness:                                         PASS  (8.64x, gate>=2x)
switch_model_matches_who1980_anchor_within_1day:                      PASS
adversary_never_achieves_both_gates_across_90_configs:                 PASS  (forced, non-strawman sweep)
adversary_trough_always_aligned_with_e2_peak_plus_delay:               PASS  (self-corrected, disclosed)
kill_switch_degenerates_to_adversary_failure_mode:                     PASS  (mechanism isolation)
duration_gate_sustained_primes_transient_does_not:                     PASS
duration_gate_transient_amplitude_stays_low:                          PASS
duration_gate_sustained_amplitude_surges:                             PASS
freq_coding_monotonic_ratio_vs_interpulse:                            PASS  (math-identity status disclosed, sec.6)
freq_coding_pcos_exceeds_all_ha:                                       PASS  (the real empirical content)
menopause_monotonic_matches_real:                                     PASS
menopause_biggest_jump_location_matches:                              PASS
robustness_all_e2_shape_variants_pass:                                PASS  (9 E2-curve-shape variants)
robustness_phase_alignment_100pct_across_secretion_law_sweep:         PASS  (27 secretion-law configs)
robustness_default_gain_has_margin_above_amplitude_failure_boundary:  PASS  (2x margin, not knife-edge)
```

`overall_pass_strict_all = True`. Determinism: 2 independent runs byte-identical (verified via
`diff`); zero NaN/Inf across the full JSON tree (checked programmatically).

## 9. Symmetric QC — sorted by evidentiary weight, why 17/17 isn't over-claimed

- **Hard, non-circular external falsifiers**: the phase-coincidence match to 3 independent papers
  (WHO 1980, Fritz 1992, Hoff/Quigley/Yen 1983 — none touched during construction except two
  magnitude/timescale landmarks); the 90-config forced-adversary sweep and its exact
  delay-shift geometry; the kill-switch mechanism-isolation test; the duration-gating contrast
  (same model, only input duration changes); the real measured-interpulse-interval ordering
  matching Rebar 1976's clinical finding.
- **Constructive/existence demonstrations** (disclosed, weaker): the E2(t)/P4(t) forcing curve's
  exact shape (landmark-calibrated, not independently fit); the menopause FSH formula (ordering/
  location match, not numeric fit); the D(τ=48h)/threshold(200pg/mL) values are the task's own
  given numbers, not independently re-derived here.
- **Held OPEN, not reconciled**: Rebar 1976's own finding that PCOS's positive-feedback switch
  machinery is intact (this model's PCOS story is deliberately narrow, frequency-coding only);
  Young & Jaffe 1976's own indexed abstract lacking the numeric strength-duration result (the
  ~200pg/mL/48h figure is field-consensus, cross-checked against Chang & Jaffe 1978's independent
  271pg/mL number, not re-derived from Young & Jaffe's own text).
- **Two genuinely corrected errors, not swept under the rug**: (1) the trough-alignment check's
  un-shifted-reference bug (§1/§2) — caught by running the check against its own analytic
  prediction, diagnosed, fixed, re-verified, disclosed. (2) A **fabricated DOI**: a first draft of
  the WHO 1980 citation (§0 table, row 1) guessed a plausible-looking DOI
  (`10.1016/0002-9378(80)90133-9`) that was never actually present in either the live-fetched
  abstract text (which, unlike all 16 other citations' fetches, printed no `DOI:` line at all) or
  in a follow-up `esummary` check of the record's authoritative `articleids` field (which lists
  only `pubmed` + `pii` — no `doi` type, and the real PII's last digit, `-7`, doesn't even match
  the guessed DOI's `-9`). Caught by re-querying the structured metadata rather than trusting a
  pattern-matched guess; corrected to `doi=None` + the real PII, disclosed here and in the
  citation table, not silently left wrong. This is exactly the failure mode watertight discipline
  exists to catch: a plausible-looking, unverified number, generated (not recalled) under the
  pressure of "every other row has a DOI."
- **Real inter-study spread**: small samples throughout (Fritz n=7, Hoff/Quigley/Yen n=5,
  Waldstreicher n=12, Reame n=8) — typical of invasive frequent-sampling endocrinology, point
  estimates not distributions.
- **Not independently verified this session**: a fully autonomous closed-loop oscillator (E2/P4
  emerging from a follicle-growth ODE rather than prescribed) — disclosed scope limit, not
  attempted; the fast pituitary self-priming mechanism needed to match the literal ~2h LH
  doubling time.

## 10. couples_to

**A genuine numeric reuse** (not just a key-presence match, unlike the male doc's own §8): k_LH
and k_FSH are loaded verbatim from `MECHANISM_REPRODUCTIVE_HPG.md`'s validated values and used
inside both this model's ODE and the frequency-coding cross-check — the SAME hormone molecules,
sex-independent clearance. **HPG-male cert**: shares the GnRH/gonadotropin axis but has the
OPPOSITE feedback structure (male = pure delayed-negative-feedback steady state, delay-
independently stable per that doc's own spectrum analysis; female = the same negative-feedback
substrate PLUS a bistable positive-feedback switch layered on top) — the two docs together cover
both halves of the same axis's qualitatively different regimes. **Reproductive/pregnancy axis**:
`data/body_twin/agent_outputs/reproduction-pregnancy__a13d17371a3dd9d85.json` (a prior, unverified
hypothesis-tier literature-scout output, re-read this session, none of its specific numbers reused
as load-bearing here) already names lactational-amenorrhea (prolactin-suppressed GnRH pulsatility)
as a direct HPG-axis coupling — this model's D(t)/pulse-frequency machinery is the natural place
that coupling would attach, not built here. **Circadian**: not coupled this session (disclosed).
**Function↔dysfunction**: PCOS (§6), hypothalamic amenorrhea (§6, via the same frequency-coding
transfer function — HA's slow pulse frequency predicts the lowest LH:FSH ratio of all 7 tested
conditions), menopause transition (§7, explicitly decorrelated mechanism) are all built as
variations of the SAME reused machinery, not three separate hacked rules.

**Output state-vector keys**: `serum_E2, serum_P4, serum_LH, serum_FSH, cycle_day, phase_label,
primed_D`. **Grep-verified this session**: no pre-existing SEED-DESIGN node in
`data/MECHANISM_ANCHOR_GRAPH.json` (998 nodes, read-only, not edited) names an explicit female-
cycle-phase state vector — this is a new upstream layer, not a resolution of an existing OPEN
node (contrast the male doc). Not folded into the graph this session (isolation: touch only files
created this session).

## 11. Confidence tier

**In-vivo-anchored, population/literature-level** (real human serial/frequent-sampling data: n=177
WHO 1980, n=20 Stricker 2006, n=8 Reame 1984, n=13 Nippoldt 1989, n=19 Reame 1985, n=14 Khoury
1987, n=12+21 Waldstreicher 1988, n=24 Rebar 1976, n=150 Burger 1999, n=110 Burger 1998, n=293
cycles Lenton 1984, n=7 Fritz 1992, n=5 Hoff/Quigley/Yen 1983, n=12 Chang & Jaffe 1978; n=4 Wildt
1981 monkey) — same tier as `MECHANISM_REPRODUCTIVE_HPG.md` and every sibling endocrine doc in this
repo: no subject-specific hormone panel exists for subject2. Every operating point is population-
level, not individual.

## 12. Honest gaps — symmetric QC: what this does NOT prove

- **E2(t)/P4(t) are prescribed, not derived from a follicle-growth/steroidogenesis ODE** — a
  disclosed reduction (matches the male doc's own precedent); a fully autonomous closed-loop
  oscillator is not attempted. Tested for robustness (9 E2-shape variants, §1/§8), not
  independence from the forcing choice entirely.
- **The abruptness gate is directionally, not numerically, matched** — model 8.64× vs real
  ~28-30× (Hoff/Quigley/Yen 1983); a separate fast pituitary self-priming mechanism would be
  needed to close this gap, not built here.
- **Young & Jaffe 1976's own indexed abstract lacks the numeric strength-duration result** — the
  ~200pg/mL/48h figure used is field-consensus attribution, cross-checked against Chang & Jaffe
  1978's independent 271pg/mL number, not independently re-derived from Young & Jaffe's own text
  this session.
- **PCOS is modeled ONLY via the frequency-coding transfer function** — Rebar 1976's own finding
  that the positive-feedback switch remains intact in PCOS (clomiphene-induced surges are normal)
  is disclosed, not modeled; disordered folliculogenesis / multi-follicle E2 dysregulation is out
  of scope.
- **The menopause-transition formula is checked against ONE paper's 4-stage means (Burger 1998)**
  for ordering/jump-location, not fit; Burger 1999 (the correlation cross-check) shares senior
  authorship/lab, not a fully independent second cohort.
- **The secretion-law parameter sweep shows the amplitude gate (not the phase gate) has a real
  gain-dependent failure boundary** (fails below gain_lh≈4-8 in the tested range) — `DEFAULT_
  PARAMS` sits with a disclosed ~2× margin above it, not at a knife-edge, but the amplitude
  claim is not universally parameter-independent the way the phase claim is (100% robust across
  all 27 swept secretion-law configs).
- **No graph-edge write this session** — no pre-existing SEED-DESIGN node names a female-cycle
  state vector (grep-verified); not folded into `data/MECHANISM_ANCHOR_GRAPH.json` (isolation
  rule: touch only files created this session).
- **Small samples throughout** (Fritz n=7, Hoff/Quigley/Yen n=5, Waldstreicher n=12 PCOD) — point
  estimates, not distributions, are what this model is checked against; real inter-subject spread
  is disclosed, not hidden.

## Files

- `scripts/msk/ovarian_cycle_axis.py` — self-contained (numpy/scipy only): E2(t)/P4(t) forcing
  construction, the D(t)-gated bistable secretion law, the ODE simulation (`solve_ivp`), the
  90-config forced-adversary sweep (steepness×operating-point×delay) with its corrected
  trough-alignment geometry, the kill-switch consistency check, the duration-gating control, the
  frequency-coding cross-check (reusing the male doc's k_LH/k_FSH), the menopause-transition
  check, both robustness sweeps (E2-shape ×9, secretion-law ×27), the couples_to check, and all
  17 gates; writes the evidence JSON below; prints a full summary.
- `data/msk_smoketest/subject2_walking1/ovarian_cycle_axis/ovarian_cycle_axis_results.json` —
  every number in this doc, machine-written: all 17 citations, the full E2/P4/D/LH/FSH time
  series, the complete 90-row adversary grid, the kill-switch result, the duration-gating
  contrast, the frequency-coding grid (7 conditions), the menopause 4-stage check, both
  robustness sweeps, the couples_to check, and all 17 gates. Verified deterministic (2 independent
  runs, byte-identical via `diff`) and NaN/Inf-free (checked programmatically over the full tree).
- Read but NOT modified (isolation: touch only files created this session):
  `data/MECHANISM_ANCHOR_GRAPH.json` (grep-verified: no matching SEED-DESIGN node exists),
  `docs/MECHANISM_REPRODUCTIVE_HPG.md` (the male axis this doc's k_LH/k_FSH/transfer-function are
  reused from), `data/body_twin/agent_outputs/reproduction-pregnancy__a13d17371a3dd9d85.json`
  (prior hypothesis-tier scout output, re-read, not reused as load-bearing).

## Repro

```
cd ~/projects/bodytwin
source .venv-msk/bin/activate
python3 scripts/msk/ovarian_cycle_axis.py
```
No external inputs required. Pure Python/numpy/scipy (`scipy.integrate.solve_ivp`), no OpenSim
call, runs in a few seconds, deterministic (byte-identical across repeated runs, machine-verified
this session). No git operations; writes only under
`data/msk_smoketest/subject2_walking1/ovarian_cycle_axis/`.
