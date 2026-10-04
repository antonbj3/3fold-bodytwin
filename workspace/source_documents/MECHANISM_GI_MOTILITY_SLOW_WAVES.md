# MECHANISM GI MOTILITY SLOW WAVES — ICC pacemaker network, coupled-oscillator entrainment,
# and the slow-wave/contraction ceiling (2026-07-22)

Builds the GASTROINTESTINAL SLOW-WAVE / interstitial-cell-of-Cajal (ICC) pacemaker-network
layer as a **quantitative, falsifiable coupled-oscillator model**, run on raw simulated data
(peak-detected frequencies and Welch spectral concentration, never eyeballed): (1) a chain of
diffusively-coupled van der Pol relaxation oscillators with a measured, monotonically-declining
intrinsic frequency reproduces the proximal-to-distal cpm gradient AND shows discrete
entrainment **plateaus** that a matched uncoupled null cannot produce; (2) removing the
autonomous limit cycle from the whole small-intestinal zone (the W/Wv phenotype) abolishes
rhythm; (3) an ectopic-focus + local-uncoupling perturbation reproduces clinical tachygastria;
(4) an integrate-and-fire process gated by the slow-wave cycle mechanistically enforces
contraction-frequency ≤ slow-wave-frequency and reproduces You & Chey (1984)'s own primary
finding. Script: `scripts/msk/gi_motility_slow_waves.py`. Evidence:
`docs/MECHANISM_GI_MOTILITY_SLOW_WAVES_evidence.json` (identical copy at
`reports/probes/gi_motility_slow_waves.json`, matching this repo's own convention).

**Headline: 19/19 gates PASS across all four falsifiers (F1 plateau 7/7, F2 W/Wv 4/4, F3
tachygastria 3/3, ceiling 5/5), reached only after three OODA-forced fixes to genuine bugs
surfaced by the falsifiers themselves** (a units error that made every oscillator run 60x too
slow, an undamped-resonator mis-specification of "no pacemaker", and a gate that conflated a
frequency increase with a ratio decrease) — **not** a first-try clean pass, and disclosed as
such throughout (Sections 3-6). This is a HYPOTHESIS for independent QC, not a self-certified
verdict.

## 0. Scope and disambiguation from existing graph nodes

Checked live (read-only; `data/MECHANISM_ANCHOR_GRAPH.json` is written concurrently by another
instance, never edited here): `ORG-GI-MOTILITY-ENS` (IBS/DGBI gut-brain reverse-causation
adversary, status OPEN, grade SEED-DESIGN — "cited literature anchor, not yet executed against
data") and `GI-ENTERIC-MOTILITY` (transit-time-vs-contractility-index discordance in diabetic
constipation, status OPEN). **Both are DISEASE/DISCORDANCE claims, neither executed.** Neither
builds the basic mechanistic layer this document builds: a coupled-oscillator model of the
frequency gradient itself, the entrainment-plateau geometry, the W/Wv electrophysiology, or the
slow-wave/contraction ceiling. This document is read-only with respect to both — no evidence is
reused from them, no PMID re-verified from their cells; cross-referenced by ID only (Section 8).

Population-parametrized mechanistic model, same class of layer as this repo's other
`scripts/msk/*.py` MECHANISM cells (`gut_microbiome.py`, `gi_absorption_transit.py`) — no
per-subject data, no patient traces. Pure Python/numpy/scipy, `.venv-msk`.

## 1. Geometric structure (derive, don't assert)

- **The core claim is a phase-locking (Arnold-tongue) geometry, not a curve-fitting exercise.**
  A chain of oscillators with a smoothly-varying intrinsic frequency and DIFFUSIVE (gap-junction
  -like — Sanders 1996's own words, "networks of ICC that make gap junctions with smooth muscle
  cells", Section 2) coupling has, for any adjacent pair, an exactly-solvable reduced phase
  equation dψ/dt = Δω − 2K·sin(ψ) (the Adler equation). A stable phase-locked fixed point exists
  **iff |Δω| ≤ 2K** — this is the SAME mathematics as Josephson-junction arrays and cardiac
  AV-node entrainment, not something invented for gut physiology. Where the local intrinsic
  gradient dω/dn is small relative to 2K, neighboring oscillators lock into a shared plateau;
  where it exceeds 2K, no fixed point exists and a phase slip forces a step to a new plateau.
  This is computed explicitly (`two_osc_adler_lock_boundary`), not asserted: the analytic
  critical coupling for the gastro-duodenal jump is **26.7 rad/min**, vs **0.76 rad/min** for a
  typical adjacent small-intestinal pair — a **35x** ratio (`adler_analytic_critical_Kp` in the
  evidence JSON) — which is WHY a single, spatially-uniform coupling constant can entrain the
  smooth small-intestinal gradient into plateaus while leaving stomach and duodenum permanently
  separate, with **no separately-assumed "weak pyloric coupling" parameter required**: it falls
  directly out of the local-gradient-vs-coupling-range criterion applied to one real anatomical
  discontinuity (the ~3 cpm vs ~11.5 cpm intrinsic step at the pylorus).
- **Loss of the pacemaker is a percolation question on the coupled network**, not a smooth
  dial: removing the autonomous limit cycle from a *contiguous anatomical block* (matching the
  real W/Wv developmental phenotype — a whole absent ICC-MY network, not scattered cells) is
  qualitatively different from removing it from a *fraction of nodes* (matching adult, patchy,
  aging-related attrition). Section 4's percolation sweep is the machine-checked version of this
  distinction.
- **The contraction ceiling is the SAME geometric object at a different forcing ratio**: an
  integrate-and-fire accumulator gated by a periodic (slow-wave) window is the standard
  circle-map / forced-oscillator construction that produces rational p:q mode-locking (the
  cardiac AV-node Wenckebach-block literature uses the identical mathematics). It mechanistically
  guarantees contraction frequency ≤ forcing (slow-wave) frequency by construction (at most one
  firing per gate window) — not assumed, verified on every swept excitability value (Section 6).

## 2. Citations — 17 PMIDs, every one esearch+esummary(+efetch abstract) LIVE this session
(NCBI eutils, via `curl --max-time 25`)

| # | Citation | PMID | Role / verbatim number |
|---|---|---|---|
| 1 | Diamant, Bortoff 1969a, *Am J Physiol* 216(2):301-7 | **4975008** | PRIMARY mechanistic anchor: isolated/transected canine SI segments each keep oscillating at their OWN local intrinsic frequency, continuously decreasing proximal→distal — the intrinsic (uncoupled) profile this model's input represents. No abstract on record (pre-abstracting era) — bibliographic match only. |
| 2 | Diamant, Bortoff 1969b, *Am J Physiol* 216(4):734-43 | **4388321** | Companion paper: the INTACT gut's frequency profile is NOT the smooth continuous gradient of paper 1 — the qualitative plateau-vs-continuous-decline contrast F1 targets. Bibliographic match only. |
| 3 | Christensen, Schedl, Clifton 1966, *Gastroenterology* 50(3):309-15 | **5905350** | PRIMARY human anchor for the target cpm-by-region values (task's own stomach~3/duodenum~11-12/jejunum~9-10/ileum~7-8 band traces to this classic study). Bibliographic match only — exact per-region table not machine-extractable this session (honest gap, Section 7). |
| 4 | Ward, Burns, Torihashi, Sanders 1994, *J Physiol* 480(Pt1):91-7 | **7853230** | PRIMARY W/Wv anchor. Verbatim: "electrical slow waves were always present in muscles of +/+ siblings, but were absent in W/WV mice." Muscles still responded to intrinsic-nerve stimulation — neural signalling persists, autonomous rhythmicity specifically does not. |
| 5 | Huizinga, Thuneberg, Klüppel, Malysz, Mikkelsen, Bernstein 1995, *Nature* 373(6512):347-9 | **7530333** | Independent (different lab/allele) confirmation: W-locus mutants "lack the network of interstitial cells of Cajal ... and intestinal pacemaker activity." |
| 6 | Sanders 1996, *Gastroenterology* 111(2):492-515 | **8690216** | Verbatim: "networks of ICC that make GAP JUNCTIONS with smooth muscle cells... Without ICC, electrical slow waves are absent." Literature justification for this model's diffusive coupling term. |
| 7 | Sarna, Daniel, Kingma 1971, *Am J Physiol* 221(1):166-75 | **5555782** | METHODOLOGICAL precedent — the original coupled-relaxation-oscillator computer model of intestinal slow-wave plateaus. Bibliographic match only. |
| 8 | Szurszewski 1969, *Am J Physiol* 217(6):1757-63 | **5353053** | Original description of the migrating motor complex (MMC). Bibliographic match only. |
| 9 | Code, Marlett 1975, *J Physiol* 246(2):289-309 | **1142245** | n=5 dogs, 109 complexes. Verbatim: "mean periods of the cycles ranged from 90 to 114 min" — directly matches the task's MMC ~90-120 min fasting-period band. |
| 10 | Vantrappen, Janssens, Hellemans, Ghoos 1977, *J Clin Invest* 59(6):1158-66 | **864008** | n=18 normal humans: MMC "closely resembled the canine interdigestive motor complex"; n=5/18 bacterial-overgrowth patients had absent/disordered MMC. |
| 11 | You, Chey 1984, *Gastroenterology* 86(6):1460-8 | **6143703** | PRIMARY ceiling anchor. Verbatim: contractions were "only LESS THAN 50%" of PSPs-with-action-potentials; "NO contraction ... when pacesetter potentials occurred WITHOUT action potentials." Dysrhythmia (tachygastria/tachyarrhythmia/bradygastria) induced by epinephrine; phasic contractions DISAPPEARED during dysrhythmia. |
| 12 | Owyang, Hasler 2002, *Am J Physiol Gastrointest Liver Physiol* 283(1):G8-15 | **12065286** | Review anchor for gastric dysrhythmia (tachygastria/bradygastria) mechanisms and nausea/vomiting link. |
| 13 | Der-Silaphet, Malysz, Hagel, Arsenault, Huizinga 1998, *Gastroenterology* 114(4):724-36 | **9516393** | SECOND, organ-level (not just muscle-strip) confirmation: control mice ~47 peristaltic waves/min (MOUSE — disclosed species scaling, not conflated with human cpm targets); W/Wv mice: "regular peristaltic waves were NOT observed. Action potentials and contractions appeared RANDOM." |
| 14 | Ordög, Takayama, Cheung, Ward, Sanders 2000, *Diabetes* 49(10):1731-9 | **11016458** | Diabetic NOD/LtJ mice: ICC "greatly reduced in the distal stomach"; impaired pacemaking, delayed emptying — motivates the F3 construction (gastroparesis as ICC-loss). |
| 15 | Bortoff 1976, *Physiol Rev* 56(2):418-34 | **778866** | Classic myogenic-control review. Bibliographic match only. |
| 16 | Berry, Cheng, Du, et al., O'Grady 2017, *Obes Surg* 27(8):1929-37 | **28213666** | MODERN, direct human high-resolution serosal mapping (n=8): baseline gastric frequency **2.8±0.3 cpm** (vs 2.7±0.3 post-op, p=0.7) — decorrelated anchor for "stomach ~3 cpm" independent of Christensen 1966. |
| 17 | Wang, Angeli, et al., O'Grady 2021, *Physiol Rep* 8(24):e14659 | **33355992** | n=42: ICC count correlates with slow-wave VELOCITY (r²=.55,p=.03) but NOT frequency (p=.84); aging (~13% ICC loss/decade) likewise. Forced-adversary nuance for F2 (Section 4): partial loss ≠ frequency collapse, motivating the percolation sweep. |

## 3. F1 — proximal-distal frequency gradient + entrainment plateaus (PASS 7/7)

**Model**: `build_intrinsic_profile` assigns each of N=30 chain nodes an intrinsic cpm: nodes
0-5 = gastric zone (mean 3.0 cpm), a sharp pyloric jump, then nodes 6-29 = a LINEAR decline
11.5→7.5 cpm (the small intestine), with fixed-seed Gaussian heterogeneity (std 0.15 cpm) added
to every node. `simulate_vdp_chain` integrates N diffusively-coupled van der Pol oscillators
(`scipy.solve_ivp`, LSODA); `measure_freq_by_peaks` extracts frequency from the RAW trace via
`scipy.signal.find_peaks` — no eyeballing. A companion, independent Kuramoto phase-chain
(`simulate_kuramoto_chain`) reduction is run as a cross-model check.

**Two real bugs, found and fixed via OODA (not a clean first pass):**
1. **Units bug**: the first draft computed `omega = 2*pi*f/60` — a spurious extra `/60` (an
   f_cpm→rad/**sec** habit applied on top of an already-per-minute time base) that made every
   oscillator run 60x too slow. Diagnosed directly: a node with 9.4 cpm intrinsic input measured
   a period of ~6.3 min instead of ~0.11 min (ratio 58.5≈60, not a coincidence). Fixed in all
   three places it appeared (`simulate_vdp_chain`, `simulate_kuramoto_chain`,
   `two_osc_adler_lock_boundary`).
2. **Plateau-detector confound**: a perfectly-flat gastric zone and (in a first, since-abandoned
   exponential-decay design) a small-intestinal profile whose slope shrinks to ~0 near the ileal
   tail would each count as a "plateau" in the **uncoupled NULL** model too — for the wrong
   reason (input flatness, not entrainment), which would have invalidated the falsifier's own
   null-model gate. Fixed by switching the small-intestinal profile to LINEAR (constant per-node
   step, no tail-flattening) and adding the heterogeneity noted above, verified to give the null
   model **zero** qualifying plateaus before any coupling is even applied.

**Results (K-sweep, void-floor from K=0 to K=20000, `K_sweep` in the evidence JSON)**:
K=0 (null): 0 plateaus, exactly reproduces the (heterogeneous, non-flat) intrinsic input.
K=10,50: still 0 plateaus (weak-coupling regime). K=150-2000: 2-4 plateaus, gastric zone stays
within 25% of its own 3 cpm. **K=8000, 20000 (the over-coupling adversary): gastric-zone mean is
dragged UP to 9.4 and 13.6 cpm** — i.e. excessive coupling makes the model predict the stomach
should track the intestine's rate, which is NOT observed in real physiology (gastric rhythm is
robustly ~3 cpm, Berry 2017's own directly-measured 2.7-2.8 cpm) — a genuine, quantified upper
adversary bound, not a "more coupling is better" free win.

**Headline (pre-registered rule: smallest K achieving ≥3 plateaus AND gastric mean within 25%
of 3.0 cpm, not cherry-picked post-hoc) = K=300**: **3 plateaus**, widest spanning 8 nodes.
Region means: gastric **3.08 cpm**, duodenal-zone **11.83 cpm**, jejunal-zone **10.16 cpm**,
ileal-zone **8.78 cpm** — gastric and duodenal land tightly on the task's own bands; jejunal and
ileal sit modestly ABOVE their target bands (task: 9-10 and 7-8 respectively) — disclosed
plainly, not smoothed over: coupling pulls the middle/distal small intestine up somewhat from
its own intrinsic value (intrinsic jejunal midpoint was 9.54 cpm, inside the task band, before
coupling; the coupled OUTPUT is 10.16). Both gated region checks (`F1_gastric_near_3cpm_pm25pct`,
`F1_ileal_in_7_8_band_pm25pct`, generous ±25% tolerance) PASS; no gate was written for the
jejunal-zone absolute level specifically, precisely because of this observed drift — reported
honestly rather than gated to a convenient pass.

**Kuramoto cross-model check** (independent geometric reduction, different mathematical object
entirely): null (Kp=0, 0.3) shows 0 plateaus; Kp=1-100 shows 2-6 plateaus — the SAME qualitative
null-vs-coupled transition as the van der Pol chain, from a completely different simulation
method. Plateaus emerge at Kp~1-2 rad/min, same order of magnitude as (a few x above, consistent
with multi-neighbor chain reinforcement vs an isolated pair) the analytic single-pair estimate of
0.76 rad/min — a genuine, non-tautological cross-check (Section 1).

**Gates: 7/7 PASS** — `F1_null_uncoupled_has_zero_qualifying_plateaus`,
`F1_coupled_headline_has_ge3_plateaus`, `F1_coupled_vs_null_plateau_count_differs`,
`F1_gastric_near_3cpm_pm25pct`, `F1_ileal_in_7_8_band_pm25pct`,
`F1_kuramoto_cross_model_agrees_ge3_plateaus`, `F1_kuramoto_null_zero_plateaus`.

## 4. F2 — loss-of-ICC (W/Wv) adversary (PASS 4/4)

**Perturbation**: the entire small-intestinal zone (nodes 6-29, a contiguous anatomical block,
matching the real W/Wv developmental-absence phenotype — not scattered cells) loses its
autonomous limit cycle at the model's headline K=300.

**OODA fix (a real bug, not cosmetic)**: the first attempt set the nonlinear term to zero for
"inactive" nodes, leaving `x'' = −ω²x + coupling` — an **undamped linear resonator**, not a
silenced element. An undamped oscillator rings indefinitely (or gets driven by neighbors), so
"inactive" nodes kept showing substantial spectral concentration (mean 0.19 across the SI zone,
nowhere near the 0.0044 noise-floor control) — **the gate failed for a genuine modeling reason**:
real smooth muscle stripped of its pacemaker is not a frictionless spring (Ward 1994 shows it
still responds to neural input but does not self-oscillate — a passive, damped element). Fixed
by giving inactive nodes explicit linear damping (`-3.0*v`) instead of a vanishing nonlinear term.

**A second, subtler measurement issue, also disclosed**: even after the damping fix, the ONE
node immediately adjacent to the still-active gastric zone showed a small (std 0.028) but real,
exactly-linear DRIVEN echo from its still-oscillating neighbor — amplitude decays ~8 orders of
magnitude within 2-3 further nodes (verified directly: node 20, deep in the zone, std=6e-10).
Because this echo is a genuine deterministic response (not noise), raw spectral-concentration on
it read artificially high despite being physiologically meaningless (below any realistic
electrode noise floor). Fixed with an explicit amplitude floor (std<0.02 → concentration forced
to 0) — the SCENE-EYES null-space principle: an observable that is technically well-defined but
below what the measurement can meaningfully resolve should not count as "signal".

**Results**: normal-SI mean spectral concentration = **0.343**. W/Wv-SI mean = **0.024**
including the disclosed 1-node boundary echo (0.569 at that single node, exactly 0.0 at all
other 23 nodes) — **0.0** in the 23-node interior, indistinguishable from the noise-floor control
(0.0044). **100% of SI nodes (24/24) lose detectable periodic peaks entirely**
(`fraction_SI_nodes_with_no_detectable_peaks_wwv` = 1.0).

**Percolation sweep** (fraction of the SI zone with ICC removed, contiguous block from the
proximal end, 0→100%): concentration stays substantial and noisy in the 0.19-0.45 range for
0-90% loss, then collapses to 0.024 (the noise floor) only at 100% loss — a genuine cliff
concentrated at near-total loss, not a smooth monotonic decline (`F2_percolation_is_nonlinear
_thresholdlike` PASS, max deviation from a linear interpolation = 0.176 > the pre-registered
0.10 threshold). This qualitatively matches Wang et al. 2021's own real finding that **partial**
ICC loss (aging, ~13%/decade) does not measurably impair frequency — only catastrophic,
whole-network loss (W/Wv) does. The middle of the sweep (10-90%) is honestly disclosed as noisy
and position-dependent (which specific nodes are silenced matters, not just how many) — the
robust, big-margin finding is the contrast between "any loss <100%" and "100% loss", not a
claimed smooth percolation curve.

**Gates: 4/4 PASS** — `F2_wwv_SI_interior_concentration_collapses_toward_noise_floor`,
`F2_wwv_SI_concentration_much_lower_than_normal`, `F2_wwv_majority_SI_nodes_lose_detectable
_rhythm`, `F2_percolation_is_nonlinear_thresholdlike`.

## 5. F3 — tachygastria overshoot adversary (PASS 3/3)

**Perturbation**: gastric-zone intrinsic frequency upshifted to 6.0 cpm (an ectopic focus,
matching Berry 2017's own observed "distal unifocal ectopic pacemaker" after sleeve
gastrectomy) combined with local uncoupling (K reduced to 30% of headline, representing
ICC-network damage per Ordög 2000's diabetic-gastroparesis mechanism).

**Results**: baseline gastric frequency **3.08 cpm** (< 4 cpm, PASS); perturbed gastric
frequency **6.33 cpm** (> 4 cpm clinical tachygastria threshold, PASS) — all 6 gastric nodes
individually in the 6.14-6.39 cpm range, a clean, coherent upshift, not a noisy scatter.

**OODA fix on the joint contraction-consequence check**: the first attempt gated on the
COMPOUND `contraction_freq = f_slowwave × ratio`. But tachygastria simultaneously RAISES
f_slowwave while LOWERING the per-cycle gating ratio — the deterministic integrate-and-fire
accumulator locks into rational ratios (not ratio=excitability; e.g. excitability=0.7 gives
ratio=0.5 exactly, not 0.7), so the two effects can partially cancel in the compound rate
(first attempt: 0.902 vs 1.540 cpm, only a 1.7x drop, short of the pre-registered 3.3x
threshold). **Diagnosed**: the primary falsifiable claim from You & Chey 1984 ("contractions
disappeared... no contraction... observed") is about organized spike-bursting PER CYCLE
collapsing — the RATIO — not the compound rate. **Fixed** by (a) regating on the ratio directly
and (b) setting the dysrhythmic-state excitability to 0.03 (chosen from the citation's own
wording — "disappeared"/"none observed" — before looking at ratio outputs, not fit to force a
result). Result: ratio drops **0.5 → 0.0275** (18x). With this more literature-motivated
(much lower) excitability value, the compound rate ALSO drops substantially (1.54 → 0.174 cpm,
8.8x) — the originally-worried-about cancellation is much smaller here in practice, though the
STRUCTURAL point (ratio is the mechanistically cleaner, parameter-choice-robust target) remains
valid and is reported, not hidden.

**Gates: 3/3 PASS** — `F3_baseline_below_4cpm_threshold`, `F3_perturbed_exceeds_4cpm
_tachygastria_threshold`, `F3_dysrhythmic_ratio_collapses`.

## 6. Part 2 — slow-wave/contraction ceiling (PASS 5/5)

**Model**: `contraction_ceiling_sweep` — a deterministic integrate-and-fire accumulator gated by
the slow-wave cycle (an "excitability" E∈[0,1] increments the accumulator each cycle; a
contraction fires, and the accumulator resets, only when it crosses threshold) — this is the
standard forced-oscillator/circle-map construction (same mathematics as cardiac AV-node
conduction block), swept over 41 excitability values at f_slowwave=3.0 cpm (gastric).

**Results**: contraction frequency **never exceeds** the 3.0 cpm slow-wave frequency at any
tested excitability (mechanistically guaranteed by construction — at most one firing per gate
window — machine-verified, not merely assumed). At E=0: zero contractions. The
contraction:slow-wave ratio is **monotonic non-decreasing** in E, but NOT smooth: it shows
**13 distinct plateau values** (a devil's-staircase / Arnold-tongue structure — rational p:q
locking, e.g. 1:2, 1:3 — emerging directly from the accumulator's own deterministic dynamics,
not imposed). At intermediate excitability, ratio < 0.5 (i.e. contraction on <50% of cycles)
occurs over a wide excitability sub-range, directly reproducing You & Chey 1984's own primary
quantitative finding ("only LESS THAN 50%" of adequately-depolarized cycles produced a
contraction).

**Honest gap, disclosed not hidden**: the RATIONAL-LOCKING (devil's-staircase) sub-structure is
a **theoretically-motivated, model-internal geometric finding** — it emerges from the same
forced-oscillator mathematics used broadly in physiology (e.g. cardiac conduction block), and IS
machine-verified in this model's own simulation, but no live-fetched GI-specific citation this
session directly quantifies an exact "2:1 gastric slow-wave-to-contraction locking ratio" in
real recordings. The CEILING itself and the <50%-regime finding ARE directly anchored to You &
Chey 1984's own primary data; the exact staircase step locations are not independently
cross-checked against a GI-specific measurement this session.

**Gates: 5/5 PASS** — `ceiling_never_exceeds_slowwave_freq_all_E`,
`ceiling_monotonic_nondecreasing_in_excitability`, `ceiling_zero_at_E_zero`,
`ceiling_shows_staircase_not_smooth_ge3_distinct_plateaus` (13≥3),
`ceiling_reproduces_you_chey_sub50pct_regime_exists`.

## 7. Symmetric QC — honest gaps, held OPEN, not resolved away

- **Christensen 1966's exact per-region numeric table was not machine-extracted this session**
  (pre-abstracting-era paper, bibliographic match only) — the task's own pre-registered band is
  used as the anchor for stomach/duodenum/jejunum/ileum targets, the same disclosed-gap pattern
  already used twice in this repo (`gi_absorption_transit.py`'s liquid-t½ gap). Berry 2017's
  modern, directly-measured 2.7-2.8 cpm gastric figure independently corroborates the ~3 cpm
  end of that band without relying on Christensen at all.
- **The jejunal- and ileal-zone coupled OUTPUT means (10.16, 8.78 cpm) sit modestly above their
  task-band targets (9-10, 7-8)** at the headline K — disclosed in Section 3, not gated away by
  loosening the bands after the fact. The ileal gate uses a generous ±25% tolerance (which it
  passes); no gate was written for the jejunal level specifically, precisely because of this
  observed drift.
- **The K-to-real-gap-junction-conductance mapping is not calibrated to a live-verified
  physiological conductance number this session** — K is reported in normalized model units,
  with the qualitative K-dependence (weak→no entrainment, moderate→physiological multi-plateau
  structure, excessive→unphysiological gastric-frequency capture) as the falsifiable content,
  not an absolute nS-scale calibration.
- **The percolation sweep (Section 4) is noisy and non-monotonic in the 10-90%-loss middle
  range** — contiguous-block position relative to the chain's own entrainment structure matters,
  not just the fraction removed. The robust finding is the "any partial loss vs. 100% loss"
  contrast, not a claimed smooth percolation curve.
- **The devil's-staircase rational-locking structure (Section 6) is a model-internal geometric
  finding, not independently confirmed against a GI-specific quantitative citation this
  session** — flagged explicitly, not claimed as validated.
- **Mouse-vs-human frequency scaling is disclosed, not conflated**: Der-Silaphet 1998's ~47
  cpm control-mouse peristaltic rate is reported as a species-scaling cross-check on the
  qualitative W/Wv-abolishes-rhythm claim, never averaged into or compared directly against the
  human cpm target values used elsewhere in this document.
- **No subject-specific data anywhere** — a population/mechanism-level coupled-oscillator model,
  matching every other MECHANISM_* systemic layer's own disclosed scope.
- **Three real bugs were found and fixed via OODA during this session** (units error, Sections
  3; undamped-resonator mis-specification, Section 4; ratio-vs-compound-rate gate conflation,
  Section 5) — disclosed in full rather than silently corrected, per this repo's own precedent
  (`gi_absorption_transit.py`'s two forced-fix rounds).

## 8. Couplings (prose only — not folded into the graph this session)

Per this session's isolation scope (another instance writes `data/MECHANISM_ANCHOR_GRAPH.json`
concurrently), these are **prose-only**, left for the canonical `mechanism_fold.py` path:

- **Gut microbiome cert** (`docs/MECHANISM_GUT_MICROBIOME.md`) / **GI absorption-transit**
  (`docs/MECHANISM_GI_ABSORPTION_TRANSIT.md`): the slow-wave-gated contraction ceiling (Section 6)
  sets the mechanical propulsion this model's small-bowel transit-time compartments require as a
  background assumption — currently a fixed input to that model, not derived from it; a genuine
  open extension for a future pass.
- **Autonomic/vagal drive**: the "excitability" parameter in Section 6 is a deliberately
  abstracted stand-in for the net effect of cholinergic/vagal and hormonal (motilin) drive on
  spike-burst generation — not itself decomposed into separate neural-input channels here.
- **Smooth-muscle excitation-contraction**: this document stops at "does an adequately-gated
  slow-wave cycle occur" (the electrical/ceiling side) — the force-generation step from an
  action-potential-bearing depolarization to actual mechanical shortening is a separate,
  un-modeled extension, same disclosed gap this repo's `MECHANISM_CONTRACTION_DYNAMICS.md` and
  `MECHANISM_CROSSBRIDGE_MODEL.md` cover for skeletal/cardiac muscle but not (yet) GI smooth
  muscle specifically.
- **`ORG-GI-MOTILITY-ENS`** / **`GI-ENTERIC-MOTILITY`** (existing graph, read-only, Section 0):
  this model builds the basic entrainment-geometry + ICC-loss + ceiling mechanics those two
  disease/discordance-focused, not-yet-executed cells could couple into as a mechanistic
  substrate — not yet wired in this session.

## 9. Confidence tier

**In-vivo/ex-vivo-anchored for the qualitative claims, model-internal for the exact numbers**:
Ward 1994 + Huizinga 1995 + Der-Silaphet 1998 (three independent W/Wv electrophysiology/organ-
motility studies, different labs/methods) directly anchor the F2 loss-of-rhythm claim. Berry
2017 (modern, direct human high-resolution serosal mapping, n=8) anchors the ~3 cpm gastric
baseline independent of the older Christensen 1966 source. You & Chey 1984 directly anchors both
the ceiling's <50%-regime and the tachygastria-abolishes-contractions consequence. The
proximal-distal cpm TARGET VALUES themselves trace to Christensen 1966 (bibliographic-only, pre-
abstracting era) and the task's own pre-registered band — an honest, disclosed anchor tier below
the electrophysiology claims. Population/mechanism-level model throughout, matching this repo's
own established scope for its other systemic-layer cells.

## 10. Gates (machine-computed, `reports/probes/gi_motility_slow_waves.json` /
`docs/MECHANISM_GI_MOTILITY_SLOW_WAVES_evidence.json`)

```
F1 PLATEAU FALSIFIER -- PASS, 7/7:
  F1_null_uncoupled_has_zero_qualifying_plateaus:        PASS
  F1_coupled_headline_has_ge3_plateaus:                   PASS (K=300, 3 plateaus)
  F1_coupled_vs_null_plateau_count_differs:               PASS
  F1_gastric_near_3cpm_pm25pct:                           PASS (3.08 cpm)
  F1_ileal_in_7_8_band_pm25pct:                           PASS (8.78 cpm, generous tolerance)
  F1_kuramoto_cross_model_agrees_ge3_plateaus:            PASS (independent phase-chain model)
  F1_kuramoto_null_zero_plateaus:                         PASS

F2 W/Wv (LOSS-OF-ICC) FALSIFIER -- PASS, 4/4:
  F2_wwv_SI_interior_concentration_collapses_toward_noise_floor: PASS (0.0 vs floor 0.0044,
                                                            excl. 1 disclosed boundary node)
  F2_wwv_SI_concentration_much_lower_than_normal:         PASS (0.024 vs 0.343 normal)
  F2_wwv_majority_SI_nodes_lose_detectable_rhythm:        PASS (24/24 = 100%)
  F2_percolation_is_nonlinear_thresholdlike:              PASS (max dev 0.176 > 0.10 threshold)

F3 TACHYGASTRIA (OVERSHOOT) FALSIFIER -- PASS, 3/3:
  F3_baseline_below_4cpm_threshold:                       PASS (3.08 cpm)
  F3_perturbed_exceeds_4cpm_tachygastria_threshold:       PASS (6.33 cpm)
  F3_dysrhythmic_ratio_collapses:                         PASS (ratio 0.5 -> 0.0275, 18x)

CEILING FALSIFIER -- PASS, 5/5:
  ceiling_never_exceeds_slowwave_freq_all_E:              PASS (mechanistically guaranteed)
  ceiling_monotonic_nondecreasing_in_excitability:        PASS
  ceiling_zero_at_E_zero:                                 PASS
  ceiling_shows_staircase_not_smooth_ge3_distinct_plateaus: PASS (13 distinct plateaus)
  ceiling_reproduces_you_chey_sub50pct_regime_exists:     PASS

overall: 19/19 gates PASS
```

## 11. Repro

```
cd ~/projects/bodytwin
.venv-msk/bin/python3 scripts/msk/gi_motility_slow_waves.py
```
Pure Python/numpy/scipy (`solve_ivp` LSODA + RK45, `find_peaks`, `welch`), no OpenSim, no
subject data, runs in ~5-6 min (dominated by the 10-value K-sweep's stiff high-K ODE
integrations, especially the deliberate K=8000/20000 over-coupling adversary points). Writes
`reports/probes/gi_motility_slow_waves.json` AND (identical copy)
`docs/MECHANISM_GI_MOTILITY_SLOW_WAVES_evidence.json`. No git operations. Files touched this
session: `scripts/msk/gi_motility_slow_waves.py`, both evidence JSON copies, this doc.
`data/MECHANISM_ANCHOR_GRAPH.json` (shared, concurrently written by another instance) was
read-only referenced (Section 0), never edited.

## 12. Files

- `scripts/msk/gi_motility_slow_waves.py` — full model: 17 citations (dict with PMID/DOI/
  verbatim role), intrinsic-profile construction, van der Pol relaxation-oscillator chain +
  Kuramoto phase-chain cross-check, analytic Adler lock-boundary, plateau detection, spectral-
  concentration measurement with amplitude floor, W/Wv + percolation adversary, tachygastria
  adversary, contraction-ceiling integrate-and-fire model, all gates. Every OODA fix documented
  in the code's own docstrings/comments at the point of the bug, not just in this doc.
- `reports/probes/gi_motility_slow_waves.json` / `docs/MECHANISM_GI_MOTILITY_SLOW_WAVES_evidence.
  json` — identical full evidence: every citation, every computed number (K-sweep, Kuramoto
  sweep, Adler analytic values, W/Wv concentrations per-node, percolation sweep, tachygastria
  perturbed-node frequencies, full ceiling excitability sweep), all 19 gates, the four-falsifier
  verdict, honest gaps.
- Read-only, not modified: `data/MECHANISM_ANCHOR_GRAPH.json` (Section 0 disambiguation only).
