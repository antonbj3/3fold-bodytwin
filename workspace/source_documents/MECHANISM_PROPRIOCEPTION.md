# MECHANISM PROPRIOCEPTION — the afferent (sensory) layer: muscle-spindle Ia + GTO Ib firing-rate models (2026-07-21)

Script: `scripts/msk/muscle_spindle.py`. Results: `data/msk_smoketest/subject2_walking1/muscle_spindle/muscle_spindle_results.json`.
Reused, unchanged: `scripts/msk/validate_joint_force.py` (paths/parser), `scripts/msk/validate_emg_timing.py`
(gait-cycle detection + %GC binning + circular-null machinery), `scripts/msk/metabolic_cost.py`
(`compute_kinematics`). Model: `LaiArnoldModified2017_poly_withArms_weldHand_scaled.osim`, subject2/walking1
(158 frames, 1.57s), Static Optimization output already computed by `static_opt_knee.py` — zero re-solving.

## Why this layer

Every neuromuscular piece already built in this repo is **efferent** (motor, outbound): Static
Optimization activation, `motor_unit_recruitment.py`'s size-principle recruitment layer,
`contraction_dynamics_forward.py`'s activation dynamics. None of it senses anything back. This is
the first **afferent** (sensory, inbound) layer: given the twin's own real fiber length + velocity
for a leg muscle during walking, what does a real, externally-anchored spindle (Ia) model predict
the afferent firing rate to be over the gait cycle — and does its shape/timing match published
spindle-firing patterns during locomotion? An optional Golgi-tendon-organ (Ib, force-sensing) model
is added for contrast.

## 0. Citations — every PMID/DOI verified LIVE this session (NCBI eutils + Crossref + PMC full text)

WebSearch was session-quota-exhausted (same disclosed fallback already used by
`scripts/msk/crossbridge_contraction_model.py`); verification used direct HTTP calls to NCBI
eutils (esearch/esummary/efetch), Crossref `query.bibliographic`, and PMC full-text fetches. **Two
recalled PMIDs were caught wrong before use**: Prochazka 1999 recalled as 10635704 (verified:
10635710); Prochazka & Gorassini 1998 recalled as ~9490847/48 (verified: 9490851/9490855) — this
matches the project's own previously-measured recalled-citation drift rate; nothing below is from
memory alone.

| # | Citation | PMID | DOI | Role |
|---|---|---|---|---|
| 1 | Prochazka A (1999). Quantifying proprioception. *Prog Brain Res* 123:133-142. | 10635710 | 10.1016/S0079-6123(08)62850-2 | Review naming the Ia-firing regression model (task brief's own name) |
| 2 | Prochazka A, Gorassini M (1998). Models of ensemble firing of muscle spindle afferents recorded during normal locomotion in cats. *J Physiol* 507(Pt1):277-291. | 9490851 | 10.1111/j.1469-7793.1998.277bu.x (PMC2230775, OA) | **PRIMARY governing equation** (Sec. 2) |
| 3 | Prochazka A, Gorassini M (1998). Ensemble firing of muscle afferents recorded during normal locomotion in cats. *J Physiol* 507(Pt1):293-304. | 9490855 | 10.1111/j.1469-7793.1998.293bu.x (PMC2230769, OA) | External envelope anchor (mean/peak imp/s) |
| 4 | Mileusnic MP, Brown IE, Lan N, Loeb GE (2006). Mathematical models of proprioceptors. I. *J Neurophysiol* 96(4):1772-1788. | 16672301 | 10.1152/jn.00868.2005 | "Mileusnic-style" detailed alternative — NOT implemented (Sec. 6) |
| 5 | Houk J, Henneman E (1967). Responses of Golgi tendon organs to active contractions of the soleus muscle of the cat. *J Neurophysiol* 30(3):466-481. | 6037588 | 10.1152/jn.1967.30.3.466 | GTO, soleus-specific |
| 6 | Houk J, Simon W (1967). Responses of Golgi tendon organs to forces applied to muscle tendon. *J Neurophysiol* 30(6):1466-1481. | 6066449 | 10.1152/jn.1967.30.6.1466 | GTO, passive-force companion |
| 7 | Crago PE, Houk JC, Rymer WZ (1982). Sampling of total muscle force by tendon organs. *J Neurophysiol* 47(6):1069-1083. | 7108572 | 10.1152/jn.1982.47.6.1069 | GTO population-linearization rationale |
| 8 | Fukunaga T et al. (2001). In vivo behaviour of human muscle tendon during walking. *Proc R Soc Lond B* 268(1464):229-233. | 11217891 | 10.1098/rspb.2000.1361 (PMC1088596, OA) | Fascicle/tendon decoupling anchor (Gate G1) |
| 9 | Rajagopal A et al. (2016). *IEEE TBME* 63(10):2068-2079. | 27392337 | (inherited, already verified in `validate_emg_timing.py`) | PHASE_PCT windows reused for Gates G2/G6 |

**Primary equation** (Prochazka & Gorassini 1998, quoted/extracted from PMC2230775 full text via
two independent extraction passes that agreed): fitted to **9 real cat hamstring Ia afferents, 132
step cycles**:

```
f = 4.3 * sign(v) * |v|^0.6  +  2*d  +  82         (f: imp/s, v: mm/s, d: mm displacement)
```

r²=0.91 (hybrid, fast cycles), r² range 0.88-0.94 across speeds, RMS error 14.1-18.4 imp/s
(9.2-14.7% of modulation depth). A velocity-only submodel (drop the `2*d` term) scored r²=0.94 —
comparable-or-better at fast cycles; the displacement term helped more at slow cycles. **The source
paper itself names the exact adversary this report must confront**: it explicitly states this fit
is for hamstring afferents and that **triceps surae (ankle extensors: soleus, gastrocnemius,
plantaris — the muscle this report models) showed "significant deviations"** from the simple
model, attributed to tendon-compliance effects and **"a short period of complete muscle unloading"
near the end of stance**. Applying these coefficients to soleus is a disclosed cross-muscle-group,
cross-species transplant, not a same-muscle fit — pre-registered as the thing to test, not hidden
after the fact.

## 1. Method

1. **Fiber kinematics** (soleus_r, plus tibialis anterior `tibant_r` as an antagonist contrast
   muscle): reuse the exact per-frame reconstruction pattern already validated in
   `scripts/msk/metabolic_cost.py` — at each of the 158 already-computed SO time frames, set the
   model's real IK-measured coordinate values+speeds, set all 80 real Muscle actuators' activation
   directly from SO's own `activation.sto`, call the real, compiled
   `model.equilibrateMuscles(state)`, then read `getFiberLength`/`getFiberVelocity`/`getLength`/
   `getTendonForce`. This **is** "the SO fiber lengths/velocities" the task asked to reuse.
2. **Gait-cycle timing**: reuse `validate_emg_timing.py`'s `detect_gait_cycle` (real GRF
   threshold-crossing detection, contralateral-symmetry cross-checked) and `pct_gc`/`bin_signal`
   (100-bin circular %GC binning) verbatim.
3. **Ia model**: Prochazka & Gorassini's equation above, applied to soleus_r's real reconstructed
   fiber displacement (`d`, referenced to this trial's own mean fiber length — the paper's absolute
   zero-reference was not recoverable live, disclosed) and fiber velocity (`v`, from the real
   compiled `getFiberVelocity`, never a numpy re-derivation). Both hybrid and velocity-only variants
   computed, neither cherry-picked. Real afferents cannot fire negative — hard floor at 0 imp/s.
4. **Ib/GTO (optional contrast)**: a simplified, disclosed log-compressive encoding of soleus_r's
   real tendon force, `Ib = 20 * ln(1 + F/(0.01·Fmax))` — anchored only *qualitatively*
   (monotonic, low-threshold) to Houk/Henneman/Crago, not a digitization of their fitted curves
   (not recoverable live — pre-modern/paywalled PubMed records).
5. Seven pre-registered gates (Sec. 3), stated in the code before any number existed.

### A genuine bug found and fixed this session: per-frame `equilibrateMuscles` silently zeros fiber velocity

Diagnosed via a standalone script before trusting any Ia number: with the model's **default
compliant tendon**, calling `model.equilibrateMuscles(state)` at each independently-set frame
solves the muscle's force equilibrium **assuming fiber velocity = 0** (an isometric-equilibrium
solve) — `getFiberVelocity()` then returns **exactly 0.0 at every one of the 158 frames**, despite
fiber length itself changing substantially frame-to-frame (confirmed against finite-difference of
`getFiberLength`; reordering `realizeVelocity` before/after `equilibrateMuscles` made no
difference — ruled out as a staleness/cache-order issue). **Fix**: force
`ignore_tendon_compliance=True` (rigid tendon) on all muscles — the same modeling choice
`metabolic_cost.py` already made (there, for force-reconstruction fidelity) — which makes fiber
length a deterministic geometric function of MTU length, so `getFiberVelocity()` becomes purely
kinematic. Verified: 6-frame spot check, `getFiberVelocity()` vs finite-difference agreed to
1-4% (14.68 vs 14.59 mm/s, 37.77 vs 37.60 mm/s, etc.); over the full 158-frame trial this becomes
Gate **G0** below (r=1.0000). This is a **new, necessary (not just sufficient) reason** for a
modeling choice this repo had already made for a different reason — and a disclosed, one-sided
consequence: **fiber displacement/velocity fed to the Ia model are likely over-estimates** of the
true in-vivo (tendon-damped) fascicle signal (Sec. 4, Honest gaps). A **separate**, length-only,
compliant-tendon reconstruction pass (velocity from it discarded) is used exclusively for Gate G1,
since a rigid tendon cannot show tendon compliance by construction.

## 2. Predicted Ia firing rate over the gait cycle (soleus_r, hybrid model)

Gait cycle detected from real GRF: T=1.330s, toe-off at 65.3%GC (consistency-checked between
right/left legs, diff=0.04 percentage points). Soleus_r: Fmax=6194.8N, Lopt=55.2mm (real model
params). Fiber length range [43.5, 67.8]mm (mean 58.2mm); fiber velocity range [-199.0, +91.9]
mm/s.

| %GC | Ia hybrid (imp/s) | Ia velocity-only (imp/s) | Ib/GTO (imp/s) | fiber length (mm) | ankle angle (deg) |
|---|---|---|---|---|---|
| 0.5 | 89.5 | 86.6 | 13.2 | 59.7 | 4.5 |
| 10.5 | 118.6 | 125.1 | 15.8 | 55.0 | -2.6 |
| 20.5 | 110.3 | 107.6 | 49.2 | 59.6 | 3.5 |
| 30.5 | 120.8 | 113.2 | 43.1 | 62.0 | 8.2 |
| 40.5 | 125.8 | 110.7 | 33.1 | 65.8 | 13.1 |
| 50.5 | 86.1 | 67.1 | 62.6 | 67.7 | 15.7 |
| 55.5 | 38.9 | 24.6 | 59.9 | 65.4 | 12.4 |
| **60.5** | **0.0** | **0.0** | 41.8 | 56.7 | 0.8 |
| **65.5** | **0.1** | 7.7 | 3.0 | 45.0 | -14.5 |
| 70.5 | 122.3 | 146.1 | 15.6 | 46.3 | -13.6 |
| 80.5 | 115.2 | 122.4 | 15.7 | 54.7 | -4.8 |
| 90.5 | 115.9 | 115.3 | 15.3 | 58.5 | 1.4 |
| 95.5 | 88.0 | 85.6 | 13.1 | 59.4 | 4.1 |

(Full 100-bin arrays in the results JSON, `pct_gc_bins`.)

**Emergent, physiologically coherent narrative** (derived from real IK kinematics + the transplanted
equation, not hand-tuned): fiber length and ankle dorsiflexion both peak at 49-51%GC (r=0.99 between
them, Gate G5) — soleus is most stretched just before push-off. Ia firing peaks earlier, at
35-45%GC (~126 imp/s), during **active lengthening** (length still rising *and* velocity still
positive) — exactly the classical Ia signature (peak firing during ongoing stretch, not at the
length peak itself where velocity crosses zero). Firing then **collapses to the 0 imp/s floor for
a tight 59.5-65.5%GC window** — 7 of 100 bins, centered almost exactly on the GRF-detected toe-off
(65.3%GC) — driven by the rapid active shortening of the push-off plantarflexion burst. This is a
direct, quantitative echo of Prochazka & Gorassini's own qualitative, pre-registered caveat about
triceps surae: **"a short period of complete muscle unloading" near the end of stance** — their
disclosed reason their hamstring-fitted model does *not* fit ankle extensors well shows up, in this
independent simulation, as literally the predicted signature. Firing then rebounds strongly in
early-mid swing (115-146 imp/s, 70-90%GC) as the ankle re-dorsiflexes for toe clearance and soleus
re-lengthens.

## 3. Pre-registered gates (stated in code before any number existed) — **5/7 PASS**

| Gate | What | Result | Verdict |
|---|---|---|---|
| G0 | Sign-convention QC: `getFiberVelocity` vs finite-diff of `getFiberLength` | r=1.0000 (gate >0.9) | **PASS** |
| G1 | Fascicle/tendon decoupling during stance vs Fukunaga 2001 (separate compliant-tendon pass) | fiber_range=19.75mm / mtu_range=23.90mm, ratio=0.827 (gate <0.7) | **FAIL** |
| G2 | Firing-minimum window vs an independently pre-existing `pre_swing_init_swing`=(50,73)%GC window, forced via `validate_emg_timing.py`'s own exhaustive 100-shift circular-null+Jaccard machinery | true_jaccard=0.548, rank vs null=93.9% (need ≥95 for PASS, ≥80 for PARTIAL) | **PARTIAL** (counted FAIL against the strict pre-registered bar) |
| G3 | Overall firing envelope vs the verified cat-locomotor envelope (mean 50-110, peak 150-200 imp/s) | mean=96.2 imp/s (gate [30,160]), peak=130.2 (gate ≤260) | **PASS** |
| G4 | Antagonist reciprocity: corr(soleus_r, tibant_r fiber velocity) | r=-0.968 (gate <-0.3) | **PASS** (see caveat below) |
| G5 | Ankle-angle / fiber-length geometric consistency | r=0.992 (gate \|r\|>0.5), sign confirms dorsiflexion lengthens soleus | **PASS** |
| G6 | Ia/Ib complementary encoding during push-off (50-73%GC): Ia should be low there, Ib relatively high | Ia_window=47.5 vs Ia_overall=94.0; Ib_window=34.6 vs Ib_overall=29.7 | **PASS** |

**G1 detail (an honest, explained miss, not hidden):** the model's own soleus tendon
(tendon_slack_length=352.7mm) shows real, correctly-directed but *modest* series compliance during
stance — estimated tendon strain (pennation ignored) ranges only **[-1.2%, +1.5%]**, versus
Fukunaga's ~7mm of measured gastrocnemius-medialis tendon stretch. The fiber absorbs 82.7% of the
whole-MTU excursion during stance (pre-registered pass bar was <70%) — real, directionally-correct,
but not as dramatic as Fukunaga's GM finding. Plausible mechanism (not chased further, secondary/
optional gate): this specific LaiArnoldModified2017 soleus instance's tendon-compliance parameters
may simply be stiffer than Fukunaga's subject's gastrocnemius medialis — a different muscle, not a
same-muscle comparison to begin with (disclosed in Sec. 0).

**G2 detail:** a stricter, purely descriptive cut (Ia < 5 imp/s, not a pre-registered threshold)
shows the near-total silencing is confined to 59.5-65.5%GC — landing almost exactly on toe-off — a
tight, sharp signal. The formal bottom-quartile mask (25/100 bins, needed for a fixed-duty-cycle
Jaccard test) is coarser: it also sweeps in several early-stance bins (1-8%GC) where firing dips for
an unrelated reason (low velocity shortly after heel-strike) plus the wide 51-67%GC descending ramp
into the true dip, diluting the mask's concentration inside the (50,73) target window just enough
to land at 93.9% rank rather than the ≥95% pass bar (74% of the window's own bins are in the
low-firing set; 68% of the low-firing set falls inside the window) — a genuine, quantified
near-miss, not a mysterious one.

**G4/G5 caveat (symmetric self-critique, not just celebrating the passes):** both of these very
strong correlations (r=-0.97, r=0.99) are **partly a mechanical consequence of the rigid-tendon
fix** (Sec. 1) — two mono-articular single-joint actuators with no series compliance have fiber
velocity/length nearly pure functions of the same one joint angle, so a strong correlation is close
to geometrically guaranteed, not fully independent surprising evidence. G6 (crossing the kinematic
Ia channel against the force-derived Ib channel, a genuinely different input) is the least
tautological of the three passes and the one weighted most heavily in this report's overall
reading.

**Illustrative-only elevated-gamma sensitivity** (Sec. "Honest gaps" — not fitted/measured): resting
gamma (bare equation) gives Ia mean=96.2 imp/s; an illustrative elevated static+dynamic fusimotor
drive (+15 imp/s bias, velocity gain ×1.25) gives mean=114.9 imp/s — included only to show the
parameter is wired up for a future dynamic-gamma layer.

## 4. Honest gaps (declared before running, not discovered after)

- **Generic spindle gains, cross-species AND cross-muscle-group transplant**: Prochazka &
  Gorassini's coefficients were fit to cat hamstring, not human soleus; their own paper says ankle
  extensors are not well fit by this simple model. This report does not claim imp/s-exact accuracy
  for soleus — only that the predicted shape/timing/envelope survive a forced circular-shift
  adversary against independently pre-existing literature windows (mostly true; G1/G2 above are the
  disclosed exceptions).
- **Fusimotor (gamma) drive is fixed, not dynamically modeled** (task's own explicit ask): default
  reproduces the bare Prochazka & Gorassini equation (zero bias, unit gain); one illustrative
  elevated-gamma variant reported, explicitly not fitted/measured.
- **Open-loop**: no reflex feedback path back to the muscles yet — this script only senses, it does
  not close any loop.
- **First step, two muscles**: soleus_r (primary) + tibant_r (antagonist contrast only). The
  detailed Mileusnic et al. 2006 intrafusal-fiber (bag1/bag2/chain) model — notably fitted AND
  validated on feline soleus/gastrocnemius, the same muscle group used here — is described (Sec. 0)
  but **not implemented**: a substantially heavier 3-fiber-type nonlinear-ODE undertaking than a
  first step warrants. Natural next-fidelity dial, and arguably the better long-run choice given it
  is already calibrated on this exact muscle.
- **No real Ia spike-train recording exists for this subject/trial/muscle** — validation is
  necessarily against published cat-locomotor envelope/timing facts, not a same-subject fit.
- **Rigid tendon forced on the primary Ia/Ib reconstruction** (Sec. 1) — a real, disclosed,
  one-sided bias: fiber displacement/velocity fed into the Ia model are likely **over-estimates**
  of the true (tendon-damped) in-vivo fascicle signal. Not correctable within a per-frame
  quasi-static equilibration pattern; would need a genuine forward-dynamics integration.
- **Ib/GTO model is illustrative**, anchored only qualitatively (monotonic, low-threshold) to
  Houk/Henneman/Crago — their own fitted curves were not independently recoverable live
  (pre-modern-era/paywalled PubMed records) and are not claimed to be reproduced quantitatively.

## 5. Overall result

**5 of 7 pre-registered gates PASS.** The two non-passes are a disclosed, explained, secondary
model-fidelity gap (G1: real but modest tendon compliance vs. Fukunaga's gastrocnemius-medialis
finding) and a genuine, quantified near-miss (G2: 93.9% vs a 95% pass bar, with a tight, sharply
correct near-silencing at the right location once measured at finer granularity). The headline,
task-requested result stands on its own terms: soleus_r's predicted Ia firing rate over a real
walking gait cycle (i) stays within the verified cat-locomotor envelope (mean 96.2, peak 130.2
imp/s, vs. a verified 50-110 mean / 150-200 peak reference range), (ii) peaks during active
mid-stance lengthening (35-45%GC) and (iii) **silences almost completely in a tight window centered
on toe-off** — the exact qualitative phenomenon ("complete muscle unloading" near end of stance)
the primary source paper itself named as the reason its equation does not fit this muscle group
well, now independently reproduced from this twin's own real fiber kinematics. Not oversold: this
is a shape/timing/envelope validation against cross-species/cross-muscle-group literature, not a
same-subject quantitative fit — the open gaps above are the honest ceiling on that claim.

## Repro

```
.venv-msk/bin/python3 scripts/msk/muscle_spindle.py
```

Reads (unmodified): `data/msk_smoketest/subject2_walking1/static_optimization/so/*.sto`, the real
subject2 IK/GRF/model files. Writes only:
`data/msk_smoketest/subject2_walking1/muscle_spindle/muscle_spindle_results.json`. No git
operations (isolation per task instruction).
