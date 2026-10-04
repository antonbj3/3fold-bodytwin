# MECHANISM SPINAL CPG LOCOMOTION — the half-center/Matsuoka rhythm generator + the fictive-locomotion falsifier (2026-07-22)

Script: `scripts/msk/spinal_cpg.py`. Results: `reports/probes/spinal_cpg_results.json`
(determinism confirmed: byte-identical md5 across 2 independent process runs this session, pure
deterministic RK4 ODE integration, no unseeded stochastic step). First node in the repo on this
topic — confirmed by a repo-wide grep of `data/MECHANISM_ANCHOR_GRAPH.json` (1038 nodes) for
CPG/CENTRAL PATTERN/HALF-CENTER/MESENCEPHALIC/MLR/FICTIVE LOCOMOTION/EPIDURAL: only an unrelated
emesis-CPG node (`NEU-EMESIS-CPG-CONVERGENCE`, a different brainstem circuit) and one passing
mention inside an unrelated auto-generated deglutition node turned up; no spinal-locomotor-CPG node
exists yet.

## 0. What this models, and the falsifier it must survive

The spinal cord generates rhythmic, alternating flexor/extensor (and left/right) locomotor output
**without** descending or rhythmic sensory input — "fictive locomotion" in deafferented, curarized,
or isolated-cord preparations (Graham Brown 1911; Grillner & Zangger 1979). The historical rival
hypothesis (Sherrington's "reflex chain": each phase of the step cycle is triggered as a peripheral
reflex response to sensory feedback from the *previous* phase's movement, with no intrinsic
oscillator at all) was falsified by exactly this kind of deafferentation/curarization experiment.
This document builds a half-center/Matsuoka oscillator, derives its oscillation condition
geometrically, and then **forces** a genuine reflex-chain adversary (in two structurally distinct
forms) through the same falsifying manipulation the historical experiments used, machine-checked
against pre-registered thresholds — not narrated.

**Model.** Two mutually-inhibiting neurons (flexor F, extensor E) per half-center, each with an
adaptation/fatigue variable v (Matsuoka 1985, PMID 2996634, the continuous-variable form of Graham
Brown's 1911 half-center concept):
```
y_i = relu(x_i)
tau  * dx_i/dt = -x_i - a*y_j - b*v_i + u + feedback_i
tau2 * dv_i/dt = -v_i + y_i
```
`u` = tonic ("MLR"/descending) drive — constant, non-rhythmic. `feedback_i = k_fb * g(m_j)` is an
optional simulated proprioceptive/mechanical loop (`m` low-pass-tracks the antagonist's own output,
i.e. "the limb has moved"); `k_fb=0` is the computational analogue of deafferentation/curarization
— sensory transmission or muscle contraction is prevented, so no movement-generated signal can
reach the circuit, while the intrinsic circuitry itself is untouched (exactly the two real
manipulations Grillner & Zangger 1979 used). Two feedback shapes are tested: `g=identity` ("linear",
proportional stretch-like feedback) and `g=`saturating sigmoid ("threshold", a relay/chain-reflex
switch) — two structurally different adversary forms, not one cherry-picked shape.

## 1. Geometric derivation (Z2 symmetry of the flexor↔extensor swap)

At the symmetric tonic-drive fixed point (`x*=u/(1+a+b)`, both neurons co-active), perturbations
decompose exactly into the **symmetric** mode (`s=xF+xE`, `w=vF+vE` — "how co-active") and
**antisymmetric** mode (`d=xF-xE`, `z=vF-vE` — "which one dominates") — the ±1 eigenspaces of the
flexor↔extensor swap, a Z2 symmetry of the circuit. By hand:
```
symmetric:      trace = -(1+a)/tau - 1/tau2         det = (1+a+b)/(tau*tau2)
antisymmetric:  trace = (a-1)/tau  - 1/tau2          det = (b-a+1)/(tau*tau2)
```
The symmetric block has trace<0, det>0 **for every a,b>0** — the "how much total activity" mode is
always stable, regardless of parameters. The antisymmetric block is the whole story: if `a<=1`,
`trace<=0` and `det>=0` **simultaneously and unconditionally** (both terms of each are already
same-signed) — the winner-take-all direction is stable too, so **no rhythm is possible at all**,
for any `b`, any `u`. `a>1` is therefore a hard, geometrically-derived **necessary condition**.

**Cross-checked, not asserted**: a finite-difference numerical Jacobian at the fixed point was
built directly from the ODE code and block-diagonalized in the (s,w,d,z) basis, for 15 (a,b)
configs (a∈{0.5,1.0,1.5,2.0,3.0}, b∈{0.5,2.5,6.0}) — **all 15/15 trace/det values match the
hand-derived closed form to the finite-difference step size** (`rtol=1e-3`). Sample (a=2.0,
b=2.5): analytic `tr_anti=19.6, det_anti=12.0`; numeric `19.600000, 12.000000`.

**Numerically forced across a diverse (b,u) grid** (necessary_condition_sweep, 120 configs: a ∈
{0.5,0.8,1.0,1.5,2.0,3.0} × b ∈ {0.5,1,2,4,8} × u ∈ {1,3,6,10}, classified by a pre-registered
oscillation detector — ≥4 peaks in the settled window, inter-peak-interval CV<0.20): **a≤1 gives
0/60 oscillating, exactly as derived** (a=0.5: 0/20, a=0.8: 0/20, a=1.0: 0/20 — the boundary case
itself is non-oscillating too, consistent with a *strict* inequality). For a>1, oscillation also
needs `b` large enough (adaptation strength) — the sweep shows a=1.5 oscillates for b≥1.0 (16/20),
a=2.0 for b≥2.0 (12/20), a=3.0 for b≥4.0 (8/20). The exact analytic boundary `b=a-1` (where
`det_anti=0`) is confirmed **non-oscillating** in all three cases (a marginal, degenerate point —
correctly outside the linear-instability region); the true sufficient threshold lies somewhere
between `b=(a-1)` and `b=2(a-1)` (the next grid point, confirmed oscillating in all three cases) —
resolved to that bracket, not claimed finer than the grid supports. This matches Matsuoka's own
abstract conclusion exactly ("the adaptation of the neurons plays a very important role for the
appearance of the oscillations") — without adaptation, `a>1` alone only buys **bistability**
(a frozen winner-take-all state), not rhythm.

## 2. The falsifier: fictive locomotion vs. the forced reflex-chain adversary

**Pre-registered** (before any run): the half-center (`b>0`, adaptation-driven) must keep
oscillating at `k_fb=0`; a reflex-chain adversary (`b=0`, rhythm from feedback alone), forced to
its strongest closed-loop form across a parameter grid, must collapse to 0% oscillating in the
*same* open-loop condition — in **both** feedback variants, not one.

| test | closed-loop (feedback on) | same configs, open-loop (feedback cut) |
|---|---:|---:|
| Adversary, linear feedback | 8/10 configs oscillate | **0/8** still oscillate |
| Adversary, threshold/relay feedback | 28/30 configs oscillate | **0/28** still oscillate |
| True half-center (adaptation), 9 (a,b) pairs × 4 k_fb values | 36/36 configs oscillate | **36/36** still oscillate (100% survival) |

The threshold-variant adversary was **forced twice**, not accepted at its first (weak) attempt: an
initial k_fb grid (0.5–5) gave 0/20 closed-loop — diagnosed (not shrugged off) as an under-forced
adversary, not evidence: the sigmoid's output is *capped* at `k_fb`, but must overcome `a*u≈10` of
inhibition, whereas the linear variant's feedback is *unbounded* (scales with the antagonist's own
steady output, ≈u) — an apples-to-oranges comparison. Re-forced with `k_fb∈{6,10,15,20,30}` ×
`theta∈{1,2,3}` × `width∈{0.3,0.6}`: 28/30 configs now oscillate cleanly (xcorr(yF,yE)≈-0.85 to
-0.92, clean antiphase), i.e. the adversary was given a genuinely strong, fair shot — and *still*
collapses 0/28 when deafferented.

**Headline single-preparation demonstration** (same system, one continuous 120 s run, feedback cut
at t=60s — mirroring the actual before/after design of a real deafferentation experiment, not two
separate runs): half-center (a=2.0, b=2.5, k_fb=1.5→0) oscillates before (1.10 Hz) **and after**
(0.27 Hz — slower once the extra feedback-driven entrainment is removed, see §3) the cut. The
reflex-chain adversary (a=2.0, b=0, k_fb=2.0→0) oscillates before (1.17 Hz) and goes **silent**
after (0 peaks, fixed point) — the same qualitative before/after signature Grillner & Zangger 1979
report in real cat spinal cord.

## 3. MLR drive → frequency: an honest negative, forced via OODA, then fixed

**Pre-registered claim**: frequency increases monotonically with tonic drive `u` (Spearman
rho≥0.9, p<0.01), holding a,b,tau,tau2 fixed. **Result: rho=0.41, p=0.17 — FAILS the bar.**

Rather than report that as a bare negative, it was diagnosed. On-attractor peak amplitude and
period (measured only in the settled window, immune to the differing-transient confound of
comparing raw trajectories at a fixed time index when the initial condition isn't rescaled with
`u`) show, across a 19× range of u (1.0→19.0): **peak amplitude / u is constant to 5 significant
figures (CV=1.4e-5), and period is constant to 4 significant figures (CV=2.2e-4)**. This is an
**exact geometric symmetry** of the model, not a numerical artifact: `relu` is positive-homogeneous
of degree 1, and the `+u` drive term scales along with it, so `F(c·X, c·u) = c·F(X, u)` for any
c>0 — the whole vector field is invariant under joint rescaling of state and drive. Consequence:
**pure tonic-drive amplitude sets burst amplitude, not frequency**, in this idealized model. This
is a real, derivable fact about the linear-adaptation Matsuoka oscillator (and is why many
CPG/robotics implementations control speed via a *separate* parameter, not the tonic input alone).

**Forced fix** (physiologically motivated, disclosed as a minimal, non-fit extension, not a
verified biological measurement): couple adaptation kinetics to drive, `tau2_eff =
tau2/(1+k_speed*u)` — real depolarizing drive increases firing rate, which increases the rate of
Ca²⁺ influx driving Ca-activated-K⁺ adaptation currents, so faster drive plausibly means faster
adaptation kinetics, not just bigger amplitude. Re-swept with `k_speed=0.15`: **rho=1.0, p=0.0**,
frequency range 0.31–0.84 Hz (18–50 cycles/min) across the same u-range — landing inside the
commonly-cited order-of-magnitude range for mammalian fictive/spinal locomotor cycle frequency
(~0.3–2 Hz; Grillner 1975, PMID 1144530, title/DOI bibliographically confirmed live this session,
no abstract available pre-dating PubMed's abstracting era — this specific numeric range is
disclosed as **not** independently re-extracted from a live-fetched primary-text number this
session, a textbook-tier plausibility check, not a decisive anchor).

## 4. Interlimb coordination: coupling SIGN sets the gait mode (measured, not just cited)

Two half-center pairs (Left, Right), coupled by a signed commissural flexor-flexor term `c`
(`c>0`=inhibitory cross-coupling, `c<0`=excitatory) — mirroring the real dissociation in Talpalar
et al. 2013 (PMID 23812590): ablating inhibitory V0 commissural neurons removes left-right
alternation at low locomotor frequency; ablating excitatory V0 neurons removes left-right synchrony
("hopping"). Swept `c` and **verified the phase-locked mode is IC-independent (a genuine attractor,
not an artifact of the starting condition)** by testing each regime from an oppositely-biased
initial condition:

| c (commissural coupling) | outcome from an asymmetric IC | outcome from the OPPOSITE-biased IC |
|---:|---|---|
| +0.5, +1.5, +3.0 (inhibitory) | anti-phase/alternating (xcorr≈-0.77 to -0.91) | **still** anti-phase, from an in-phase-biased start (xcorr≈-0.79 to -0.82) |
| -0.5, -1.0 (excitatory) | in-phase/synchrony (xcorr=+1.0) | **still** in-phase, from an anti-phase-biased start (xcorr≈+1.0) |
| 0 (no coupling) | whatever the IC set (marginal) | **flips** with the IC (0.999 in-phase from an in-phase start, -0.91 from an asymmetric start) — correctly NOT a third enforced mode, the mathematically-expected neutral case for two truly uncoupled identical oscillators |

A genuine, bidirectional, IC-independent bifurcation in coupling sign — not a coincidence of one
initial condition. **Disclosed limitation, not hidden**: excitatory coupling stronger than
`c≈-1.4` makes this unsaturated-linear model numerically run away (confirmed finite for c=-1.3,
-1.4; non-finite for c≤-1.5) — a real structural gap of the simplified linear-outside-ReLU
excitatory pathway (real neurons have firing-rate ceilings preventing this), not a third gait mode.
This is a qualitative, coupling-*sign* mechanism check, explicitly **not** a quantitative
reproduction of Talpalar 2013's own frequency-dependence thresholds, and not a claim to model the
specific walk→trot→gallop kinematic sequence Shik et al. 1966 (PMID 6000625, no abstract available
— Russian-language, pre-abstracting era, confirmed via its own bibliographic record and via the
secondary review Ryczko & Dubuc 2013, PMID 23360276) describe under graded MLR stimulation
strength — that specific quantitative threshold structure is an honest, disclosed gap.

## 5. SCI dysfunction ↔ epidural-stimulation restoration (the contrast pole)

Same half-center (a=2.0, b=2.5): **intact** (u=5 throughout) oscillates (0.269 Hz). **"Complete
SCI"** (descending drive severed, u=0 from t=0): **silent** — 0 peaks, max output over the settled
window = exactly 0.0 (the trivial fixed point; no motor output at all, matching complete loss of
supraspinal drive). **SCI, then exogenous tonic drive substituted at t=60s** (modeling sustained
epidural electrical stimulation of the isolated cord below the lesion): silent before (0/0 peaks),
**oscillating after** (0.269 Hz, matching the intact baseline exactly) — rhythm resumes from the
SAME intrinsic mechanism, once a tonic-drive-equivalent is restored by any source, descending or
exogenous.

This mirrors, and is anchored by, the real literature — not a tautology gate, since the real
finding is an independently measured fact about real spinal cords, and my model's `u=0→u=5` step is
a schematic stand-in for it, not a fit to it: **Dimitrijevic, Gerasimenko & Pinter 1998** (PMID
9928325, full abstract verified live) directly state the mechanism this section models: "externally
controlled sustained electrical stimulation of the spinal cord can replace the tonic drive generated
by the brain" — measured in human paraplegics via a real 25–60 Hz, 5–9 V epidural train inducing
"rhythmic, alternating stance and swing phases." **Harkema et al. 2011** (PMID 21601270, first
human case, n=1, 16-electrode L1-S1 array, 170 prior training sessions/26 months, 29 stimulation
experiments) achieved full weight-bearing standing (4.25 min) and locomotor-like patterns.
**Angeli et al. 2018** (PMID 30247091, n=4 chronic motor-complete SCI): 2/4 achieved independent
over-ground walking after 278 sessions/85 weeks and 81 sessions/15 weeks of combined stimulation +
training; all 4 achieved standing/trunk stability. **Courtine et al. 2009** (PMID 19767747, rat,
complete transection): epidural stimulation + serotonergic agonists restored "full weight-bearing
treadmill locomotion... almost indistinguishable from voluntary stepping" within a week — explicitly
attributed to "a combination of central pattern-generating capability and the ability of these
spinal circuits to use sensory afferent input" (a hybrid mechanism, honestly more than my pure-tonic
-substitution toy model captures — real epidural stimulation works substantially through recruiting
proprioceptive afferents, not by literally injecting a tonic current into interneurons; disclosed
gap, not hidden, see §7). **van den Brand et al. 2012** (PMID 22654062, rat): showed the deeper
recovery is training- and cortex-dependent (automated treadmill-only training did *not* promote
translesional plasticity) — my model has no training/plasticity dimension at all, an honest scope
limit, not a claim about how real recovery is earned.

## 6. Citations — 18 sources, every PMID/DOI verified LIVE this session (NCBI eutils / Crossref)

| # | citation | id | verified as |
|---|---|---|---|
| 1 | Brown TG (1911) "The intrinsic factors in the act of progression in the mammal." *Proc R Soc Lond B* 84(572):308-19 | DOI 10.1098/rspb.1911.0077 (no PMID, pre-PubMed) | Crossref bibliographic match live |
| 2 | Grillner S, Zangger P (1979) "On the central generation of locomotion in the low spinal cat." *Exp Brain Res* 34(2):241-61 | PMID 421750, DOI 10.1007/BF00235671 | Full abstract fetched live — THE decisive fictive-locomotion anchor |
| 3 | Grillner S (1985) "Neurobiological bases of rhythmic motor acts in vertebrates." *Science* 228(4696):143-9 | PMID 3975635, DOI 10.1126/science.3975635 | Full abstract fetched live |
| 4 | Grillner S, Wallén P (1985) "Central pattern generators for locomotion, with special reference to vertebrates." *Annu Rev Neurosci* 8:233-61 | PMID 2984978, DOI 10.1146/annurev.ne.08.030185.001313 | Title/journal/DOI verified live; no abstract text in PubMed for this record |
| 5 | Grillner S (1975) "Locomotion in vertebrates: central mechanisms and reflex interaction." *Physiol Rev* 55(2):247-304 | PMID 1144530, DOI 10.1152/physrev.1975.55.2.247 | Title/journal/DOI verified live; no abstract (pre-abstracting era) |
| 6 | Shik ML, Severin FV, Orlovskiĭ GN (1966) "[Control of walking and running by means of electric stimulation of the midbrain]." *Biofizika* 11(4):659-66 | PMID 6000625 | Title/journal verified live; Russian, no abstract, pre-abstracting era |
| 7 | Ryczko D, Dubuc R (2013) "The multifunctional mesencephalic locomotor region." *Curr Pharm Des* 19(24):4448-70 | PMID 23360276, DOI 10.2174/1381612811319240011 | Full abstract fetched live — secondary confirmation of Shik 1966's discovery + cross-species/clinical generality |
| 8 | Matsuoka K (1985) "Sustained oscillations generated by mutually inhibiting neurons with adaptation." *Biol Cybern* 52(6):367-76 | PMID 2996634, DOI 10.1007/BF00449593 | Full abstract fetched live — directly validates this doc's derivation approach ("no stable stationary state for constant input ⟹ sustained oscillation... adaptation plays a very important role") |
| 9 | Kiehn O (2006) "Locomotor circuits in the mammalian spinal cord." *Annu Rev Neurosci* 29:279-306 | PMID 16776587, DOI 10.1146/annurev.neuro.29.051605.112910 | Full abstract fetched live |
| 10 | Kiehn O (2016) "Decoding the organization of spinal circuits that control locomotion." *Nat Rev Neurosci* 17(4):224-38 | PMID 26935168, DOI 10.1038/nrn.2016.9 | Full abstract fetched live |
| 11 | Talpalar AE et al. (2011) "Identification of minimal neuronal networks involved in flexor-extensor alternation in the mammalian spinal cord." *Neuron* 71(6):1071-84 | PMID 21943604, DOI 10.1016/j.neuron.2011.07.011 | Full abstract fetched live — real inhibitory-interneuron identity (Ia INs) behind this model's abstract "mutual inhibition" |
| 12 | Talpalar AE, Bouvier J, Borgius L, Fortin G, Pierani A, Kiehn O (2013) "Dual-mode operation of neuronal networks involved in left-right alternation." *Nature* 500(7460):85-8 | PMID 23812590, DOI 10.1038/nature12286 | Full abstract fetched live — direct anchor for §4's coupled-oscillator result |
| 13 | Courtine G et al. (2009) "Transformation of nonfunctional spinal circuits into functional states after the loss of brain input." *Nat Neurosci* 12(10):1333-42 | PMID 19767747, DOI 10.1038/nn.2401 | Full abstract fetched live |
| 14 | van den Brand R et al. (2012) "Restoring voluntary control of locomotion after paralyzing spinal cord injury." *Science* 336(6085):1182-5 | PMID 22654062, DOI 10.1126/science.1217416 | Full abstract fetched live |
| 15 | Angeli CA et al. (2018) "Recovery of Over-Ground Walking after Chronic Motor Complete Spinal Cord Injury." *N Engl J Med* 379(13):1244-50 | PMID 30247091, DOI 10.1056/NEJMoa1803588 | Full abstract fetched live |
| 16 | Harkema S et al. (2011) "Effect of epidural stimulation of the lumbosacral spinal cord on voluntary movement, standing, and assisted stepping after motor complete paraplegia: a case study." *Lancet* 377(9781):1938-47 | PMID 21601270, DOI 10.1016/S0140-6736(11)60547-3 | Full abstract fetched live |
| 17 | Dimitrijevic MR, Gerasimenko Y, Pinter MM (1998) "Evidence for a spinal central pattern generator in humans." *Ann N Y Acad Sci* 860:360-76 | PMID 9928325, DOI 10.1111/j.1749-6632.1998.tb09062.x | Full abstract fetched live — explicit tonic-drive-substitution claim + 25-60Hz/5-9V numbers |
| 18 | Ijspeert AJ (2008) "Central pattern generators for locomotion control in animals and robots: a review." *Neural Netw* 21(4):642-53 | PMID 18555958, DOI 10.1016/j.neunet.2008.03.014 | Full abstract fetched live |

**Recall-drift instances caught this session** (disclosed, not smoothed over): (1) the Talpalar
"dual-mode" left-right coordination paper was initially mis-recalled as *Neuron* 2013 — live search
corrected it to *Nature* 2013 (PMID 23812590), and separately surfaced a genuinely different,
also-relevant 2011 *Neuron* paper by the same lab (PMID 21943604) not originally sought, which
became an additional real anchor. (2) An initial esearch for Shik et al. 1966's own primary paper
returned, as its top hit, a *2013 secondary review* (Ryczko & Dubuc, PMID 23360276) rather than the
primary 1966 record — caught by checking the returned publication year against the expected 1966
date, then re-searching to find the actual primary record (PMID 6000625); the secondary review was
kept as a separate, disclosed, corroborating citation rather than discarded.

## 7. Scope, honest gaps (first-step, not final)

- **Schematic time constants**, not fit to measured spinal interneuron biophysics: `tau=50ms`
  (membrane/synaptic), `tau2=2.5s` (adaptation) are illustrative choices giving a legible period,
  disclosed as such (same discipline as this repo's other textbook-tier-constant disclosures).
  The resulting baseline frequency (0.27–1.17 Hz across configs tested) sits at the right
  order of magnitude for commonly-cited mammalian fictive-locomotion cycle rates, but this specific
  numeric range was **not** independently re-extracted from a live-fetched primary-source number
  this session — flagged, not silently assumed.
- **The frequency-vs-drive "fix" (§3) is a disclosed, unverified model extension**, not a measured
  biological mechanism: real MLR pathways almost certainly recruit both amplitude- and
  frequency-affecting mechanisms in parallel; this document demonstrates *sufficiency*
  (drive-coupled adaptation kinetics CAN restore monotonic frequency scaling) not which mechanism
  biology actually uses.
- **The SCI/epidural-stimulation model is a pure tonic-drive substitution**, not a model of the real
  electrode/dorsal-root-afferent recruitment biophysics; Courtine et al. 2009's own abstract
  explicitly attributes real recovery to a *hybrid* of intrinsic pattern-generation AND sensory
  afferent use, and van den Brand et al. 2012 show a training/cortex-dependence this model has no
  mechanism for at all (no plasticity, no learning, no cortex).
- **§4's coupled-oscillator gait-mode result is a coupling-sign qualitative mechanism check**,
  explicitly not a quantitative reproduction of Talpalar 2013's own frequency-dependence thresholds
  or of Shik et al. 1966's specific walk→trot→gallop stimulation-intensity numbers (the latter's
  primary text was not independently accessible this session — Russian, pre-abstracting era).
- **Two neurons per side is the minimum viable half-center**, not a segmental, multi-joint, or
  multi-limb (beyond the 2-oscillator L/R toy in §4) network; no real flexor/extensor motoneuron
  pools, no muscle/tendon mechanics, no ground-contact/load feedback of the kind that in real
  locomotion also entrains phase (Duysens & Pearson-type findings) — this document models the
  *intrinsic* oscillator claim specifically, which is exactly what the fictive-locomotion falsifier
  isolates, by design.
- **The excitatory-coupling runaway boundary (§4, c≤-1.5)** is a genuine structural limitation of
  the unsaturated linear-outside-ReLU excitatory pathway, not resolved by adding a saturating
  ceiling in this pass — disclosed, not patched.
- Every simulation is **deterministic** (fixed initial conditions, no RNG) — the "diverse
  instance-space" claims are diversity over swept *parameters*, not over stochastic seeds; this is
  appropriate for a claim about the geometry of one deterministic ODE, not for a claim requiring
  noise-robustness (not attempted here).

## 8. couples_to

- `docs/MECHANISM_CORTICOSPINAL_MOTOR_COMMAND.md` (+ its evidence.json) — the complementary
  descending pathway: that document models the corticospinal/M1 tract's population-vector movement
  *direction* code and conduction-delay physics (voluntary, precise commands); this document models
  the separate, phylogenetically older MLR→reticulospinal→spinal-CPG rhythm-*generation* substrate
  the corticospinal tract does not itself generate (a real anatomical and functional distinction,
  not an overlap — confirmed by reading that doc's own scope, not re-verified biologically this
  session).
- `docs/MECHANISM_STRETCH_REFLEX.md`, `docs/MECHANISM_GAIT_KINEMATICS_FIDELITY.md`,
  `docs/MECHANISM_SPINE_GAIT_VBR.md`, `docs/MECHANISM_BIARTICULAR_MUSCLE.md` — the twin's gait/
  athlete-movement axis this CPG's output would ultimately have to drive (motoneuron pools, muscle
  activation, kinematics); none of these model the rhythm-generation layer itself, a genuine,
  non-overlapping gap this document fills. Read-only, not edited this session.
- `data/MECHANISM_ANCHOR_GRAPH.json` — no pre-existing node found for this topic (§ header); this is
  the first spinal-locomotor-CPG node.
- The function↔dysfunction axis: §5's SCI/epidural-stimulation contrast pole is this document's own
  built-in dysfunction leg, not a pointer to a separate doc (none exists yet for spinal cord injury
  specifically).
- The "compute = mind" substrate framing: this is offered as a concrete, textbook example of a
  neural oscillator running on the inherited compute-boundary framing (COORDINATOR.md's own mission
  language) — a pointer, not a re-derivation of that framing.

## 9. Repro

```
cd ~/projects/bodytwin
python3 scripts/msk/spinal_cpg.py
```
Pure numpy (fixed-step RK4) + scipy (`find_peaks`, `spearmanr`); no external data, no GPU, no
OpenSim. Runtime ≈85s on one core. Writes `reports/probes/spinal_cpg_results.json` (verified
byte-identical across 2 independent runs this session, md5 `7efcd788e717e3a703902baeed6f9916`).
No git operations.
