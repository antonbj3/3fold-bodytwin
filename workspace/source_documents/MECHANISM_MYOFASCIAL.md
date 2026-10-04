# MECHANISM MYOFASCIAL FORCE TRANSMISSION — fidelity-audit-#3 gap, quantified + prototyped (2026-07-21)

Executes the fidelity-audit follow-up on `docs/MECHANISM_FIDELITY_ARCHITECTURE_AUDIT.md` §3 row
"Fascia (myofascial transmission)": **"~zero. Repo-wide search finds no lateral-force-transmission
model."** OpenSim (and every fork in this repo) treats each muscle as an INDEPENDENT line
actuator — force flows only along its own tendon path. Real muscles also transmit force
LATERALLY, through epimysium/perimysial connective tissue, to adjacent synergists
("epimuscular myofascial force transmission", Huijing/Baan/Maas, VU Amsterdam). This document
(1) quantifies that effect from the primary literature (every number **live-verified this
session**, not recalled), (2) ranks where the independent-actuator assumption is most exposed in
THIS repo's own model, and (3) builds + measures a first lumped gastrocnemius↔soleus coupling
against the existing, already-certified independent-actuator Static-Optimization solution on the
real subject2 `walking1` trial. Script: `scripts/msk/myofascial_transmission.py`. Run:
`source_repository/.venv-msk/bin/python3 scripts/msk/myofascial_transmission.py`.
Output: `data/msk_smoketest/myofascial_transmission/myofascial_transmission_results.json`.

**Verification method** (stated up front, per this repo's own discipline): every citation below
was fetched LIVE this session via Europe PMC's REST API (machine JSON against PubMed's own
indexed record — PubMed's own HTML page returned a reCAPTCHA wall to automated fetch, worked
around via the EBI mirror). The sourcing pass used a subagent; the **3 numbers that directly
calibrate this script's coupling constants** were independently re-fetched a SECOND time by the
orchestrating agent (not merely trusted) — this second pass caught a real, if minor, error in
the subagent's report: PMID 12655613's author order was mis-stated (Maas listed first; the
verified record gives Huijing as first author). Corrected here.

---

## 1. Literature — the quantitative lateral-transmission fraction

### 1a. Rat anterior crural compartment (EDL/TA/EHL) — Huijing/Baan/Maas lineage, VU Amsterdam

The classic paradigm: surgically alter one synergist's relative position/length, measure a
target muscle's force at **both** its proximal and distal tendons simultaneously (any
difference between the two ends is force that entered/left along the muscle's own length via
the surrounding connective tissue, not through either tendon).

| Study | PMID / DOI | Quantitative finding |
|---|---|---|
| Huijing & Baan 2001, *Arch Physiol Biochem* 109(2):97-109 | 11780782 / 10.1076/apab.109.2.97.4269 | Blunt EDL–TA dissection: force ↓~10% at all lengths (no shift in optimum/slack length). Full lateral fasciotomy: length range ↑~47%. |
| Huijing & Baan 2001, *Acta Physiol Scand* 173(3):297-311 | 11736692 / 10.1046/j.1365-201x.2001.00911.x | F_proximal − F_distal: 0 to **+22.7%** of F_proximal at long muscle lengths, 0 to **−24.5%** at short lengths. |
| Maas, Baan & Huijing 2001, *J Biomech* 34(7):927-40 | 11410176 / 10.1016/s0021-9290(01)00055-0 | Establishes the proximal≠distal EDL force mechanism via TA/EHL relative-length manipulation (mechanism paper, no single headline %). |
| **Huijing, Maas & Baan 2003**, *J Morphology* 256(3):306-21 — *author order corrected this session* | 12655613 / 10.1002/jmor.10097 | Intact: proximo-distal EDL force difference = **"0–14% of proximal force"**, constant ~14% at most lengths. Post-fasciotomy: falls to 0–5%. Post-full-isolation: "no longer significantly different from zero" (active force) — **passive**-force difference persists even after isolation. |
| Maas, Meijer & Huijing 2005, *Cells Tissues Organs* 181(1):38-50 | 16439817 / 10.1159/000089967 | TA+EHL lengthened +12mm: proximal EDL force **+9.5%** (+0.14N), distal EDL force **−11.8%** (−0.21N). Blunt dissection cut the distal-curve amplitude 39%. |
| Huijing 2009, *J Biomech* 42(1):9-21 (ISB Muybridge Award Lecture) | 19041975 / 10.1016/j.jbiomech.2008.09.027 | Historical/conceptual review; explicitly hedges: "even if the quantitative effects in terms of force would prove small […]" — the field's own leading author does not claim the effect is always large. |
| Maas & Sandercock 2010, *J Biomed Biotechnol* 2010:575672 (review) | 20396618 / 10.1155/2010/575672 | Synthesis: supraphysiological dissection/preps → "substantial force"; physiologically-relevant **intact** conditions → "role of this myofascial pathway is **small**." |
| Huijing 1999, *J Biomech* 32(4):329-45 | 10213024 / 10.1016/s0021-9290(98)00186-9 | Conceptual: "muscle as a collagen fiber reinforced composite" — the shear-lag/composite load-transfer picture this prototype's coupling law is geometrically modeled on. No quantitative fraction given. |

### 1b. Triceps-surae-adjacent literature — the direct analog to this build's target pair

The EDL/TA papers above are a **different compartment** (anterior crural). These two are the
anatomically closest analog to gastrocnemius↔soleus, and are the two numbers this build's
coupling constants are directly anchored to — **both independently re-verified this session**:

- **Rijkelijkhuizen, Baan, de Haan, de Ruiter & Huijing 2005**, *J Exp Biol* 208(Pt24):4715-25.
  **PMID 15601884 / DOI 10.1242/jeb.01360.** Rat medial gastrocnemius (GM) + plantaris,
  progressive dissection: with GM's surrounding connective tissue fully **intact**, **up to
  40.5±5.9%** (mean±SEM) of plantaris force is transmitted onto the **calcaneus** despite
  plantaris having no tendon attachment there — entirely via GM's fascial envelope. Once GM is
  fully dissected free (isolated except neurovascular pedicle) **and** returned to its reference
  relative position, transmission is "no [longer] relevant" — reappears only if relative
  position is again shifted. → calibrates this build's **SENSITIVITY** regime, `f_max = 0.405`.
- **Maas & Sandercock 2008**, *J Appl Physiol* 104(6):1557-67. "Are skeletal muscles independent
  actuators? Force transmission from soleus muscle in the cat." **PMID 18339889 / DOI
  10.1152/japplphysiol.01208.2007.** A **different, independent lab** (Northwestern, not VU
  Amsterdam — removes single-lab-lineage confound for exactly this pair). Method: vary **knee**
  angle (70–140°) at fixed ankle angle — changes gastrocnemius length/position relative to
  soleus (soleus does not cross the knee) while soleus's own MTU length is held constant by the
  ankle — the *same* "one synergist's relative position changes" paradigm, applied to the exact
  pair this build models. Result: soleus ankle moment was **unchanged** by knee-angle-driven
  gastrocnemius repositioning when tissue was intact. After soleus tenotomy, moment fell
  **55±16%** but did not vanish (a myofascial path exists) — **confounded** by the isolated
  soleus shortening 16.0±0.6mm under contraction vs only 1.0±0.1mm intact (a large secondary
  position shift, not a clean single-variable test). Repositioned to the intact-matched
  position, moment "approached zero." Their own conclusion, quoted directly: **"the intact cat
  soleus muscle appears to act mechanically as an independent actuator"** under physiological
  relative position. → calibrates this build's **DEFAULT** regime slack zone (`Δ_slack = 0.20`).

### 1c. Human in-vivo evidence — checked explicitly, honest gap

| Study | PMID / DOI | Finding |
|---|---|---|
| Bojsen-Møller et al. 2010, *J Appl Physiol* | 20884838 / 10.1152/japplphysiol.01381.2009 | Human, n=7, ultrasonographic **displacement** (not force) used as a loading proxy between MG/soleus/FHL. "Force **may** be transmitted"; "only limited evidence" for a triceps-surae→FHL transfer. **No % given.** |
| Yucesoy et al. 2018, *J Mech Behav Biomed Mater* | 28892760 / 10.1016/j.jmbbm.2017.08.040 | Human, but **pathological**/intraoperative (12 limbs, 7 cerebral-palsy patients) — not representative of healthy modeling assumptions. One condition: spastic semitendinosus force **+33.3%**. |
| Finni, de Brito Fontana & Maas 2023, *J Biomech* 154:111575 (review) | 37120913 / 10.1016/j.jbiomech.2023.111575 | As of 2023: "most direct evidence is from animal experiments; studies on humans also suggest functional implications" — still **no quantified human fraction**, 22+ years after the original rat work. |

### 1d. Synthesis (pre-registered characterization, anchored to the numbers above)

Under **intact/physiological relative position**, the effect is **small-to-moderate** (0–14%
typical, rat EDL/TA; **~0%**, cat soleus at true intact position — the most directly relevant
result for this build's pair) — reaching 22–40% only under either (a) length **extremes** (rat
EDL, ±22.7/−24.5% at long/short length) or (b) an intact-fascia-but-fully-exposed dissection
stage (rat GM→plantaris, 40.5%). **No direct human in-vivo quantitative fraction exists** in any
source checked. For gastrocnemius↔soleus specifically, the single most directly relevant,
independent-lab study is a genuine small-effect finding at physiological position — reported
here as a **valid result**, not downgraded or hidden, per this task's own instruction that a
small in-vivo effect is a legitimate finding.

---

## 2. Where the independent-actuator assumption is most exposed in THIS repo's model

Ranked by (real force magnitude × anatomical adjacency × literature availability), using facts
already measured elsewhere in this repo, or live-measured in §3 below:

1. **Triceps surae (gasmed_r/gaslat_r/soleus_r) — this build's target.** Highest-force
   posterior-compartment group in gait (peak Static-Optimization force 1140–1540 N per head this
   trial, §3); gastrocnemius sits directly superficial to soleus across a well-defined fascial
   plane (the anatomical site Bojsen-Møller 2010 imaged); **near-identical ankle moment arms**
   (live-measured §3: −49.9 / −50.4 / −47.5 mm, within ~5% of each other) — mechanically these
   three are already the closest thing to "one shared output" of any muscle group in this
   model; and, uniquely among this model's muscle groups, **direct literature exists for this
   exact synergist neighborhood** (§1b). Note: this repo already flags a **related but distinct**
   simplification here — `scripts/msk/tendon_elastic_energy.py` documents gasmed_r/gaslat_r/
   soleus_r as **3 independent tendons** rather than 1 confluent Achilles tendon. That is a
   distal **tendon-convergence topology** question; this build addresses a mechanistically
   different path — **mid-belly lateral shear** — which has zero representation of any kind in
   the existing model, untouched by any tendon-compliance fix.
2. **Erector spinae / multifidus block (88 muscles, `scripts/msk/add_erector_spinae.py`).** The
   single densest multi-layer paraspinal compartment in this repo's own model by muscle count,
   and human thoracolumbar fascia is a well-studied real force-transmission structure in the
   broader literature (Vleeming/Willard-type work — **not fetched/verified this session**,
   flagged as a pointer for future work, not a claim). Whether it is the largest candidate by
   force-fraction affected is **untested here** — an honest gap, not asserted.
3. **Deep posterior/tarsal-tunnel compartment** (tibpost_r/fdl_r/fhl_r), already flagged in
   `scripts/msk/foot_multisegment.py` as a shared fibro-osseous canal — a plausible secondary
   candidate, not analyzed quantitatively in this pass.

Only candidate #1 is built and measured quantitatively below, per this task's explicit scope.
(Quadriceps were considered and set aside: the 4 heads already converge on one shared patellar
tendon, a mechanically different — and largely already-captured — coupling than a mid-belly
fascial shear path between separately-tendoned synergists.)

---

## 3. The geometric model

The physical degree of freedom the connective tissue actually senses is **relative sliding**
between two adjacent muscle bellies, not either muscle's absolute force — the shear-lag/
composite load-transfer picture Huijing's own 1999 title uses ("muscle as a collagen fiber
reinforced composite"), not a heuristic. From **live** model geometry (never assumed):

```
e_m(t) = [L_MTU,m(t) - L_MTU,m(zero-pose)] / L_opt,m
```

each muscle's own excursion away from a shared anatomical reference (all coordinates = 0),
normalized by its own optimal fiber length (the standard muscle-mechanics normalization — makes
the two muscles' excursions comparable despite different absolute L_opt: gasmed_r L_opt=74.0mm,
soleus_r L_opt=55.2mm, live-extracted from the model). Gastrocnemius crosses the **knee and
ankle**; soleus crosses **only the ankle** — so knee flexion moves gastrocnemius relative to
soleus without moving soleus at all, the exact "vary one synergist's relative position" lever
Maas & Sandercock (2008) used surgically, arising here for free from real gait kinematics:

```
Delta(t) = e_gastroc(t) - e_soleus(t)          [relative excursion mismatch]
```

`gasmed_r` is used as the "gastrocnemius" representative — **machine-justified**, not assumed:
its length trajectory correlates with `gaslat_r`'s at **r = 0.999436** over the real trial
(pre-registered floor 0.99), so the two heads move in near-perfect lockstep and collapsing them
to one excursion node loses essentially no information.

A shear/spring coupling with an explicit **slack (dead) zone** — directly justified by Maas &
Sandercock's own words for physiological position ("connecting tissues remain slack or operate
in early stress-strain regions"), not a numerical convenience:

```
frac(t) = f_max * clip( sign(Delta) * max(0, |Delta(t)| - Delta_slack) / Delta_scale, -1, +1 )
```

`frac(t) > 0` (gastrocnemius relatively more lengthened than soleus, beyond slack) ⇒
gastrocnemius sheds `frac(t)` of **its own** total tendon force to soleus (the literature's own
convention: "X% of F_proximal" — a fraction of the *donor's* own force). `frac(t) < 0` ⇒ soleus
sheds `|frac(t)|` of its own force to gastrocnemius. This is symmetric, bounded (can never shed
more than the donor has), and **conserves total triceps-surae tendon force by construction**
(machine-verified below — a deliberate simplification of a true 3+-body fascial network:
Rijkelijkhuizen 2005 shows force can also leak directly onto **bone**, a third path this 2-body
lumped model cannot represent; flagged in §6).

**Two regimes, each anchored to a different, named, verified literature number** (never an
invented stiffness constant):

| Regime | Δ_slack | f_max | Anchor |
|---|---|---|---|
| **DEFAULT** (physiological) | 0.20 | 0.14 | Maas & Sandercock 2008's physiological near-zero finding (large slack zone) + Huijing/Maas/Baan 2003's "14% typical" as the ceiling if ever breached |
| **SENSITIVITY** (upper-bound) | 0.00 | 0.405 | Rijkelijkhuizen 2005's 40.5% GM→plantaris/calcaneus number, no slack — "what if this pair's fascia behaved like the most-coupled triceps-surae-adjacent case actually measured" |

`Δ_scale = 0.10` (the breach magnitude at which `frac` saturates to `f_max`) is a disclosed
modeling dial in both regimes, not a measured constant — no study measures a subject/species-
specific stiffness for this exact pair (§6).

---

## 4. Results — measured on the real subject2 `walking1` trial (158 frames, 100 Hz)

All numbers below are machine-computed by `scripts/msk/myofascial_transmission.py`, read
directly from its own JSON output — not hand-derived.

**Kinematic signal (real, not synthetic):** `Delta(t)` ranges **[−0.2757, −0.0237]**
(std 0.0685) — **always negative** this trial: soleus is consistently relatively more lengthened
than gastrocnemius vs the zero-pose reference, driven mainly by knee flexion during swing
(knee_angle_r reaches 65.7°, ankle_angle_r ranges ±16°). `corr(|Delta(t)|, F_triceps_surae_total(t))
= −0.274` — **out of phase**: the relative-mismatch signal is largest during swing, when
triceps-surae force is near its minimum, and smallest during the stance/push-off window when
force peaks (t≈0.5–0.6s, per the independent-actuator SO baseline: gasmed_r 1140N, gaslat_r
406N, soleus_r 1448N peaks).

**Machine checks (both regimes):** force-conservation invariant `max|err| = 4.5e-13 N`
(**PASS**, tolerance 1e-6 N — floating-point exact); null control (`f_max=0`) reproduces the
independent-actuator baseline **exactly** (`max|transfer| = 0.0` N, **PASS**).

### DEFAULT regime (Δ_slack=0.20, f_max=0.14 — physiological-anchored)

- Max transfer over the whole cycle: **7.97 N**, occurring at t=0.92s — **during swing**, where
  both muscles are near their force floor (F_gas_tot=48.8N, F_sol=75.2N).
- At gastrocnemius's own peak-force frame (t=0.51s, 1538N) and soleus's own peak-force frame
  (t=0.61s, 1448N): transfer is **-0.00 N (0.0%)** — the slack zone comfortably contains this
  gait cycle's mismatch during the entire high-force stance window.
- Net ankle plantarflexion moment: peak 108.54 N·m, max change **0.011 N·m = 0.01% of peak**
  — **NEGLIGIBLE** (pre-registered ceiling 2.0%).
- **Reading:** at physiological gait ROM, with the coupling calibrated to the literature's own
  most-relevant physiological finding, the twin's existing independent-actuator prediction for
  this pair is **essentially unchanged** — a genuine, forced (not lazy) small-effect finding,
  consistent with Maas & Sandercock (2008)'s own conclusion for this exact pair.

### SENSITIVITY regime (Δ_slack=0.00, f_max=0.405 — upper-bound)

- Max transfer over the whole cycle: **420.1 N** at t=0.63s (F_gas_tot=871N, F_sol=1413N —
  both substantial, near their own stance-phase levels).
- At gastrocnemius's own peak-force frame (t=0.51s): transfer = **−93.0 N (−6.0%** of that
  frame's gastrocnemius force).
- At soleus's own peak-force frame (t=0.61s, 1448N): transfer = **−411.4 N (−28.4%** of that
  frame's soleus force) — a **substantial** individual-muscle redistribution.
- Net ankle plantarflexion moment: peak 108.54 N·m, max change **1.705 N·m = 1.57% of peak** —
  still **NEGLIGIBLE** by the pre-registered 2.0% ceiling, but **>100× larger** than the DEFAULT
  regime's 0.01%.
- **Reading:** even at the literature's own largest verified triceps-surae-adjacent fraction,
  the **joint-level** mechanics (net ankle moment) stays small — because gastrocnemius and
  soleus moment arms are nearly equal (§2, within ~5%), so redistributing force **between** them
  barely changes their **sum**. But the **individual-muscle-level** effect is not small (−28.4%
  of soleus's own peak) — meaning a per-muscle downstream prediction (this repo's own
  `metabolic_cost.py` and `crossbridge_contraction_model.py` both operate on soleus
  individually) could be measurably wrong even when the whole-joint torque looks fine. This is
  the genuine, symmetric, dual-sided finding: **the independent-actuator assumption is more
  defensible at the joint level than at the single-muscle level for this pair**, and the gap
  between the two is now quantified, not asserted.

---

## 5. Was the "small effect" finding forced, or lazy?

Per this task's own discipline: a one-shot negative is not accepted without forcing the
strongest adversary first. Here, the adversary to the DEFAULT regime's near-zero finding is
"the coupling mechanism itself is too weak/broken to ever show a real effect." That adversary
is falsified by the SENSITIVITY regime: the **same** code path, same conservation invariant,
same real kinematic signal — driven with a different, still literature-anchored (not invented)
parameter choice — produces a substantial, correctly-conserved, individually-large (−28.4%)
redistribution. The mechanism works; the DEFAULT regime's near-zero result is a genuine,
literature-anchored small-effect prediction, not a degenerate model.

---

## 6. Honest gaps (full list, stated up front in the script's own docstring too)

1. **Illustrative lumped first step, not continuum FEM.** Real myofascial transmission is a
   continuously-distributed shear field over the whole muscle-belly contact area; this collapses
   it to one scalar spring between two lumped nodes.
2. **Coupling constants are literature-anchored dials, not subject-specific measurements.**
   `Δ_slack`, `Δ_scale`, `f_max` are chosen to match the closest verified numbers (§1b), not
   measured for this subject, this species, or even human tissue at all — no such measurement
   exists anywhere in the literature checked (§1c).
3. **Every quantitative number is rat or cat.** Zero direct human force-partitioned data exists
   for any species/muscle-pair combination checked (§1c) — this is the literature's own gap, not
   a search failure; explicitly re-confirmed by a 2023 review (Finni/de Brito Fontana/Maas)
   still reporting no quantified human fraction.
4. **2-body-only — cannot represent 3rd-party leakage.** Rijkelijkhuizen (2005) shows force can
   transmit directly onto **bone** (calcaneus) via a muscle's fascial envelope, not only to a
   named neighboring muscle. This model's force-conservation-between-the-pair construction
   cannot represent that third path.
5. **Sign/redistribution convention is a defensible modeling choice, not independently verified
   for this pair.** The "more-relatively-lengthened muscle sheds to the less-lengthened one"
   rule follows the shear-lag/composite picture (§3) and is consistent with the EDL proximal/
   distal literature's qualitative direction, but no study directly measures the sign of
   gastrocnemius↔soleus transfer under gait-like (non-surgical) relative-position changes.
6. **n=1 subject, one trial (`walking1`)** — the same breadth-of-validation ceiling already
   flagged repo-wide in `docs/MECHANISM_FIDELITY_ARCHITECTURE_AUDIT.md`'s scope note; not a new
   gap introduced here, but not separately closed here either.
7. **`gasmed_r` stands in for "gastrocnemius"** (both heads combined for the force ledger, but
   the excursion signal driving `Delta(t)` uses `gasmed_r` alone) — justified by the measured
   0.999436 length-trajectory correlation with `gaslat_r`, not assumed.
8. **The DEFAULT/SENSITIVITY split is a bracketing device, not a confidence interval** — the
   true value for this specific human pair is unmeasured and could in principle sit outside
   this bracket; the bracket's endpoints are real literature numbers, but neither is a direct
   measurement of gastrocnemius↔soleus in an intact human.
9. **Erector spinae / thoracolumbar fascia (§2, candidate #2) was not quantitatively assessed**
   — flagged as the largest-by-muscle-count candidate, but whether it is the largest by
   force-fraction affected is untested.

---

## Files

- `scripts/msk/myofascial_transmission.py` — the build + verification script (literature
  constants with PMID/DOI, live geometry extraction, both coupling regimes, conservation +
  null-control machine checks, net-ankle-moment analysis, JSON output).
- `data/msk_smoketest/myofascial_transmission/myofascial_transmission_results.json` —
  machine-readable results (muscle properties, correlation checks, both regimes' full
  158-frame time series, thresholds, literature dict).
- Reused, read-only: `scripts/msk/validate_joint_force.py` (`parse_mot`, model/IK trial path
  constants), `data/msk_smoketest/subject2_walking1/static_optimization/so/
  walking1_StaticOptimization_force.sto` (the existing, already-verified independent-actuator
  baseline this build redistributes — not re-run, no new Static Optimization was performed).
- Context: `docs/MECHANISM_FIDELITY_ARCHITECTURE_AUDIT.md` §3 (the gap this build addresses),
  `scripts/msk/tendon_elastic_energy.py` (the related-but-distinct 3-independent-tendons
  finding), `scripts/msk/foot_multisegment.py` (the deep-posterior-compartment candidate #3).

ISOLATION: bodytwin only; external subject2 model/IK read in place (never written to); no git
commit or push performed (per task).
