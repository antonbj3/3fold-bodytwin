# MECHANISM TENDON ELASTIC ENERGY — the body's own spring (Achilles/patellar/hamstring), subject2 walking1 (2026-07-21)

Complements the resistance-band (`docs/MECHANISM_MSK_ELASTIC_BAND.md`, `PathSpring`) and knee-ligament
(`docs/MECHANISM_KNEE_LIGAMENTS.md`, `Blankevoort1991Ligament`) elastic-element work already done with the
body's **own** elastic material: the Achilles/triceps-surae tendon as a biological spring, plus the patellar
(quadriceps) and hamstring tendon groups. Script: `scripts/msk/tendon_elastic_energy.py`
(`.venv-msk/bin/python3 scripts/msk/tendon_elastic_energy.py`, exit 0). REUSES the already-computed,
already-gate-verified Static Optimization tendon forces for subject2 `walking1`
(`data/msk_smoketest/subject2_walking1/static_optimization/so/`, the same file
`docs/MECHANISM_STATIC_OPT.md` and `validate_ankle_force.py` build on) — **no re-solve, no Ipopt re-run.** The
only new computation is a pure forward-kinematics pass (muscle-tendon-unit path length at the
already-validated IK poses), combined with the model's own live Millard tendon force-length **curve
object** (introspected, not hand-reimplemented) to invert force→strain and integrate strain→stored energy.

## Headline

| tendon group | peak stored energy (J) | when | vs Lichtwark & Wilson 2005 hopping anchor (38 J) |
|---|---:|---:|---:|
| **Achilles / triceps-surae** (gasmed_r+gaslat_r+soleus_r) | **7.76 J** | t=0.53 s (mid-trial, clean) | ratio 0.20 |
| **Patellar / quadriceps** (recfem_r+vasmed_r+vaslat_r+vasint_r) | **4.83 J** (full window) / **2.84 J** (interior-only) | t=1.57 s (recording boundary, flagged §5) / t=0.49 s | — |
| **Hamstring** (bflh_r+bfsh_r+semimem_r+semiten_r) | **1.55 J** | t=1.19 s (mid-trial, clean) | — |

**Achilles/triceps-surae elastic-vs-fiber positive-work split** (§6): of **25.9 J** total positive
muscle-tendon-unit work delivered by the triceps-surae group over the gait cycle, **7.82 J (29.8%)** comes
from **elastic tendon recoil** and **18.46 J (70.2%)** from **active muscle-fiber shortening**. Walking
stores roughly a **fifth** of what Lichtwark & Wilson (2005) directly measured in one-legged hopping (a far
more elastic-energy-demanding activity) — consistent with the task's own expectation that walking stores
far less than running/hopping, not an over-claim.

**Both watertight self-tests (synthetic, known-answer controls) PASS**, all 3 external-anchor plausibility
gates PASS, and the one internal-consistency check that initially screamed a problem (274–3644% "relative
error") was diagnosed via OODA to be a metric artifact (divide-by-near-zero at direction reversals), not a
real bug — fixed at the source, not waved through (§4.3).

---

## 0. Headline honesty flag — the task's own "Ker 1987 PMID 3808070" citation is the WRONG paper

Verified LIVE via NCBI eutils (`efetch`, abstract text fetched directly, not recalled), matching this
repo's established citation-verification discipline (`docs/MECHANISM_KNEE_LIGAMENTS.md` §0):

**PMID 3808070 = Ker RF, Bennett MB, Bibby SR, Kester RC, Alexander RM. "The spring in the arch of the
human foot." Nature. 1987;325(7000):147-9. DOI 10.1038/325147a0.**

This paper's own abstract is explicitly about the **plantar arch** of the human foot, a *different* elastic
structure from the Achilles tendon, and explicitly contrasts itself with prior tendon-specific work:

> "Among the elastic structures involved, **the tendons of distal leg muscles have been shown to be
> important** [prior work]. Here we show that **the elastic properties of the arch of the human foot are
> also important**."

Citing it for "Achilles tendon stores tens of J in running" would misattribute a foot-arch-spring finding to
the Achilles tendon. **Not silently substituted — flagged, then replaced with directly-verified sources**:

- **Lichtwark GA, Wilson AM. "In vivo mechanical properties of the human Achilles tendon during one-legged
  hopping." J Exp Biol. 2005;208(Pt 24):4715-25. PMID 16326953, DOI 10.1242/jeb.01950.** (This IS the
  task's "Lichtwark & Wilson 2005" — correctly named, independently confirmed.) Verbatim: *"An average of 38
  J of energy was recovered from the elastic recoil of the tendon, which contributes 16% of the total
  average mechanical work of the hop (254 J)."* Peak strain measured: 8.3%.
- **Fukunaga T, Kubo K, Kawakami Y, Fukashiro S, Kanehisa H, Maganaris CN. "In vivo behaviour of human
  muscle tendon during walking." Proc Biol Sci. 2001;268(1464):229-33. PMID 11217891, DOI
  10.1098/rspb.2000.1361.** A **walking**-specific (not hopping/running) anchor. Verbatim: *"the tendon
  stretched by ca. 7 mm during single support, and recoiled in push-off... the muscle contracted
  near-isometrically in the stance phase."*
- **Lichtwark GA, Bougoulias K, Wilson AM. "Muscle fascicle and series elastic element length changes along
  the length of the human gastrocnemius during walking and running." J Biomech. 2007;40(1):157-64. PMID
  16364330, DOI 10.1016/j.jbiomech.2005.10.035.** Qualitative mechanism support, verbatim: *"Muscle
  fascicles acted relatively isometrically during the stance phase during walking, however during running
  the fascicles shortened throughout the stance phase... high velocity shortening during take off in both
  walking and running achieved by recoil of the SEE [series elastic element]."* — predicts exactly the
  near-isometric-fascicle / tendon-recoil mechanism this build measures.
- **Alexander RM, Bennet-Clark HC. "Storage of elastic strain energy in muscle and other tissues." Nature.
  1977;265(5590):114-7. PMID 834252, DOI 10.1038/265114a0.** Real paper, existence/title/DOI verified live —
  but PubMed carries **no indexed abstract** for this pre-abstract-era Nature note, so a specific Joules
  figure could not be independently verified live from it this session. Cited for historical
  existence/context only; **no number is attributed to it** in this report.
- **Giddings, Beaupre, Whalen & Carter 2000** (already verified live in this repo,
  `data/msk_smoketest/subject2_walking1/ankle_force_validation/`): Achilles tendon **force** (not energy),
  walking = 390 %BW. Used below (§4.4) as an independent force-level cross-check.

---

## 1. Method — the geometric derivation (no re-solve)

Let `L_mt(t)` be a muscle-tendon unit's path length (pure forward kinematics), `L_tendon(t)` its tendon's own
length, and `F(t) ≥ 0` its tendon tension (from the existing Static Optimization `force.sto`, reused
unchanged). Define the along-tendon-axis fiber projection `L_fiber_proj(t) := L_mt(t) - L_tendon(t)` — exact
by this definition, independent of any pennation-modeling detail. Then:

```
E(t)     = Fmax · tendon_slack_length · [curve.calcIntegral(1+strain(t)) − curve.calcIntegral(1)]   (Joules)
dE/dt    = F(t) · v_tendon(t)                              (exact chain rule on E(strain(t)))
P_out(t) = −F(t) · v_mt(t)                                  (power the MTU delivers to the skeleton)
```

Since `v_mt = v_tendon + v_fiber_proj` (exact, from the length identity), `P_out(t) = [−dE/dt(t)] +
P_fiber_out(t)`, where `P_fiber_out(t) := −F(t)·v_fiber_proj(t)` is the along-tendon-axis fiber/contractile-
element power — the convention that makes the power balance close **exactly**, with no pennation-rate
residual (the standard tendon-vs-fascicle work-splitting convention in this literature, e.g. Lichtwark 2007
above). Only the tendon's own force-length curve is used — the fiber's active/passive force-length or
force-velocity curves are never touched, exactly the task's "Millard tendon force-length properties" scope.

**Curve API, confirmed live** (`.venv-msk`, opensim 4.6-2026-06-22-85aaf64): the model's muscles are
`Millard2012EquilibriumMuscle` (the OpenSim class implementing Millard, Uchida, Seth & Delp 2013 "Flexing
Computational Muscle," J Biomech Eng 135(2):021005 — the paper the "Millard2010" shorthand in the task
refers to). Each exposes `getTendonForceLengthCurve()` → a live `TendonForceLengthCurve` object whose
argument convention is **normalized tendon length `x = 1+strain`** (confirmed: `calcValue(1.0)=0`,
`calcValue(1.049)=1.000` exactly matching `strain_at_one_norm_force=0.049`, linear slope beyond =
`stiffness_at_one_norm_force_in_use=28.061 = 1.375/0.049`), plus a **`calcIntegral(x)`** method — the curve's
own antiderivative, used directly rather than hand-reimplementing the curve's shape (avoids introducing a
second, possibly-mismatched curve).

## 2. Two watertight self-tests (synthetic controls) — run FIRST, both PASS

Per the discipline that a build must force its own concept via known-answer controls before trusting real
data:

**[1] Curve-integral usage.** At 4 strains (0.01–0.07), `curve.calcIntegral()`-derived energy was
cross-checked two independent ways: (a) an independent from-scratch trapezoidal quadrature of the same
curve's `calcValue()` over 20,001 points, and (b) OpenSim's own **native** `Muscle.getTendonPotentialEnergy(state)`,
computed on a **fresh** `model.initSystem()` state per test point (reproducing — and avoiding — the exact
stale-cached-state gotcha already documented in `docs/MECHANISM_KNEE_LIGAMENTS.md` §3: reusing one state
across `setFiberLength()` calls silently freezes cached energy at its first-computed value). All 4 strains:
**0.0000% relative error both ways** (e.g. strain=0.049 → my formula = OpenSim native = 29.933442 J,
bit-for-bit).

**[2] Power-decomposition sign/algebra.** A synthetic sinusoidal `strain(t)` and an independent synthetic
along-tendon fiber-length trajectory (**both with known closed-form analytic derivatives**, no finite
differencing) were pushed through the actual production functions (table-based inversion, energy
integration, then the same Savitzky-Golay derivative pipeline used on real gait data). Results: round-trip
strain inversion 0.0000%, recovered `v_tendon` vs analytic truth 0.0000% (RMS-normalized), recovered
`P_fiber_out` vs the independently-defined analytic ground truth **0.044%** max relative error. This is the
kind of known-answer toy test that would catch an exact-sign error (this repo has a documented history of
exactly that failure mode, `docs/MECHANISM_SIGN_BUG_AUDIT.md`) — it did not find one.

## 3. Muscle tendon properties (11 muscles, live-introspected)

| muscle | Fmax (N) | tendon_slack_length (m) | e0 (strain@1-norm-force) | k (stiffness@1-norm-force) |
|---|---:|---:|---:|---:|
| gasmed_r | 3115.5 | 0.4852 | 0.049 | 28.06 |
| gaslat_r | 1575.1 | 0.4690 | 0.049 | 28.06 |
| soleus_r | 6194.8 | 0.3527 | 0.049 | 28.06 |
| recfem_r | 2191.7 | 0.5374 | 0.049 | 28.06 |
| vasmed_r | 2747.8 | 0.2516 | 0.049 | 28.06 |
| vaslat_r | 5148.8 | 0.2671 | 0.049 | 28.06 |
| vasint_r | 1697.4 | 0.2479 | 0.049 | 28.06 |
| bflh_r | 1313.2 | 0.4060 | 0.049 | 28.06 |
| bfsh_r | 557.1 | 0.1278 | 0.049 | 28.06 |
| semimem_r | 2201.0 | 0.4078 | 0.049 | 28.06 |
| semiten_r | 591.3 | 0.3032 | 0.049 | 28.06 |

**Confirmed live: all 11 muscles use the identical DEFAULT tendon-curve shape** (e0/k, hence the entire
force-strain relationship, differ ONLY in the (per-muscle) `tendon_slack_length`/`Fmax` scale — the curve's
*dimensionless shape* is model-wide generic). `ignore_tendon_compliance=False` confirmed for all 11 (elastic
tendon active, not the rigid-tendon simplification). This is the single most important honesty flag in this
report (§7.1).

## 4. Per-muscle peak stored energy + force cross-check (all 11 muscles, full 158-frame gait cycle)

| muscle | peak F (N) | peak F/Fmax | peak strain | peak stretch (mm) | peak E (J) | at t (s) |
|---|---:|---:|---:|---:|---:|---:|
| gasmed_r | 1140.0 | 36.6% | 0.0250 | 12.13 | 5.646 | 0.51 |
| gaslat_r | 405.7 | 25.8% | 0.0198 | 9.30 | 1.570 | 0.53 |
| soleus_r | 1447.8 | 23.4% | 0.0186 | 6.55 | 3.965 | 0.61 |
| recfem_r | 843.8 | 38.5% | 0.0259 | 13.90 | 4.781 | 1.57 (boundary, §5) |
| vasmed_r | 185.4 | 6.7% | 0.0075 | 1.90 | 0.140 | 0.06 |
| vaslat_r | 694.2 | 13.5% | 0.0124 | 3.31 | 0.948 | 0.06 |
| vasint_r | 70.6 | 4.2% | 0.0054 | 1.33 | 0.035 | 0.06 |
| bflh_r | 104.7 | 8.0% | 0.0085 | 3.45 | 0.146 | 1.25 |
| bfsh_r | 60.4 | 10.8% | 0.0106 | 1.36 | 0.034 | 1.14 |
| semimem_r | 472.1 | 21.4% | 0.0175 | 7.13 | 1.410 | 1.19 |
| semiten_r | 31.8 | 5.4% | 0.0064 | 1.94 | 0.024 | 1.25 |

**4.1 Round-trip inversion gate — PASS, all 11 < 0.001%.** Forward-evaluating the curve at each recovered
strain reproduces the target SO-derived force to ≤0.0008% for every muscle at every frame.

**4.2 Peak strains are tiny (0.5–2.6%), well inside the curve's nonlinear "toe" region** (which ends at
strain=4.9%) — walking never drives any of these 11 tendons into the stiff linear region. This is the
direct, geometric reason the strain-inversion is numerically delicate at low force (§4.3): the curve's own
derivative `dF/dstrain → 0` as strain→0 (confirmed live, `calcDerivative(1.0,1)=0.0000`), so near slack the
*inverse* map has `dstrain/dF → ∞`, amplifying any small force fluctuation into a proportionally larger
strain fluctuation.

**4.3 An OODA-forced diagnosis, not a wave-through.** An internal consistency check — `dE/dt` (computed via
Savitzky-Golay differentiation of `E(t)` directly, what production actually uses) vs `F(t)·v_tendon(t)` (the
exact chain-rule expression, via differentiating the strain-inverted tendon length) — initially reported
**274–3644% pointwise relative error** on real data, despite the identical check passing at 0.0007% on the
smooth synthetic control (§2). This was NOT accepted as "some numerical noise" and moved past. Diagnosis:
printing the raw signals showed `dE/dt` and `F·v_tendon` actually agree closely in **scale** (RMS 15.22 vs
15.25 W, peak range −42.1..+47.5 vs −42.3..+47.7 W) with a small **absolute** difference (RMS 0.26 W, max
1.11 W) — the enormous *relative* error occurred exactly at samples where `dE/dt≈0` (a direction-reversal
point of `E(t)`, an inevitable near-zero denominator for ANY trial with more than one local extremum) — a
"fixed threshold on a statistic that is regime-blind near a zero-crossing" artifact, not a physics/algebra
bug. **Fixed at the source**: replaced the pointwise-relative-error gate with an RMS/peak-**scale**-normalized
one (both for this real-data check and, symmetrically, inside the self-test itself). Re-measured:
**RMS-normalized relative error 1.7–2.7%, peak-normalized 2.3–3.8%** across the 3 triceps-surae muscles —
gate (<15%) **PASS**. A cheap, independent 3rd differentiation method (`numpy.gradient`, plain central
difference) confirms `dE/dt` itself is not noise-dominated (RMS deviation from Savitzky-Golay's own estimate
is 3.6% of signal scale).

**4.4 Over-determination against an independently-built sibling script.** This script's peak per-muscle
tensions (gasmed_r 1139.972 N, gaslat_r 405.675 N, soleus_r 1447.798 N) match
`data/msk_smoketest/subject2_walking1/ankle_force_validation/ankle_force_validation_results.json`'s own
`per_muscle_peak_tension_N_context` **bit-for-bit (0.0000% difference)** — both scripts independently read
the same underlying SO output correctly. That sibling script's own **time-aligned** scalar-sum peak (320.26
%BW at t=0.59s) was already gated **PASS** there against the Giddings et al. 2000 published walking Achilles
force anchor (390 %BW, ratio 0.82, within its pre-registered (0.30,3.00) plausibility band) — meaning the
force trajectory this energy calculation is built on is *itself* already literature-anchored at the force
level, one layer beneath the energy-specific anchors below (§6).

## 5. Peak GROUP-level ("per tendon") stored energy — the headline number

Peak of the **summed** group energy trajectory (not the sum of each muscle's own individual peak, which
generally occurs at a different instant per muscle — both reported, for full disclosure):

| group | peak (summed-then-peaked) | sum of individual muscle peaks | same instant? |
|---|---:|---:|---|
| Achilles/triceps-surae | **7.760 J** @ t=0.53 s | 11.180 J | different (0.51/0.53/0.61 s) |
| Patellar/quadriceps | **4.828 J** @ t=1.57 s | 5.904 J | different (1.57/0.06/0.06/0.06 s) |
| Hamstring | **1.551 J** @ t=1.19 s | 1.613 J | different (1.25/1.14/1.19/1.25 s) |

**Boundary-effect flag on the patellar/quadriceps group (forced check, not assumed clean):** its peak sits
at **t=1.57 s — the very last recorded IK/SO sample**, driven by `recfem_r` (rectus femoris), whose force is
**still rising steeply** at that instant (26 N at t=1.46s → 844 N at t=1.57s, no sign of turning over). This
repo's own `docs/MECHANISM_STATIC_OPT.md` documents a real prior failure mode at exactly this kind of
boundary (spline/extrapolation edge effects) — so this was checked, not trusted. Two things distinguish this
from that prior bug: (a) t=1.57s is the genuine data ceiling, not an extrapolated point beyond it (already
fixed upstream, in `static_opt_knee.py`'s own construction of this SO run); (b) `recfem_r` shows the
**identical steep-rise signature** at the *other* boundary too (25.0→25.0→42.7→95.1 N over t=0.07–0.09s,
comfortably mid-trial with data on both sides) — i.e. this is very likely a genuine, repeatable
early-stance/loading-response recruitment burst for this biarticular hip-flexor/knee-extensor (this walking
trial spans **two** right-foot stance phases, `docs/MECHANISM_STATIC_OPT.md` §5), not a numerical artifact —
but the trial happens to **end mid-burst** for the second occurrence, so the true peak of that specific event
is plausibly **higher** than the 4.83 J captured (right-censored, not a fabricated number). The conservative,
non-truncated alternative — the interior-only peak (excluding the outer 5 frames each side, same convention
`docs/MECHANISM_STATIC_OPT.md` already established) — is **2.837 J at t=0.49s**. Neither the Achilles nor the
hamstring group's peak is boundary-adjacent (0.53s and 1.19s respectively, both comfortably interior).

## 6. Achilles/triceps-surae: elastic recoil vs. active fiber, positive-work fraction

Per-muscle decomposition (interior window, excluding the outer 5 frames each side for derivative validity):

| muscle | positive MTU-out work (J) | elastic-recoil positive work (J) | fiber-active positive work (J) |
|---|---:|---:|---:|
| gasmed_r | 8.711 | 5.684 | 7.329 |
| gaslat_r | 3.869 | 1.608 | 2.621 |
| soleus_r | 13.435 | 4.901 | 12.395 |
| **GROUP TOTAL** | **25.905** | **7.822** | **18.462** |

Reconciliation check `|W_out − (W_elastic+W_fiber)| / W_out` = **1.47%** (expected to be nonzero in general —
`max(0,a+b) ≠ max(0,a)+max(0,b)` pointwise when signs can differ at some instants — small here, disclosed
not hidden).

**FRACTION ELASTIC = W_elastic / (W_elastic + W_fiber) = 0.2976 → 29.8%.**

**Scope precision (important, do not conflate with Lichtwark & Wilson's own 16%):** this 29.8% is the
fraction of the triceps-surae group's **own** positive muscle-tendon-unit work. Lichtwark & Wilson's 16% is
the fraction of the **entire hop's total mechanical work (254 J)** — a much broader, whole-body/whole-leg
denominator. These are **not the same percentage scale** and must not be read as "walking (29.8%) is more
elastic than hopping (16%)" — they answer different questions (MTU-internal split vs. whole-activity
energetics). The **absolute** Joules comparison (§0, §5: 7.76 J walking vs. 38 J hopping, ratio 0.20) is the
valid cross-activity comparison; the **percentage** is only valid compared within this same MTU-scoped
definition.

**Context cross-check (decorrelated method, not a strict identity):** total positive `ankle_angle_r` JOINT
work over the same window, computed independently via the already-verified-clean-for-ankle
`InverseDynamicsTool` moment (`docs/MECHANISM_MSK_ELASTIC_BAND.md` §4 confirms ankle, unlike knee/hip, is NOT
corrupted by the documented OpenSim 4.6 bug) times the IK-derived angular velocity = **21.87 J**. The
triceps-surae MTU total (25.90 J) is **118%** of this — plausible, not alarming: gastrocnemius is biarticular
and also crosses the knee, so its own MTU-level power reflects work exchanged at *both* ends of its path,
not exclusively the ankle DOF; a >100% ratio from a strict-subset-of-crossing-muscles comparison in a
different reference frame is expected, not a red flag.

## 7. Honest gaps (full list — read before trusting any number quantitatively)

1. **DOMINANT UNCERTAINTY (flagged per the task's own explicit ask):** the tendon force-length curve *shape*
   (`strain_at_one_norm_force=0.049`, `stiffness_at_one_norm_force=28.06`) is the OpenSim
   `Millard2012EquilibriumMuscle` **default**, confirmed identical across all 11 target muscles (and, as
   shipped, essentially model-wide) — not ultrasound/MRI-measured, subject-specific tendon compliance. A
   stiffer or more compliant real tendon changes stored energy substantially (energy scales roughly with
   1/stiffness for a given force, given the curve's toe-then-linear shape) — every Joule number in this
   report inherits this uncertainty.
2. **Achilles tendon is modeled as THREE INDEPENDENT tendons** (gasmed_r/gaslat_r/soleus_r each with their
   own `tendon_slack_length`/curve instance), not the single confluent structure they form anatomically.
   "Achilles tendon energy" here = the sum of three separately-modeled springs, not a biofidelic
   single-free-body tendon. The same caveat applies to the **patellar/quadriceps** group (4 anatomically-
   converging muscles, modeled as 4 independent tendons).
3. The **hamstring** group is not one anatomical tendon even in reality (biceps femoris, semitendinosus via
   pes anserinus, and semimembranosus have 3 distinct distal insertions) — reported as a group sum by scope
   convention, not because the anatomy or model treats it as a single structure.
4. **Zero tendon hysteresis** in this model as used here (`Millard2012EquilibriumMuscle`'s tendon curve is a
   purely elastic spring — everything stored is, by construction, fully returned). Real tendon has measured
   hysteresis (17–35% inter-quartile range, per Lichtwark & Wilson 2005's own hopping data) — this likely
   makes the 29.8% "fraction elastic" an **upper bound** on true recoil efficiency, not a hysteresis-corrected
   estimate.
5. **Patellar/quadriceps peak is right-censored** (§5): the full-window peak (4.83 J) occurs at the last
   recorded frame during an apparently still-rising, genuine recruitment burst — likely an undercount of
   that burst's true peak. The interior-only alternative (2.84 J) is the conservative, non-truncated number.
6. **Along-tendon-axis fiber power convention**, not the fiber's own true fiber-frame power (these differ by
   a pennation-angle projection factor) — the standard convention in this literature for an exactly-closed
   power balance (§1), disclosed not hidden.
7. **The "fraction elastic" denominator scope differs from Lichtwark & Wilson's own 16%** (§6) — not
   directly comparable percentages; only the absolute-Joules comparison (§0) is cross-activity-valid.
8. **Right leg only, single trial** (subject2 `walking1`) — same scope caveat as every other cert in this
   family; no claim of generality across subjects/speeds/trials.
9. **Static Optimization's own frame-by-frame force estimate is not measured EMG** — inherited, unchanged
   structural caveat from `docs/MECHANISM_STATIC_OPT.md`; SO force is an activation-minimizing solution.
10. The ankle-joint-work context cross-check (§6) combines an ID-Tool moment with an IK-derived angular
    velocity — a genuinely different computational path than the MTU-force-based method, so it is a
    plausibility cross-check, not a strict energy-conservation identity.

## 8. Gates (pre-registered, machine-checked)

| gate | threshold | measured | verdict |
|---|---|---:|---|
| Self-test 1: curve-integral usage (trapz + native OpenSim PE) | <0.1% relerr | 0.0000% | PASS |
| Self-test 2: power-decomposition sign/algebra (synthetic) | <3% relerr | 0.044% (fiber power), 0.0007% (dE/dt) | PASS |
| Sibling SO convergence/activation gate (independently re-verified) | frames/NaN/bounds | 158/158, 0 NaN, in bounds | PASS |
| Round-trip force→strain→force inversion, all 11 muscles | <1% | max 0.0008% | PASS |
| dE/dt (savgol) vs F·v_tendon (chain rule), RMS-normalized | <15% | 1.7–2.7% (3 muscles) | PASS |
| Achilles/triceps-surae peak energy vs Lichtwark&Wilson 2005 hopping ceiling | 0.02–38 J | 7.76 J (ratio 0.20) | PASS |
| gasmed_r/gaslat_r peak tendon stretch vs Fukunaga 2001 walking anchor (~7mm) | 1–30 mm band | 12.13 / 9.30 mm | PASS |

**Overall: PASS** (script exit code 0).

## Files

- `scripts/msk/tendon_elastic_energy.py` — the full pipeline (self-contained, re-runnable; imports
  `validate_joint_force.py` for proven `parse_mot`/Savitzky-Golay code and `static_opt_knee.py` for the
  proven SO-gate re-verification code; two synthetic self-tests run first and hard-fail the script if either
  does not pass).
- `data/msk_smoketest/subject2_walking1/tendon_elastic_energy/tendon_elastic_energy_results.json` — every
  number in this document, machine-written (self-test detail, muscle properties, per-muscle and per-group
  energy trajectory peaks, the full elastic-vs-fiber decomposition, all anchors with verbatim quotes/PMIDs/
  DOIs, all gates).
- Reused, unchanged: `data/msk_smoketest/subject2_walking1/static_optimization/so/` (SO force/activation),
  `data/msk_smoketest/subject2_walking1/joint_force_validation/walking1_id.sto` (ankle moment context,
  §6), the scaled `LaiArnoldModified2017_poly_withArms_weldHand_scaled.osim` model (read-only, mtime
  unchanged).

## Next step

1. If a subject-specific tendon compliance estimate is ever available (ultrasound-measured
   strain-at-a-given-force for THIS subject's Achilles/patellar tendons), re-run with a per-muscle
   `strain_at_one_norm_force`/`stiffness_at_one_norm_force` override — the dominant uncertainty (§7.1) would
   shrink from "generic default" to "subject-measured."
2. A running or hopping trial for subject2 (if one exists in the LabValidation corpus) would let this exact
   script report a same-subject, same-model, same-method walking-vs-running/hopping comparison — a much
   tighter test of the "walking stores far less" claim than comparing against a different population's
   published number.
3. Left leg / other trials, if generality beyond this single trial is later needed (same scope caveat as
   the rest of this cert family).
