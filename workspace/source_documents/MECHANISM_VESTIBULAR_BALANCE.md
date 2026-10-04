# MECHANISM VESTIBULAR/BALANCE — VOR canal dynamics + inverted-pendulum postural sway (2026-07-22)

Resolves the vestibular/postural-balance layer named in `docs/MECHANISM_STATE.md` (`ORG-POSTURAL-BALANCE`,
"fall-risk") and `bt_memory/SESSION_HANDOFF.md`'s remaining-hardening list, coupling the sensory system
to the twin's existing MSK postural work. Two decorrelated, mechanistic (geometry-derived, not curve-fit)
models: (1) the vestibulo-ocular reflex (VOR) as a canal torsion-pendulum transfer function, and (2)
quiet-standing postural sway as a single-link inverted pendulum with delayed active feedback. Scripts:
`scripts/msk/vestibular_vor.py`, `scripts/msk/postural_sway_pendulum.py`. Raw results:
`reports/probes/vestibular_vor_results.json`, `reports/probes/postural_sway_pendulum_results.json`.
Consolidated citation ledger: `docs/MECHANISM_VESTIBULAR_BALANCE_evidence.json`.

**Scope note before anything else**: this repo's graph (`data/MECHANISM_ANCHOR_GRAPH.json`) already
contains an **executed, MEASURED-B PASS** cell, `ORG-POSTURAL-BALANCE` (2026-07-18), that analyzed REAL
PhysioNet HBEDB force-plate data (163 subjects, byte-verified via SHA256 against the dataset's own
manifest) and found the Romberg ratio to be **1.155 (firm surface) / 1.304 (foam surface)** — not the
task brief's textbook "~1.5-2×" prior. That existing cell is a **statistical analysis of real COP data**;
it is not touched or re-run here (isolation discipline: re-verify, build a companion, don't edit). This
doc's postural script instead builds the **mechanistic model** that existing cell didn't have, and uses
its own already-measured, real numbers as the external falsifier anchor — a stronger anchor than a fresh
literature pull. Two further existing nodes (`SNS-VESTIBULAR-BALANCE`, `SENS-VESTIBULAR-BALANCE`) are
SEED-DESIGN-only clinical/epidemiological hypotheses (fall-risk vs. lesion laterality, DHI handicap
correlation) — genuinely different scope from the basic-physiology transfer-function model built here;
confirmed by direct read, not assumed.

## Part 1 — VOR: canal torsion-pendulum transfer function

### Method, in one paragraph

The semicircular canal's cupula-endolymph system is a heavily overdamped torsion pendulum (Fernandez &
Goldberg 1971, PMID 5000363). Because its mechanical time constant (~ms) is far below the physiological
band, the canal-only transfer function in the 0.001-2 Hz range collapses to a single-pole **high-pass
filter** of head velocity: `H(s;tau) = tau*s/(1+tau*s)`, `gain = wt/sqrt(1+wt^2)`, `phase = 90-atan(wt)`
(deg), `wt = 2*pi*f*tau`. Central "velocity storage" (vestibular nuclei/nodulus-uvula loop) re-uses the
identical shape with an EXTENDED tau (Karmali 2019, PMID 31239137: peripheral tau=5.7s, extended
centrally to 6-30s) — pushing the low-frequency corner down and rescuing low-frequency gain. Both regimes
are evaluated; the closed-form gain/phase is cross-checked against an independent `scipy.signal` LTI
evaluation (two different computational paths, not one formula trusted on its own).

### Citations — verified LIVE this session (NCBI eutils/PubMed/DOI fetch, not recalled)

| source | verified as | what it anchors |
|---|---|---|
| Fernandez C, Goldberg JM (1971). *J Neurophysiol* 34(4):661-75. **PMID 5000363** | Bibliographic match confirmed (title/journal/vol/pages); no abstract available (pre-abstracting era) | PRIMARY source, torsion-pendulum canal sinusoidal-response dynamics. |
| Karmali F (2019). "The velocity storage time constant: Balancing between accuracy and precision." *Prog Brain Res*. **PMID 31239137**, DOI 10.1016/bs.pbr.2019.04.038 | Abstract fetched live | **Canal tau=5.7s; central velocity-storage tau range 6-30s.** Refines the task brief's generic "4-6s/15-20s" prior with a specific, live-sourced number. |
| McGarvie LA et al (2015). "The Video Head Impulse Test (vHIT)... Age-Dependent Normative Values of VOR Gain." *Front Neurol*. DOI 10.3389/fneur.2015.00154 | **Full open-access article fetched live** (doi.org -> frontiersin.org redirect chain followed) | Horizontal-canal VOR gain 95% CI "include or are very close to 1.0" from age 10-19y through 80-89y — anchors this model's gain-plateau=1.0 normalization. Regime caveat: vHIT is a high-acceleration impulse, not sinusoidal rotary chair. |
| Bouveresse A et al (1998). *Acta Otorhinolaryngol Belg* 52(3):207-14. **PMID 9810455** | Abstract fetched live | n=52 healthy subjects, 0.01-0.64 Hz pseudorandom rotary-chair stimulation — confirms a **linear transfer-function** description of the human VOR is empirically supported over exactly this band. Exact per-frequency numbers paywalled (honest gap, disclosed). |
| Baloh RW et al (1988). "Ultralow vestibulo-ocular reflex time constants." *Ann Neurol* 23(1):32-7. **PMID 3422799**, DOI 10.1002/ana.410230107 | Abstract fetched live | DISEASE population (Chiari I / brainstem-cerebellar atrophy, n=3): tau<2s; "no response... below 0.2Hz but normal gain...above 0.4Hz" — an independent, different-instance-space (diseased, not healthy) confirmation of the corner-frequency structure. |

### Headline results (all numbers from `vestibular_vor.py`'s own JSON output, not hand-typed)

| check | result |
|---|---|
| Machine crosscheck (closed-form vs. independent scipy LTI evaluation) | max error **2.2e-16** (peripheral), **1.1e-16** (central), **1.1e-16** (disease) — PASS |
| **Falsifier 1** (pre-registered: gain in [0.90,1.00] at 0.1/0.5/1.0 Hz) | peripheral-only (tau=5.7s): **0.9632, 0.9984, 0.9996**. Central (tau=17.5s): **0.9959, 0.9998, 1.0000**. **PASS_all = True** |
| Void-floor, frequency sweep (0.0001-3.16 Hz at tau=17.5s) | gain spans **0.011 to 1.000** — a real, non-degenerate, ~90-fold dynamic range, not a disguised constant |
| Void-floor, tau sweep (1-30s at f=0.03 Hz) | gain ratio (tau=30 / tau=1) = **5.32×**, strictly monotonic — PASS |
| Disease dissociation (Baloh 1988, structural/qualitative) | low-to-high-band gain drop: disease (tau=1.5s) **0.340** vs. healthy (tau=17.5s) **0.0072** — a ~47× larger low-frequency-specific deficit in the disease case, matching the clinical dissociation's DIRECTION and corner-frequency STRUCTURE |

**Low-vs-high frequency split — held OPEN by design, per task instruction**, not force-collapsed to one
verdict: at 0.1-1 Hz, peripheral-only and central-extended models converge (both ≈0.96-1.0, matching the
"~0.9-1.0" falsifier). Below ~0.05 Hz (the caloric-equivalent regime) they diverge sharply — e.g. at 0.01
Hz, peripheral-only gain is 0.336 vs. central (tau=17.5s) 0.686 vs. central-high (tau=30s) even higher.
This divergence, not a single number, IS the physics: which time constant governs is exactly what
distinguishes an intact-central-integrator subject from a Baloh-1988-type lesion.

### Honest gaps
- Bouveresse 1998's exact per-frequency numbers were not extractable live (paywalled abstract) — the
  0.9-1.0 mid-band falsifier is checked against this script's own closed-form model, corroborated (not
  duplicated) by McGarvie 2015's plateau finding and Bouveresse's frequency-range/linearity confirmation.
- vHIT (McGarvie 2015) anchors the gain PLATEAU only, not the frequency-dependent SHAPE — a different
  regime (amplitude/velocity nonlinearity) than this linear model's frequency response.
- The Baloh 1988 reproduction is qualitative/structural (direction + corner-effect), not a quantitative
  refit — real disease likely also carries a central GAIN deficit this single-tau model doesn't represent.
- Central gain-plateau calibration (the adaptive neural gain-control loop) is not mechanistically modeled,
  only the frequency-dependent shape (the high-pass corner) is.

## Part 2 — Postural sway: inverted pendulum + delayed feedback

### Method, in one paragraph

Quiet standing = a single-link inverted pendulum pivoting at the ankle (point mass m at COM height h,
I=m·h²). Gravity supplies a destabilizing torque mgh·θ (small angle). Newton's law about the ankle:
`I·θ'' = mgh·θ - T_A(t) + ξ(t)`, where the corrective ankle torque comes from a **delayed** feedback law
on a noisy position estimate — `T_A(t) = Kp·θ_hat(t-Δ) + Kd·θ'(t-Δ)`, `θ_hat=θ+η` — because Loram & Lakie
(2002, PMID 12482906) directly measured passive ankle stiffness at only **91±23% of the critical value**
(insufficient alone; active modulation required), and Peterka (2002, PMID 12205132) found feedback
**delay**, not the damping coefficient, is the primary mechanism producing effective damping. Center of
pressure follows directly from taking moments about the ankle: `x_COP(t) = T_A(t)/(mg)`. Mass (78.2 kg)
is read live from this twin's own `metabolic_cost_results.json` (subject2) for cross-cell consistency;
COM height uses Winter's ~0.56×height anthropometric fraction (disclosed population estimate, subject2's
own stature not found recorded).

### A real bug the machine crosscheck caught — not swept under the rug

The first attempt (damping ratio ζ=0.35, sized for the undelayed system) **failed** the exact-vs-Monte-
Carlo crosscheck outright: the Monte Carlo simulation exploded to θ~1e21 rad (z-score 5.8-5.9σ against
the frequency-domain prediction). Diagnosis (Orient step, not skipped): a numerical right-half-plane
root search of the characteristic equation `D(s)=I·s²-mgh+(Kp+Kd·s)e^(-sΔ)` (grid scan + Nelder-Mead
polish, `max_real_part_root()` in the script) found a genuine unstable root at **Re(s)=+0.44** once the
150ms delay was included — the parameter combination is **actually unstable**, not a simulation artifact;
the frequency-domain PSD-integral formula had been silently returning a finite-but-meaningless number
because that formula doesn't itself check for closed-loop stability. Fix: raise damping to **ζ=1.0**
(Kd=180.9 N·m·s/rad), numerically re-verified stable (margin **-0.452** at 150ms delay). After the fix,
the crosscheck **passes**: z=0.21 (θ), z=1.40 (θ̇), both within 4σ. Every variance this script reports is
now gated by this explicit stability check — an unstable configuration is flagged and excluded, not
silently reported.

### External falsifier anchor: this twin's OWN already-measured `ORG-POSTURAL-BALANCE` (not fresh literature)

`ORG-POSTURAL-BALANCE` (`data/MECHANISM_ANCHOR_GRAPH.json`, `bt_memory/LEDGER.jsonl`, executed 2026-07-18,
MEASURED-B PASS) analyzed real HBEDB force-plate data and found: **Romberg ratio (EC/EO) = 1.155 (firm) /
1.304 (foam)** on COP velocity/path-length, but **NO significant effect on RMS/area** (p=0.58/0.22, an
honest null) — confirmed as a genuine PhysioNet dataset (163 subjects, EO/EC×firm/foam, 100Hz, 10Hz
low-pass) via a direct live fetch of `physionet.org/content/hbedb/` this session, matching the internal
claim exactly (not fabricated). **This metric-specific pattern — velocity carries the signal, RMS/area
does not — is the real target**, not just "sway increases eyes-closed."

### Forced adversary #1 (OODA): does naive noise-amplitude scaling explain the pattern? No — refuted analytically AND numerically

The obvious first guess — eyes-closed just increases sensory-noise variance, same shape, bigger amplitude
— is refuted by a general LTI/spectral fact, stated before running anything: scaling a single noise
source's variance (shape unchanged) scales **every** linear output's variance, position and velocity
alike, by the **identical** factor. Numerically confirmed to 4 significant figures: for k=1.5, 2.0, 3.0×
noise, `RMS_ratio` **exactly equals** `velocity_ratio` (dissociation_ratio = **1.000** in all three cases).
This mechanism **cannot**, even in principle, reproduce the twin's own measured dissociation. Clean,
decisive, forced negative — not a near-miss.

### System-shape mechanisms tested (Peterka 2002 motivates both: delay and stiffness are what actually change with sensory reweighting)

| mechanism | RMS ratio | velocity ratio | dissociation ratio | reproduces target pattern? |
|---|---:|---:|---:|---|
| delay +10/20/30 ms | 1.10 / 1.25 / 1.47 | 1.11 / 1.26 / 1.50 | ~1.05-1.06 | No |
| stiffness ×0.97/0.94/0.90 | 1.06 / 1.18 / 1.64 | 0.93 / 0.87 / 0.79 | **-1.17 to -0.32** (wrong direction) | No |
| combined (delay+stiffness+noise) | 1.32 / 1.64 | 1.16 / 1.31 | 0.48-0.52 | No |

**Honest negative on the precise quantitative target** (dissociation_ratio > 2.0, i.e. velocity-ratio
excess at least double the RMS-ratio excess): none of the tested single mechanisms clear this
pre-registered bar within the safely-stable parameter region. This is reported plainly, not hidden behind
the headline PASSes above.

### The OODA loop is not abandoned here — a genuine, already-computed structural clue

Rather than stopping at the negative, the already-computed void-floor stability-margin sweep (no
additional simulation — pure arithmetic on numbers already on disk) was checked for a second signature:
does the **velocity/position variance ratio itself** shift as the delay-only system approaches its
stability margin? It does, monotonically: **1.65 (50ms) → 1.81 (100ms) → 1.94 (150ms) → 1.99 (175ms) →
2.03 (200ms, at the edge of the stability boundary)**. Operating closer to the margin measurably shifts
spectral content toward velocity relative to position — the right *direction* for the target dissociation,
though not yet a full quantitative reproduction of it. **Left explicitly OPEN**: the natural next step is
a mechanism that combines the EO→EC condition change with near-margin operation (physiologically
plausible — real postural control is known to operate close to, not far from, its stability limit), not
claimed as resolved here.

### Void-floor / non-degeneracy: the stability margin genuinely governs the response

| delay | stability margin | var(θ) | note |
|---:|---:|---:|---|
| 50 ms | -1.061 | 2.41e-3 | far from boundary |
| 100 ms | -0.805 | 3.53e-3 | |
| 150 ms | -0.452 | 6.67e-3 | (baseline EO) |
| 175 ms | -0.251 | 1.21e-2 | |
| 200 ms | -0.045 | 6.64e-2 | **27.6× the 50ms value** — right at the edge |
| 210-300 ms | +0.036 to +0.649 | **UNSTABLE — excluded** | crosses the margin, as it must |

A real, monotonic, machine-verified margin-governs-amplitude structure (the same structural role as
σ_min in the reduced-representation work this project's other lanes verified, applied here to a genuinely
different geometry — a linear-control stability margin, not a matrix singular value; an analogy, not an
identity, stated as such). The instability onset (~205ms) sits just above this script's disclosed 100-200ms
physiological delay estimate — i.e. real postural gains plausibly operate close to, but on the stable
side of, this exact margin, consistent with Loram & Lakie's and Peterka's emphasis on precisely-tuned
active control rather than a comfortably overdamped passive system.

### Independent plausibility check (non-tautological, not fitted to match)

The model's own predicted sway spectrum peaks at **0.20 Hz** with **99.9% of power below 1 Hz** — matching
the well-documented real-world fact that quiet-stance COP sway concentrates below ~1-2 Hz, a check this
script's parameters were not tuned to hit.

### Honest gaps
- Absolute noise scale is arbitrary; only EO/EC **ratios** are falsifier-relevant (not claimed otherwise).
- `height_m=1.73` is a disclosed population estimate, not subject2's own measured stature.
- `delay=150ms` is a disclosed round estimate (100-200ms literature range); a live search this session
  for one specific pinned PMID for the exact millisecond figure did not converge (two eutils queries,
  zero hits) — disclosed as an honest gap, not silently asserted as literature-extracted.
- `zeta=1.0` (Kd) is **not** a free aesthetic choice — it is the minimum damping this script's own
  root-search found *necessary* for stability at 150ms delay; ζ=0.35 was tried first and is genuinely
  unstable there (see the crosscheck-failure story above).
- The 10Hz noise band-limit is matched to HBEDB's own stated low-pass filter (real, dataset-derived) but
  is fundamentally a *regularization* need (true white noise gives a divergent velocity-variance integral
  through a proportional-gain path) — disclosed as such, not hidden as a physiological measurement.
- 1-D (single-axis) pendulum; real HBEDB COP is 2-D (AP+ML) — the EC/EO ratio's axis-invariance is
  argued, not independently verified against the real dataset's own AP/ML breakdown.
- `max_real_part_root()` is a practical grid+Nelder-Mead search, not a rigorous Nyquist/argument-principle
  stability proof — adequate to catch the gross instability found this session and gate every reported
  number, but a smaller residual root elsewhere in the complex plane cannot be fully excluded.

## Symmetric QC — what is and is not claimed

**Claimed (PASS, machine-verified, externally anchored)**:
1. The VOR canal-dynamics transfer function reproduces measured mid-band gain (0.96-1.0 at 0.1-1 Hz vs.
   measured ~0.9-1.0), is non-degenerate (90-fold gain range across frequency, 5.3× across tau), and its
   low/high-frequency structure is qualitatively corroborated by an independent disease population.
2. The postural inverted-pendulum-with-delay model is internally consistent (two independent
   computational paths agree to <1.4σ after a real instability bug was found and fixed) and exhibits a
   genuine, monotonic stability-margin-governs-sway-amplitude structure, plus an independently-plausible
   (not fitted) low-frequency-dominated sway spectrum.

**NOT claimed (explicitly OPEN, honest negative)**:
3. A full quantitative mechanistic reproduction of the twin's own measured metric-specific Romberg
   dissociation (velocity ratio 1.155-1.304 significant, RMS/area null) — the naive mechanism is
   analytically refuted, three system-shape mechanisms tested fall short of the pre-registered
   dissociation-ratio threshold, and a directionally-suggestive (not closed) structural clue is disclosed
   as the next OODA iteration.

`couples_to`: `ORG-POSTURAL-BALANCE` (real-data anchor, unedited), `SNS-VESTIBULAR-BALANCE` /
`SENS-VESTIBULAR-BALANCE` (clinical/epidemiological siblings, different scope, confirmed non-overlapping
by direct read), the MSK postural/spine-gait cells (`docs/MECHANISM_SPINE_GAIT_VBR.md`,
`docs/MECHANISM_MZ_PELVIS_MOMENT.md`) via the shared ankle-torque/COP mechanical framing, and this twin's
own `metabolic_cost_results.json` (subject2 mass, reused for cross-cell consistency).

## Repro

```
cd ~/projects/bodytwin
source .venv-msk/bin/activate
python scripts/msk/vestibular_vor.py           # writes reports/probes/vestibular_vor_results.json
python scripts/msk/postural_sway_pendulum.py   # writes reports/probes/postural_sway_pendulum_results.json
```
Both are pure numpy/scipy (no OpenSim, no external data download). The postural script reads
`data/msk_smoketest/subject2_walking1/metabolic_cost/metabolic_cost_results.json` read-only if present
(falls back to a disclosed default mass otherwise). Total runtime: VOR script <2s; postural script ~1-2
minutes (dominated by the Monte Carlo crosscheck + root-search stability gates across ~20 parameter
points). No git operations; no writes outside `reports/probes/` and this doc pair.
