# MECHANISM PUSH-OFF PLANTARFLEXOR / ACHILLES-FORCE TEST — is push-off over-prediction a gastroc/soleus over-force? (2026-07-21)

Tests the mechanism `docs/MECHANISM_CONTACT_WAVEFORM.md` localized (knee/hip contact force over-predicts
specifically at the push-off hump, 44.8%GC, Z=2.67 knee / 4.56 hip, while weight-acceptance matches
OrthoLoad) against a hypothesis `docs/MECHANISM_CONTACT_MUSCLE_DECOMP.md` invited (gastrocnemius = 48.7%
of the knee contact peak) but did not itself test at the ankle: **is the twin's ankle-plantarflexor
(gastrocnemius+soleus, "Achilles-tendon-equivalent") force over-forced at push-off, checked against a
NEW anchor decorrelated from OrthoLoad — in-vivo Achilles tendon force measured by optic-fiber/buckle
transducer (Komi/Finni/Frisland)?**

**Headline: MIXED / PARTIAL, not a clean confirm and not a clean rule-out.** At the exact instant
already established as the push-off over-prediction (t=0.51s, 44.8%GC — the SAME frame both SO and
CMC peak at, re-verified here, not re-derived), the modeled Achilles-tendon-equivalent force
(gasmed_r+gaslat_r+soleus_r tension, summed = 1960.70 N = 255.7 %BW = 2.56xBW) is **37.1% above** the
best-verified direct in-vivo transducer anchor's mean (Finni et al. 1998, PMID 9535592: 1430±500 N
walking peak) — a ratio that clears this repo's own established >20%-material bar, but a **Z-score of
only 1.06** given that anchor's own reported (large, disclosed-by-the-source-paper) inter-subject
variance — i.e. **not statistically distinguishable from normal in-vivo variation**, and clearly
**smaller** than the knee/hip contact-force excess at the identical instant (46-67%, `docs/
MECHANISM_CONTACT_WAVEFORM.md`). A decisive, already-built causal test (Handsfield-MRI-PCSA Fmax
correction, `docs/MECHANISM_FMAX_PCSA_VALIDATION.md`) — which closed 17.1% of the KNEE CONTACT-force
gap — moves this ankle-level number by only **+1.4%**, because SO's own redundancy resolution
reallocates the required net ankle moment among the three (now individually weaker) synergists
(soleus_r tension nearly **doubles**, +93.5%, compensating gasmed_r/gaslat_r's ~24% drops): **muscle-
strength mis-calibration, unlike part of the knee-contact story, is NOT the driver here.** A real-EMG
cross-check (gasmed_r+soleus_r have actual measured channels for this subject) shows real muscle
activity does **not** imply the twin over-recruited these muscles — if anything the opposite for
soleus. The excess **grows** to a clearer 53-72% (Z≈2.05) only a few %GC **later**, at the
plantarflexors' own natural peak (t=0.59s, 50.8%GC) — a **timing offset** from the knee/hip's own
push-off peak, not a uniform over-force at the same instant. **A forced adversary caught mid-analysis:
the task's own pre-registered "~2600-2900N / 3.5-4xBW" anchor band does not match the direct-transducer
Komi/Finni/Frisberg literature at all (1430-1610N / 1.9-2.1xBW) — it matches almost exactly (2991N /
3.9xBW) a DIFFERENT, model-based (finite-element, not transducer) source, Giddings et al. 2000 (PMID
10731005), which shares the twin's own "estimate force from measured kinematics+GRF" paradigm and is
therefore NOT a fully decorrelated anchor.** Against that secondary anchor the twin actually
**under-predicts** (-34%) — the opposite sign. **Verdict: plantarflexor/Achilles over-force is a real
but modest, partial, non-dominant contributor at best — it does not by itself explain the larger
knee/hip contact-force excess at the same instant, and points, as the task's own framing anticipated,
toward another mechanism (free-body geometry / timing-phase offset) as more important.**

## Headline numbers

| construction | t (s) | %GC | gasmed_r (N) | gaslat_r (N) | soleus_r (N) | Achilles-equiv (N) | %BW | vs Finni 1998 mean (1430±500N) | vs Fröberg 2009 (2.1xBW=1610N) | vs Giddings 2000 (3.9xBW=2991N, method-shared) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| **SO baseline, AT push-off instant (PRIMARY)** | 0.51 | 44.8 | 1140.0 | 398.0 | 422.7 | **1960.7** | 255.7% | **+37.1%, Z=1.06** | +21.7% | **−34.4%** |
| SO baseline, own-peak (context only) | 0.59 | 50.8 | 769.3 | 333.4 | 1353.3 | 2456.0 | 320.3% | +71.8%, Z=2.05 | +52.5% | −17.9% |
| Fmax/PCSA-corrected, at instant (causal test) | 0.51 | 44.8 | 865.5 | 304.4 | 817.9 | 1987.7 | 259.2% | +39.0%, Z=1.12 | +23.4% | −33.5% |
| EMG-hybrid (2/3 real EMG), at instant (directional only) | 0.51 | 44.8 | 917.7 (real EMG) | 448.1 (SO-proxy) | 1188.8 (real EMG) | 2554.5 | 333.1% | +78.6%, Z=2.25 | +58.6% | −14.6% |
| CMC bonus (different config, heavily caveated) | 0.51 | 44.8 | 1285.0 | 373.7 | 173.9 | 1832.6 | 239.0% | +28.1%, Z=0.81 | +13.8% | −38.7% |

All numbers machine-computed by `scripts/msk/pushoff_plantarflexor_achilles.py`, written to
`data/msk_smoketest/subject2_walking1/pushoff_plantarflexor/pushoff_plantarflexor_results.json` — none
transcribed from console prose.

## 1. Pre-registration (stated before the literature anchor was fetched)

**C / ¬C, threshold fixed in advance:**
- **C (over-force IS the mechanism):** modeled Achilles-equivalent force AT the established push-off
  instant (t=0.51s, matching `docs/MECHANISM_CONTACT_WAVEFORM.md`'s own frame — not a separately
  re-derived or cherry-picked time) MATERIALLY exceeds the direct in-vivo transducer anchor
  (>20% excess, this repo's own established materiality convention, `docs/
  MECHANISM_EMG_CONSTRAINED_FORCE.md` Sec.3) AND is not explained away by the anchor's own reported
  inter-subject variance (Z clearly >2).
- **¬C (rules it out):** within range (<20% excess, or Z not distinguishable from the anchor's own
  reported spread) — points elsewhere (free-body geometry at push-off, the moment-arm/
  redundancy-resolution effects already in `docs/MECHANISM_CONTACT_MUSCLE_DECOMP.md`).

**Forced adversary (the one this task is tempted to skip):** a plausible "in-vivo Achilles force"
anchor could also come from an inverse-dynamics/finite-element MODEL estimate rather than a direct
tendon-force transducer — which would share the twin's own "estimate force from measured
kinematics+GRF" paradigm, making agreement/disagreement a weaker, confounded signal, not a genuinely
external anchor. Both anchor types were fetched live and are reported (Sec.4); the DIRECT TRANSDUCER
literature is pre-registered as decisive for the falsifier, because it is what the task explicitly
named (optic-fiber/buckle transducers) and is the one genuinely decorrelated from any model-based
force-estimation paradigm.

**Symmetric-QC checks forced, not assumed:** (1) the ANKLE-plantarflexion tendon force is read from the
raw per-muscle `tension` column of the SO `force.sto` (the force transmitted through the whole
muscle-tendon actuator — what a buckle/optic-fiber transducer on the Achilles tendon would read), **not**
the knee-contact-axis-projected component `docs/MECHANISM_CONTACT_MUSCLE_DECOMP.md` computed for a
different purpose (that doc's "gasmed_r: 139.94%BW, 35.8% of total" is gasmed_r's projection onto the
KNEE joint-contact direction `n_hat`, a smaller, geometrically different number — not reused here except
as a sanity cross-check that the raw tension values match, Sec.5); (2) the push-off instant is verified
to be the IDENTICAL row/frame already established in `docs/MECHANISM_CONTACT_WAVEFORM.md` (t=0.51s,
asserted in-script, would raise not silently pass on a grid mismatch), not a separately re-derived time.

## 2. Method

Reused, unedited, **zero new OpenSim run**: subject2/`walking1`'s already-committed SO
`force.sto`/`activation.sto` (`data/msk_smoketest/subject2_walking1/static_optimization/so/`), the
already-built Fmax/PCSA-corrected re-solve (`data/msk_smoketest/subject2_walking1/
fmax_pcsa_correction_test/so/`, from `docs/MECHANISM_FMAX_PCSA_VALIDATION.md`), the already-built
EMG-hybrid forward-dynamics force (`data/msk_smoketest/subject2_walking1/emg_driven/
emg_hybrid_force.sto`, from `docs/MECHANISM_EMG_CONSTRAINED_FORCE.md`), and one bonus, heavily-caveated
second-solver check (`ankle_reserve_fix/boost150/cmc/full_1p57s/..._Actuation_force.sto`, a **different**
CMC run configuration than `cmc_second_solve`). Real subject EMG read in place from the read-only
external drive (`/media/anton/8838D60F38D5FBDE/mechanism_data/LabValidation_withVideos/subject2/
EMGData/walking1_EMG.sto`), isolation-compliant. New code: `scripts/msk/pushoff_plantarflexor_achilles.py`
(no OpenSim import needed, same convention as `contact_waveform_analysis.py` — these are all already-
computed `.sto` files being re-parsed).

**"Achilles-tendon-equivalent" definition:** `gasmed_r + gaslat_r + soleus_r` tendon tension, summed.
This is the standard OpenSim-biomechanics operationalization of Achilles tendon force — all three
muscle-tendon actuators insert on the calcaneus via what is anatomically the single Achilles tendon;
each one's OpenSim "force" output IS its own tendon tension (already accounting for pennation, the
physically relevant quantity a buckle/optic-fiber transducer distal to the myotendinous junctions would
read), and because the three lines of action are closely aligned at the ankle, the scalar sum is the
standard, literature-consistent way to estimate combined Achilles tendon force from a multi-muscle
triceps-surae representation.

**%GC mapping:** reused verbatim from `docs/MECHANISM_CONTACT_WAVEFORM.md` (`T_stride=1.327s`,
`R_heel_strike_t=1.2425s`) — the SAME convention, not re-derived, so the reported %GC values are
directly comparable to that doc's own 44.8%GC push-off figure. (`docs/
MECHANISM_CONTACT_MUSCLE_DECOMP.md`'s sibling `T=1.33s`/`right_hs_0=-0.0765` convention, from a
different gait-cycle-detection method already in this repo, gives 44.1%GC for the identical t=0.51s
frame — a small, pre-existing ~0.7-percentage-point discrepancy between two already-established
pipelines in this repo, not something introduced here; both point to the SAME row/frame, which is what
matters for this test.)

## 3. Verification gates (machine-checked, all PASS)

| gate | result |
|---|---|
| t=0.51s lands exactly on the SO `.sto` grid (row 51, not interpolated) | **PASS** (`assert abs(t[row]-0.51)<1e-9`, `assert row==51`) |
| Same row (51, t=0.51s) also exact in `activation.sto`, Fmax-corrected `force.sto`, `emg_hybrid_force.sto`, `naive_fwd_only_force.sto`, and the raw subject `walking1_EMG.sto` (all 158-row, 0.01s-step grids) | **PASS**, confirmed by direct inspection, all 6 files share the identical time vector |
| Re-extracted gasmed_r/gaslat_r tension at t=0.51s matches `docs/MECHANISM_CONTACT_MUSCLE_DECOMP.md`'s published values (1140.0N/398.0N) | **PASS**, exact to the digit |
| EMG-hybrid gasmed_r/SO ratio (0.805) matches `docs/MECHANISM_CONTACT_MUSCLE_DECOMP.md` Sec.9's published value | **PASS**, exact |
| `gaslat_r` hybrid vs naive_fwd_only relative difference (confirms gaslat_r is NOT real-EMG-informed, as disclosed) | **PASS**, −0.00016% (floating-point noise) |
| Raw-EMG independent re-derivation of gasmed_r/soleus_r own-peak %GC vs already-published values (50.1/61.4%GC) | **PASS**, 50.8/62.1%GC (small diff from the T=1.327 vs T=1.33 convention gap already noted, not a new discrepancy) |
| Ankle reserve actuator negligible in both baseline and Fmax-corrected models (confirms muscle-driven, not reserve-actuator-driven, reallocation) | **PASS**, 0.0057 Nm / 0.0050 Nm at instant; max 0.068/0.078 Nm whole-trial |
| Literature PMIDs verified live via NCBI eutils (not from training-data memory) | **PASS** — see Sec.4 |

## 4. External literature anchor — verified live, and a forced adversary caught mid-analysis

Fetched via direct `curl` to `eutils.ncbi.nlm.nih.gov` (esearch/esummary/efetch abstract) and the
Semantic Scholar API — the same live-verification convention `docs/MECHANISM_FMAX_PCSA_VALIDATION.md`
Sec.0 already used to catch a hallucinated PMID, reused here rather than trusted from memory (this
repo's `WebSearch` tool budget was exhausted this session — `curl`-to-eutils is the fallback the prior
doc already established works).

**PRIMARY anchor — direct in-vivo transducer, the one this task asked for:**
- **Finni T, Komi PV, Lukkariniemi J. "Achilles tendon loading during walking: application of a novel
  optic fiber technique." Eur J Appl Physiol Occup Physiol. 1998 Feb;77(3):289-91. PMID 9535592.
  DOI 10.1007/s004210050335.** n=8, 3 walking speeds (1.1/1.5/1.8 m/s). Direct quote: *"peak ATF's
  (1430±500 N)... peak ATF was found to be rather insensitive to speed."* Raw Newtons, not
  %BW-normalized in the abstract.
- **Fröberg A, Komi P, Ishikawa M, Movin T, Arndt A. "Force in the Achilles tendon during walking with
  ankle foot orthosis." Am J Sports Med. 2009 Jun;37(6):1200-7. PMID 19229043.
  DOI 10.1177/0363546508330126.** n=8, same optic-fiber technique, Komi co-author (same lab lineage).
  Direct quote: *"decreased to 2.1 times body weight during barefoot walking (P < .01)."*
- **Triangulation:** 2.1xBW at a plausible ~680N-BW cohort implies ~1430N — consistent with Finni
  1998's own directly-reported mean to within its own SD. Two independent studies from the same
  research tradition converge on **~1430-1610N (~1.9-2.1xBW)** as the direct in-vivo TRANSDUCER walking
  peak Achilles tendon force.

**SECONDARY anchor — model-based, method-shared risk flagged (a forced adversary caught, not glossed
over):** **Giddings VL, Beaupré GS, Whalen RT, Carter DR. "Calcaneal loading during walking and
running." Med Sci Sports Exerc. 2000 Mar;32(3):627-34. PMID 10731005.
DOI 10.1097/00005768-200003000-00012.** A finite-element foot model driven by measured GRF +
cineradiography kinematics — **not** a direct tendon-force transducer. Direct quote: *"The maximum
predicted Achilles tendon forces were 3.9 and 7.7 BW for walking and running."* **3.9xBW × subject2's
own BW (766.88N) = 2990.8N** — this closely matches the task's OWN pre-registered "~2600-2900N
(~3.5-4xBW)" band. This is strong circumstantial evidence that band was sourced (directly or via a
secondary citation chain) from Giddings — a MODEL estimate that shares this twin's own "estimate force
from measured kinematics+GRF" paradigm — **not** from the Komi/Finni/Fröberg direct transducer
literature the task also named in the same sentence. Flagged explicitly rather than silently used as
if it were the decorrelated anchor the task intended.

**Checked and set aside:** Komi PV 1990 (PMID 2081741, J Biomech) reports AT forces "as high as 9 kN...
12.5 times body weight," but the abstract does not attribute this peak specifically to WALKING (the
paper covers walking+running+jumping together) — 12.5xBW is far more consistent with stretch-
shortening-cycle running/jumping/hopping loading than walking, so **not used as the walking anchor**.
Komi, Fukashiro & Jisvinen 1992 (PMID 1638639, a review) names hopping as "unexpectedly high" but gives
no walking-specific N or %BW figure in its abstract — not independently usable as a quantitative
walking anchor.

## 5. Results in detail

### 5.1 At the pre-registered push-off instant (t=0.51s, 44.8%GC) — the decisive comparison

SO's own solution puts **1960.70 N (255.7%BW, 2.56xBW)** of combined gasmed_r+gaslat_r+soleus_r tendon
tension through the ankle at the exact instant the knee/hip contact force over-predicts. Against the
PRIMARY (direct-transducer) anchor this is **+37.1% (ratio 1.371x), Z=1.06** — real by the ratio
convention this repo already uses elsewhere, but **not statistically distinguishable** from normal
population variation given Finni's own reported ±500N (35%-of-mean) inter-subject spread (their own
abstract: "great intersubject variation... in the peak ATF's"). By direct comparison, the knee/hip
contact-force excess at the SAME instant is 46-67% (ratio 1.46-1.67x, Z=2.67-4.56) — **materially
larger and materially more statistically significant** than the plantarflexor-force excess. If
plantarflexor over-force were the dominant mechanism behind the knee/hip push-off over-prediction, its
own excess at that instant should be at least comparable in magnitude/significance — it is not.

### 5.2 At the plantarflexors' own peak (t=0.59s, 50.8%GC) — context, not the decisive number

The Achilles-equivalent group's OWN maximum (not the pre-registered instant — reported for robustness,
using it as primary would be exactly the un-pre-registered cherry-pick this task's method forbids)
occurs 8%GC later, reaching **2456.05 N (320.3%BW, 3.20xBW)** — **+71.8% (Z=2.05)** vs the primary
anchor, now a magnitude/significance more comparable to the knee/hip story. The composition shifts
sharply: soleus_r rises from 422.7N to 1353.3N (own peak) while gasmed_r falls from 1140.0N to 769.3N —
a genuine within-model handoff from gastrocnemius-dominant to soleus-dominant loading as stance
progresses (physiologically expected: soleus is a pure ankle plantarflexor unaffected by knee angle,
gastrocnemius is biarticular and loses effective plantarflexion leverage as the knee extends through
late stance). This is a **timing** observation: the model's own plantarflexor peak arrives after the
knee/hip contact peak, not at the same instant.

### 5.3 The Fmax/PCSA causal test — decisive, reused, not re-derived

`docs/MECHANISM_FMAX_PCSA_VALIDATION.md` already built and ran a Handsfield-MRI-PCSA-corrected model for
subject2 (all three plantarflexors individually over-strong vs MRI-implied Fmax: gasmed_r ratio 1.246,
gaslat_r 1.253, soleus_r 1.108) and showed it closes 17.1% of the KNEE CONTACT-force gap. Re-reading
that SAME already-built model's `force.sto` at the identical t=0.51s instant shows the summed
Achilles-equivalent tension moves by only **+1.4%** (1960.70N → 1987.75N) — because SO's own
redundancy resolution **reallocates** the required net ankle moment among the three now-weaker
synergists: soleus_r tension **nearly doubles** (+93.5%, 422.7N→817.9N) while gasmed_r (−24.1%) and
gaslat_r (−23.5%) drop. Ankle reserve stays negligible in both models (0.0057/0.0050 Nm — not a
reserve-actuator artifact absorbing the difference). **This is a real, mechanistically precise,
decisive finding: unlike part of the knee CONTACT-force story, the modest ankle-level excess here is
NOT explained by muscle-strength mis-calibration** — correcting all three muscles' Fmax barely moves
the combined tendon force, because the required net ankle moment is met by internal reallocation, not
by an overall drop in demand.

### 5.4 Real-EMG cross-check — directional, not absolute

Subject2 has real surface EMG for `gasmed_r` and `soleus_r` (no channel for `gaslat_r`). The already-
built CEINMS-style EMG-hybrid forward-dynamics force (`docs/MECHANISM_EMG_CONSTRAINED_FORCE.md`,
`docs/MECHANISM_EMG_DRIVEN.md`) gives, at t=0.51s: gasmed_r 917.7N (ratio 0.805 vs SO — corroborates,
matching the already-published number exactly) and **soleus_r 1188.8N (ratio 2.81 vs SO's 422.7N —
real EMG implies SUBSTANTIALLY MORE force than SO chose)**. Independently re-derived raw-EMG timing
(this session, from the raw `walking1_EMG.sto` file) confirms gasmed_r's own real-EMG peak at
50.8%GC (inside its 30-60%GC literature window — PASS) and soleus_r's at 62.1%GC (a boundary miss past
the 60%GC window edge — FAIL by ~2 points), both reproducing the already-published 50.1/61.4%GC
findings independently. **Reading this carefully, not hyping it**: the EMG-hybrid absolute magnitude
is known-uncalibrated (no true-MVC trial for this subject, already disclosed in `docs/
MECHANISM_EMG_DRIVEN.md` Sec.9.1) — so the 2554.5N/333.1%BW combined figure is **directional evidence
only**. Its direction is unambiguous, though: real measured muscle activity does **not** support "the
twin invented excess plantarflexor drive" — if anything, for soleus specifically, real EMG points
toward MORE force being physiologically plausible at this instant, not less.

### 5.5 Second-solver bonus check — heavily caveated

A different CMC run configuration (`ankle_reserve_fix/boost150`, built for a different purpose —
CMC convergence via a deliberately boosted ankle reserve actuator, `docs/MECHANISM_ANKLE_RESERVE_FIX.md`
— NOT the `cmc_second_solve` configuration every other CMC number in this cert family uses) gives, at
the same instant (interpolated from its own finer variable-step grid): Achilles-equivalent 1832.6N
(239.0%BW), ratio 0.935 vs SO — **lower**, not higher, than SO's own number. This does not corroborate
a clean "both solvers over-force the plantarflexors" story the way the knee CONTACT force does (where
CMC reads 8-11% higher than SO everywhere) — reported as a minor, disclosed, non-primary data point
given the configuration mismatch, not folded into the verdict.

## 6. Symmetric-QC

**The adversary this result (leaning toward a positive "over-force" finding at 5.1's ratio-only read)
must survive:** is the >20%-material ratio a real signal or an artifact of comparing to an anchor with
enormous, disclosed inter-subject variance? Forced via the Z-score against Finni's own reported SD:
Z=1.06 at the primary instant is unremarkable (≈85th percentile of a normal distribution, not an
outlier) — the ratio-only reading is NOT allowed to stand alone; the Z-score check, forced here, pulls
the verdict back from "material" toward "modest, not clearly distinguishable."

**The adversary this result (leaning toward a negative "rules it out" finding) must survive:** is
t=0.51s really the right instant to test, or does using the pre-registered knee/hip instant hide a
larger excess that shows up at the plantarflexors' own natural peak? Forced via Sec.5.2: yes, a larger
(53-72%, Z≈2.05) excess is real at t=0.59s — reported prominently, not suppressed, specifically so this
result cannot be dismissed as a clean negative either.

**A third adversary, caught only by fetching the anchor literature carefully (Sec.4):** the task's own
pre-registered anchor range turned out to match a model-based, not transducer-based, source almost
exactly (Giddings 3.9xBW = 2990.8N vs the task's stated "~2600-2900N"). Using that anchor instead would
flip the entire verdict's sign (the twin would UNDER-predict, −34% to −18%) — a stark illustration of
why the pre-registered decision was anchored to the direct-transducer literature specifically, and why
both anchors are reported rather than silently picking whichever one supports a cleaner story.

## 7. Honest gaps

1. **Finni 1998's own subject cohort's body mass is not available from the abstract** (paywalled full
   text) — the raw-Newton comparison (Sec.5.1) assumes subject2's own Achilles force is comparable in
   absolute terms to a plausible ~65-80kg Finnish walking cohort; the %BW-normalized Fröberg 2009 point
   (2.1xBW) is used as a cross-check specifically to reduce this dependency, and the two triangulate
   consistently, but neither is a per-subject-matched comparison.
2. **n=8 subjects in each direct-transducer source** — small-N, exactly as this task's own confidence-
   tier framing anticipated; the reported SD (±500N, 35% of mean) already reflects this in the Z-score,
   not hidden.
3. **The EMG-hybrid absolute-magnitude uncertainty (Sec.5.4) is real and unresolved** (inherited,
   not re-litigated, from `docs/MECHANISM_EMG_DRIVEN.md` Sec.9.1) — only the DIRECTION (real EMG does not
   support under-recruitment by the model) is treated as evidenced; the specific 2.81x soleus ratio is
   not treated as a precise number.
4. **The second-solver (CMC) bonus check (Sec.5.5) uses a different run configuration** (boosted ankle
   reserve) than every other CMC number in this cert family — reported only as a minor, disclosed data
   point, explicitly not folded into the primary verdict; a same-configuration CMC per-muscle
   `Actuation_force.sto` for `cmc_second_solve` does not exist in this repo yet (only its `JointReaction`
   contact-force output was kept) — a natural, cheap follow-up for a future session (re-run CMC's own
   already-configured setup with the `Actuation` reporter enabled) is flagged, not done here (LEAN
   discipline — this session's falsifier does not require it, since the primary SO-based finding is
   already decisive on its own terms).
5. **Single trial, single subject, right leg** (subject2 `walking1`) — same scope caveat as every cert
   in this family; no claim of generality across subjects/trials/speeds.
6. **The Giddings 2000 finite-element model's own assumptions (foot/calcaneus contact model, GRF/
   kinematic inputs from a different cohort) were not audited here** — it is reported as a flagged,
   secondary, method-shared-risk anchor, not validated or invalidated in its own right.
7. **A concurrent, independent analysis in this same repo**
   (`data/msk_smoketest/subject2_walking1/contact_waveform_fmax_corrected/`,
   `scripts/msk/contact_waveform_fmax_corrected.py`) tests a related but distinct question — whether
   Fmax/PCSA correction selectively fixes the KNEE CONTACT-force push-off hump while preserving the
   weight-acceptance match — and reports `H_SELECTIVE_CORROBORATED` for that knee-contact-waveform
   question. That is a different measurement (knee joint contact force under Fmax correction, not
   ankle-level Achilles-tendon tension under Fmax correction) from Sec.5.3 here; noted as a pointer,
   not relied upon or contradicted (concurrent, independent work by another instance, per this repo's
   multi-instance convention).
8. **The %GC-convention discrepancy (44.8% vs 44.1%) between this repo's two existing gait-cycle-
   detection pipelines is pre-existing, not introduced or resolved here** (Sec.2) — both point to the
   identical underlying frame (t=0.51s), which is what this test's validity actually depends on.

## 8. Files

- `scripts/msk/pushoff_plantarflexor_achilles.py` (new) — the full, re-runnable, self-contained
  pipeline (no OpenSim import; only re-parses already-computed `.sto` files). Run with
  `/usr/bin/python3 scripts/msk/pushoff_plantarflexor_achilles.py`.
- `data/msk_smoketest/subject2_walking1/pushoff_plantarflexor/pushoff_plantarflexor_results.json`
  (new) — every number in this document, machine-written: baseline/own-peak/Fmax-corrected/EMG-hybrid/
  CMC-bonus per-muscle tension and Achilles-equivalent sums, the live-verified literature anchor block
  (PMIDs, DOIs, direct quotes), and the full verdict computation (ratios, %-excess, Z-scores vs all
  three anchors, for all four constructions).
- Reused, unedited, zero new OpenSim run: `data/msk_smoketest/subject2_walking1/static_optimization/
  so/walking1_StaticOptimization_{force,activation}.sto`; `data/msk_smoketest/subject2_walking1/
  fmax_pcsa_correction_test/so/walking1_StaticOptimization_force.sto`; `data/msk_smoketest/
  subject2_walking1/emg_driven/{emg_hybrid_force.sto,naive_fwd_only_force.sto}`; `data/msk_smoketest/
  subject2_walking1/ankle_reserve_fix/boost150/cmc/full_1p57s/
  subject2_walking1_cmc_full_1p57s_Actuation_force.sto`.
- Read in place, never modified (isolation-compliant, read-only external drive):
  `/media/anton/8838D60F38D5FBDE/mechanism_data/LabValidation_withVideos/subject2/EMGData/
  walking1_EMG.sto`.
- External literature verified live via `curl` to `eutils.ncbi.nlm.nih.gov` and
  `api.semanticscholar.org` (no credentials, read-only): PMID 9535592 (Finni et al. 1998), PMID
  19229043 (Fröberg et al. 2009), PMID 10731005 (Giddings et al. 2000), PMID 2081741 (Komi 1990,
  checked and set aside), PMID 1638639 (Komi/Fukashiro/Jisvinen 1992, checked and set aside).
- **Not touched**: `scripts/msk/contact_waveform_analysis.py`, `scripts/msk/contact_muscle_decomp.py`,
  `scripts/msk/fmax_pcsa_validation.py`, `scripts/msk/fmax_pcsa_correction_test_subject2.py`,
  `scripts/msk/emg_constrained_force.py`, `scripts/msk/emg_driven.py`, and the concurrent
  `scripts/msk/contact_waveform_fmax_corrected.py` (independent, in-progress work by another instance).

No git commit, no git push performed (isolation respected, per `COORDINATOR.md` §1 and the task). All new
files are untracked, for the coordinator to commit.
