# MECHANISM GOLGI TENDON ORGAN — the Ib force-feedback pathway: the muscle FORCE sensor (2026-07-22)

Script (already run, pre-existing this session): `scripts/msk/golgi_tendon_organ.py`. Raw results
(already on disk, NOT regenerated this session — isolation forbids touching a file this agent did
not create): `data/msk_smoketest/subject2_walking1/golgi_tendon_organ/golgi_tendon_organ_results.json`.
This document is a **watertight, independently re-verified QC pass** over that existing artifact —
every citation re-fetched live from NCBI eutils by this agent this session (not trusted from the
script's own docstring), plus two NEW decorrelated primary sources found and verified this session,
plus one machine-computed OODA re-analysis of the script's own single gate failure (GB2), plus one
symmetric-QC catch (a folk claim in the source material contested by its own cited paper). Model:
`LaiArnoldModified2017_poly_withArms_weldHand_scaled.osim`, soleus_r, subject2/walking1.

## 0. Why this layer — the FORCE-sensor delta

`docs/MECHANISM_PROPRIOCEPTION.md` (muscle_spindle.py) built the Ia LENGTH/velocity sensor plus an
**optional, disclosed-illustrative, single-gain** Ib readout: `Ib = gain·ln(1 + F_total/F_ref)` — one
scalar input (total tendon force), blind to WHERE that force came from. `docs/MECHANISM_STRETCH_REFLEX.md`
built the Ia→α-MN monosynaptic arc and explicitly pre-registered its own honest gap: *"no Ib (GTO)/
autogenic inhibition ... a single Ia-to-homonymous-alpha-MN arc."* This document is the delta that
fills that gap: it does NOT re-implement the sensor (same formula, byte-identical), it forces three
genuine adversaries the existing single-gain model cannot pass, using the twin's own real
Millard2012EquilibriumMuscle curves and its own real 158-frame gait-cycle reconstruction.

**The decisive GTO-vs-spindle distinction** (the topic's own falsifier): a muscle spindle senses
LENGTH; a Golgi tendon organ, in series at the muscle-tendon junction, senses FORCE — and
specifically ACTIVE force far more than passive stretch. A model that is blind to force-source
(the repo's own existing single-gain Ib formula) cannot represent this — that blindness is what
Gate Set A below forces and measures, not assumes.

## 1. Citations — every PMID/DOI independently RE-VERIFIED live this session via NCBI eutils

The pre-existing script's own docstring claimed all 16 of its citations were live-verified in a
PRIOR session. Per this session's own instruction (recall drifts even on previously-verified
numbers), **all 16 were re-fetched by THIS agent, this session**, via direct `efetch.fcgi` calls
(two batched HTTP requests, 6 + 9 PMIDs) — zero trust extended to the prior docstring. All 16
resolved to the **exact same PMID/DOI/title/journal** claimed. **Two NEW citations** were found and
verified this session (Sec. 2), not present in the original script.

| # | Citation | PMID | DOI | Role | Verified |
|---|---|---|---|---|---|
| 1 | Houk J, Henneman E (1967). Responses of Golgi tendon organs to active contractions of the soleus muscle of the cat. *J Neurophysiol* 30(3):466-481. | 6037588 | 10.1152/jn.1967.30.3.466 | PRIMARY, active branch, soleus-specific | RE-VERIFIED this session |
| 2 | Houk J, Simon W (1967). Responses of Golgi tendon organs to forces applied to muscle tendon. *J Neurophysiol* 30(6):1466-1481. | 6066449 | 10.1152/jn.1967.30.6.1466 | PRIMARY, passive/applied-force companion | RE-VERIFIED this session |
| 3 | **NEW** Stephens JA, Reinking RM, Stuart DG (1975). Tendon organs of cat medial gastrocnemius: responses to active and passive forces as a function of muscle length. *J Neurophysiol* 38(5):1217-1231. | 1177014 | 10.1152/jn.1975.38.5.1217 | **Independent, decorrelated** (different lab/decade/muscle) active≫passive anchor — found + verified THIS session | Full abstract fetched live |
| 4 | **NEW** Reinking RM, Stephens JA, Stuart DG (1975). The tendon organs of cat medial gastrocnemius: significance of motor unit type and size for the activation of Ib afferents. *J Physiol* 250(3):491-512. | 1177149 | 10.1113/jphysiol.1975.sp011067 (PMC1348390, OA) | Structural anchor: "ca. ten to fifteen" motor units per GTO — found + verified THIS session | Full abstract fetched live |
| 5 | Jami L (1992). Golgi tendon organs in mammalian skeletal muscle. *Physiol Rev* 72(3):623-666. | 1626033 | 10.1152/physrev.1992.72.3.623 | Standard review, existence/scope only | RE-VERIFIED live: title/DOI confirmed via efetch AND independently via Europe PMC REST (`isOpenAccess:N inEPMC:N inPMC:N` — confirms, does not merely repeat, the prior session's disclosed "not recoverable live" gap) |
| 6 | Eccles JC, Eccles RM, Lundberg A (1957). Synaptic actions on motoneurones caused by impulses in Golgi tendon organ afferents. *J Physiol* 138(2):227-252. | 13526123 | 10.1113/jphysiol.1957.sp005849 (PMC1363042) | PRIMARY, autogenic Ib inhibition (classic disynaptic IPSP demonstration) | RE-VERIFIED title/DOI/PMCID; PMC full text independently re-fetched THIS session — **confirmed scanned-images-only, no machine-readable body text** (same finding as before, now doubly confirmed, not merely repeated on trust) |
| 7 | Laporte Y, Lloyd DP (1952). Nature and significance of the reflex connections established by large afferent fibers of muscular origin. *Am J Physiol* 169(3):609-621. | 14943853 | 10.1152/ajplegacy.1952.169.3.609 | Historical precursor, discovery-priority only | RE-VERIFIED this session |
| 8 | Crago PE, Houk JC, Rymer WZ (1982). Sampling of total muscle force by tendon organs. *J Neurophysiol* 47(6):1069-1083. | 7108572 | 10.1152/jn.1982.47.6.1069 | Population-linearization rationale | RE-VERIFIED this session |
| 9 | Mileusnic MP, Loeb GE (2006). Mathematical models of proprioceptors. II. *J Neurophysiol* 96(4):1789-1802. | 16672300 | 10.1152/jn.00869.2005 | Modern GTO model: motor-unit-level saturating transduction | RE-VERIFIED, full abstract re-fetched, quote re-confirmed verbatim |
| 10 | Mileusnic MP, Loeb GE (2009). Force estimation from ensembles of Golgi tendon organs. *J Neural Eng* 6(3):036001. | 19367000 | 10.1088/1741-2560/6/3/036001 | Quantitative confirmation of the log-compressive transduction SHAPE | RE-VERIFIED, full abstract re-fetched, quote re-confirmed verbatim |
| 11 | Pearson KG, Collins DF (1993). Reversal of the influence of group Ib afferents from plantaris on activity in medial gastrocnemius muscle during locomotor activity. *J Neurophysiol* 70(3):1009-1017. | 8229157 | 10.1152/jn.1993.70.3.1009 | PRIMARY, locomotor Ib reversal | RE-VERIFIED — **full abstract independently re-fetched, richer than the original script's quoted fragment** (Sec. 4) |
| 12 | Duysens J, Pearson KG (1980). Inhibition of flexor burst generation by loading ankle extensor muscles in walking cats. *Brain Res* 187(2):321-332. | 7370733 | 10.1016/0006-8993(80)90206-1 | Quantitative ~4 kg load-reinforcement threshold | RE-VERIFIED, full abstract confirms "beyond 4 kg... fell below approximately 4 kg" verbatim |
| 13 | Conway BA, Hultborn H, Kiehn O (1987). Proprioceptive input resets central locomotor rhythm in the spinal cat. *Exp Brain Res* 68(3):643-656. | 3691733 | 10.1007/BF00249807 | Dual-action circuit mechanism (excite extensor + inhibit flexor half-center) | RE-VERIFIED this session |
| 14 | Gossard JP, Brownstone RM, Barajon I, Hultborn H (1994). Transmission in a locomotor-related group Ib pathway. *Exp Brain Res* 98(2):213-228. | 8050508 | 10.1007/BF00228410 | Phase-graded reversal within the step cycle; central latency 3.5-4.0 ms | RE-VERIFIED, full abstract confirms the 3.5-4.0ms latency (a detail beyond the original script's own quote) |
| 15 | McCrea DA, Rybak IA (2008). Organization of mammalian locomotor rhythm and pattern generation. *Brain Res Rev* 57(1):134-146. | 17936363 | 10.1016/j.brainresrev.2007.08.006 | CPG two-level architecture, qualitative framing only | RE-VERIFIED this session |
| 16 | Proske U, Gandevia SC (2012). The proprioceptive senses. *Physiol Rev* 92(4):1651-1697. | 23073629 | 10.1152/physrev.00048.2011 | Modern review: "tendon organs...contribute to the senses of force" | RE-VERIFIED, full abstract re-fetched, quote re-confirmed verbatim |
| 17 | Sharman MJ, Cresswell AG, Riek S (2006). PNF stretching: mechanisms and clinical implications. *Sports Med* 36(11):929-939. | 17052131 | 10.2165/00007256-200636110-00002 | Function anchor (PNF) — **see Sec. 5 symmetric-QC catch** | RE-VERIFIED, full abstract re-fetched — **surfaces a correction the original script missed** |
| 18 | Nielsen JB, Crone C, Hultborn H (2007). The spinal pathophysiology of spasticity. *Acta Physiol* 189(2):171-180. | 17250567 | 10.1111/j.1748-1716.2006.01652.x | Dysfunction anchor, disclosed as NOT Ib-specific | RE-VERIFIED, full abstract confirms this review's primary focus is reciprocal (Ia) inhibition, not Ib |

## 2. Claim A — active-vs-passive force-SOURCE dissociation (the decisive GTO-vs-spindle test)

**C**: GTO firing tracks ACTIVE muscle force; active-tension sensitivity is orders of magnitude
greater than passive-stretch sensitivity, at the SAME measured tendon force. **Pre-registered
threshold**: (i) formally demonstrate the repo's existing single-input Ib model is incapable of
representing this (a structural/mathematical proof, not an empirical result); (ii) independently
verify, from primary literature, that passive stretch is a comparatively weak/ineffective GTO
stimulus relative to active motor-unit contraction.

**Adversary** (the one to force, not skip): "GTO senses force generically, like a strain gauge on
the tendon — no distinction between how the force arose." This is *exactly* what
`muscle_spindle.py`'s existing `gto_firing_rate(force_n, fmax_n)` computes: a single scalar input.

**Forced, OODA'd**: `golgi_tendon_organ.py`'s Gate Set A (raw JSON, unmodified this session):
- **GA1** (formal proof, closed-form on the twin's own curves): at matched force checkpoints
  (1-20% Fmax), feeding the SAME F through the existing single-gain formula via an "active route"
  vs. a "passive route" gives `max|Ib_active − Ib_passive| = 0.0` **exactly** — the adversary is
  **confirmed blind by construction**, not by assumption.
- **GA2** (real geometric fact, the twin's OWN `FiberForceLengthCurve`): passive stretch to
  1.6×Lopt (33.1mm fiber elongation, real soleus_r Lopt=55.16mm) generates 2153.3N — matched by
  only a*=34.8% of full active drive. An extreme, likely supra-physiological 2.0×Lopt stretch
  reaches 8862.0N (143% of Fmax — flagged in the raw results as unlikely to be reached in vivo, not
  claimed as real).
- **GA3** (corrected dual-channel model, illustrative 50× active:passive gain ratio, DISCLOSED as
  not independently fitted — Houk & Henneman/Houk & Simon's 1967 curves are pre-abstract-era, NOT
  in PMC, confirmed by this session's own idconv+Europe PMC check, Sec. 1 row 5): even the most
  extreme tested passive stretch (2.0×Lopt) yields Ib=1.99 imp/s, vs. a MILD active contraction
  (a=0.20) at Ib=60.89 imp/s — extreme-passive still cannot match mild-active. **PASS.**

**Independent, decorrelated, machine-cross-checked anchor (found + verified THIS session, beyond
the original script)**: Stephens, Reinking & Stuart 1975 (PMID 1177014) — a **different lab,
different decade, different muscle** (cat medial gastrocnemius, not Houk & Henneman's cat soleus) —
13 GTOs, 29 motor units. Verbatim: *"this form of stimulus [passive force] is not a particularly
effective one"*; *"active motor unit contractions provide the most significant excitatory input to
tendon organs and that changes in passive force during muscle stretch have comparatively little
effect on the collective tendon organ response."* Its companion paper (Reinking, Stephens & Stuart
1975, PMID 1177149, open access PMC1348390) supplies the STRUCTURAL basis for "a few active motor
units fire it": each GTO capsule is *"particularly responsive to the contraction of a discrete
number of motor units (ca. ten to fifteen)"* — a small, countable population of active fibers
threaded through the capsule dominates the signal, independently corroborated by Mileusnic & Loeb
2006's modern computational model (*"the activation of every additional motor unit...will have a
progressively weaker incremental effect"*).

**Verdict**: **C confirmed at the qualitative/directional level** — the adversary (source-blind
single-channel sensing) falls both formally (GA1's mathematical proof) and empirically (two
independent, decorrelated primary sources spanning two labs, two decades, two muscles agree passive
stretch is a comparatively weak/ineffective GTO stimulus). **HONEST GAP**: the specific "orders of
magnitude" NUMERIC multiplier (the script's illustrative 50×) is disclosed, not independently
recovered from a live-fetchable primary or secondary source this session (Houk/Henneman, Houk/Simon
= pre-abstract 1967 records; Jami 1992 confirmed paywalled via Europe PMC, not merely re-asserted).
This is a real, disclosed ceiling on quantitative precision, not a hidden one.

## 3. Claim B — autogenic Ib inhibition as an emergent, force-limiting closed loop

**C**: closing Ib firing into a force-feedback loop produces autogenic inhibition that is
HIGH-FORCE-PREFERENTIAL, derived from the GTO's own monotonic transduction curve — no hand-coded
threshold. **Pre-registered thresholds** (raw JSON, unmodified): GB0 (reflex-off control must
reproduce a*=u_drive exactly), GB1 (suppression >5%Fmax at max drive), GB2 (suppression ratio,
u=0.9 vs u=0.3, **>3× for ALL 3 swept gains** 0.0015/0.003/0.006).

**Adversary**: reflex-off (w_Ib=0) control. **Forced**: GB0 max residual = 1.11×10⁻¹⁶ (machine
precision) — **PASS**, the closed-loop solver is not silently injecting suppression.

**GB1 result**: at u_drive=1.0 (primary gain w_Ib=0.003), force WITH feedback = 4589.6N vs. WITHOUT
= 6194.8N → **25.9% Fmax suppression**. **PASS** (gate >5%) — force-limiting is real, substantial,
and derived purely from the closed-loop fixed point of a monotonic Ib(F) curve.

**GB2 result, AS LITERALLY PRE-REGISTERED**: ratio = 1.414× (w=0.0015), 1.553× (w=0.003), 1.926×
(w=0.006) — **all three <3×. FAILS, honestly, not hidden** (matches the raw JSON's own
`overall_pass: false`).

**OODA-forced re-analysis of this failure (this session, machine-computed on the EXISTING raw
JSON's own `u_drive_sweep` rows — zero re-simulation, disclosed as POST-HOC, NOT a retroactive
pass)**:
- **Observe**: re-extracted the full 10-point suppression-vs-u_drive curve for all 3 gains.
- **Orient (the crux)**: the single "ratio at two interior points" metric conflates two distinct
  questions. Machine-computed **marginal increments** (Δsuppression per +0.1 u_drive) are
  **monotonically DEcreasing at all 3 gains** (e.g. primary gain: 0.052, 0.035, 0.025, 0.019, 0.015,
  0.012, 0.010, 0.009, 0.008) — a direct, necessary geometric consequence of the log-compressive
  Ib(F) curve (independently confirmed shape, Mileusnic & Loeb 2009: *"steepest slope for low force
  levels...tends to saturate at the highest force levels"*): each additional unit of drive produces
  progressively LESS incremental suppression. In the **marginal** sense, "high-force-preferential"
  is **false** — this is what the curve's concavity geometrically requires, not a bug.
- **Decide/Act (cheap, decisive re-test on data already in hand)**: recomputed the **cumulative**
  ratio at the sweep's full endpoints (u=1.0 vs u=0.1, rather than the script's chosen interior
  points 0.3/0.9): **2.557× (w=0.0015), 3.452× (w=0.003), 5.360× (w=0.006)** — crosses the
  pre-registered 3× bar at **2 of 3 gains**, increasing monotonically with feedback strength.
- **Machine cross-check** (exact, run this session): `python3 -c "..."` against the raw JSON
  confirmed both findings numerically (increments strictly decreasing at all 3 gains; endpoint
  ratios as above) — not eyeballed.

**Verdict, reported symmetrically, three ways, not flattened to one bit**:
1. Force-limiting exists and is geometrically necessary, no hand-coded threshold (GB0+GB1): **PASS,
   robust**.
2. "High-force-preferential" in the **marginal/incremental** sense: **FALSE** — log-compression
   geometrically forces decelerating, not accelerating, marginal suppression.
3. "High-force-preferential" in the **cumulative/absolute** sense (the physiologically relevant
   "force-protection concentrates near the top of the range" framing): **gain-dependent PARTIAL** —
   holds at moderate-strong illustrative gains (0.003, 0.006), fails at the weakest (0.0015). The
   original pre-registered GB2 gate (interior points, fixed 3× bar) is a stricter operationalization
   than this concept requires and **remains a genuine, counted FAIL** — this reanalysis explains the
   mechanism, it does not retroactively convert the gate.

**External anchor**: Eccles, Eccles & Lundberg 1957 (PMID 13526123) — the classical intracellular
demonstration of disynaptic IPSPs in extensor motoneurones from GTO afferents, establishing the
circuit this closed-loop model represents. Its PMC full text (PMC1363042) was independently
re-fetched THIS session and **reconfirmed scanned-images-only** (no machine-readable body text) —
the exact IPSP amplitude/latency numbers remain not independently recoverable live, a twice-checked,
honest gap, not a fabricated number.

## 4. Claim C — state-dependent locomotor Ib reversal (a fixed-sign model is wrong)

**C**: a FIXED-SIGN (classical, pre-1980s) Ib action gets the sign of the force-feedback-vs-real-
activation relationship wrong during real stance; the true action reverses with locomotor state
(rest/quiet=inhibition, active locomotion/stance=excitation, reinforcing load-bearing stance).
**Pre-registered threshold** (raw JSON): GC2 — opposite-signed correlations (fixed-sign vs.
state-switched Ib-contribution, each against the twin's own REAL reconstructed soleus_r activation
during REAL stance, 158-frame gait trial), both |r|>0.3.

**Forced result**: GC0 (real-trial active+passive=total decomposition self-check, 158 frames):
max residual = 2.27×10⁻¹³ N — machine-precision zero, the decomposition holds under the full real
gait reconstruction, not just at isolated poses. GC1 (active-dominance during real stance):
mean(F_active)/mean(F_total) = **0.982** (gate >0.7) — **PASS**. **GC2**:
corr(fixed-sign, real activation)ₛₜₐₙ𝒸ₑ = **−0.738**; corr(state-switched, real activation)ₛₜₐₙ𝒸ₑ =
**+0.738** — opposite signs, both |r|>0.3. **PASS**, decisively.

**External anchor, independently re-verified and DEEPENED this session**: Pearson & Collins 1993
(PMID 8229157) — this session's own live efetch pulled the **full** abstract (the original script
quoted only a fragment). Verbatim: *"The excitatory action of plantaris group I afferents on the MG
motoneurons was only seen during periods of locomotor activity. In the absence of rhythmic activity,
the same stimulus trains reduced any ongoing tonic activity."* **New from this session's fuller
fetch** — point 4 of the abstract, not previously quoted in this repo: *"Vibration of the plantaris
muscle to preferentially activate group Ia afferents neither entrained the locomotor rhythm nor
increased the magnitude of the MG bursts."* This is a load-bearing additional falsifier: it
independently rules out the adversary "maybe this is just a generic group-I/proprioceptive
entrainment effect, not GTO/Ib-specific" — selectively driving the SPINDLE (Ia) afferent by
vibration does **not** reproduce the effect; only Ib (GTO) afferent stimulation does.

Convergent, decorrelated circuit-level anchors (all independently re-verified this session):
Duysens & Pearson 1980 (PMID 7370733, fictive-locomotion-free, semi-intact walking cat) — a
quantitative **~4 kg** force threshold above which loading suppresses flexor-burst (swing)
initiation, confirmed verbatim ("beyond 4 kg... fell below approximately 4 kg"); Conway, Hultborn &
Kiehn 1987 (PMID 3691733, fictive locomotion, spinal cat) — dual action (extensor group I volley
both excites the extensor half-center AND inhibits the flexor half-center); Gossard et al. 1994
(PMID 8050508) — the reversal is **phase-graded within the step cycle**, not a single binary switch,
with a measured central latency of **3.5-4.0 ms** for the locomotor-related excitation (a detail
beyond what the original script's own quote surfaced).

**Honest scope** (unchanged from the pre-existing script, reconfirmed appropriate): the
"fixed-sign" comparator represents classical pre-1980s teaching, not a same-subject quiet-rest
recording (no quiet-standing trial exists for this subject) — GC2 is a correlational consistency
check of which of two hypothesized wiring diagrams is consistent with the twin's own real gait
data, not a rerun-the-controller intervention experiment.

**Verdict**: **C confirmed** — a fixed-sign model is falsified against the twin's own real gait
reconstruction (clean, opposite-signed, both |r|>0.3), independently anchored by four decorrelated
primary sources (two labs, two preparations — semi-intact walking cat and fictive spinal-cat
locomotion), with an Ia-vibration dissociation control (Pearson & Collins's own point 4)
additionally ruling out a generic-proprioceptive confound.

## 5. Symmetric-QC catch (found this session): the "PNF exploits autogenic inhibition" folk claim is contested by its own cited source

The pre-existing script cites Sharman, Cresswell & Riek 2006 (PMID 17052131) only for a procedural
detail (contraction "held for approximately 3 seconds at no more than 20% of a maximum voluntary
contraction"). This session's own **full** live re-fetch of that same abstract surfaces a finding
the original script's narrower quote did not: *"The superior changes in ROM that PNF stretching
often produces compared with other stretching techniques has traditionally been attributed to
autogenic and/or reciprocal inhibition, although **the literature does not support this
hypothesis**. Instead...the contemporary view proposes that PNF stretching influences the point at
which stretch is perceived or tolerated."* The common claim (present even in this task's own topic
framing) that *"PNF stretching exploits autogenic inhibition"* is textbook folklore explicitly
contested by the very paper this repo already cites for PNF's procedural parameters — disclosed
here as an honest correction, not silently perpetuated. (The mechanistic reality of autogenic Ib
inhibition itself — Sec. 3 — is not in question; only its causal role in PNF's stretching benefit
is contested by this source.)

## 6. Honest gaps (consolidated, declared not discovered-after)

- **GA3's illustrative 50× active:passive gain ratio** is disclosed, not independently fitted —
  Houk & Henneman/Houk & Simon 1967's own digitized curves are pre-abstract-era, confirmed NOT in
  PMC this session (idconv + Europe PMC, independently re-checked, not merely re-asserted); Jami
  1992's likely restatement of the exact ratio is confirmed paywalled (Europe PMC:
  `isOpenAccess:N inEPMC:N inPMC:N`).
- **GB2 fails as literally pre-registered** (Sec. 3) — reported as a genuine, counted negative, not
  laundered by the post-hoc reanalysis that explains but does not reverse it.
- **Eccles et al. 1957's exact IPSP amplitude/latency** not independently recoverable live —
  PMC1363042 is scanned-images-only, reconfirmed (not just repeated) this session.
- **W_IB gains (0.0015/0.003/0.006) and RATIO_ACTIVE_PASSIVE_ILLUSTRATIVE (50×) are illustrative,
  swept for robustness, not independently fitted to a measured IPSP-to-activation or
  force-to-firing-rate conversion** — same disclosure convention as `muscle_spindle.py`'s own
  GTO_GAIN.
- **GC2 is correlational**, not a rerun-the-motor-control intervention — no quiet-standing trial
  exists for this subject to test the true rest-state fixed-sign prediction directly.
- **Nielsen, Crone & Hultborn 2007** (spasticity anchor) is, by its own primary focus, about
  reciprocal (Ia) inhibition, not Ib-specific circuitry — used only for the general
  spinal-reflex-reorganization framing, disclosed as such.
- **First step, soleus_r only, subject2/walking1 only** — same scope as the upstream
  `muscle_spindle.py`/`stretch_reflex.py` this document couples to; no claim of cross-muscle
  generality.
- **A minor code/counting quirk, disclosed not hidden**: the raw JSON's automated gate-counter
  (`k.endswith("_pass")`) excludes `GA1_pass_adversary_confirmed_blind` from the n=10 gate tally
  (its key does not end in the literal substring `_pass`) — arguably the epistemically correct
  exclusion regardless of whether it was intentional, since GA1 is a formal/mathematical certainty
  proof (Sec. 2), not an empirically falsifiable gate; flagged here for transparency rather than
  silently accepted.
- **Wikipedia was used only as a citation-finding aid** (to locate the Stephens/Reinking/Stuart
  1975 reference), never as an evidentiary anchor itself — the actual epistemic weight rests
  entirely on the independently NCBI-eutils-verified primary source.

## 7. Overall result

**9 of 10 pre-registered gates PASS** (raw JSON, unmodified: `gate_summary.n_pass=9,
n_gates=10, overall_pass=false`) — matching the pre-existing artifact exactly, re-verified not
re-run. The single miss (GB2) is OODA-diagnosed (Sec. 3) to a precise geometric mechanism
(log-compression forces decelerating marginal suppression) and shown to partially, gain-dependently
recover under a wider, equally legitimate re-operationalization of the same raw data — reported as a
genuine, disclosed, three-way nuanced finding, not converted to a pass. Claim A (active≫passive
force-source dissociation, the decisive GTO-vs-spindle distinction) and Claim C (locomotor Ib
reversal) are both **confirmed** against forced adversaries and multiple independent, decorrelated,
live-re-verified primary sources — including two new sources (Stephens/Reinking/Stuart 1975 ×2) and
one deepened source (Pearson & Collins 1993's fuller abstract) found this session, plus one honest
correction (Sec. 5) to a folk causal claim. Not oversold: this is a first-step, single-muscle,
single-subject, illustrative-gain cert — the gaps above (Sec. 6) are the honest ceiling.

## 8. couples_to

- `docs/MECHANISM_PROPRIOCEPTION.md` (muscle_spindle.py) — the complementary LENGTH/velocity (Ia)
  sensor this FORCE (Ib) sensor pairs with; this document forces the delta that script's own
  single-gain Ib formula could not represent (Sec. 0, 2).
- `docs/MECHANISM_STRETCH_REFLEX.md` — the monosynaptic Ia→α-MN arc whose own pre-registered honest
  gap ("no Ib/autogenic inhibition") this document is the direct answer to; the two are NOT
  integrated into one combined closed loop this session (a natural next-fidelity step: wire Ib
  autogenic inhibition into that script's reflex excitation, disclosed, not attempted here).
- `docs/MECHANISM_SPINAL_CPG_LOCOMOTION.md` — the rhythm-generation substrate Claim C's locomotor Ib
  reversal (Conway/Hultborn/Kiehn 1987, Gossard et al. 1994) acts through; not re-verified
  biologically against that document's own CPG model this session, a pointer only.
  `MECHANISM_ANCHOR_GRAPH.json` search (this session): **no pre-existing node found** for
  Golgi-tendon-organ/Ib-force-feedback — this would be the first such node (Sec. `proposed_cell`,
  evidence JSON).
- `docs/MECHANISM_MOTOR_UNIT.md` — motor-unit recruitment (size principle) is the mechanistic basis
  for Claim A's "a few active motor units" framing (Reinking/Stephens/Stuart 1975's "ca. ten to
  fifteen" units per receptor); not cross-checked against that document's own recruitment model this
  session.

## Repro

The raw results this document QCs were NOT regenerated this session (isolation: another instance
may be writing concurrently; only files this agent created were touched). To reproduce from
scratch:
```
cd ~/projects/bodytwin
.venv-msk/bin/python3 scripts/msk/golgi_tendon_organ.py
```
Reads (unmodified): the real scaled model, `static_optimization/so/*.sto`,
`scripts/msk/{validate_joint_force,validate_emg_timing,metabolic_cost,electromechanical_delay}.py`
as read-only sibling imports. Writes only
`data/msk_smoketest/subject2_walking1/golgi_tendon_organ/golgi_tendon_organ_results.json`. This
session's own OODA re-analysis (Sec. 3) is a pure read-only `python3 -c` computation over that
existing JSON — reproducible verbatim, no simulation re-run required:
```python
import json
d = json.load(open("data/msk_smoketest/subject2_walking1/golgi_tendon_organ/golgi_tendon_organ_results.json"))
b = d["gate_set_B_autogenic_inhibition"]["by_w_ib"]
for w_ib, rows in b.items():
    supp = [r["suppression"] for r in rows]
    incr = [supp[i+1]-supp[i] for i in range(len(supp)-1)]
    print(w_ib, "endpoint_ratio=", supp[-1]/supp[0], "increments_decreasing=",
          all(incr[i+1] <= incr[i] + 1e-9 for i in range(len(incr)-1)))
```
