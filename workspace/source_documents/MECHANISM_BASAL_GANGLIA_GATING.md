# MECHANISM BASAL GANGLIA GATING — direct/indirect push-pull action-selection gate, PD/HD double dissociation, beta-oscillation spectral falsifier (2026-07-22)

Builds and MEASURES the basal ganglia direct(D1/Go)-indirect(D2/NoGo) circuit as a
movement GAIN-CONTROL gate on the cortically-selected action, with dopamine (SNc)
setting the D1/D2 balance. Script: `scripts/msk/basal_ganglia_gating.py`. Evidence:
`data/basal_ganglia_gating/basal_ganglia_gating_results.json`.

**Provenance note, stated up front (per this repo's own discipline — nothing hidden until
the reader stumbles on it):** a repo-wide check found this exact topic already has TWO
pre-existing SEED-DESIGN nodes in `data/MECHANISM_ANCHOR_GRAPH.json`
(`NEU-BASAL-GANGLIA-GATING`, cluster `brain-basal-ganglia`; and a near-duplicate
`AUTO-BG-ACTION-GATING-CELL-IMPLEMENT-DIRECT-D`, cluster `autonomous-expansion`, carrying
a spurious unrelated `GOAL-CANCER-IMMUNE-WARGAME-TWIN` couples_to entry — apparent
graph-building noise, flagged here, **not edited**, per isolation). Both nodes cite the
same scout output (`data/body_twin/agent_outputs/basal-ganglia-gating__a93bab9e997bceb06.json`,
re-read this session and treated as an **unverified hypothesis**, not trusted) and are
`status:OPEN`, `mechanism_grade:SEED-DESIGN` ("cited literature anchor, not yet executed
against data") — i.e. the node has been through ACQUIRE+DESIGN
(`MECHANISM_HARDENED_CONVENTIONS.md` §4b) but not MEASURE/CERTIFY. This build advances that
gap: it independently re-derives and **measures** two of the SEED-DESIGN node's three named
anchors — (b) Kravitz optogenetic causal dissociation (full text fetched, exact % numbers
used) and (c) Reiner HD pathway-loss (re-verified, extended with the Albin 1992
presymptomatic confirmation) — and **adds** the beta-oscillation electrophysiological
signature this session's own task explicitly named, which the SEED-DESIGN node did not
carry at all. It does **not** execute anchor (a) (Frank et al ON/OFF-medication **Go/NoGo
learning-bias** magnitude — a cognitive/RL claim, decorrelated from but adjacent to this
build's movement-**gating** claim) or the habit/procedural-memory (DMS→DLS) coupling that
node also names — a disclosed **partial** contribution, not a resolution, same precedent as
`MECHANISM_CIRCADIAN_RHYTHM.md`'s disclosed partial coverage of its own 5 pre-existing nodes.
No graph-edge write this session (isolation: touch only files created this session).

**A caught citation-drift instance** (measured, not asserted): the pre-existing scout output
claims Fearnley & Lees 1991 (PMID 1933245) gives a "logistic-regression ~31%" SNc-neuron-loss
estimate at symptom onset. This session's own **live re-fetch** of that abstract (verbatim,
§1 below) contains **no such number and no logistic regression** — the paper's own words are
"a 68% cell loss in the lateral ventral tier and a 48% loss in the caudal nigra as a whole",
matching (not contradicting) the sibling `MECHANISM_DOPAMINE_KINETICS.md` cert's own
independently-verified figures. The scout's 31% figure is **not used anywhere in this
build** — a live, concrete instance of this repo's own measured ~62–67%
citation-drift-from-memory rate, caught by direct primary-source verification rather than
trusted.

## The two falsifiers, verdicts stated up front (nothing hidden)

> **Falsifier 1 (the double dissociation)**: does a circuit with two independently-lesionable,
> oppositely-signed channels (a dopamine-tone/receptor-pharmacology axis; a
> cell-population/survival axis) converging on one output node (GPi/SNr) reproduce PD
> (SNc dopamine loss → GPi output **up** → bradykinesia) **and** HD (striatal indirect-MSN
> loss → GPi output **down** → chorea) — **while a collapsed single-channel ("ungated")
> null, forced to its strongest fair form, cannot get both signs right simultaneously**?

**YES for the real model, and the forced adversary genuinely falls.** Pre-registered bands:
PD dopamine-loss ∈ **[48%, 90%]** (reused **read-only** from the sibling
`dopamine_kinetics.py`'s own verified Fearnley-&-Lees/Kordower/Cheng-2010 band — a live
structural coupling, not a re-typed duplicate: the script `json.load`s that file's own
`compensatory_reserve_toy` numbers `[48,68]` + `80` at runtime); HD indirect(D2)-MSN
survival ∈ **[10%, 50%]**, direct(D1) survival ∈ **[85%, 100%]** ("relatively spared",
Reiner 1988/Albin 1992's own qualitative language — no exact abstract-level percentage
exists for either, disclosed). Across the **entire** anchored PD band (9-point sweep),
ΔGPi is **always positive** (+0.398 to +0.738); across the **entire** anchored HD grid
(20-point sweep, D2-survival × D1-survival), ΔGPi is **always negative** (−0.162 to
−0.324) — **opposite sign, robust across the whole diverse instance-space, not a knife-edge
point.** The forced adversary (§3below) — the single most natural, non-strawman
"collapsed" model, both lesions reducing one aggregate striatal-output channel, swept
across **both** possible sign conventions of its one free parameter — produces the
**same** sign for PD and HD under **both** branches (w=+1: both positive; w=−1: both
negative) — it cannot be rescued by any choice of its own free parameter. This is a rank
argument, not a tuning contest: a single output channel fed by two lesions that both enter
with the same sign can never span an opposite-sign target pair.

> **Falsifier 2 (beta oscillations, decorrelated/dynamical)**: does a minimal STN-GPe
> delay-loop model show a genuine spectral (Hopf) bifurcation into the beta band
> (13–30 Hz) as loop gain rises with dopamine loss — **while an open-loop (broken
> feedback) null structurally cannot oscillate at any gain**?

**YES.** The closed-loop characteristic equation `(τ₁s+1)(τ₂s+1) + G·e^(−sD) = 0` has an
exact marginal-stability boundary (solved via `s=iω`) at **f\* = 27.68 Hz, G\* = 4.025**
for an illustrative-but-disclosed literature-plausible primary point (τ₁=τ₂=10 ms,
round-trip delay D=6 ms) — **inside the pre-registered [13,30] Hz band**, not tuned to hit
it (the (τ,D) grid was fixed from physiological plausibility *before* computing outcomes;
**13/16 = 81%** of the 4×4 sensitivity grid independently lands in-band). Continuing the
dominant root by Newton iteration across a gain sweep shows Re(s) crosses zero
**monotonically** exactly at G≈G\*, with the emergent frequency staying inside a widened
[11,32] Hz band for **100%** of the post-threshold sweep. Driving G by the **same**
D2/indirect-drive function used in Falsifier 1 (no second free mechanism), the loop
crosses G\* at **exactly 50% DA loss** — inside the anchored PD band. The open-loop null
(GPe→STN connection severed, G≡0) reduces the characteristic equation to a plain quadratic
whose two roots are **exactly** −1/τ₁, −1/τ₂ = **−100, −100** (purely real, zero imaginary
part) — **for any value** of the remaining single connection's weight, since the product
defining G is identically zero. It cannot oscillate at any gain — a structural
impossibility, not a parameter-fishing result.

## 1. Citations — 23 PMIDs, every one verified LIVE this session (NCBI eutils
esearch/esummary/efetch, verbatim text; 2 PMC full-text pulls: Kravitz 2010, Reiner 1988)

| # | Citation | PMID / DOI | Role |
|---|---|---|---|
| 1 | Albin RL, Young AB, Penney JB (1989). "The functional anatomy of basal ganglia disorders." *Trends Neurosci* 12(10):366-75. | **2479133**, DOI 10.1016/0166-2236(89)90074-x | Founding direct/indirect rate-balance model; abstract itself states the hyperkinetic=indirect-loss / hypokinetic=output-increase double dissociation. |
| 2 | DeLong MR (1990). "Primate models of movement disorders of basal ganglia origin." *Trends Neurosci* 13(7):281-5. | **1695404**, DOI 10.1016/0166-2236(90)90110-v | Companion rate model; groups Huntington's + hemiballismus as the hyperkinetic pole, PD as the hypokinetic pole of one circuit-level spectrum. |
| 3 | Mink JW (1996). "The basal ganglia: focused selection and inhibition of competing motor programs." *Prog Neurobiol* 50(4):381-425. | **9004351**, DOI 10.1016/s0301-0082(96)00042-1 | Winner-take-all/surround-suppression framing: BG broadly inhibit competing programs, focally disinhibit the selected one. |
| 4 | Kravitz AV, Freeze BS, Parker PR, Kay K, Thwin MT, Deisseroth K, Kreitzer AC (2010). "Regulation of parkinsonian motor behaviours by optogenetic control of basal ganglia circuitry." *Nature* 466(7306):622-6. | **20613723**, DOI 10.1038/nature09159, PMC3552484 | **Decisive causal anchor, full text fetched.** SNr firing: D1-stim 8.6±3.0% of baseline (n=8); D2-stim 162±19% (n=4). Bidirectional freezing/locomotion/6-OHDA-rescue. |
| 5 | Reiner A, Albin RL, Anderson KD, D'Amato CJ, Penney JB, Young AB (1988). "Differential loss of striatal projection neurons in Huntington disease." *PNAS* 85(15):5733-7. | **2456581**, DOI 10.1073/pnas.85.15.5733, PMC281835 | HD indirect(enkephalin+)-preferential loss, early/middle stages; n=17 specimens. PMC record is scanned/pre-digital — no exact grade-by-grade % extractable this session. |
| 6 | Albin RL, Reiner A, Anderson KD, Dure LS 4th, Handelin B, Balfour R, Whetsell WO Jr, Penney JB, Young AB (1992). "Preferential loss of striato-external pallidal projection neurons in presymptomatic Huntington's disease." *Ann Neurol* 31(4):425-30. | **1375014**, DOI 10.1002/ana.410310412 | Confirms indirect-pathway loss precedes SYMPTOM ONSET (presymptomatic allele carriers) — the push-pull imbalance predates overt chorea. |
| 7 | Bergman H, Wichmann T, DeLong MR (1990). "Reversal of experimental parkinsonism by lesions of the subthalamic nucleus." *Science* 249(4975):1436-8. | **2402638**, DOI 10.1126/science.2402638 | Causal reversal: STN lesion in MPTP monkeys reduces akinesia/rigidity/tremor — direct test of "excess STN drive causes parkinsonism." |
| 8 | Bergman H, Wichmann T, Karmon B, DeLong MR (1994). "The primate subthalamic nucleus. II. Neuronal activity in the MPTP model of parkinsonism." *J Neurophysiol* 72(2):507-20. | **7983515**, DOI 10.1152/jn.1994.72.2.507 | STN spontaneous rate 19±10→26±15 spikes/s after MPTP (n=220→326 cells); 4-8Hz oscillatory-cell fraction 2%→16% (STN), 0.6%→25% (GPi). |
| 9 | Filion M, Tremblay L (1991). "Abnormal spontaneous activity of globus pallidus neurons in monkeys with MPTP-induced parkinsonism." *Brain Res* 547(1):142-51. | **1677607** | Independent lab/species (cynomolgus vs Bergman's African green): GPi rate INCREASES, GPe rate DECREASES in MPTP parkinsonism — direction-only (abstract gives no exact Hz), corroborates Bergman 1994. |
| 10 | Hamada I, DeLong MR (1992). "Excitotoxic acid lesions of the primate subthalamic nucleus result in transient dyskinesias of the contralateral limbs." *J Neurophysiol* 68(5):1850-8. | **1479448**, DOI 10.1152/jn.1992.68.5.1850 | STN lesion → dyskinesias emerge ~20min post-lesion, peak 60-80min, resolve <4h — behavioral companion to #11. |
| 11 | Hamada I, DeLong MR (1992). "Excitotoxic acid lesions of the primate subthalamic nucleus result in reduced pallidal neuronal activity during active holding." *J Neurophysiol* 68(5):1859-66. | **1479449**, DOI 10.1152/jn.1992.68.5.1859 | **Exact matched pair**: GPi 69.8±21.6→47.4±22.6 spikes/s (n=169→180, P<0.001); GPe 63.6±25.1→41.0±18.1 (n=218→208) after STN lesion. Same-node-different-lesion-site analogy to HD. |
| 12 | Brown P, Oliviero A, Mazzone P, Insola A, Tonali P, Di Lazzaro V (2001). "Dopamine dependency of oscillations between subthalamic nucleus and pallidum in Parkinson's disease." *J Neurosci* 21(3):1033-8. | **11157088**, DOI 10.1523/JNEUROSCI.21-03-01033.2001, PMC6762327 | GPi-STN coherence dominated by <30Hz OFF medication; L-DOPA collapses it, new ~70Hz peak emerges — direct dopamine-dependence of the oscillation. n=4 patients/6 recordings. |
| 13 | Kühn AA, Kupsch A, Schneider GH, Brown P (2006). "Reduction in subthalamic 8-35 Hz oscillatory activity correlates with clinical improvement in Parkinson's disease." *Eur J Neurosci* 23(7):1956-60. | **16623853**, DOI 10.1111/j.1460-9568.2006.04717.x | L-DOPA-induced beta-power reduction correlates with UPDRS improvement: r=0.811 (total, P<0.001), r=0.835 (akinesia-rigidity, P<0.001), NOT tremor. n=9, 17 sides. |
| 14 | Kühn AA, Kempf F, Brücke C, et al., Brown P (2008). "High-frequency stimulation of the subthalamic nucleus suppresses oscillatory beta activity in patients with Parkinson's disease in parallel with improvement in motor performance." *J Neurosci* 28(24):6165-73. | **18550758**, PMC6670522 | DBS/HFS suppresses beta for 12s post-stim, in parallel with motor gain; explicit "beta (13-30 Hz)" band definition matching this build's own. n=11. Decorrelated (stimulation-based, not pharmacological) confirmation of #13. |
| 15 | Krack P, Batir A, Van Blercom N, et al. (2003). "Five-year follow-up of bilateral stimulation of the subthalamic nucleus in advanced Parkinson's disease." *NEJM* 349(20):1925-34. | **14614167** | STN-DBS: off-medication UPDRS motor improved 54% at 5yr (P<0.001), ADL 49%. n=49 (42 completers). |
| 16 | Limousin P, Krack P, Pollak P, Benazzouz A, Ardouin C, Hoffmann D, Benabid AL (1998). "Electrical stimulation of the subthalamic nucleus in advanced Parkinson's disease." *NEJM* 339(16):1105-11. | **9770557** | Seminal STN-DBS efficacy: off-medication UPDRS-III improved 60% at 1yr (P<0.001); L-DOPA dose halved. n=24. |
| 17 | Little S, Pogosyan A, Neal S, et al., Brown P (2013). "Adaptive deep brain stimulation in advanced Parkinson disease." *Ann Neurol* 74(3):449-57. | **23852650**, PMC3886292 | Closed-loop beta-triggered aDBS: UPDRS improved 66%(unblinded)/50%(blinded) vs no-stim, 29%/27% BETTER than continuous DBS, with 56% LESS stimulation time — beta is a usable causal control variable, not epiphenomenal. n=8. |
| 18 | Holgado AJ [Nevado Holgado AJ], Terry JR, Bogacz R (2010). "Conditions for the generation of beta oscillations in the subthalamic nucleus-globus pallidus network." *J Neurosci* 30(37):12340-52. | **20844130**, PMC6633459 | The established STN-GPe-loop mechanism class this build's Part F implements; own 3 conditions (loop-gain, delay-dominance, cortical/striatal-drive-balance) match this build's geometry. Full parameter table copyright-blocked even via its own open PMC record (`pmc-prop-open-access:no`) — confirmed structurally, not abandoned after one try. PubMed indexes the first author as "Holgado AJ" (drops "Nevado-"), which is *why* author-based searches for "Nevado-Holgado AJ[Author]" initially returned zero hits — resolved via title-only search. |
| 19 | Richfield EK, Penney JB, Young AB (1989). "Anatomical and affinity state comparisons between dopamine D1 and D2 receptors in the rat central nervous system." *Neuroscience* 30(3):767-77. | **2528080**, DOI 10.1016/0306-4522(89)90168-1 | Mechanistic crux: D1 receptor primarily LOW-affinity-agonist-state (RH=21±6%), D2 primarily HIGH-affinity (RH=77±3%) — the ~3.67x ratio sets this build's EC50 asymmetry. |
| 20 | Frank MJ (2005). "Dynamic dopamine modulation in the basal ganglia: a neurocomputational account of cognitive deficits in medicated and nonmedicated Parkinsonism." *J Cogn Neurosci* 17(1):51-72. | **15701239**, DOI 10.1162/0898929052880093 | Prior-art computational precedent (independent architecture from this build): phasic DA dynamically sets the BG go/no-go threshold on a cortical command. |
| 21 | Frank MJ, Seeberger LC, O'Reilly RC (2004). "By carrot or by stick: cognitive reinforcement learning in parkinsonism." *Science* 306(5703):1940-3. | **15528409**, DOI 10.1126/science.1102941 | The SEED-DESIGN node's named anchor (a): PD OFF-medication learn better from negative than positive outcomes; L-DOPA reverses the bias. NOT executed this session (a learning-bias claim, decorrelated from this build's movement-gating claim) — disclosed scope boundary. |
| 22 | Schultz W, Dayan P, Montague PR (1997). "A neural substrate of prediction and reward." *Science* 275(5306):1593-9. | **9054347**, DOI 10.1126/science.275.5306.1593 | Classic RPE/dopamine-teaching-signal context (couples to reward/dopamine node); not itself re-derived here. |
| 23 | Fearnley JM, Lees AJ (1991). "Ageing and Parkinson's disease: substantia nigra regional selectivity." *Brain* 114(Pt 5):2283-301. | **1933245**, DOI 10.1093/brain/114.5.2283 | **Re-verified live this session specifically to adjudicate a scout-vs-sibling discrepancy** (see provenance note above). Verbatim: "68% cell loss in the lateral ventral tier and a 48% loss in the caudal nigra as a whole" at symptom onset; 45%/decade PD loss-rate vs 4.7%/decade aging. No "31%"/logistic-regression figure anywhere in the abstract. |

**Read-only reused, not independently re-verified this session** (already live-verified by
the sibling `MECHANISM_DOPAMINE_KINETICS.md` cert): Kordower et al 2013 (PMID 23884810,
50-90% terminal DA-marker loss) and Cheng et al 2010 (PMID 20517933, ~80% striatal DA
loss) — both feed this build's PD severity band only via the sibling's own JSON, not
re-fetched here.

## 2. Method — geometric, machine-checked (not narrated)

**Part A/B (static push-pull gate).** Two Hill-function channels of dopamine tone DA:
`d1(DA) = DA^n/(EC50_D1^n+DA^n)` (increasing — D1 is Gs-coupled/facilitating, but
LOW-affinity, so needs high/phasic DA to engage) and `d2(DA) = 1 − DA^n/(EC50_D2^n+DA^n)`
(decreasing — D2 is Gi-coupled/inhibitory on the indirect MSN, but HIGH-affinity, so it is
substantially engaged even at low tonic DA, and the indirect MSN is disinhibited as DA
falls). The EC50 **ratio** (EC50_D1/EC50_D2 = 3.67) is fixed by Richfield 1989's own
affinity-state numbers; the **absolute** scale (EC50_D2=0.75) is a disclosed free choice
picked so the PD and HD effect magnitudes land in a comparable order of magnitude (not
tuned toward any sign outcome — the ratio, which the sign/direction actually depends on,
is the only citation-anchored piece). GPi output is `GPi_norm = 1 − k1·(d1−d1_ref) +
k2·(d2−d2_ref)` with **k1=k2=1** (symmetric unit weights — any asymmetry in the OUTCOME
must come from the drive changes, not a hand-tuned gain), and movement gain (the thalamic
disinhibition reaching cortex) is `gain = 1/GPi_norm`. PD lesions the DA-tone axis; HD
lesions the survival axis of each Hill channel independently (direct "relatively spared",
indirect "much more affected", per Reiner 1988's own words). **This is a rank argument**:
the double dissociation requires the (DA-lesion, survival-lesion) → (ΔGPi_PD, ΔGPi_HD) map
to be rank-2 (two channels, oppositely signed); collapsing to rank-1 (§3) provably cannot
span an opposite-sign target.

**Part F (STN-GPe delay loop, the spectral/geometric falsifier).** Linearized STN (x₁) and
GPe (x₂) populations, GPe→STN inhibitory (weight w_ge, delay d₂) and STN→GPe excitatory
(weight w_eg, delay d₁), each with its own integration time constant τ. Substituting
exponential trial solutions gives the closed-loop characteristic equation
`(τ₁s+1)(τ₂s+1) + G·e^(−sD) = 0` where `G=w_eg·w_ge` (the loop-gain product) and
`D=d₁+d₂` (round-trip delay) — exactly the form Nevado-Holgado/Terry/Bogacz 2010's own
three stated conditions describe (loop gain; delay-dominance; drive balance). The marginal
(Hopf) boundary is solved EXACTLY via `s=iω`, and the dominant root is CONTINUED by Newton
iteration across a gain sweep spanning the boundary — a genuine eigenvalue/spectral
derivation, not a heuristic. G is driven by the SAME `d2(DA)` function from Part A/B (no
second free mechanism introduced for this falsifier).

## 3. The forced adversary — single-channel null, swept across both sign branches

The most natural, non-strawman "collapsed" alternative: BOTH lesions reduce ONE aggregate
striatal-output variable multiplicatively (`striatal_output = DA_tone × cell_survival`),
driving GPi with one signed weight w (`GPi_norm = 1 − w·(striatal_output − 1)`) — the
tempting, reasonable-sounding intuition that "dopamine loss" and "neuron loss" are both
just instances of "the striatum does less." Steelmanned by testing **both** algebraic
branches of its one free parameter (w=+1, the anatomically-correct sign since
striatal projections are GABAergic/inhibitory; w=−1, checked for completeness):

| w | sign(ΔGPi_PD) | sign(ΔGPi_HD) | same sign (structural failure)? |
|---|---|---|---|
| +1.0 | +1 | +1 | **YES — fails** |
| −1.0 | −1 | −1 | **YES — fails** |

**Both branches fail identically** — no choice of the null's one free parameter rescues
it, because both lesions enter the SAME channel with the SAME sign by construction. This
is the forced, machine-checked kill: `gate_adversary_forced_and_falls_both_sign_branches =
True`.

## 4. Gates — 15/15 pre-registered PASS (two carry an honesty asterisk)

```
ab_pd_sign_robust_positive:                    PASS (ΔGPi∈[+0.398,+0.738] across full 9-pt PD band)
ab_hd_sign_robust_negative:                    PASS (ΔGPi∈[-0.324,-0.162] across full 20-pt HD grid)
ab_double_dissociation:                        PASS (opposite sign, robust, not a single point)
c_adversary_forced_and_falls_both_sign_branches: PASS (w=+1 AND w=-1 both same-sign-fail)
d_kravitz_sign_match:                          PASS (D1-stim<baseline<D2-stim, matches Kravitz sign)
d_kravitz_ratio_order_of_magnitude:             PASS† (model 7.10x, measured 18.84x -- both >=3x)
e_stn_direction_matches_bergman1994:            PASS (model D2-drive up + Bergman 19->26Hz, +36.8%)
e_gpi_pd_direction_matches_filion_tremblay1991:  PASS (both GPi-up in PD/MPTP)
e_gpi_hd_direction_matches_hamada_delong1992:    PASS‡ (model GPi-down + Hamada GPi -32.1%/GPe -35.5%)
f_primary_in_beta_band:                        PASS (27.68 Hz in [13,30])
f_sensitivity_grid_majority_in_beta:            PASS (13/16 = 81% of physiologically-motivated grid)
f_crosses_within_pd_band:                      PASS (crosses G* at exactly 50% DA loss)
f_re_s_monotonic_and_crosses_zero:              PASS (Re(s): -34.5 -> +28.4 monotonically over sweep)
f_open_loop_never_oscillates:                  PASS (roots exactly -100,-100, zero imaginary part)
f_frequency_stable_post_threshold:              PASS (100% of post-threshold sweep stays in [11,32]Hz)
```

**† Kravitz ratio asterisk**: the SIGN and the qualitative bidirectional-swing structure are
decisively matched (D2-stim model=154.0% vs measured=162% — remarkably close; D1-stim
model=21.7% vs measured=8.6% — same order of magnitude, undershoots the suppression). The
model **ratio** (7.10x) undershoots the measured ratio (18.84x) by ~2.6x — the coarse
2-channel abstraction captures the bidirectional STRUCTURE, not a quantitative fit to the
exact optogenetic magnitude. Treat the sign+order-of-magnitude as decisive, the exact
ratio as a hypothesis, per this repo's own asterisk convention.

**‡ Hamada&DeLong cross-check asterisk**: Hamada & DeLong 1992 lesioned the STN directly (a
hemiballismus/hyperkinesia model), **not** striatal indirect-MSNs (the actual HD lesion
site) — this is a same-node-different-lesion-site STRUCTURAL ANALOGY (both collapse net
excitatory drive reaching GPi via the indirect/STN relay), not a claim that HD and
STN-lesion hemiballismus are the same disease. DeLong 1990's own abstract textually groups
them as the two hyperkinetic-pole examples, supporting the shared-mechanism-class reading.

## 5. couples_to — this task's named targets, cross-checked against the graph node's own list

This task named: the **athlete-movement axis** (motor-unit/weightlifting/sprint), the
**dopamine-kinetics cert**, the **motor-cortex/corticospinal cert** (in flight), the
**compute=mind substrate**, and the **function↔dysfunction axis**. The pre-existing
SEED-DESIGN graph node independently named: reward/dopamine cell, gait cell, motor-control
cell, procedural/habit-memory cell — these substantially **agree** (dopamine-kinetics ~
reward/dopamine cell; motor-cortex/corticospinal ~ motor-control cell), a convergent
signal that this is a correctly-scoped coupling set, not an invented one.

- **dopamine-kinetics** (`docs/MECHANISM_DOPAMINE_KINETICS.md`): a REAL, live, on-disk
  structural coupling, not a prose pointer — `basal_ganglia_gating.py` `json.load`s
  `data/dopamine_kinetics/dopamine_kinetics_results.json` at runtime and reuses its own
  verified `[48,68]%`/`80%` PD-severity numbers as this build's own DA-loss axis (falls
  back to a disclosed literature band if that file is absent).
- **athlete-movement axis (motor-unit/weightlifting/sprint)**: this build's own computed
  movement-gain multiplier — **0.575–0.716× normal in PD** (bradykinesia-direction) vs
  **1.194–1.479× normal in HD** (chorea-direction) — is a genuine NEW number, directly
  usable as a multiplicative scale factor on whatever peak force/recruitment a downstream
  corticospinal-command-to-motor-unit model computes (`scripts/msk/motor_unit_recruitment.py`,
  `scripts/msk/weightlifting_triple_extension_model.py`), not yet wired in (that integration
  is a follow-on step, disclosed, not performed this session — no edit to those files).
- **motor-cortex/corticospinal cert (in flight)**: read-only checked
  (`data/body_twin/agent_outputs/motor-cortex-corticospinal__ac1205fe1a04450a1.json`, a
  scout output, itself names "basal ganglia (thalamocortical loops gating SMA/M1 action
  selection & initiation timing)" as one of its own couplings) — mutually confirms the
  coupling direction; that cert is not yet built, so this remains a **prose-level**
  coupling until it lands.
- **function↔dysfunction axis**: this build's own PD/HD falsifier IS the
  function↔dysfunction pair — the gated model reproduces both dysfunction poles from the
  SAME functional circuit under two different lesions; the forced adversary (§3) cannot.
- **compute=mind substrate**: conceptual pointer only, no computed coupling attempted
  (action-selection is named in the task as a canonical neural computation on the
  inherited compute boundary; this build supplies one concrete instance of it).

## 6. Symmetric QC — honest gaps, held OPEN, not smoothed over

- **The hyperdirect pathway (cortex→STN, Nambu's fast global "stop") is NOT separately
  modeled as a third channel** in the Part A/B steady-state gate — it appears only
  implicitly as a constant background drive term in Part F's cortical-input condition. A
  full 3-channel model, and specifically the hyperdirect pathway's temporal-precedence role
  in fast global movement-stopping (stop-signal-task framework), is out of scope this
  session — disclosed, not silently substituted. The classic Albin/DeLong rate model itself
  does not require a hyperdirect channel for the double dissociation either, so this is a
  correctly-scoped, not an under-scoped, omission for THIS falsifier.
- **The absolute EC50 scale (not the ratio) is a disclosed free/illustrative choice** — see
  §2. Without it the PD:HD effect-magnitude ratio can range from ~1.6x to ~24x depending on
  the arbitrary absolute placement (tested directly, both extremes still preserve the sign
  structure — the double-dissociation SIGN result is robust to this choice; the specific
  magnitudes reported are not).
- **The HD survival bands ([10,50]% indirect, [85,100]% direct) are illustrative**, not
  extracted from an exact published percentage — Reiner 1988's own PMC record is a
  pre-digital SCANNED article (full text not machine-readable this session); Albin 1992's
  abstract is qualitative-only ("evident loss" vs "relatively spared"). Same disclosed
  limitation the sibling dopamine-kinetics cert already hit for equally old founding papers.
- **The STN-GPe τ/delay values (Part F) are illustrative, not extracted** from
  Nevado-Holgado/Terry/Bogacz 2010's own parameter table — confirmed structurally
  copyright-blocked this session (`pmc-prop-open-access:no` on the direct efetch, despite a
  PMC id existing). The GEOMETRY (Hopf-boundary derivation, root continuation) is exact;
  only the plugged-in τ/delay magnitudes are a disclosed free choice, restricted to a
  physiologically-motivated neighborhood fixed BEFORE computing the sensitivity-grid
  outcome (81% in-band), not shopped after seeing results.
- **The DA-loss→loop-gain G(DA) linkage (Part F) is this build's own illustrative
  construction**, not an independent measurement — it reuses the SAME D2-drive function as
  Part A/B for internal consistency (one mechanism, not two free-floating ones), but the
  specific functional FORM linking indirect-pathway tone to STN-GPe loop gain is not itself
  citation-sourced.
- **Kravitz 2010's exact ratio is undershot by ~2.6x** (§4†) — sign and qualitative
  bidirectional swing are decisive, the magnitude is not.
- **Part G (late-HD flips back toward PD-like) is an internal structural prediction, NOT
  externally verified this session** — Albin et al 1990 Ann Neurol's "rigid HD"
  characterization (the paper that would anchor this) was not fetched. In the specific
  illustrative parameters tested (D1 survival 35%, D2 survival 15%), GPi_norm moves from
  0.754 toward 0.770 (attenuates the hyperkinetic deviation) but does **not** cross back
  above 1.0 with these EC50 parameters — reported exactly, not rounded up to "flips sign."
- **Frank et al 2004's ON/OFF-medication Go/NoGo learning-bias anchor (the SEED-DESIGN
  node's anchor (a)) is verified as a citation but NOT executed as a falsifier this
  session** — it is a cognitive/reinforcement-learning claim, decorrelated from (not the
  same measurement as) this build's movement-gating claim; a genuine follow-on, not
  smuggled in as if already covered.
- **Filion & Tremblay 1991's exact Hz numbers are not extractable** from its abstract (only
  qualitative direction: GPi up, GPe down) — used as a direction-only corroboration of
  Bergman 1994, not a second magnitude.
- **This is a population-level, not subject2-specific, build** — no subject2 EMG/kinematic
  panel exists to certify the movement-gain multiplier against; same disclosed scope as
  every sibling neuro-circuit layer in this repo.

## 7. Confidence tier

**Mixed-tier, per component, stated precisely (not smoothed into one label):**
electrophysiology-anchored for the firing-rate direction checks (in-vivo MPTP-primate and
STN-lesion primate single-unit recordings: Bergman 1994, Filion & Tremblay 1991, Hamada &
DeLong 1992); optogenetic-causal-anchored for the Kravitz sign/order-of-magnitude check
(in-vivo transgenic mouse); human-electrophysiology-and-clinical-anchored for the beta/DBS
chain (Brown 2001, Kühn 2006/2008, Little 2013, Krack 2003, Limousin 1998); and
illustrative-toy for the specific EC50 absolute scale, HD survival percentages, and STN-GPe
τ/delay magnitudes (disclosed throughout, same discipline as `dopamine_kinetics.py`'s own
Part 2/3 toys).

## Files

- `scripts/msk/basal_ganglia_gating.py` — self-contained (numpy/scipy only), builds the
  static push-pull gate, the forced adversary (both sign branches), the Kravitz check, the
  3 firing-rate direction cross-checks, the late-HD disclosed-not-verified prediction, and
  the full STN-GPe delay-loop spectral falsifier (Hopf boundary + root continuation + open-
  loop null); writes the JSON below; prints all 15 gates + overall_pass.
- `data/basal_ganglia_gating/basal_ganglia_gating_results.json` — every number in this doc,
  machine-written. Verified deterministic (2 independent runs, byte-identical via `diff`)
  and NaN/Inf-free (checked programmatically over the full JSON tree, 0 found).
- Read-only input: `data/dopamine_kinetics/dopamine_kinetics_results.json` (PD severity
  band, live structural coupling, not modified).
- Read but NOT modified (isolation): `data/MECHANISM_ANCHOR_GRAPH.json` (the 2 pre-existing
  SEED-DESIGN nodes discussed above), `data/body_twin/agent_outputs/basal-ganglia-gating__a93bab9e997bceb06.json`
  (the prior scout output — its 31% Fearnley-&-Lees claim was checked live and found not to
  exist in the source; its Mink/Kravitz/Reiner/Frank/Schultz citations were independently
  re-verified rather than trusted, and all confirmed accurate),
  `data/body_twin/agent_outputs/motor-cortex-corticospinal__ac1205fe1a04450a1.json` (in-flight
  coupling check).

## Repro

```
cd ~/projects/bodytwin
source .venv-msk/bin/activate
python3 scripts/msk/basal_ganglia_gating.py
```
Pure numpy/scipy (root-finding + Newton continuation over the delay-loop characteristic
equation), no external data dependency beyond an optional read-only
`dopamine_kinetics_results.json` (degrades to a disclosed literature-band fallback if
absent). <5s wall time, deterministic, no git operations. Writes only under
`data/basal_ganglia_gating/`.
