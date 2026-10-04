# MECHANISM AMYGDALA FEAR-CONDITIONING CIRCUIT — dual-pathway timing, the associative-learning gate, lesion necessity, extinction-not-erasure, PTSD dysfunction pole (2026-07-22)

Builds and MEASURES a quantitative model of the amygdala threat/fear-conditioning
circuit: the dual-pathway timing (thalamic "low road" direct MG→LA vs cortical
"high road" indirect MG→AuCx→LA, LeDoux lineage) unified with the associative-
learning gate (an LTP-consistent delta rule at LA synapses) via a single
geometric mechanism (NMDA receptor slow kinetics as a temporal coincidence
window); the lesion-necessity falsifier (forcing a "conditioning is amygdala-
independent" adversary to its strongest fair form); the extinction-is-new-
learning-not-erasure falsifier (renewal/reinstatement/spontaneous-recovery/
savings); and the PTSD dysfunction pole (impaired extinction *recall*, not
acquisition — matching what was actually measured, not the common
misconception). Script: `scripts/msk/amygdala_fear_circuit.py`. Evidence:
`data/amygdala_fear_circuit/amygdala_fear_circuit_results.json`.

**Relation to a pre-existing graph node (not edited this session — isolation
rule "touch only files you create"; `data/MECHANISM_ANCHOR_GRAPH.json`, node
`NEU-AMYGDALA-THREAT`)**: that node is **SEED-DESIGN** status — a rich cited
literature synthesis (LeDoux 2000, Rogan/Staubli/LeDoux 1997, Milad & Quirk
2002, Ciocchi et al 2010, Namburi et al 2015, the S.M. patient literature,
LeDoux & Brown's higher-order-feeling theory), explicitly marked "cited, not
executed against data" by its own `cert_design.verify.note`, and its own
`proposed_cell` explicitly calls for building a reduced-order state-space cell
with a falsifier test that was **not built**. This document independently
re-verifies (not re-uses/trusts) 4 PMIDs that node also cites (LeDoux 2000,
Rogan 1997, Milad & Quirk 2002, Ciocchi 2010 — all re-fetched live this
session), adds **19 new PMIDs** that node does not cite (the NMDA-block
falsifier, the LA/CeA lesion-specificity controls, the cued-vs-contextual
double dissociation, the human triple dissociation, the thalamic-vs-cortical
dual-route lesion dissociation, the <15ms latency + NMDA-kinetics chain, the
extinction reinstatement/renewal literature, and the PTSD imaging anchors) —
and, crucially, **actually builds and runs** the missing executed model that
node's own proposed_cell calls for (via a different, independently-anchored
falsifier design — that node's own proposed test requires digitizing an fMRI
BOLD-SCR co-decay curve from a specific figure, which this session did not
attempt, consistent with the "never eyeball a figure, measure into
physics-chosen numbers" discipline). No graph-edge write performed this
session (isolation).

## The falsifiers, verdicts stated up front (nothing hidden)

> **Falsifier 1 (dual-pathway timing + geometric unification)**: does the
> circuit reproduce the measured thalamic-fast/cortical-slow timing
> asymmetry, and does the thalamic-vs-cortical LESION dissociation emerge
> from pure graph topology (not a fitted parameter)?

**YES**, with one number pinned to a primary source and one deliberately
disclosed as unpinned. Quirk, Repa & LeDoux 1995 (n, freely-moving rats,
parallel LA single-unit recording): fear conditioning's effect on tone-evoked
LA responses is "greatest on the shortest latency (**less than 15 ms**)
components" — landing inside the task's own named 12–15 ms band. Li, Stutzmann
& LeDoux 1996 independently confirms the *ordering* (thalamic "rapid but
relatively impoverished," cortical "slower but richer") via a different
method (in vivo intracellular/extracellular recording), but its own abstract
does not state the cortical route's exact latency in ms — disclosed, not
fabricated. The thalamic-vs-cortical **lesion** dissociation (cortex-lesion
spares, thalamus/MG-lesion abolishes, either pathway alone is sufficient,
combined lesion disrupts) is reproduced from **pure graph reachability** on a
3-node DAG (MG is the shared upstream source feeding BOTH the direct MG→LA
edge and the indirect MG→AuCx→LA edge) — zero free parameters — and a forced
**adversary topology** (cortex has its own independent, non-MG-relayed input)
is shown to **fail** against the real data (it would incorrectly predict that
an MG lesion spares conditioning). Machine-verified: `gate_true_topology_
matches_real_pattern=True`, `gate_adversary_topology_fails_to_match_real_
data=True`, `gate_equipotentiality_each_pathway_alone_sufficient=True`.

> **Falsifier 2 (the associative-learning gate)**: does CR magnitude grow
> with CS–US pairings (with a genuine ZERO floor for unpaired presentations),
> and does NMDA blockade reproduce Miserendino et al 1990's own double
> dissociation (blocks acquisition, not expression)?

**YES.** A delta-rule (Rescorla–Wagner-form) LTP gate run for 20 trials
produces a monotonic, **saturating/negatively-accelerating** curve (both
properties machine-checked on the raw trace, not eyeballed) reaching 95.1% of
ceiling by trial 7 — a free/uncalibrated rate constant, reported as a
qualitative shape match to Quirk/Repa/LeDoux 1995's own finding that
conditioning-induced change occurs "often within the first several trials,"
**not** a fitted numeric prediction. The **unpaired void-floor** control (20
unpaired CS/US presentations) produces **exactly zero** growth — machine-
matching Rogan, Stäubli & LeDoux 1997's own verbatim finding that
potentiation "do[es] not occur if the CS and US remain unpaired." **NMDA
block** (η=0): applied to a naive animal, acquisition is fully abolished
(W stays at 0.0); applied to an *already-conditioned* animal, the existing
weight is untouched and CR is unaffected — machine-reproducing Miserendino et
al 1990's own two-condition double dissociation ("blocks acquisition, but not
expression") exactly, not just a single directional match.

> **Falsifier 3 (lesion necessity + the forced "amygdala-independent"
> adversary)**: does LA/CeA lesion abolish CR while adjacent-tissue and
> other-amygdala-nuclei "sham" lesions do not — and does the STRONGEST fair
> form of an amygdala-independent adversary (a different medial-temporal
> structure substituting) still fail?

**YES**, on both the anatomical-specificity control and the forced
adversary. LeDoux, Cicchetti, Xagoraris & Romanski 1990: LA lesion abolishes
conditioning; lesions of the **striatum above or cortex adjacent to** LA had
**no effect**. Nader, Majidishad, Amorapanth & LeDoux 2001 sharpens this
within the amygdala itself: LA **or** CeA lesion abolishes; lesions of the
**basal, accessory-basal, or medial nuclei had little effect**. The forced
adversary is **not a strawman**: Phillips & LeDoux 1992 shows hippocampal
lesion **spares CUED conditioning** (while abolishing contextual) — a real
candidate for "maybe a different structure substitutes." Machine-checked
against the model: hippocampal "lesion" leaves cued CR at ceiling (=1.0,
matching control) — the adversary **fails** for the cued case specifically,
because amygdala lesion (not hippocampal lesion) is the one that abolishes
it. Bechara et al 1995's human triple dissociation (amygdala-lesion patient:
conditioning lost, declarative facts intact; hippocampus-lesion patient: the
reverse; both-lesioned patient: neither) is reproduced as an exact structural
isomorphism, cross-checked field-by-field against the model. **Scope
boundary, disclosed not hidden**: a genuinely DIFFERENT, amygdala-independent
threat response exists for INTEROCEPTIVE stimuli (35% CO₂-panic surviving
amygdala loss, per this repo's own prior `NEU-AMYGDALA-THREAT` node) — a
different stimulus category, not a counterexample to the claim as scoped
here (exteroceptive cued/contextual Pavlovian conditioning).

> **Falsifier 4 (extinction = new learning, NOT erasure)**: does an
> "erasure" adversary (extinction directly decrements the original trace)
> fail to reproduce renewal, reinstatement, spontaneous recovery, AND
> re-extinction savings, while a "new inhibitory learning" model reproduces
> all four from ONE mechanism?

**YES — and this is the decisive, most heavily forced part of this build**,
including a real bug caught and fixed by the model's own robustness sweep,
not smoothed over. The erasure adversary (a real, historically serious
position — Quirk 2002's own abstract frames the debate exactly this way: "it
has been hypothesized that extinction does not erase conditioning, but forms
a new memory") is forced computationally: extinction trials directly
decrement W_LA toward zero. Under this model, ALL FOUR test manipulations
(context switch/renewal, US-alone/reinstatement, time-passage/spontaneous-
recovery, re-extinction/savings) correctly produce **no effect** — there is
nothing left to renew, reinstate, or recover, and re-extinguishing a
re-acquired trace takes exactly as long as the original (machine-verified:
`gate_erasure_adversary_falls_all_four_absent=True`). The "new inhibitory
learning" model (a separate, context-bound trace W_ext, gated off outside the
extinction context and overridden by a US-alone presentation, decaying faster
than W_LA over time) reproduces **all four**
(`gate_new_learning_all_four_present=True`), externally anchored to Bouton &
Bolles 1979 (reinstatement, context-gated), Bouton & King 1983 (ABA renewal),
and Quirk 2002's own directly measured numbers: conditioned freezing
"gradually recovered with time to reach **100% by day 10**" AND, despite
complete behavioral recovery, "rats showed **savings** in their rate of
re-extinction" — proving the extinction trace persists *underneath* a
behaviorally-complete-looking return of fear. Milad & Quirk 2002 supplies the
causal mechanism: infralimbic (mPFC) neurons fire to the CS *only* during
successful extinction recall, and *electrically stimulating* IL alone (no
extinction training) is *sufficient* to lower freezing.

**The forced bug** (OODA, not hidden): the first version of the savings check
used only an integer "trials-to-criterion" count, which **failed** at 2 of 7
head-start fractions swept (0.05, 0.1) — not because savings was absent, but
because a coarse integer discretization can mask a small continuous effect.
Diagnosed: the underlying recurrence `w_ext(t+1) = (1−η)·w_ext(t) + η·W_LA` is
an **affine contraction** — for any two starting points, their difference at
time t is `(1−η)^t` times their initial difference, which **never changes
sign** for η∈(0,1). Ordering is therefore preserved *forever*, a provable
geometric property, not a curve-fit. Fixed with a continuous monotonicity
measure (is the head-started trajectory strictly ahead of the from-zero
trajectory at *every* trial along the horizon) — robust across **all 7/7**
sweep points (`gate_savings_robust_by_continuous_monotonicity=True`), while
the cruder integer read is honestly reported alongside as
discretization-sensitive, not deleted.

> **Falsifier 5 (PTSD dysfunction pole)**: does a PTSD-like parameter shift
> reproduce the measured day-1-normal / day-2-impaired pattern (Milad et al
> 2009), without falling into the common misconception that PTSD = simply
> "more/faster conditioning" (which the real data itself refutes)?

**YES**, and getting this right required NOT doing the naive thing. Milad et
al 2009 (n=16 PTSD, n=15 trauma-exposed non-PTSD controls, 2-day fMRI+SCR
protocol): "SCR data revealed **no significant differences between groups**
during acquisition and extinction" on day 1 — so the model deliberately holds
the acquisition rate η **identical** across groups
(`gate_no_day1_acquisition_extinction_difference=True`, exactly 0 difference
by construction — a disclosed *modeling choice* reflecting the real finding,
not a discovery). On day 2 (extinction recall), the model applies a
*reduced* mPFC/vmPFC extinction-recall gate efficacy for the PTSD-like case
(0.35 vs 0.85, illustrative contrast, disclosed as not independently fit to
an exact reported effect size) — producing net residual fear of **0.65 vs
0.15** (PTSD-like vs normal), directionally matching Milad et al 2009's own
measured pattern: "PTSD subjects showed impaired recall of extinction
memory," with greater amygdala activation during day-1 *learning* and lesser
hippocampus/vmPFC + greater dACC activation during day-2 *recall*. Rauch,
Shin & Phelps 2006's review-level model corroborates the general shape
("exaggerated amygdala responses... and deficient frontal cortical
function... deficient hippocampal function").

**22/22 machine-computed gates PASS** (`overall_pass=True`) — §5. 2
independent runs produce byte-identical JSON (verified via `diff`); zero
NaN/Inf anywhere in the output tree (checked programmatically, inside the
script's own `main()`). **Not all 22 gates carry equal evidentiary weight —
see the symmetric self-grading in §4** before reading the gate list as "22
independent discoveries."

## 1. Citations — every PMID below verified LIVE this session (NCBI
E-utilities: esearch → esummary/efetch, direct curl, raw MEDLINE abstract
text fetched and read — not recalled from memory, not WebFetch-summarized;
this repo's own prior finding across sibling docs is a measured ~62–67%
citation-drift rate from memory alone). 27 PMIDs machine-cross-checked
against the raw fetched records; 1 classic/textbook-tier reference confirmed
to have **0 PubMed hits** live (a book chapter, not journal-indexed — the
same disclosed treatment this repo gives Bernheimer 1973 / Nagatsu 1964).

| # | Citation | PMID | Role |
|---|---|---|---|
| 1 | LeDoux JE (2000). Emotion circuits in the brain. *Annu Rev Neurosci* 23:155-84. | **10845062** | Circuit topology review anchor. |
| 2 | Romanski LM, LeDoux JE (1992). Equipotentiality of thalamo-amygdala and thalamo-cortico-amygdala circuits in auditory fear conditioning. *J Neurosci* 12(11):4501-9. | **1331362** | THE dual-route equipotentiality/redundancy anchor: "destruction of either pathway alone had no effect... combined lesions... disrupted conditioning." |
| 3 | LeDoux JE, Sakaguchi A, Reis DJ (1984). Subcortical efferent projections of the medial geniculate nucleus mediate emotional responses conditioned to acoustic stimuli. *J Neurosci* 4(3):683-98. | **6707732** | MG(thalamus) lesion abolishes both autonomic+behavioral CR; auditory-cortex lesion has NO effect; inferior-colliculus lesion replicates MG-lesion effect. |
| 4 | Li XF, Stutzmann GE, LeDoux JE (1996). Convergent but temporally separated inputs to lateral amygdala neurons from the auditory thalamus and auditory cortex use different postsynaptic receptors. *Learn Mem* 3(2-3):229-42. | **10456093** | Dual-pathway timing (qualitative: thalamic "rapid... impoverished," cortical "slower... richer"); NMDA specifically involved in thalamic (not cortical) transmission; proposes the coincidence-integration hypothesis this doc's Part 1 builds on. |
| 5 | Quirk GJ, Repa C, LeDoux JE (1995). Fear conditioning enhances short-latency auditory responses of lateral amygdala neurons. *Neuron* 15(5):1029-39. | **7576647** | THE <15ms fast-component latency anchor; conditioning's plasticity effect concentrated specifically there. |
| 6 | Weisskopf MG, LeDoux JE (1999). Distinct populations of NMDA receptors at subcortical and cortical inputs to principal cells of the lateral amygdala. *J Neurophysiol* 81(2):930-4. | **10036290** | Independent (slice-physiology) confirmation of pathway-specific NMDA receptor populations. |
| 7 | Weisskopf MG, Bauer EP, LeDoux JE (1999). L-type voltage-gated calcium channels mediate NMDA-independent associative LTP at thalamic input synapses to the amygdala. *J Neurosci* 19(23):10512-9. | **10575047** | DISCLOSED NUANCE: thalamic-input LTP itself is NMDA-**independent** in vitro (Ca-channel-dependent instead). |
| 8 | Lester RA, Clements JD, Westbrook GL, Jahr CE (1990). Channel kinetics determine the time course of NMDA receptor-mediated synaptic currents. *Nature* 346(6284):565-7. | **1974037** | NMDA EPSC decay = "several hundred milliseconds" — the coincidence-window anchor (hippocampal culture, disclosed cross-region use). |
| 9 | Miserendino MJD, Sananes CB, Melia KR, Davis M (1990). Blocking of acquisition but not expression of conditioned fear-potentiated startle by NMDA antagonists in the amygdala. *Nature* 345(6277):716-8. | **1972778** | THE NMDA-block-acquisition-not-expression falsifier. |
| 10 | Rogan MT, Stäubli UV, LeDoux JE (1997). Fear conditioning induces associative long-term potentiation in the amygdala. *Nature* 390(6660):604-7. | **9403688** | LA-LTP behavioral correlate; explicit verbatim null for unpaired CS/US. |
| 11 | McKernan MG, Shinnick-Gallagher P (1997). Fear conditioning induces a lasting potentiation of synaptic currents in vitro. *Nature* 390(6660):607-11. | **9403689** | In vitro companion; presynaptic AMPA facilitation. |
| 12 | LeDoux JE, Cicchetti P, Xagoraris A, Romanski LM (1990). The lateral amygdaloid nucleus: sensory interface of the amygdala in fear conditioning. *J Neurosci* 10(4):1062-9. | **2329367** | LA lesion abolishes; adjacent striatum/cortex lesion spares (anatomical-specificity control). |
| 13 | Nader K, Majidishad P, Amorapanth P, LeDoux JE (2001). Damage to the lateral and central, but not other, amygdaloid nuclei prevents the acquisition of auditory fear conditioning. *Learn Mem* 8(3):156-63. | **11390635** | LA+CeA necessary; basal/accessory-basal/medial nuclei spared — within-amygdala specificity. |
| 14 | Ciocchi S, Herry C, Grenier F, et al (2010). Encoding of conditioned fear in central amygdala inhibitory circuits. *Nature* 468(7321):277-82. | **21068837** | CEl required for acquisition; CEm output neurons drive conditioned freezing (independently re-verified live this session). |
| 15 | Bechara A, Tranel D, Damasio H, Adolphs R, Rockland C, Damasio AR (1995). Double dissociation of conditioning and declarative knowledge relative to the amygdala and hippocampus in humans. *Science* 269(5227):1115-8. | **7652558** | Human triple dissociation (amygdala-only / hippocampus-only / both-lesioned patients). |
| 16 | LaBar KS, LeDoux JE, Spencer DD, Phelps EA (1995). Impaired fear conditioning following unilateral temporal lobectomy in humans. *J Neurosci* 15(10):6846-55. | **7472442** | Corroborating human lesion cohort, larger than the n=1 case-study lineage. |
| 17 | Phillips RG, LeDoux JE (1992). Differential contribution of amygdala and hippocampus to cued and contextual fear conditioning. *Behav Neurosci* 106(2):274-85. | **1590953** | Amygdala necessary for BOTH cued+contextual; hippocampus only contextual — forces the amygdala-independent/hippocampal adversary. |
| 18 | Campeau S, Davis M (1995). Involvement of subcortical and cortical afferents to the lateral nucleus of the amygdala in fear conditioning measured with fear-potentiated startle. *J Neurosci* 15(3):2312-27. | **7891169** | Modality-specific double dissociation (auditory thalamus necessary, cortex not); a retraining-compensation nuance disclosed in §6. |
| 19 | Bouton ME, Bolles RC (1979). Role of conditioned contextual stimuli in reinstatement of extinguished fear. *J Exp Psychol Anim Behav Process* 5(4):368-78. | **528893** | Reinstatement, context-gated. |
| 20 | Bouton ME, King DA (1983). Contextual control of the extinction of conditioned fear. *J Exp Psychol Anim Behav Process* 9(3):248-65. | **6886630** | Renewal (ABA renews, Expts 1/3; reversed-order ABC shows NO renewal, Expt 4 — disclosed asymmetry). |
| 21 | Rescorla RA (2004). Spontaneous recovery. *Learn Mem* 11(5):501-9. | **15466300** | Review; honestly discloses that "the empirical literature on its determinants is relatively sparse and quite mixed." |
| 22 | Quirk GJ (2002). Memory for extinction of conditioned fear is long-lasting and persists following spontaneous recovery. *Learn Mem* 9(6):402-7. | **12464700** | THE decisive extinction-not-erasure anchor: 7 CS-US pairings/20 extinction trials protocol; 100% spontaneous recovery by day 10; savings in re-extinction rate. |
| 23 | Milad MR, Quirk GJ (2002). Neurons in medial prefrontal cortex signal memory for fear extinction. *Nature* 420(6911):70-4. | **12422216** | Infralimbic neurons fire during extinction recall only; IL stimulation alone is sufficient to lower freezing (causal, not just correlational). |
| 24 | Milad MR, Pitman RK, Ellis CB, et al (2009). Neurobiological basis of failure to recall extinction memory in PTSD. *Biol Psychiatry* 66(12):1075-82. | **19748076** | n=16 PTSD/15 control: NO day-1 SCR difference; day-2 impaired recall + amygdala hyperactivation (learning) + hippocampus/vmPFC hypoactivation + dACC hyperactivation (recall). |
| 25 | Rauch SL, Shin LM, Phelps EA (2006). Neurocircuitry models of posttraumatic stress disorder and extinction. *Biol Psychiatry* 60(4):376-82. | **16919525** | Review-level PTSD model: exaggerated amygdala + deficient frontal/hippocampal function. |
| 26 | Bordi F, LeDoux JE (1994). Response properties of single units in areas of rat auditory thalamus that project to the amygdala. I. *Exp Brain Res* 98(2):261-74. | **8050512** | Thalamic response-latency/tuning properties. |
| 27 | Bordi F, LeDoux JE (1994). ...II. Cells receiving convergent auditory and somatosensory inputs... *Exp Brain Res* 98(2):275-86. | **8050513** | CS(auditory)+US(somatosensory) convergence onto the SAME thalamic cells projecting to the amygdala. |
| 28 | Rescorla RA, Wagner AR (1972). A theory of Pavlovian conditioning. In: Black AH, Prokasy WF (eds), *Classical Conditioning II*. Appleton-Century-Crofts. | — (book chapter, confirmed 0 PubMed hits live) | Classic/theoretical tier — the canonical delta-rule formalism Part 2 instantiates. |

## 2. Geometric structure — dual-pathway timing and the associative-learning
gate as ONE mechanism, not two

**The graph argument (thalamic-vs-cortical lesion dissociation)**: model 3
nodes {MG (auditory thalamus), AuCx (auditory cortex), LA} as a DAG with edges
`MG→LA` (direct) and `MG→AuCx→LA` (indirect) — MG is the **shared upstream
source** feeding both. Pure BFS reachability from `input` to `LA` under 4
lesion conditions reproduces the real pattern **exactly**, with **zero free
parameters**:

| condition | LA reachable? | real data |
|---|---|---|
| no lesion | YES | spared |
| cortex (AuCx) lesion | YES (direct edge intact) | **spared** (LeDoux/Sakaguchi/Reis 1984; Campeau & Davis 1995) |
| thalamus (MG) lesion | **NO** (both edges lose their shared source) | **abolished** (same 2 sources) |
| combined lesion | NO | abolished |

This is not "thalamus > cortex" in any functional-priority sense — it is a
**topology** fact: MG sits upstream of *both* branches, so removing the
shared source node severs both, while removing only the downstream cortical
branch leaves the other branch's own independent source untouched. A forced
**adversary topology** (cortex has its own independent, non-MG-relayed
input) is shown to fail: under that topology, an MG lesion would incorrectly
leave LA reachable via cortex — contradicting the real data
(`gate_adversary_topology_fails_to_match_real_data=True`). Romanski & LeDoux
1992's more surgically refined (hemidecussated) design isolates each pathway
with its *own* intact upstream source and finds each *individually*
sufficient (equipotentiality) — reproduced here as two independent
single-pathway reachability checks, both True.

**The NMDA-coincidence-window unification**: Li, Stutzmann & LeDoux 1996
propose that "the slow time course of NMDA currents could provide LA cells
with a mechanism to integrate the inputs arriving rapidly from the thalamus
and somewhat later from the cortex." Lester, Clements, Westbrook & Jahr 1990
give the quantitative anchor: NMDA-receptor-mediated EPSCs decay over
"**several hundred milliseconds**" (conservatively used here as 300 ms) —**100
times longer** than the channel's own mean open time, and vastly longer than
AMPA's rapid onset/decay. Since the thalamic component is independently
pinned at **<15 ms** (Quirk/Repa/LeDoux 1995) and any plausible cortical-route
gap must be sub-second (both are components of one short-latency tone-evoked
response), the NMDA window (300 ms) comfortably exceeds even a deliberately
generous, disclosed-as-unpinned 200 ms gap ceiling (margin: **1.5×**,
machine-computed, not a huge margin but real and conservative-in-the-
adversarial-direction). A **void-floor check** confirms this is not a
vacuous test: if NMDA decayed like AMPA (~2 ms), the same test would
correctly **FAIL** (2 < 200) — the gate genuinely discriminates on the real
physiological difference between receptor kinetics, not on an
always-true tautology. **Disclosed nuance, not hidden**: Weisskopf, Bauer &
LeDoux 1999 found the thalamic-input synapse's own in-vitro LTP is
NMDA-**independent** (L-type Ca-channel-dependent) — so this coincidence-
window story is a proposed, plausible mechanism (the hypothesis's own
authors' words), not a fully closed, uncontested account.

## 3. The associative learning gate — delta-rule LTP, the void floor, and the
NMDA-block double dissociation

A single delta rule, `ΔW = η·(W_max − W)` on paired (CS+US) trials, `ΔW=0` on
unpaired trials — the canonical Rescorla–Wagner form (classic/theoretical
tier, §1 #28), instantiating the LTP mechanism Rogan/Staubli/LeDoux 1997 and
McKernan/Shinnick-Gallagher 1997 measured directly (auditory-evoked field
potentiation in vivo; presynaptic AMPA facilitation in vitro). Machine-checked
properties of the raw 20-trial trace (η=0.35, disclosed free/uncalibrated):
zero at trial 0, strictly monotonic increasing, strictly concave/negatively-
accelerating (the discrete second difference of the trace is ≤0 at every
step) — reaching 95.1% of ceiling by trial 7, a qualitative-shape match
(not a fitted number) to Quirk/Repa/LeDoux 1995's own finding that
conditioning's effect occurs "often within the first several trials" and to
Quirk 2002's use of 7 pairings as a standard protocol.

**The void floor**: 20 UNPAIRED CS/US presentations (the fair, non-degenerate
baseline) produce **exactly zero** growth — machine-matching Rogan et al
1997's own verbatim finding. Paired growth (≈0.9998 of ceiling) vastly
exceeds this true zero, not merely a smaller number.

**The NMDA-block double dissociation** (Miserendino et al 1990): applying η=0
to a **naive** animal during acquisition trials leaves W at exactly 0.0 — new
learning is fully abolished. Applying the SAME block to an **already-
conditioned** animal (W built up via normal η first) leaves the existing W,
and therefore CR, completely unaffected — expression of the prior memory is
spared. This is not a single directional finding; it is the SAME two-line
model reproducing BOTH halves of Miserendino's own reported dissociation.

## 4. Symmetric self-grading — not all 22 gates carry equal weight

Per this task's own "kills are auditable, an unjustified pass is as much a
failure as an unjustified kill" discipline, here is an honest breakdown
(full detail: `docs/MECHANISM_AMYGDALA_FEAR_CIRCUIT_evidence.json`
→ `gates_summary.symmetric_self_grading_of_gate_strength_not_all_equal_weight`):

- **Tier A — decisive, non-trivial falsifiers** (a genuine alternative
  structure was built and shown to fail, or a real bug was caught and fixed):
  the adversary-topology graph test (§2); the NMDA-coincidence void-floor
  check (§2); the NMDA-block double dissociation (§3); the hippocampal
  forced adversary (§ falsifier 3); all four extinction structural
  predictions plus the caught-and-fixed savings-monotonicity proof (§
  falsifier 4, the strongest result in this build).
- **Tier B — formalization/consistency checks**: most lesion-specificity
  and dissociation tables (they verify the model was built correctly to
  MATCH a real reported pattern, not an independent derivation of it); the
  acquisition curve's shape properties (mathematically guaranteed by the
  delta-rule equation for η∈(0,1)); the PTSD day-1-invariance gate (true
  by a deliberate, disclosed modeling choice, not a discovery).
- **Tier C — transcription-fidelity only**: the "15ms ≤ 15ms" check (confirms
  `PREREG` faithfully carries the cited number, nothing more).

None of this weakens the falsifiers in §Falsifier 1–5 above — it is a
disclosure of WHICH gates are doing real evidentiary work versus which are
housekeeping, so that "22/22 PASS" is read correctly rather than as 22
independent discoveries of equal strength.

## 5. Pre-registered gates — 22/22 PASS, machine-computed (not narrated)

```
part1_true_topology_matches_real_lesion_pattern:            PASS
part1_adversary_topology_fails_on_real_data:                 PASS
part1_equipotentiality_reproduced:                           PASS
part1_thalamic_latency_matches_primary_source_le15ms:        PASS (<=15ms)
part1_nmda_coincidence_window_exceeds_gap:                   PASS (300ms > 200ms, 1.5x margin)
part2_acquisition_zero_at_start:                             PASS
part2_acquisition_monotonic:                                 PASS
part2_acquisition_saturating_concave:                        PASS
part2_unpaired_void_floor_exact_zero:                        PASS
part2_paired_greatly_exceeds_unpaired:                       PASS
part2_nmda_block_abolishes_acquisition:                      PASS
part2_nmda_block_spares_expression:                          PASS
part3_lesion_specificity_matches_real_data:                  PASS
part3_hippocampal_adversary_fails_for_cued:                  PASS
part3_full_rat_double_dissociation_matches:                  PASS
part3_human_triple_dissociation_matches:                     PASS
part4_erasure_adversary_falls_all_four_absent:               PASS
part4_new_learning_all_four_present:                         PASS
part4_savings_robust_by_continuous_monotonicity:             PASS (7/7; integer-discretized version: 5/7, diagnosed+disclosed)
part5_no_day1_group_difference_by_construction:              PASS (0 diff, disclosed as deliberate choice)
part5_ptsd_impaired_extinction_recall:                       PASS (0.650 vs 0.150)
couples_hpa_timescale_computed:                              PASS
```

`overall_pass=True`. Determinism: 2 independent runs produce byte-identical
JSON (verified via `diff`); zero NaN/Inf anywhere in the output tree (checked
programmatically inside the script's own `main()`).

## 6. couples_to — HPA-cortisol (concrete re-derived number), compute=mind
substrate, dopamine/serotonin, function↔dysfunction

**HPA-cortisol** (a concrete, machine-computed number, not a prose pointer —
reads `hpa_cortisol_axis_results.json` READ-ONLY, not modified): the fast
subcortical amygdala threat-output route (thalamic input <15ms; the full
LA→CeA→PAG/hypothalamus effector chain padded to a disclosed, unpinned ≤100ms
ceiling) precedes the HPA axis's own independently-measured fastest
cortisol-detectable component (15 min = 900,000 ms, from that sibling doc's
own already-certified derived-latency chain) by a computed **60,000×** (bare
thalamic-input figure) to **9,000×** (padded full-circuit ceiling) — both
endpoints independently sourced, 3–4+ orders of magnitude, matching the
textbook picture that fast behavioral/autonomic defense precedes slower
endocrine mobilization.

**compute=mind substrate** (conceptual, not numeric): the delta-rule
associative-learning gate this doc builds and machine-verifies is a
canonical coincidence-detection + error-correction weight-update algorithm —
the same algorithmic class the broader "compute=mind" framing points to
generically. Here it is shown to be literally implemented in wetware via
NMDA-receptor voltage+ligand-gated coincidence detection (§2) — a concrete
instance of the general thesis, not just an assertion of it.

**dopamine/serotonin neuromodulation** — disclosed as **prose/citation tier
only**, not computed: no numeric coupling was built between amygdala
plasticity and `scripts/msk/dopamine_kinetics.py` / `serotonin_system.py`
this session. A real, well-known additional regulatory layer (dopaminergic
salience/prediction-error modulation of LA plasticity; serotonergic
anxiolytic modulation of amygdala reactivity) exists but is not quantified
here — an honest scope limit, matching `MECHANISM_HPA_CORTISOL.md`'s own
disclosed treatment of its immune coupling.

**function↔dysfunction axis** — BUILT, not just pointed to: §5/Falsifier 5's
PTSD dysfunction-pole model IS this coupling (reduced extinction-recall gate
efficacy, acquisition rate deliberately held constant to match Milad et al
2009's own no-day-1-difference finding).

## 7. Confidence tier

Per this task's own pre-registration and matching this repo's established
precedent (`MECHANISM_HPA_CORTISOL.md`, `MECHANISM_DOPAMINE_KINETICS.md`):
**rodent-lesion/electrophysiology + human-lesion-case + human-fMRI-imaging
anchored** (LeDoux/Davis/Quirk/Milad lineage — real rat single-unit
recording, n's per-study in the tens; real rat lesion cohorts; a small
(n=1–3) human lesion case lineage; a real n=16/15 human fMRI+SCR cohort) —
one tier below a subject-specific in-vivo measurement (no subject2
fear-conditioning/lesion/imaging panel exists in this repo to certify
against — the same disclosed scope every sibling neuro/endocrine doc already
carries).

## 8. Honest gaps — symmetric QC: what this does NOT prove

- **The cortical-route exact latency (ms) is not verbatim-pinned.** Li,
  Stutzmann & LeDoux 1996 gives only the qualitative ordering; a live `elink`
  check this session confirmed that paper has no PMC full-text deposit. The
  200ms gap ceiling used in §2 is explicitly an unpinned, deliberately
  generous (conservative-against-the-hypothesis) construction, not an
  independently-sourced number.
- **A genuine literature tension, forced and disclosed**: Weisskopf, Bauer &
  LeDoux 1999 found thalamic-input LTP itself is NMDA-**independent** in
  vitro, despite NMDA receptors being present and pathway-specific there
  (Weisskopf & LeDoux 1999). The "NMDA = coincidence detector" story is a
  proposed mechanism (its own originating authors' hypothesis), not fully
  closed.
- **The savings-robustness sweep caught a real artifact** in this script's
  own first version (integer trials-to-criterion failed at 2/7 head-start
  fractions) — diagnosed and fixed via a provable monotonicity argument, not
  smoothed over; both the flawed and the fixed metric are reported in the
  raw JSON.
- **Renewal is not simply "any context switch."** Bouton & King 1983's own
  Experiment 4 found no renewal in a specific reversed-context-order
  arrangement — a real, order-dependent structure this doc's toy renewal
  gate does not capture.
- **Human lesion evidence is small-n.** Bechara et al 1995's triple
  dissociation is 3 patients; the S.M. lineage is congenital
  (Urbach-Wiethe), not an acute adult lesion — not equivalent, not smoothed
  into a population-level claim.
- **The PTSD gate-efficacy parameters (0.85 vs 0.35) are illustrative**,
  reproducing the DIRECTION of Milad et al 2009's finding, not independently
  fit to a reported effect size (that paper's abstract does not give an
  exact SCR-ratio number).
- **η=0.35 (acquisition rate) is free/uncalibrated**, disclosed as such — no
  primary source with an explicit number-of-pairings-vs-percent-CR empirical
  curve was pinned down live this session (a genuine, searched-for, disclosed
  gap; the closest anchors found — Quirk 2002's 7-pairing protocol,
  Quirk/Repa/LeDoux 1995's "first several trials" — are qualitative/
  protocol-level, not a digitized curve).
- **Rescorla & Wagner 1972 is a book chapter**, confirmed live this session
  to have 0 PubMed hits — cited as classic/theoretical tier, the same
  disclosed treatment this repo gives Bernheimer 1973 / Nagatsu 1964.
- **No graph-edge write this session** — folding into `NEU-AMYGDALA-THREAT`
  or any other graph node requires the separate `mechanism_fold →
  fold_gate_v2` path, not performed here (isolation: touch only files
  created this session).
- **Single, population-level build** — no subject2 fear-conditioning panel
  exists to certify against; every operating point here is generic/
  literature-level, the same disclosed scope every sibling MECHANISM
  neuro/endocrine doc in this repo already carries.

## Files

- `scripts/msk/amygdala_fear_circuit.py` — self-contained (numpy only, no
  OpenSim); builds the dual-pathway graph-reachability model + NMDA-
  coincidence-window check, the delta-rule associative-learning gate + void
  floor + NMDA-block dissociation, the lesion-necessity + forced-adversary
  cross-tabulations, the extinction-not-erasure structural model (with the
  caught-and-fixed savings-robustness sweep), the PTSD dysfunction-pole
  model, and the HPA-timescale couples_to computation; writes the evidence
  JSON below; prints all 22 gates + overall_pass.
- `data/amygdala_fear_circuit/amygdala_fear_circuit_results.json` — every
  number in this doc, machine-written: all 27 citations + 1 classic
  reference, every falsifier computation, the full savings-robustness sweep
  (both the flawed integer metric and the fixed continuous metric), and all
  22 gates. Verified deterministic (2 independent runs, byte-identical JSON
  via `diff`) and NaN/Inf-free (checked programmatically).
- `docs/MECHANISM_AMYGDALA_FEAR_CIRCUIT_evidence.json` — the structured
  evidence summary (task, confidence tier, gates summary with symmetric
  self-grading, couples_to, distinct-from-existing-graph-node analysis,
  citations, honest gaps).
- Raw fetch logs (verbatim NCBI E-utilities responses, this session):
  `data/raw_fetch/amygdala_efetch_batchA_latency_nmda_ltp.txt` through
  `..._batchG_nmda_kinetics_nuance.txt`, plus
  `data/raw_fetch/amygdala_efetch_ciocchi2010_indep_reverify.txt`.
- Input read (READ-ONLY, not modified):
  `data/msk_smoketest/subject2_walking1/hpa_cortisol_axis/hpa_cortisol_axis_results.json`
  (for the couples_to HPA-timescale computation, §6); `data/MECHANISM_ANCHOR_
  GRAPH.json` (the pre-existing `NEU-AMYGDALA-THREAT` node this doc relates
  to — read, not edited).

## Repro

```
cd ~/projects/bodytwin
source .venv-msk/bin/activate
python3 scripts/msk/amygdala_fear_circuit.py
```

Pure Python/numpy, no OpenSim call, no upstream MECHANISM JSON required (the
HPA couples_to computation degrades gracefully — `couples_to.hpa_cortisol_
timescale.available=False` — if that sibling file is absent, affecting only
§6, not the 21 falsifier gates). Runs in under 2 seconds, deterministic. No
git operations; writes only under `data/amygdala_fear_circuit/`.
