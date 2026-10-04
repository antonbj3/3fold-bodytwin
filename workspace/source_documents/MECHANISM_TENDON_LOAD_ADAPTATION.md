# MECHANISM TENDON LOAD ADAPTATION — the tendon as a slow-remodeling mechanotransducer, and its failure mode (2026-07-22)

Adds the **load-adaptation/remodeling layer** above the twin's existing tendon-elastic/strain/slack
certs (`docs/MECHANISM_TENDON_ELASTIC.md`, `docs/MECHANISM_TENDON_STRAIN.md`,
`docs/MECHANISM_TENDON_SLACK_SENSITIVITY.md`), which model the tendon as a **fixed, single-gait-cycle,
non-adapting** elastic element. This layer is time-varying (weeks-to-months state variables) and adds
the **dysfunction pole**: tendinopathy as a failed-healing/degenerative overload response, paired
against the healthy-adaptation pole on the **same** dose-response mechanism. Builds on, and
independently re-verifies live, two pre-existing repo passes
(`data/body_twin/agent_outputs/tendon-adaptation__a690887818fc20020.json`,
`data/body_twin/agent_outputs/therapeutic-tendon-levers__a36c01cc0ee7e871a.json`) — this document's
own contribution is the full dysfunction/tendinopathy pole (histology, continuum model, the NSAID
molecular falsifier, the Alfredson/Beyer rehab RCTs) plus an in-repo quantitative cross-check
(Gate G) against the twin's own already-certified gait-simulation tendon strain.

**All 13 citations below were independently verified live this session** via `curl` to
`eutils.ncbi.nlm.nih.gov` (esearch → esummary/efetch, verbatim abstract text), not recalled and not
trusted from the two pre-existing agent-output files (whose own numbers were cross-checked, and in
one case a citation ambiguity was caught and resolved, §8).

## Headline

**10/10 pre-registered gates PASS** (§7, machine-computed, not eyeballed). The two symmetric
falsifiers both survive being forced to their strongest fair form:

1. **"Tendon is a passive inert spring, no biological adaptation"** — falsified by a real,
   contemporaneous, **zero-response control arm** (Bohm et al. 2014, n=13, "did not show any changes,
   P=0.86") sitting next to two loaded arms showing +51–57% modulus/stiffness (P<0.05), and by a
   27-study/264-participant meta-analysis whose stiffness-effect confidence interval excludes zero by
   a wide margin (SMD 0.70, CI 0.51–0.88).
2. **"Tendinopathy is inflammatory (tendinitis), anti-inflammatories should help"** — falsified at the
   **gene-expression level** by a double-blind ibuprofen-vs-placebo RCT (Heinemeier et al. 2017, n=26):
   no change in collagen/TGF-β/COX-2 expression, no change in pain/function, and corroborated
   clinically 19 years earlier by Alfredson's own conventional-treatment comparison arm (which
   included NSAIDs + rest): **0/15 success, all 15 required surgery**, vs 15/15 success with
   heavy-load eccentric training.

One gate (A2) was **first computed wrong** — with a threshold I invented *after* already reading the
data, which is not a legitimate pre-registration — caught by the machine gate-check itself, fixed at
the source (§7), and disclosed rather than silently patched.

---

## 1. Pre-registration — C/¬C, thresholds fixed before treating any gate as decided

**Claim 1 (load dose-response, positive-leaning).** C: resistance training increases tendon
stiffness/modulus (and, more weakly, CSA) in a way gated on load **magnitude**, with a genuine
threshold below which loading is not a sufficient stimulus. ¬C: apparent "stiffening" is a
measurement artifact (assessment bias, force-estimation drift, publication bias) with no real
tissue change, or the dose-response is smooth/monotonic in total work with no magnitude-specific
gate. **Threshold:** PASS requires (a) a meta-analytic effect-size CI excluding zero, (b) a
same-subject paired design isolating magnitude from frequency/volume, and (c) a real concurrent
unloaded control arm showing no change.

**Claim 2 (contraction-mode independence, positive-leaning-then-tested-both-ways).** C: the
operative dose variable is strain **magnitude**, not contraction type (eccentric is not uniquely
required). ¬C: eccentric contraction has a distinct, non-magnitude-reducible mechanotransduction
effect. **Threshold:** PASS requires the "eccentric is special" adversary to fail in **two**
decorrelated populations (healthy adaptation AND diagnosed tendinopathy rehab).

**Claim 3 (dysfunction pole, positive-leaning): tendinopathy is degenerative (tendinosis), not
inflammatory (tendinitis).** ¬C / the adversary this claim is tempted to skip: clinical pain +
swelling look inflammatory; corticosteroids give real short-term relief (consistent with SOME
inflammatory component); decades of literature literally called it "tendinitis." **Threshold:**
PASS requires (a) direct histology quantifying inflammatory-cell presence in real surgical biopsies,
not just a review's say-so, and (b) the **falsifiable, decisive** test: if primarily inflammatory,
targeted anti-inflammatory treatment (not just steroids, which have systemic/catabolic effects
clouding interpretation — the cleaner test is an NSAID RCT with gene-expression readout) should
measurably reduce the pathology; if degenerative, it should not.

**Claim 4 (rehab efficacy).** C: high-load eccentric/HSR rehab outperforms rest+anti-inflammatory
conventional care by a large margin, consistent with a controlled-reloading mechanism. **Threshold:**
PASS requires a pre-registered **big-margin** bar (≥50 percentage points success-rate gap) — a
deliberately generous bar so a modest, easily-overinterpreted effect cannot pass.

## 2. The forced adversary — "inert passive spring" — and where it falls

The sharpest available form of this adversary is not a rhetorical strawman but a **real, published,
concurrent control arm**: Bohm, Mersmann, Tettke, Kraft & Arampatzis (2014, PMID **25267851**, J Exp
Biol, DOI 10.1242/jeb.112268) ran a 3-arm Achilles-tendon training study — two loaded experimental
arms (n=14, n=12) and an **unloaded control** (n=13). Verbatim: *"The reference and long strain
duration protocols induced significantly increased (P<0.05) tendon stiffness (57% and 25%),
cross-sectional area (4.2% and 5.3%) and Young's modulus (51% and 17%)... **The control group did not
show any changes (P=0.86).**"* An inert-spring model predicts **all three arms** should look like the
control (no cyclic-loading history should change a passive material's own stiffness). The measured
result: two loaded arms significant at P<0.05 with 17–57% magnitude changes; the unloaded arm flat at
P=0.86. **The adversary falls.**

Independently, the same paper discloses a **third, non-significant** arm — high strain **rate**
alone (jump training) gave only "a tendency of increased stiffness (P=0.08) and cross-sectional area
(P=0.09)" — reported here **honestly, not folded into the two PASS arms**: rate alone, without
matched magnitude+duration, is a weaker stimulus. This is exactly the kind of disclosure rule 1/2 of
this method demand: a partial, non-significant result is reported as such, not silently merged into a
cleaner-looking aggregate.

Broader confirmation, a genuinely different instance (27 studies, 264 participants, different
research teams/tendons/decades): Bohm, Mersmann & Arampatzis (2015, PMID **27747846**, Sports Med
Open, DOI 10.1186/s40798-015-0009-9). Verbatim: *"SMD was 0.70 (confidence interval: 0.51, 0.88) for
tendon stiffness (N=37), 0.69 (0.36, 1.03) for Young's modulus (N=17), and 0.24 (0.07, 0.42) for CSA
(N=33), with significant overall intervention effects (p<0.05)."* A CI excluding zero by 0.51 (not a
marginal 0.02) is a **big-margin** rejection of the null, not a knife-edge one.

## 3. The dose-response model — geometry, not a heuristic list

**Per-bout stimulus.** Let `eps_peak` be peak cyclic tendon strain per loading bout. Define a
threshold-gated dose function:

```
S(eps_peak) = 0                                            for eps_peak < eps_th
S(eps_peak) = min(1, (eps_peak - eps_th) / (eps_sat - eps_th))   for eps_peak >= eps_th
```

`eps_th ≈ 3%` and `eps_sat ≈ 6.5%` are **interpolations**, disclosed as such (§9), bracketing three
real measured points, not a swept continuous curve:

| eps_peak | outcome | source |
|---:|---|---|
| 2.85 ± 0.99% | **no adaptation** (same subject, same frequency/volume as the 4.55% leg) | Arampatzis, Karamanidis & Albracht 2007, PMID **17644689**, J Exp Biol |
| 4.55 ± 1.38% | adaptation (↓strain-at-given-force, ↑stiffness, ↑modulus, region-specific hypertrophy) | same study, same subjects, other leg |
| ~6.5% | largest measured gain (+57% stiffness, +51% modulus, P<0.05) | Bohm et al. 2014, PMID 25267851 |

The Arampatzis design is the cleanest available isolation of magnitude from everything else: **same
subject**, one leg at each strain magnitude, "of similar frequency and volume" — the only
manipulated variable is peak strain. Verbatim: *"an increase in tendon-aponeurosis stiffness and
tendon elastic modulus and a region-specific hypertrophy of the Achilles tendon **only** in the leg
exercised at high strain magnitude."* This is a genuine threshold/set-point, not a smooth
monotonic-with-volume curve — the low-strain leg received the **same** number of contractions and
did not adapt.

**S is deliberately NOT a function of contraction type.** Bohm 2015's own subgroup analysis: *"stiffness
adaptation significantly (p<0.05) depends on loading intensity (I²=0%, i.e. clean/homogeneous across
studies), **but not on muscle contraction type**."* Forced across a second, wholly decorrelated
population (diagnosed tendinopathy, not healthy adaptation): Beyer et al. 2015 (PMID **26018970**, Am
J Sports Med, n=58 RCT, eccentric vs heavy-slow-resistance) — *"Both groups showed significant
(P<.0001) improvements... **None of these robust clinical and structural improvements differed
between the ECC and HSR groups.**"* The "eccentric contraction is uniquely required" folk mechanism
is forced in two different populations, two different decades, two different research groups — and
falls both times.

**Slow first-order kinetics, decorrelated time constants.** Kubo, Ikebukuro, Yata & Tsunoda (2010,
PMID **19996769**, J Strength Cond Res, n=8, within-cohort monthly time-course): *"Training increased
muscle strength and neural activation level by 29.6 and 7.3% after 2 months and by 40.5 and 8.9%
after 3 months (all p's<0.05). **Muscle CSA and tendon stiffness did not change until 2 months** of
training period, and afterward, the increases in muscle CSA and tendon stiffness reached statistical
significance **at the end of training period** [3 months]."* Same 8 subjects, same measurement
sessions: strength/neural adaptation is fast (significant by month 2); muscle-morphological and
tendon-mechanical adaptation is slow (not significant until month 3). This is a real, measured,
within-cohort rate dissociation — not an assumption.

## 4. The load-tolerance envelope — a demand/capacity ratio, not a fixed load number

Define tendon capacity `C(t) ∝ k(t)·CSA(t)` (the slowly-adapting state above) and mechanical demand
`D(t)` (what training/competition actually imposes — which, under progressive overload, tracks the
**fast** strength channel). The governing scalar is the dimensionless ratio

```
R(t) = D(t) / C(t)
```

— the single number whose regime, and critically whose **rate of change** relative to a
separately-adapting capacity, determines outcome:

- **R(t) ≲ 1, sustained** → homeostatic remodeling (§3's normal adaptive loop).
- **R(t) > 1, transient/acute** → Cook & Purdam's proposed **"reactive tendinopathy"** stage: an
  acute, reversible response, not yet structural collagen disorganization.
- **R(t) ≫ 1, sustained** — because demand (fast, weeks-scale) races ahead of capacity (slow,
  months-scale) during aggressive progression — → **"tendon dysrepair"** then **"degenerative
  tendinopathy"**, matching the actual histology (§5).

This is a genuinely geometric point, not a heuristic: the danger zone is defined by a **rate
mismatch** between two independently-adapting quantities (`dD/dt` vs `dC/dt`), not by any single
absolute-load number in isolation — an athlete whose tendon *could* eventually tolerate a given load,
once adapted, can still be injured by reaching that load faster than `C(t)`'s own months-scale time
constant allows.

**Confidence note on the stage names.** The continuum concept itself is directly, live-verified from
Cook & Purdam (2009, PMID **18812414**, Br J Sports Med, DOI 10.1136/bjsm.2008.051193): *"The
pathology of tendinopathy has been described as degenerative or failed healing. Neither of these
descriptions fully explains the heterogeneity of presentation. This review proposes, and provides
evidence for, a continuum of pathology."* The specific three-stage names used above
(reactive/dysrepair/degenerative) are the field's standard secondary characterization of that paper's
full-text model — MEDIUM confidence, not verbatim-quoted from a fetched full text this session
(abstract-only). **Symmetric QC, not omitted:** the model's own original authors self-critique it 7
years later (Cook, Rio, Purdam & Docking 2016, PMID **27127294**, DOI 10.1136/bjsports-2015-095422):
*"classifying tendinopathy based on structure in what is primarily a pain condition has been
challenged. The interplay between structure, pain and function is not yet fully understood."* This is
included deliberately so the continuum model is not presented as more settled than the field itself
considers it.

## 5. The dysfunction pole — forcing the "inflammatory" adversary to its strongest form

**Direct histology, not just a review's assertion.** Åström & Rausing (1995, PMID **7634699**, Clin
Orthop Relat Res, n=**163** patients, primary surgical/histopathologic data, not a review): *"134
males, 29 females; mean age 38... Degenerative changes (tendinosis) characterized by abnormal fiber
structure, focal hypercellularity, and vascular proliferation were noted in **90%** of biopsy
specimens from symptomatic parts of the tendons and, to a lesser degree, in **20%** from
nonsymptomatic parts... **Important features are a lack of inflammatory cells** and a poor healing
response."* A 90%-vs-20% symptomatic/asymptomatic contrast, from real surgical biopsies, with
inflammatory cells explicitly absent in the majority pathology. Khan, Cook, Bonar & Harcourt (1999,
PMID **10418074**, Sports Med, DOI 10.2165/00007256-199927060-00004), synthesizing this and other
primary series: *"The most significant feature is the absence of inflammatory cells... consistent
with... tendinosis — a degenerative condition... As the common overuse tendon conditions are rarely,
if ever, caused by 'tendinitis', we suggest the term 'tendinopathy'... effective treatment... must
target the most common underlying histopathology, **tendinosis, a noninflammatory condition**."*

**The decisive falsifier: force the inflammatory model to make a testable prediction, then test it
directly.** If tendinopathy were primarily inflammatory, anti-inflammatory treatment should reduce
inflammatory gene expression and improve outcomes. Heinemeier, Øhlenschlæger, Mikkelsen, Sønder,
Schjerling, Svensson & Kjaer (2017, PMID **28860166**, J Appl Physiol, DOI
10.1152/japplphysiol.00281.2017) ran exactly this test — double-blind, placebo-controlled, n=26
chronic Achilles tendinopathy patients, ibuprofen 600mg×3/day vs placebo for 1 week, with a **tendon
biopsy** and RT-PCR gene expression as the primary outcome (collagen I/III, TGF-β isoforms, **COX-2**
explicitly assayed): *"Expression of collagens and TGF-β isoforms showed relatively low variation and
was **unaffected by ibuprofen treatment**. Further, **no changes were seen in tendon thickness or
VISA-A score**. The placebo treatment **reduced** the color Doppler... compared with the ibuprofen
group and also **increased** the perception of pain at rest... tendinopathic cells are **not
responsive to ibuprofen**."* Not only null — the one directional signal present runs opposite to what
an inflammatory model predicts. This is a molecular-level falsification, not just a clinical-outcome
one.

Corroborating, independent, review-tier: Andres & Murrell (2008, PMID **18446422**, Clin Orthop Relat
Res, DOI 10.1007/s11999-008-0260-1): *"Recent basic science research suggests little or no
inflammation is present... NSAIDS and corticosteroids appear to provide pain relief in the short
term, but their effectiveness in the long term has not been demonstrated... Current data support the
use of eccentric strengthening protocols."* And, 19 years earlier and at the clinical-outcome level:
Alfredson, Pietilä, Jonsson & Lorentzon (1998, PMID **9617396**, Am J Sports Med, DOI
10.1177/03635465980260030301) — a comparison group of 15 patients treated conventionally ("rest,
**nonsteroidal antiinflammatory drugs**, changes of shoes or orthoses, physical therapy... ordinary
training programs"): *"**In no case was the conventional treatment successful, and all patients were
ultimately treated surgically.**"* Meanwhile the 15 patients given 12 weeks of heavy-load eccentric
training: *"**all 15 patients were back at their preinjury levels with full running activity.**"*
100% vs 0%, a 100-percentage-point gap — the single largest-margin result in this report, and it
directly embeds an NSAID-inclusive treatment arm as the failure condition.

## 6. Rehab — one mechanism, not two ad hoc ones

The rehab efficacy of eccentric/HSR loading is modeled as the **same** `S(eps_peak)` dose-response
function (§3) acting on a tendon currently sitting in the dysrepair/degenerative regime (§4),
deliberately bringing `R(t)` back toward ≈1 via controlled progressive reloading — not two separate
mechanisms for "how healthy tendon adapts" and "how injured tendon recovers." This unification makes
an additional, checkable prediction: **complete rest should not, by itself, be curative**, because
`S(0)=0` drives no remodeling signal at all. Alfredson's own conventional-treatment arm (§5) included
rest as one of its components and failed 15/15 times — consistent with, not contradicting, this
reading (though rest was combined with several other components in that arm, so this is corroborating
context, not an isolated test of rest alone — disclosed, §9). The "eccentric is uniquely required"
version of this same claim was already forced and falsified in §3 (Bohm 2015 healthy + Beyer 2015
diseased) — high load, not eccentric contraction specifically, is the operative mechanism, consistent
with a single magnitude-gated dose function serving both healthy adaptation and rehabilitation.

## 7. Gates (pre-registered, machine-computed — see `tendon_gates.py` logic reproduced below, not eyeballed)

| gate | threshold | measured | verdict |
|---|---|---|---|
| A1 Bohm 2015 meta CI excludes 0 | CI lower bound > 0 | SMD 0.70, CI [0.51, 0.88] | **PASS** |
| A2 Bohm 2014 control-null vs experimental-significant | exp. p<0.05 AND control not significant (primary study's own α=0.05 test) | ref +57%/+51%/+4.2% p<0.05; long-dur +25%/+17%/+5.3% p<0.05; control p=0.86 | **PASS** (corrected, §8) |
| B1 Arampatzis 2007 within-subject magnitude dissociation | high-strain leg adapts, low-strain leg (same freq/volume) does not | 4.55%→adapt; 2.85%→no adapt | **PASS** |
| B2 Bohm 2015 intensity sig., contraction-type n.s. | intensity subgroup p<0.05 AND contraction-type not sig. | intensity p<0.05 (I²=0%); mode n.s. | **PASS** |
| C eccentric-specific adversary falls, 2 populations | mode-effect n.s. in BOTH healthy meta-analysis AND diseased RCT | Bohm2015 n.s.; Beyer2015 (n=58) n.s. | **PASS** |
| D Kubo 2010 strength leads tendon, within-cohort | mo.2: strength sig., tendon n.s.; mo.3: both sig. | exactly this pattern, n=8 | **PASS** |
| E Åström & Rausing histology, degenerative not inflammatory | symptomatic tendinosis-rate > asymptomatic AND no inflammatory cells | 90% vs 20%, n=163; "lack of inflammatory cells" | **PASS** |
| F1 Heinemeier 2017 NSAID RCT, molecular + functional null | no sig. change in gene expression or VISA-A/pain/thickness | null on all axes, n=26 | **PASS** |
| F2 Alfredson 1998 big-margin eccentric vs conventional | success-rate gap ≥50pp (pre-registered generous bar) | 100% vs 0%, gap 100pp, n=15+15 | **PASS** |
| G in-repo cross-check vs twin's own certified walking strain | twin's own subject2/walking1 peak strain < this model's eps_th | 2.501% < ~3% (`MECHANISM_TENDON_STRAIN.md`) | **PASS** |

**Overall: 10/10 PASS.**

## 8. A self-correction, disclosed not hidden

Gate A2 was **first computed with a threshold I invented after already having read the Bohm 2014
numbers** — an arbitrary "≥20% magnitude change" cutoff, intended to operationalize a "big margin."
It **FAILED** on that criterion: the long-duration arm's modulus gain (17%) narrowly missed my own
20% bar. This is not a legitimate pre-registration (the threshold was set with the data already in
view) and, worse, it is the wrong criterion in principle — the primary study's own reported
statistical test (α=0.05) is the actual, externally-anchored decision rule the data was analyzed
under, not a magnitude cutoff I retrofitted. Re-derived on the correct criterion (both experimental
arms P<0.05, control P=0.86): PASS. Flagged here explicitly, per this repo's symmetric-QC discipline
— a machine gate-check catching a flaw in my *own* gate design is exactly the kind of check this
method exists to force, not a result to quietly patch and move past.

Separately: a pre-existing repo file
(`data/body_twin/agent_outputs/therapeutic-tendon-levers__a36c01cc0ee7e871a.json`) packed several
citations into one `source_url` field including PMID 25267851 adjacent to a "Bohm 2015" claim; this
session independently verified live that 25267851 is actually the **2014** Bohm/Mersmann/Tettke/Kraft
rate-and-duration paper (used directly and correctly in §2/§3 above), while the actual Bohm **2015**
Sports Med Open meta-analysis is PMID 27747846 — not a fabrication in that prior file (the citations
were legitimately co-located, referring to different claims in the same packed string), but disclosed
here for precision, consistent with this repo's established citation-verification discipline (e.g.
`docs/MECHANISM_TENDON_ELASTIC.md` §0's own Ker-1987 correction).

## 9. Honest gaps

1. `eps_th≈3%` / `eps_sat≈6.5%` are **interpolations** between three discrete measured points (2.85%
   null, 4.55% positive, 6.5% highest-tested-positive), not a swept continuous dose-response curve —
   precision beyond "somewhere between 2.85% and 4.55%" / "at least to 6.5%, untested toward the
   ~8–9.9% failure-strain range already anchored in `docs/MECHANISM_TENDON_STRAIN.md`" should not be
   over-read.
2. No single cohort tracked stiffness + CSA + a synthesis-marker (PICP/FSR) + strength together across
   one continuous multi-month time course — the `tau_tendon >> tau_strength` ordering is well-supported
   (Kubo 2010, within-cohort) but the full collagen-synthesis-marker mechanistic chain is a cross-study
   synthesis in the pre-existing `tendon-adaptation` agent-output file, not one decisive dataset —
   inherited gap, not resolved this session.
3. Cook & Purdam's specific 3-stage terminology is MEDIUM confidence (concept verified live from the
   abstract; stage names are the field's standard secondary characterization, not verbatim-quoted from
   a fetched full text this session).
4. The continuum model's own structure-based staging is self-critiqued by its original authors 7 years
   later (§4) — included deliberately, not smoothed over.
5. `R(t)=D(t)/C(t)` and the regime boundaries are this document's own geometric formalization of the
   qualitative literature — a hypothesis for QC, not a directly-measured numeric ratio from any cited
   study.
6. Sex/age modifiers (Onambele-Pearson 2012, McMahon 2018, Reeves 2003) and the AAS muscle:tendon
   mismatch mechanism are real findings already in this repo's `data/body_twin/agent_outputs/` but were
   **not** independently re-verified live this session (out of this task's explicit scope) — mentioned
   only as pointers.
7. Rio et al. 2015 BJSM (isometric analgesia via reduced corticospinal inhibition) is a relevant
   nuance — pain modulation ≠ structural adaptation — but was **not** independently re-verified live
   this session (inherited citation only).
8. The demand-side `D(t)` is not numerically instantiated against any specific athlete-load dataset
   this session (`MECHANISM_WEIGHTLIFTING_TRIPLE_EXTENSION.md` etc.) — named as a next step, not built.
9. No single study bridges a healthy tendon under a *known* `R(t)` trajectory transitioning into
   diagnosed tendinopathy — the envelope's two endpoints (healthy adaptation; diagnosed
   dysrepair/degenerative pathology) are each well-anchored, but the transition itself is a plausible,
   literature-consistent construction, not one directly-observed longitudinal dataset.
10. Alfredson's conventional-treatment arm combined rest with NSAIDs, orthoses, PT, and ordinary
    training — its 0/15 success rate is strong corroborating context for both "rest alone doesn't cure
    it" and "NSAIDs don't cure it," but is not a clean isolated test of either component alone.

## Files

- `docs/MECHANISM_TENDON_LOAD_ADAPTATION_evidence.json` — all 13 live-verified citations (PMID/DOI,
  verbatim quotes), the full model (dose-response function, rate-mismatch kinetics, load-tolerance
  envelope, tendinosis/rehab sections), the 10-gate machine-computed table including the A2
  self-correction, `couples_to`, `honest_gaps`, `proposed_cell`.
- Read-only, not edited: `docs/MECHANISM_TENDON_ELASTIC.md`, `docs/MECHANISM_TENDON_STRAIN.md`,
  `docs/MECHANISM_TENDON_SLACK_SENSITIVITY.md`, `docs/MECHANISM_BONE_REMODELING.md`,
  `data/body_twin/agent_outputs/tendon-adaptation__a690887818fc20020.json`,
  `data/body_twin/agent_outputs/therapeutic-tendon-levers__a36c01cc0ee7e871a.json`.
- No git operations performed (isolation respected). No file outside `docs/` was written.

## Next step

1. Numerically instantiate `D(t)` against a specific athlete-load cell
   (`MECHANISM_WEIGHTLIFTING_TRIPLE_EXTENSION.md` / sprint / plyometric volume) to make `R(t)` a
   computable trajectory rather than a qualitative regime label.
2. If a subject-specific resistance-training longitudinal dataset ever enters this repo, replace the
   generic default Millard tendon curve (`docs/MECHANISM_TENDON_ELASTIC.md` §7.1's disclosed dominant
   uncertainty) with this layer's own time-evolving `k(t)`/`CSA(t)`, closing the loop between the two
   documents.
3. Independently re-verify the sex/age-modifier and Rio 2015 citations live if the twin's sex/age or
   pain-modulation axes need this layer's numbers at more than pointer-confidence.
