# MECHANISM MUSCLE HYPERTROPHY — the mechanical-load → mTORC1 → protein-accretion training adaptation (2026-07-22)

**Status: HYPOTHESIS (literature-synthesis + machine-verified-arithmetic layer, `status:OPEN` per this
repo's node schema) — awaiting an independent QC pass.** Builds the resistance-training hypertrophy
chain — mechanical tension (primary) + metabolic stress + damage → mechanotransduction → mTORC1
activation → ↑myofibrillar muscle-protein-synthesis (MPS) exceeding breakdown → net protein accretion →
fiber-CSA/PCSA increase → ↑force capacity — as a falsifier-tested model, feeding this twin's existing
`docs/MECHANISM_FMAX_PCSA_VALIDATION.md` / `docs/MECHANISM_SPECIFIC_TENSION.md` / `docs/
MECHANISM_FORCE_VELOCITY_POWER_RFD.md` static+dynamic strength substrate. All arithmetic below is direct,
reproducible Python (shown inline, not hand-narrated) — this task's scope is two files + a summary line,
not a full repo-integrated script+PREREG pipeline, so no standalone script was built; every number is
still machine-computed, not eyeballed. ISOLATION: bodytwin only, read-only on every file not created this
session, no git commit/push.

**Pre-existing graph context (read-only, not edited this session — isolation discipline, another instance
may write concurrently):** `data/MECHANISM_ANCHOR_GRAPH.json` (1038 nodes) already carries four `status:OPEN`
relatives this document couples to rather than duplicates: `MSK-MUSCLE-HYPERTROPHY-TIMECOURSE` (the
apparent-CSA-vs-true-growth edema dissociation + session-count framework, Damas et al.), `MOL-MTORC1-
NUTRIENT-SIGNALING` (mTORC1-necessity via a *lifespan* readout, NIA-ITP), `MSK-SATELLITE-CELL-RESERVE`
(already tracks the McCarthy-vs-Egner depletion dispute this doc's §6 independently re-lands on), and the
already-**executed** `MSK-SARCOPENIA-COMPOSITE` (NHANES dynapenia-dissociation cert, n=2810). This
document's job is the two falsifiers the task specifies (mTORC1 necessity; the systemic-hormone-hypothesis
refutation) plus the volume dose-response and the geometric PCSA/Fmax coupling — none of which the four
existing nodes already cover — and to state the couplings concretely (§7), not just as prose pointers.

## The two falsifiers, verdicts stated up front

> **Falsifier 1 (mTORC1 necessity)**: does pharmacological mTORC1 blockade (rapamycin) BLOCK the
> load-induced MPS increase (human) and hypertrophy (animal, causal loss-of-function AND gain-of-function),
> forcing the strongest in-field rival (calcineurin/NFAT) to fall in the same experiments?

**YES — PASS**, on a pre-registered threshold of "rapamycin blocks, not merely attenuates, in ≥1 human RCT
AND ≥1 animal causal study spanning multiple hypertrophy-induction models." Bodine et al. 2001 (rat, *Nat
Cell Biol*): rapamycin **blocked hypertrophy in all models tested**; calcineurin inhibitors (cyclosporin A,
FK506) — the forced adversary — **did not blunt hypertrophy** in the same models, and calcineurin/NFAT was
not even activated during hypertrophy in vivo. Genetic activation of Akt/mTOR was **sufficient** to cause
hypertrophy and prevent atrophy; genetic blockade **blocked** hypertrophy. Drummond et al. 2009 (human, *J
Physiol*): rapamycin blocked the early (1–2 h) ~40% contraction-induced MPS rise, blocked the 6-fold S6K1
phosphorylation rise, and blocked the ~25% eEF2 dephosphorylation (translation elongation) — see §2.

> **Falsifier 2 (systemic-hormone hypothesis)**: forced to its strongest fair (CAUSAL, not merely
> correlational) form, does the transient post-exercise testosterone/GH/IGF-1 surge drive hypertrophy?

**ADVERSARY FALLS on its pre-registered causal threshold — PASS for the refutation — but a weaker,
disclosed, NOT-fully-closed residual survives.** West et al. 2010 (*J Appl Physiol*, the actual controlled
manipulation: identical local arm loading, one arm's contralateral leg-work spikes systemic hormones,
15 wk) found CSA +12% (low-hormone) vs +10% (high-hormone) — **nominally reversed** from the hypothesis'
predicted direction, interaction P = 0.25. Acute MPS/signalling (West 2009) was likewise null (P = 0.72,
high-hormone arm numerically *lower*). **Forced independent-lab adversary** (Ahtiainen et al. 2003,
Finland — fully decorrelated from the Canada/McMaster group behind every other hormone citation here):
found CORRELATIONAL basal-T/ΔT-vs-outcome relationships in small subgroups (n=8 each) — this survives
*partially*, as a trait-level/confounded claim, but does not overturn the causal test (no local-stimulus
control, small n, plausible multiple-comparisons noise — see §3.4). Symmetric verdict, not smoothed over.

---

## 0. Citations — 14 fresh PMIDs live-verified this session (NCBI eutils esearch→esummary→efetch, raw
`curl --max-time 25`, never WebFetch-summarized), + 5 read-only citations already live-verified by sibling
graph nodes this session did not re-derive (provenance disclosed per row)

**Recall-drift caught live, not trusted** (this repo's own measured ~62% citation-drift-from-memory rate):
Bodine 2001 *Nat Cell Biol* recalled as PMID `11715018` — **wrong**; live esearch/esummary gives
**`11715023`**. West & Phillips 2012 recalled as `22105450` — **wrong**; real PMID **`22105707`**. West &
Phillips 2010's review recalled as `20959701` (off by one digit) — **wrong**; real PMID **`20959702`**
(the `20959701` record is an unrelated patellar-mobility-test paper). Egner et al.'s satellite-cell paper
recalled as "2013" — **wrong year**; the real paper is **2016** (PMID `27531949`; `24004949`, my first
guess, is an unrelated plant-biology paper). All four corrections shown, not silently fixed.

| # | Citation | PMID / DOI | n | Role |
|---|---|---|---|---|
| 1 | Phillips SM, Tipton KD, Aarsland A, Wolf SE, Wolfe RR (1997). *Mixed muscle protein synthesis and breakdown after resistance exercise in humans.* Am J Physiol 273(1 Pt1):E99-107. | **9252485**, DOI `10.1152/ajpendo.1997.273.1.E99` | 8 (4M/4F) | THE MPS timecourse anchor: FSR +112%(3h)/+65%(24h)/+34%(48h) above rest, P<0.01; net balance significant at all points, persisting to 48h. |
| 2 | Bodine SC, Stitt TN, Gonzalez M, et al. (2001). *Akt/mTOR pathway is a crucial regulator of skeletal muscle hypertrophy and can prevent muscle atrophy in vivo.* Nat Cell Biol 3(11):1014-9. | **11715023**, DOI `10.1038/ncb1101-1014` | NR (rat, multiple in vivo models) | Falsifier 1 decisive anchor: rapamycin blocks hypertrophy in ALL models; calcineurin/NFAT adversary FALLS; genetic activation sufficient, blockade necessary. |
| 3 | Drummond MJ, Fry CS, Glynn EL, et al. (2009). *Rapamycin administration in humans blocks the contraction-induced increase in skeletal muscle protein synthesis.* J Physiol 587(Pt7):1535-46. | **19188252**, DOI `10.1113/jphysiol.2008.163816`, PMC2678224 | NR in abstract | Falsifier 1, human leg: rapamycin blocks early (1-2h) ~40% MPS rise, 6-fold S6K1-P rise, ~25% eEF2 dephosphorylation. |
| 4 | Schoenfeld BJ, Ogborn D, Krieger JW (2017). *Dose-response relationship between weekly resistance training volume and increases in muscle mass: a systematic review and meta-analysis.* J Sports Sci 35(11):1073-82. | **27433992**, DOI `10.1080/02640414.2016.1210197` | 34 treatment groups / 15 studies | Volume dose-response meta: continuous-sets P=0.002 (+0.37%/set); 2-level categorical P=0.03 (3.9% diff); 3-level categorical only a TREND, P=0.074 (disclosed, not hidden). |
| 5 | West DW, Kujbida GW, Moore DR, et al. (2009). *Resistance exercise-induced increases in putative anabolic hormones do not enhance muscle protein synthesis or intracellular signalling in young men.* J Physiol 587(Pt21):5239-47. | **19736298**, DOI `10.1113/jphysiol.2009.177220`, PMC2790261 | 8 | Falsifier 2, acute leg: MPS rest=0.040, LH=0.071, HH=0.064 %/h; hormone effect P=0.72 (HH numerically LOWER than LH despite huge hormone gap). |
| 6 | West DW, Burd NA, Tang JE, et al. (2010). *Elevations in ostensibly anabolic hormones with resistance exercise enhance neither training-induced muscle hypertrophy nor strength of the elbow flexors.* J Appl Physiol 108(1):60-7. | **19910330**, DOI `10.1152/japplphysiol.01147.2009`, PMC2885075 | 12 | Falsifier 2, DECISIVE causal leg: 15wk training, CSA +12%(LH)/+10%(HH), interaction P=0.25; strength gain equal, no HH advantage. |
| 7 | West DW, Phillips SM (2012). *Associations of exercise-induced hormone profiles and gains in strength and hypertrophy in a large cohort after weight training.* Eur J Appl Physiol 112(7):2693-702. | **22105707**, DOI `10.1007/s00421-011-2246-z`, PMC3371329 | 56 | Falsifier 2, correlational leg: T/IGF-1 AUC r≈0 vs LBM/strength; GH AUC r=0.36(P=0.006,typeI)/0.28(P=0.04,typeII) fiber CSA only; cortisol AUC r=0.29(P<0.05,LBM)/0.35(P<0.01,typeII CSA) — weak, disclosed, not zero. |
| 8 | West DW, Phillips SM (2010). *Anabolic processes in human skeletal muscle: restoring the identities of growth hormone and testosterone.* Phys Sportsmed 38(3):97-104. | **20959702**, DOI `10.3810/psm.2010.10.1814` | review | Synthesis: "GH...is not anabolic toward...myofibrillar muscle tissue in healthy individuals"; "local mechanisms...are predominant." |
| 9 | Cuthbertson D, Smith K, Babraj J, et al. (2005). *Anabolic signaling deficits underlie amino acid resistance of wasting, aging muscle.* FASEB J 19(3):422-4. | **15596483**, DOI `10.1096/fj.04-2640fje` | 44 (young+old) | Dysfunction pole: basal MPS indistinguishable by age; EAA-stimulated MPS blunted in old via reduced mTOR/p70S6k/4E-BP1/eIF2B activation, insulin-independent (clamped); NFκB up. |
| 10 | McCarthy JJ, Mula J, Miyazaki M, et al. (2011). *Effective fiber hypertrophy in satellite cell-depleted skeletal muscle.* Development 138(17):3657-66. | **21828094**, DOI `10.1242/dev.068858`, PMC3152923 | mouse (Pax7-DTA, >90% SC ablation) | Satellite-cell leg (pro): 2wk overload gives ~2-fold mass/CSA gain WITHOUT satellite cells; myonuclear-domain expansion substitutes for new nuclei; hyperplasia/regeneration (small component) blunted. |
| 11 | Egner IM, Bruusgaard JC, Gundersen K (2016). *Satellite cell depletion prevents fiber hypertrophy in skeletal muscle.* Development 143(16):2898-906. | **27531949**, DOI `10.1242/dev.134411` | mouse (same Pax7-based strategy) | Satellite-cell leg (contra, FORCED ADVERSARY): overload hypertrophy PREVENTED in both plantaris+EDL; attributes McCarthy's contrary finding to mass-as-proxy + regenerating-fiber inclusion artifacts. Genuine, disclosed, unresolved field dispute. |
| 12 | Burd NA, Holwerda AM, Selby KC, et al. (2010). *Resistance exercise volume affects myofibrillar protein synthesis and anabolic signalling molecule phosphorylation in young men.* J Physiol 588(Pt16):3119-30. | **20581041**, DOI `10.1113/jphysiol.2010.192856`, PMC2956949 | 8 | Mechanistic dose-response bridge: 1SET MPS 2.3-fold@5h→baseline@29h; 3SET 3.1-fold@5h→2.3-fold@29h (bigger amplitude AND longer duration); p70S6K-MPS correlation r=0.338, P=0.033. |
| 13 | Ahtiainen JP, Pakarinen A, Alen M, Kraemer WJ, Häkkinen K (2003). *Muscle hypertrophy, hormonal adaptations and strength development during strength training in strength-trained and untrained men.* Eur J Appl Physiol 89(6):555-63. | **12734759**, DOI `10.1007/s00421-003-0833-3` | 16 (8+8) | FORCED independent-lab adversary for Falsifier 2: basal-T-mean vs Δisometric-force r=0.84(P<0.01) in strength-athletes; acute-ΔT vs ΔCSA r=0.76(P<0.05) in non-athletes — correlational, small-n, opposite-direction CSA/force response between subgroups (§3.4). |
| 14 | Bodine SC, Latres E, Baumhueter S, et al. (2001). *Identification of ubiquitin ligases required for skeletal muscle atrophy.* Science 294(5547):1704-8. | **11679633**, DOI `10.1126/science.1065874` | rat (transcript profiling) + mouse (KO) | Dysfunction-pole molecular mirror: MuRF1/MAFbx are the universal atrophy-program E3 ligases across atrophy models; MAFbx overexpression sufficient for atrophy; MAFbx-KO/MuRF1-KO mice resistant. |

| # | Citation (read-only, verified by a sibling graph node this session, NOT independently re-verified here) | PMID | Sibling node |
|---|---|---|---|
| 15 | Damas F, et al. (2016). *Early- and later-phases muscle protein synthesis... hypertrophy.* Eur J Appl Physiol. | 26280652 | `MSK-MUSCLE-HYPERTROPHY-TIMECOURSE` |
| 16 | Damas F, et al. (2016). *Resistance training-induced changes in integrated myofibrillar protein synthesis are related to hypertrophy only after attenuation of muscle damage.* J Physiol. | 27219125 | `MSK-MUSCLE-HYPERTROPHY-TIMECOURSE` |
| 17 | Damas F, et al. (2018). Review reframing hypertrophy onset by training-session count. | 29282529 | `MSK-MUSCLE-HYPERTROPHY-TIMECOURSE` |
| 18 | Harrison DE, et al. (2009). Rapamycin extends lifespan, NIA-ITP. Nature. | 19587680 | `MOL-MTORC1-NUTRIENT-SIGNALING` |
| 19 | Miller RA, et al. (2011). Rapamycin lifespan replication, NIA-ITP. J Gerontol A. | 20974732 | `MOL-MTORC1-NUTRIENT-SIGNALING` |

---

## 1. Geometric structure — force = stress × area, and hypertrophy is a parallel-vs-series growth-axis choice

**The governing continuum-mechanics identity** (this repo's own convention, `docs/
MECHANISM_FMAX_PCSA_VALIDATION.md`/`MECHANISM_SPECIFIC_TENSION.md`, unchanged): `PCSA = volume /
optimal_fiber_length (Lopt)`, `Fmax = PCSA × specific_tension`. This is a genuine **area × stress**
geometric relationship, not an empirical curve-fit. Resistance-training hypertrophy is predominantly
**radial** growth — sarcomeres added **in parallel** (raising physiological cross-sectional area, PCSA)
— rather than **longitudinal** growth (sarcomeres added **in series**, raising Lopt without raising PCSA,
a different-axis adaptation more associated with specific eccentric/long-length training paradigms).
**This parallel-vs-series distinction is the geometric reason hypertrophy training predominantly shifts
the F-V-P curve's force axis (F0/Fmax) rather than its velocity axis (Vmax, which scales with fiber
length/sarcomere-series-count per `docs/MECHANISM_FORCE_VELOCITY_POWER_RFD.md`'s own Hill-hyperbola
framework)** — a concrete, falsifiable, machine-checkable asymmetry prediction for a future cert (§9),
not yet tested here. The `Lopt`-unchanged assumption underlying this is a standard textbook premise,
**not independently re-verified live this session** — disclosed, not hidden (§8).

**Concrete re-computed coupling number** (not a prose pointer — machine-computed inline this session):
taking this twin's own already-published `vaslat_r` (vastus lateralis) row from `docs/
MECHANISM_FMAX_PCSA_VALIDATION.md` (`Fmax_model@ST=60 N/cm² = 4262.8 N`, read-only, unchanged) and West et
al. 2010's own measured 15-week training CSA gain (+10% HH / +12% LH, PMID 19910330):

```
Fmax_base = 4262.8 N  (this twin's own vaslat_r, ST=60N/cm^2 convention)
+10% (HH) -> 4689.1 N   (+426.3 N)
+12% (LH) -> 4774.3 N   (+511.5 N)
```

This is a **first-order geometric prediction** (constant specific tension + constant Lopt ⇒ ΔFmax/Fmax ≈
ΔPCSA/PCSA ≈ ΔCSA/CSA), **not** a measured-vs-predicted Fmax validation — West 2010's abstract reports the
CSA gain but not a matching dynamometry Fmax percentage, so the loop is not closed this session (§8). It
demonstrates the coupling is a concrete number, not a hand-wave, and sets up the exact comparison a future
cert would need real per-subject dynamometry to close.

**Pulse-superposition geometry (dose-response mechanism, machine-checked, no fabricated AUC):** each bout
produces a time-localized elevated-MPS pulse (Phillips 1997: +112%→+65%→+34% decaying over 48h). Burd et
al. 2010 show more within-bout volume (1 vs 3 sets) increases **both** the pulse's amplitude (3.1-fold vs
2.3-fold at 5h) **and** its duration (still 2.3-fold at 29h for 3SET vs already back to baseline for 1SET)
— a genuinely measured 2×2 monotonic ordinal fact, not an interpolated curve. **Deliberately not
attempted**: fitting a continuous decay curve and integrating a numeric AUC from only 2-3 abstract-level
points per paper — this would fabricate false precision from too few points with no independent anchor to
check the integral against (the generative-narrative-substitution trap this repo's own memory flags), so
instead the machine-checkable claim is restricted to **sign/direction convergence across two independent
measurement scales** (fine-grained single-bout tracer kinetics vs coarse population dose-response
meta-regression):

```
Burd amplitude (3SET > 1SET @ 5h):  True   (3.1 vs 2.3)
Burd duration  (3SET > 1SET @ 29h): True   (2.3 vs 1.0)
Schoenfeld continuous-sets slope positive & P<0.05: True (P=0.002)
Schoenfeld 2-level categorical P<0.05:              True (P=0.03)
Schoenfeld 3-level categorical P<0.05 (disclosed miss): False (P=0.074)
CROSS-LEVEL CONVERGENCE PASS: True
```

Two decorrelated levels of measurement (single-bout mechanistic tracer kinetics; 15-study population
meta-regression) agree on the **sign** of the dose-response — this is the falsifiable, machine-checkable
geometric claim, not a narrated "more volume is better."

---

## 2. Falsifier 1 in depth — mTORC1 necessity, forced adversary (calcineurin/NFAT) falls

Bodine et al. 2001 (#2) tested **both** major in-vitro-implicated hypertrophy pathways **in the same in
vivo models**: Akt/mTOR was up during hypertrophy, down during atrophy; rapamycin (selective mTOR blocker)
**"blocked hypertrophy in all models tested, without causing atrophy in control muscles."** The forced
adversary, calcineurin/NFAT: **"the calcineurin pathway was not activated during hypertrophy in vivo, and
inhibitors of calcineurin, cyclosporin A and FK506 did not blunt hypertrophy."** Sufficiency+necessity
closes both directions: **"genetic activation of the Akt/mTOR pathway was sufficient to cause hypertrophy
and prevent atrophy in vivo, whereas genetic blockade of this pathway blocked hypertrophy in vivo."** This
is a clean instance-space sweep within one paper (multiple hypertrophy models + multiple atrophy models +
pharmacological blockade + genetic gain-of-function + genetic loss-of-function, one rival pathway directly
falsified in parallel).

Drummond et al. 2009 (#3) replicates the causal-necessity direction in **humans** (a decorrelated species
and a decorrelated modality — pharmacology, not genetics): rapamycin pretreatment **blocked the early (1-2
h) ~40% contraction-induced MPS increase**; **S6K1 phosphorylation rose 6-fold in controls, unchanged with
rapamycin**; **eEF2 phosphorylation (which gates translation elongation) fell ~25% post-exercise in
controls but was unaltered with rapamycin** — i.e., rapamycin blocked the de-repression of elongation, a
downstream mechanistic lock-step with the upstream mTORC1 block, not just a correlated coincidence.
Interesting disclosed nuance the abstract itself flags: rapamycin **also** blocked contraction-induced
ERK1/2 phosphorylation (unexpected — ERK1/2 is canonically MAPK-pathway, not mTORC1-pathway), while eIF4E
(Ser209, a known MNK1 target) phosphorylation was **unaffected** — the authors read this as evidence
rapamycin does not directly inhibit MAPK signalling, i.e., an indirect/crosstalk effect, not a clean single-
node story. Disclosed, not smoothed over.

**Decorrelated corroboration (read-only, sibling node `MOL-MTORC1-NUTRIENT-SIGNALING`, #18-19):** mTORC1
necessity is not a hypertrophy-specific artifact — the **same causal node**, blocked the **same way**
(rapamycin), produces a measurable causal effect on a **completely different organismal phenotype**
(lifespan extension in genetically heterogeneous mice, independently replicated across 3 NIA-ITP test
sites in two separate published cohorts). Different phenotype, different consortium, same drug/target —
this is a genuine additional decorrelation axis for "mTORC1 blockade has real, reproducible, causally
necessary biological consequences," not a re-statement of the hypertrophy finding.

**Verdict: PASS.** Two species (rat, human), two modalities (pharmacological, genetic gain+loss-of-
function), multiple hypertrophy/atrophy induction models, one directly-tested rival pathway falsified in
the same experiments, plus a decorrelated-phenotype corroboration. This is the diverse instance-space the
watertight method requires before accepting C.

---

## 3. Falsifier 2 in depth — the systemic-hormone hypothesis, forced to its strongest (causal) form

### 3.1 The naive hypothesis and why it is intuitive but wrong to test correlationally alone

The folk/bodybuilding-literature hypothesis: the post-exercise testosterone/GH/IGF-1 spike is *why* heavy
compound lifting builds more muscle than isolation work. West & Phillips 2010 (#8) name this directly:
**"a persistent belief (both in scientific literature and among recreational weightlifters)"** — precisely
the kind of intuitive-but-untested claim this method exists to force.

### 3.2 Acute leg — West et al. 2009 (#5)

Within-subject design: one arm's local-only (LH) elbow-flexor exercise vs the other arm's identical
exercise **plus** a superimposed high-volume leg session specifically to spike systemic hormones (HH).
Confirmed manipulation worked: **T, GH significant P<0.001; IGF-1 P<0.05, elevated only in HH.** Outcome:
MPS rose with exercise (P<0.05) but **"with no effect of elevated hormones (P=0.72)"** — and HH's MPS
(0.064 %/h) was **numerically lower**, not higher, than LH's (0.071 %/h). p70S6K phosphorylation rose with
exercise, no between-condition difference. Verbatim: **"Local mechanisms are likely to be of predominant
importance for the post-exercise increase in MPS."**

### 3.3 The decisive causal leg — West et al. 2010 (#6), a genuine controlled manipulation, not a correlation

This is the **pre-registered decisive test**: the same LH/HH split, sustained for **15 weeks of actual
training**, outcome = measured hypertrophy and strength, not just acute MPS. This is a true **factorial
orthogonal contrast** — local mechanical stimulus is held IDENTICAL between arms; only the systemic-hormone
axis is varied — so a hormone effect, if real, has nothing else to hide behind. Result: **CSA +12% (LH)
vs +10% (HH), P<0.001 for the training effect, condition×training interaction P=0.25** — not significant,
and the point estimate runs **opposite** the hypothesis' predicted direction. Type I and type II fiber CSA
both increased with training, no hormone effect on either. Strength increased in both arms, no between-
condition difference. Verbatim conclusion: **"exposure of loaded muscle to acute exercise-induced
elevations in endogenous anabolic hormones enhances neither muscle hypertrophy nor strength."**

### 3.4 Forced independent-lab adversary — Ahtiainen et al. 2003 (#13), and why it does not overturn §3.3

Every citation in §3.2-3.3 (and the #7, #8 syntheses) comes from **one lab** (McMaster University,
Phillips group) — a genuine common-mode risk this document does not paper over. Searching specifically for
a decorrelated replication (different country, different decade-adjacent era, different design) surfaces
Ahtiainen et al. 2003 (Jyväskylä, Finland — fully independent of the McMaster group): in strength-athletes
(n=8), **basal-mean testosterone correlated with strength gain, r=0.84, P<0.01**; in non-athletes (n=8),
**acute training-induced ΔT correlated with CSA gain, r=0.76, P<0.05**. This is a real, live-verified,
disclosed finding that runs counter to §3.2-3.3 — forced to its strongest form, not swept aside.

**Why this adversary survives only partially, not fully:**
1. **Correlational, not causal** — individual differences in T could reflect a third factor (androgen-
   receptor density/sensitivity, training history, genetics) that independently drives both T and
   hypertrophy, rather than T causing growth. West 2010 (#6) is the only design here that actually
   *manipulates* the hormone axis while holding the mechanical stimulus fixed — the correct way to
   separate cause from confound — and it is null.
2. **Small n, plausible multiple-comparisons noise**: the abstract's own text implies up to ~3 hormones ×
   2 outcomes × 2 subgroups ≈ 12 plausible comparisons, of which 2 are reported significant — consistent
   with, but not proof of, chance-level false discovery at n=8/subgroup. The original paper's full
   correlation table was not fetched this session (abstract-only depth, disclosed in §8), so this concern
   is not resolved, only flagged.
3. **Directionally inconsistent between subgroups**: strength-athletes gained *less* CSA (-1.8%,
   nonsignificant direction) and less strength (+3.9%) than non-athletes (+5.6% CSA, +20.9% strength) over
   the *same* 21-week protocol — a training-status interaction that complicates any single T→hypertrophy
   story, since the correlation only appears in different subgroups for different outcomes.
4. West & Phillips 2012's (#7) own large cohort (n=56, McMaster) is **not** a clean zero either — it finds
   weak GH (r=0.28-0.36) and cortisol (r=0.29-0.35) correlations with fiber-CSA gains specifically (not
   lean mass consistently, and **not strength at all**), while T and IGF-1 are flatly null. Cortisol being
   *positively* associated with hypertrophy is itself a reversal of the "catabolic hormone" folk framing —
   a genuinely odd, disclosed residual, not force-fit into either a clean pro- or anti-hormone story.

**Verdict**: the adversary, forced to its strongest CAUSAL, pre-registered form, **FALLS** — §3.3 is a true
manipulated experiment and is a clean, well-powered-for-its-claimed-effect null. A weaker, correlational,
GH/cortisol/trait-level residual **survives partially and is left OPEN**, not force-closed either direction
— this is the symmetric-QC discipline: a rejection of the strong hormone hypothesis does not license
pretending the correlational literature is unanimous, and it isn't.

---

## 4. Volume dose-response (§1's cross-level convergence table is the machine-checked core; narrative context here)

Schoenfeld et al. 2017 (#4, meta-analysis, 34 treatment groups from 15 studies) found a **graded**
dose-response: each additional weekly set per muscle associates with +0.37% muscle-size gain (continuous
model, P=0.002); a coarser lower-vs-higher-volume split shows a 3.9-percentage-point gap (P=0.03); the
finer <5/5-9/10+ sets-per-week 3-level split is only a **trend** (P=0.074) — disclosed as a genuine
non-significant finer-grained miss, not hidden to keep the headline clean. Burd et al. 2010 (#12) supplies
the proximate mechanism at the single-bout level: 3 sets vs 1 set produce a **larger-amplitude AND
longer-duration** MPS elevation (§1), with p70S6K phosphorylation tracking the MPS response (r=0.338,
P=0.033) — a modest but real coupling between the signalling proxy and the outcome it is supposed to
predict.

---

## 5. Dysfunction pole — anabolic resistance (aging) and the mirror-image atrophy program

Cuthbertson et al. 2005 (#9, n=44 young+old) locates the aging deficit precisely: **basal MPS was
indistinguishable by age** — this is not a resting-state deficiency — but **EAA-stimulated MPS was
blunted in the elderly**, tied to **reduced expression/activation (phosphorylation) of mTOR, p70S6k,
4E-BP1, and eIF2B** after amino acid ingestion, an effect **independent of insulin signalling** (plasma
insulin was clamped at basal values, ruling out an insulin-resistance confound) and **associated with
increased NFκB** (an inflammation-linked transcription factor). This is the mechanistic floor under
"anabolic resistance": the same mTORC1 node Falsifier 1 shows is *necessary* for the hypertrophic response
is the node whose *responsiveness* — not its baseline activity — degrades with age.

Bodine et al. 2001's companion Science paper (#14, same lab/year as #2, a disclosed same-lab pairing) is
the **mirror-image molecular machinery**: transcript profiling across multiple rat atrophy models (disuse,
denervation, glucocorticoid-relevant catabolic states) isolated **MuRF1 and MAFbx (atrogin-1)** as the
small, universal subset of E3 ubiquitin ligases activated across *all* atrophy models tested. Overexpressing
MAFbx in myotubes was **sufficient** to produce atrophy; mice genetically lacking either MAFbx or MuRF1 were
**resistant** to atrophy — the same necessity+sufficiency causal structure Falsifier 1 established for the
hypertrophy-building direction, now for the tearing-down direction. This is the molecular substrate for the
"glucocorticoid/cachexia atrophy" pole the task names; a dedicated glucocorticoid-specific human RCT anchor
was not pulled this session (honest gap, §8).

**Couples concretely** to the already-**executed** `MSK-SARCOPENIA-COMPOSITE` node (NHANES n=2810, grip×ALMI
dynapenia dissociation, strength declining faster than mass — a real, machine-verified downstream OUTCOME):
Cuthbertson's blunted-signalling finding is the plausible upstream MECHANISM feeding that composite's
already-measured strength-mass dissociation — a mechanism→outcome pairing, not a duplicate claim. Also
couples to `AUTO-SARCOPENIA-CELL-EWGSOP2-3-DOMAIN-GATE-ST` (SEED-DESIGN, case-definition layer).

---

## 6. Satellite cells — donation for LARGE hypertrophy is a genuinely open dispute, not settled fact

The task's framing ("satellite cells donate myonuclei for large hypertrophy") is **directionally
reasonable but contested at the mechanistic-necessity level** — forcing it surfaced a real, unresolved,
disclosed field dispute, independently re-landed-on here via a fresh NCBI search (not by first reading the
sibling node that already tracks it, below):

- McCarthy et al. 2011 (#10): Pax7-DTA mice, >90% satellite-cell ablation, 2-week synergist-ablation
  overload — **~2-fold mass/CSA gain, indistinguishable from vehicle controls**, WITHOUT satellite cells.
  Myonuclear number did not rise (confirmed via reduced BrdU+ myonuclei) — the myonuclear *domain*
  expanded instead. Fiber hyperplasia/regeneration (a small component) *was* blunted.
- Egner et al. 2016 (#11), forced adversary, same general Pax7-based genetic strategy: **overload
  hypertrophy was PREVENTED** in satellite-cell-deficient mice, in both plantaris and EDL. They attribute
  McCarthy's contrary result to methodological artifacts — reliance on whole-muscle **mass** as a fiber-
  hypertrophy proxy (mass can rise via non-fiber tissue) and inclusion of **regenerating fibers**
  (naturally larger, a different biological process) in the analyzed pool. Verbatim: **"there is currently
  no model in which functional, sustainable hypertrophy has been unequivocally demonstrated in the absence
  of satellite cells."**

**Not re-litigated here**: `data/MECHANISM_ANCHOR_GRAPH.json`'s `MSK-SATELLITE-CELL-RESERVE` node (status
OPEN, read-only, not edited this session) already carries this exact McCarthy-vs-Egner dispute in more
depth, plus Sousa-Victor et al. 2014's causal-rescue geriatric-senescence finding and human correlational
legs (Verdijk 2007, Petrella 2008) this document does not re-derive. This section's role is narrow: confirm
the task's satellite-cell claim needs the "for LARGE/sustained hypertrophy specifically, contested" hedge,
not accept it as settled, and note the independent convergence (this session's fresh search and the
pre-existing node land on the identical two decisive papers).

---

## 7. Couples_to — concrete, read-only, not edited this session

- **`docs/MECHANISM_FMAX_PCSA_VALIDATION.md` / `docs/MECHANISM_SPECIFIC_TENSION.md`** — supplies the
  `Fmax = PCSA × specific_tension` identity and the `vaslat_r` 4262.8 N baseline §1's coupling number
  re-uses directly.
- **`docs/MECHANISM_FORCE_VELOCITY_POWER_RFD.md`** — the parallel-vs-series geometric argument (§1) predicts
  hypertrophy training shifts that doc's Hill-hyperbola F0/Fmax axis specifically, not Vmax; that doc's own
  Aagaard 2002/Tillin 2010/Klass 2008 neural-vs-contractile RFD decomposition is the reason this document
  does not claim a 1:1 CSA-to-strength mapping (§8).
- **`docs/MECHANISM_MUSCLE_DAMAGE.md`** — the eccentric damage-repair-adaptation cycle this task names as
  "IN FLIGHT"; not found as a separate in-progress file this session (only the already-published damage doc
  exists), so coupled to the existing doc only, read-only.
- **`MSK-MUSCLE-HYPERTROPHY-TIMECOURSE`** (graph node, OPEN) — the apparent-CSA-vs-true-myofibrillar-growth
  edema dissociation and training-session-count reframing (Damas et al.); directly qualifies this
  document's §1 pulse-superposition picture — a single bout's MPS pulse must repeat across the ~10-18+
  sessions that node's own Damas 2018 citation identifies before apparent CSA reflects true growth, not
  just swelling.
- **`MOL-MTORC1-NUTRIENT-SIGNALING`** (graph node, OPEN) — the decorrelated lifespan-phenotype corroboration
  used in §2.
- **`MSK-SATELLITE-CELL-RESERVE`** (graph node, OPEN) — the deeper satellite-cell treatment §6 points to.
- **`MSK-SARCOPENIA-COMPOSITE`** (graph node, **executed**, ASSUMED/verified) — the downstream dynapenia
  outcome §5's mechanism feeds.
- **`AUTO-SARCOPENIA-CELL-EWGSOP2-3-DOMAIN-GATE-ST`** (graph node, SEED-DESIGN) — case-definition layer,
  dysfunction pole.

---

## 8. Honest gaps (disclosed, not silently dropped)

1. **Abstract-only depth throughout.** No full-text PDF/PMC fetch was performed this session (time/scope
   budget; the task specifies two files + a summary line, not a full-text deep-dive); every number quoted
   is exactly what each paper's own PubMed abstract states — no fabricated timepoints or interpolated
   values anywhere in this document.
2. **No MPS-pulse AUC was computed** — deliberately. Only 2-4 timepoints per paper are available at
   abstract depth; fitting and integrating a curve through that few points, with no independent measurement
   to check the integral against, would manufacture false precision. The cross-level sign/direction check
   (§1, §4) is the honest substitute.
3. **Ahtiainen 2003's multiple-comparisons concern (§3.4) is flagged, not resolved** — the full paper's
   complete correlation table was not fetched this session.
4. **The satellite-cell McCarthy-vs-Egner dispute (§6) is presented as genuinely open** — this document
   does not adjudicate it, and neither, as far as could be determined without re-deriving it, does the
   field as of the most recent citation found.
5. **The `Lopt`-unchanged assumption (§1) underlying the PCSA/Fmax coupling recompute is a standard
   textbook premise, not independently re-verified live this session.**
6. **The §1 Fmax coupling recompute is one-sided** — West 2010's abstract gives the CSA percentage gain but
   not a matching measured dynamometry Fmax/strength percentage, so this document computes the
   geometrically-predicted Fmax change but does not close the loop against a measured Fmax change from the
   same cohort.
7. **No dedicated glucocorticoid-specific human RCT anchor was pulled** for the "glucocorticoid/cachexia
   atrophy" half of the dysfunction pole (§5) — the MuRF1/MAFbx molecular machinery (#14) is general-atrophy,
   not glucocorticoid-specific; a glucocorticoid-dose-response human study is a natural next fetch.
8. **No standalone script/PREREG/data-pipeline was built** this session, unlike several sibling MECHANISM
   docs — matching this task's explicit two-file+summary-line scope; every number shown is nonetheless
   directly, deterministically reproducible from the inline Python shown in §1.
9. **The fold into `data/MECHANISM_ANCHOR_GRAPH.json` is a natural next step, not performed here** —
   isolation rule ("touch only files you create"); another instance may write concurrently.

---

## 9. Proposed cell (for future MEASURE/CERTIFY promotion)

`id: MSK-HYPERTROPHY-MTOR-VS-HORMONE-GATE` — ingest this document's two falsifiers as a standing
machine-checkable template: (a) any future signalling-pathway claim in this twin must show blockade
BLOCKS the phenotype in the SAME experiment a rival pathway's blockade does NOT, mirroring Bodine 2001's
mTOR-vs-calcineurin structure; (b) any future "circulating factor X drives adaptation Y" claim must be
tested with a LOCAL-STIMULUS-HELD-CONSTANT factorial manipulation (West 2010's design) before a
correlational finding (Ahtiainen 2003's design) is accepted as causal. Natural first measurement: pull
subject2/subject3/subject4's own longitudinal training-load history (if any exists in this twin's data) and
test the §1 Fmax-scales-with-ΔCSA / Vmax-does-not asymmetry prediction against real dynamometry, closing
the one-sided coupling flagged in §8.6.
