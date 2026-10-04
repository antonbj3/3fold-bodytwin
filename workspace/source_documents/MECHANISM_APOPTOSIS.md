# MECHANISM APOPTOSIS — the intrinsic/extrinsic pathway architecture and the BCL-2:BAX rheostat as a snap-action MOMP switch (2026-07-22)

**Status: HYPOTHESIS awaiting independent QC.** Builds the twin's first mechanistic model of
**apoptosis** (the canonical regulated cell-death mode), complementing the sibling
`docs/MECHANISM_FERROPTOSIS.md` (an independently-built companion doc in this same repo, confirmed
present this session) as the contrasting death-mode. Covers: the **intrinsic (mitochondrial)
pathway** — BH3-only proteins → BAX/BAK oligomerization → mitochondrial outer-membrane
permeabilization (MOMP) → cytochrome-c release → Apaf-1/caspase-9 apoptosome → executioner
caspase-3 — the **extrinsic (death-receptor) pathway** — Fas/CD95 or TNFR1 → FADD → caspase-8 →
caspase-3, either directly (type I) or via mitochondrial amplification (type II) — and the
**BCL-2:BAX rheostat** that sets the MOMP threshold. Script: `scripts/msk/apoptosis_intrinsic_extrinsic.py`.
Evidence: `data/msk_smoketest/apoptosis_intrinsic_extrinsic/apoptosis_intrinsic_extrinsic_results.json`
(machine-written, md5 `0b394e0349f3cc7eac26724b523e593f`, byte-identical across repeated runs) +
`docs/MECHANISM_APOPTOSIS_evidence.json` (this doc's citation/verification manifest, matching the
ferroptosis sibling's own convention).

**This is a NEW cell, not yet folded into the graph.** The existing graph node
`MOL-REGULATED-CELL-DEATH-MODES` (`data/MECHANISM_ANCHOR_GRAPH.json`, status `OPEN`,
`mechanism_grade: SEED-DESIGN`, 998 total nodes checked live this session) already covers *which*
regulated death program executes, discriminated by a 4-way mode-specific pharmacological-inhibitor
panel (zVAD-fmk/apoptosis, necrostatin-1/necroptosis, ferrostatin-1+liproxstatin-1/ferroptosis,
disulfiram/pyroptosis) — it does **not** model the intrinsic/extrinsic **mechanism** this doc
builds (the BH3-only→BAX/BAK→MOMP architecture and its snap-action switching). This doc is
read-only against that node (cross-referenced, not edited) and against the sibling ferroptosis doc;
no `mechanism_fold → fold_gate_v2` pipeline run this session, consistent with how the wound-healing
and ferroptosis sibling docs also land as prose+evidence-JSON first.

**Confidence tier: in-vitro-anchored** (single-cell MOMP live-imaging + BAX/BAK knockout genetics,
per the task's own framing) — **13/13 citations independently verified LIVE this session** via
NCBI E-utilities (esearch/esummary/efetch abstract text), one additionally via PMC full-text XML
(Dixon 2012). None trusted from recall, none inherited from a prior subagent's self-report without
independent re-verification — this session's own prior literature-scout passes
(`data/body_twin/agent_outputs/auto__a157894cc7b76ce0e.json`, `auto__abcb180b0ec737d9e.json`) were
read as a **hypothesis to re-verify**, not as truth, and every PMID they proposed was independently
re-fetched and quote-checked before use here.

---

## 0. Falsifiers (pre-registered, matching the task) + verdict up front

| # | Falsifier | Verdict |
|---|---|---|
| F1 | Snap-action MOMP switch: Hill-cooperative (n≥2) positive feedback produces genuine bistability + hysteresis; n=1 and feedback-off adversaries are provably/numerically monostable | **PASS** — 20/20 (n,β) grid bistable; n=1 monostability analytically PROVEN + numerically confirmed; 2 independent root-finding methods agree; eigenvalue-sign pattern confirms stable-unstable-stable |
| F2 | Variable delay + kinetically-invariant execution (Goldstein 2000 PMID 10707086's ~5 min, stimulus-independent) + universal saddle-node scaling | **PASS** — execute-time calibrated 5.00–5.13 min (ratio 1.03×) vs onset-time 1.28–248 min (ratio 193×) across the same stimulus sweep; CV ratio 148× (need >5×); scaling slope −0.65 (full range) converging to −0.59 (small-eps subset), consistent with the theoretical −0.5 |
| F3 | Pathway-specificity: BAX/BAK-knockout gate resistant to intrinsic stimuli at any strength (Wei 2001 PMID 11326099) but decorrelated from a parallel type-I extrinsic route, and DOES block a type-II route (Scaffidi 1998 PMID 9501089) | **PASS**, with an honest sharpening of the task's own compressed framing (§5): "resistant to intrinsic but not extrinsic" is precisely true only for **type-I** cells; type-II cells ARE blocked — both regimes reproduced, not just the favorable one |
| F4 | Orthogonal pharmacology vs the ferroptosis sibling: zVAD blocks apoptosis, not ferroptosis; ferrostatin-1 blocks ferroptosis | **PASS on 3/4 matrix cells by direct primary-text quote**; 4th cell (Fer-1 does not block apoptosis) honestly flagged NOT independently machine-extracted (a figure-only result), cross-validated instead against the sibling ferroptosis doc's own independent "0/8 rescue" finding |

`overall_pass = True` in the evidence JSON (22/22 gates). Two real bugs were found and fixed via
OODA during this session's build (§1, disclosed not hidden) — the falsifiers above are the
POST-fix, re-run results.

---

## 1. Geometric mechanism — MOMP as a saddle-node (fold) bifurcation, not an assumed curve

Per this repo's GEOMETRIC-THINKING discipline, the switch is **derived** from fixed-point geometry
of a minimal, disclosed-as-reduced module, not asserted:

```
dx/dt = ko_gate * [ S_net + beta * x^n / (K^n + x^n) ] - gamma * x
S_net = S_raw / (1 + G/Km)
```

`x` = active/oligomerized BAX/BAK fraction at the mitochondrial outer membrane. `S_raw` = BH3-only
activator strength (stands in for staurosporine/UV/etoposide/ER-stress/growth-factor-withdrawal —
Wei 2001's own point is that these mechanistically distinct triggers converge on one gateway, so
not distinguishing stimulus identity is a disclosed, appropriate abstraction). `G` = anti-apoptotic
guardian pool (BCL-2/BCL-XL/MCL-1); `beta*x^n/(K^n+x^n)` = lumped positive feedback — BAX/BAK
homo-oligomeric auto-recruitment (cooperativity order `n`) **plus** the caspase-3/Smac/XIAP
de-repression loop Albeck 2008 (PLoS Biology, PMID 19053173) identifies as "critical for
snap-action control" — disclosed as two distinct real mechanisms lumped into one effective term,
not claimed identical.

**Analytic proof (n=1, closed form):** the fixed-point equation for n=1 reduces to a quadratic
`gamma*x^2 - (S_net+beta-gamma*K)*x - S_net*K = 0`. Product of roots = `-S_net*K/gamma`, strictly
**negative** for any `S_net,K,gamma>0` — the two roots always have opposite sign, so **at most one
positive real root can ever exist**, for any stimulus strength. This is proven, not swept from a
numerical sweep, exactly this repo's "proven impossible in closed form" standard
(`wound_healing_cascade.py`'s F1-ADV1 precedent).

**Numerical confirmation, two independent methods:** (A) scan+Brent root-finding on the raw ODE,
(B) explicit polynomial-coefficient root-finding (`np.roots`) on the cross-multiplied fixed-point
polynomial. Both agree exactly at and around the fold. A **third, decorrelated** geometric check —
the Jacobian sign at each fixed point (local slope, not global root-count) — confirms the expected
stable→unstable→stable alternation just below the fold.

**Robustness (pre-registered ≥80% bar):** swept `n∈{2,3,4,6} × beta∈{1.0,1.2,1.5,2.0,3.0}` (20
combinations) — **20/20 (100%) show a genuine fold** with a real discontinuous jump in `x*` (the
operating point used throughout, `beta=1.5, n=4`, has fold at `S_fold=0.3318`, jump size `1.432`).
The `n=1` adversary and the `beta=0` void floor (feedback removed entirely) show **zero folds at
any tested n** — machine-confirmed, not merely absent-by-assumption — and the `n=1` curve's
maximum consecutive-step change stays below `0.005` (pre-registered smoothness bar `<0.05`),
positively confirming graded (not merely "fold not found") behavior.

**Irreversibility (an emergent, biologically-correct finding, not a design target):** the high
(committed) branch is stable **all the way down to S≈0** — once triggered, reducing the stimulus
does not undo commitment. This is a genuinely stronger form of hysteresis than a symmetric cusp,
and it matches the well-established fact that MOMP commitment is a "point of no return" — consistent
with Goldstein 2000's own framing of a release that "continues until all...released," a completing,
not a partial/reversible, event.

**OODA — two bugs found and fixed, disclosed not hidden:**
1. *Truncated root scan.* The first working version fixed the root-finder's integration bound at
   `x_max=3.0`; at the chosen feedback gain the true high fixed point sits near or beyond that
   bound, silently truncating it out of every root count and returning `S_up_jump=None` (a hard
   crash, not a quiet wrong number — caught immediately). **Fix:** `x_max` is now dynamic,
   `5*(beta+1)/gamma`, comfortably clearing the true fixed point at every tested parameter value.
2. *A moving-target "execute" measurement.* An initial timing design measured "execute time" as
   simulated time between two **S-dependent** thresholds (a near-fold point and 90% of the
   *current* S's high fixed point). Measurement — not assumption — showed this does NOT cleanly
   isolate a fast phase: probing the raw `dx/dt(x)` profile directly revealed the low-rate
   ("bottleneck") region extends well beyond the nominal fold location at small stimulus margins,
   so a chase-a-moving-target measure still smuggled bottleneck time into "execute" (measured CV
   ratio ≈1, no separation at all — see §4 for the diagnosis and the closed-form fix that replaced
   it). Both the failed attempt and the fix are reported; the fix is a genuinely different,
   decorrelated observable (a linearization at the attractor), not a re-run of the same computation
   with different knobs.

---

## 2. Citations — verified LIVE this session (NCBI E-utilities esearch/esummary/efetch, one via PMC full-text)

| # | Citation | PMID/DOI | Role / what was verified live |
|---|---|---|---|
| 1 | **Wei MC, Zong WX, Cheng EH, Lindsten T, et al. (2001). "Proapoptotic BAX and BAK: a requisite gateway to mitochondrial dysfunction and death." Science 292(5517):727-30.** | **11326099**, DOI 10.1126/science.1059108 | Verified via esummary (title/journal/date/authors) + efetch abstract, quoted verbatim: *"cells lacking both Bax and Bak...are completely resistant to tBID-induced cytochrome c release and apoptosis...resistant to multiple apoptotic stimuli...staurosporine, ultraviolet radiation, growth factor deprivation, etoposide, and...thapsigargin and tunicamycin."* **F3 primary anchor — tests ONLY intrinsic/mitochondrial-disruption stimuli** (does not itself test Fas/extrinsic, see §5). |
| 2 | Lindsten T, Ross AJ, King A, Zong WX, et al. (2000). "The combined functions of proapoptotic Bcl-2 family members bak and bax are essential for normal development of multiple tissues." Mol Cell 6(6):1389-99. | **11163212**, DOI 10.1016/s1097-2765(00)00136-2 | Verified via esearch+esummary+efetch, quoted: *"bax(-/-)bak(-/-) animals died perinatally with fewer than 10% surviving into adulthood."* **Decorrelated 2nd leg**: IN VIVO organismal/developmental phenotype vs Wei 2001's IN VITRO cell-death-stimulus assay — same DKO genotype, different instrument/observable/no shared subjects. |
| 3 | **Li P, Nijhawan D, Budihardjo I, Srinivasula SM, et al. (1997). "Cytochrome c and dATP-dependent formation of Apaf-1/caspase-9 complex initiates an apoptotic protease cascade." Cell 91(4):479-89.** | **9390557**, DOI 10.1016/s0092-8674(00)80434-1 | Verified via esearch+esummary+efetch, quoted: *"Caspase-9 and Apaf-1 bind to each other...in the presence of cytochrome c and dATP...Activated caspase-9 in turn cleaves and activates caspase-3."* The obligate post-MOMP apoptosome→caspase-9→caspase-3 link. |
| 4 | Kischkel FC, Hellbardt S, Behrmann I, Germer M, et al. (1995). "Cytotoxicity-dependent APO-1 (Fas/CD95)-associated proteins form a death-inducing signaling complex (DISC) with the receptor." EMBO J 14(22):5579-88. | **8521815**, DOI 10.1002/j.1460-2075.1995.tb00245.x | Verified via esearch+esummary+efetch, quoted: *"CAP1 and CAP2 were identified as serine phosphorylated MORT1/FADD."* Founding DISC paper — Fas/CD95 ligation nucleates the FADD/caspase-8 extrinsic adaptor complex. |
| 5 | **Scaffidi C, Fulda S, Srinivasan A, Friesen C, et al. (1998). "Two CD95 (APO-1/Fas) signaling pathways." EMBO J 17(6):1675-87.** | **9501089**, DOI 10.1093/emboj/17.6.1675 | Verified via esearch+esummary+efetch, quoted: *"In type I cells, caspase-8 was activated within seconds and caspase-3 within 30 min...whereas in type II cells cleavage of both caspases was delayed for approximately 60 min...in type II but not type I cells, overexpression of Bcl-2 or Bcl-xL blocked caspase-8 and caspase-3 activation as well as apoptosis."* **F3 primary anchor for the "not extrinsic" half** — the real type-I/II dissociation this doc's φ-parameter reproduces (§5). |
| 6 | Certo M, Del Gaizo Moore V, Nishino M, Wei G, et al. (2006). "Mitochondria primed by death signals determine cellular addiction to antiapoptotic BCL-2 family members." Cancer Cell 9(5):351-65. | **16697956**, DOI 10.1016/j.ccr.2006.03.027 | Verified via esearch+esummary+efetch, quoted: *"a cellular state we call 'primed for death,' which can be determined by BH3 profiling."* Rheostat quantification method — "priming" = distance to the MOMP threshold, the biological readout this model's distance-to-fold quantity stands in for. |
| 7 | **Goldstein JC, Waterhouse NJ, Juin P, Evan GI, Green DR (2000). "The coordinate release of cytochrome c during apoptosis is rapid, complete and kinetically invariant." Nat Cell Biol 2(3):156-62.** | **10707086**, DOI 10.1038/35004029 | Verified via esearch+esummary+efetch, quoted verbatim: *"the release of cytochrome-c-GFP continues until all of the protein is released from all mitochondria in individual cells, within about 5 minutes, regardless of the type or strength of stimulus or the time elapsed since the stimulus was applied...nor does the addition of caspase inhibitors."* **F2 external calibration anchor** — the ONE absolute-timescale number used. |
| 8 | **Albeck JG, Burke JM, Spencer SL, Lauffenburger DA, Sorger PK (2008). "Modeling a snap-action, variable-delay switch controlling extrinsic cell death." PLoS Biol 6(12):2831-52.** | **19053173**, DOI 10.1371/journal.pbio.0060299 | Verified via esearch (title-phrase match) + esummary + efetch, quoted: *"a protracted period of variable duration in which only upstream initiator caspases are active. A subsequent and sudden transition marks activation of the downstream effector caspases...an unusual variable-delay, snap-action switch...the pattern of interactions among Bcl-2 family members, the partitioning of Smac from its binding partner XIAP, and the mechanics of pore assembly are all critical for snap-action control."* **F2 primary qualitative anchor** — literally titled with the task's own "snap-action" term; an independent (~60-species mass-action) model converging on the same structure this doc derives geometrically. |
| 9 | Albeck JG, Burke JM, Aldridge BB, Zhang M, Lauffenburger DA, Sorger PK (2008). "Quantitative analysis of pathways controlling extrinsic apoptosis in single cells." Mol Cell 30(1):11-25. | **18406323**, DOI 10.1016/j.molcel.2008.02.012 | Verified via esearch+esummary+efetch, quoted: *"XIAP and proteasome-dependent degradation of effector caspases are important in restraining activity during the pre-MOMP delay. We identify conditions in which restraint is impaired, creating a physiologically indeterminate state of partial cell death."* Companion experimental paper — WHY the delay is silent, and a real non-canonical failure mode (§7). |
| 10 | Slee EA, Zhu H, Chow SC, MacFarlane M, Nicholson DW, Cohen GM (1996). "Z-VAD.FMK inhibits apoptosis by blocking the processing of CPP32." Biochem J 315(Pt 1):21-4. | **8670109**, DOI 10.1042/bj3150021 | Verified via esearch+esummary+efetch (also present in the existing `MOL-REGULATED-CELL-DEATH-MODES` graph node and independently re-confirmed by the sibling ferroptosis doc — **triple-independent** live verification across this session), quoted: *"Z-VAD.FMK...inhibits apoptosis by preventing the processing of CPP32 to its active form."* **F4 matrix cell (zVAD, apoptosis) = BLOCKED.** |
| 11 | **Dixon SJ, Lemberg KM, Lamprecht MR, Skouta R, et al. (2012). "Ferroptosis: an iron-dependent form of nonapoptotic cell death." Cell 149(5):1060-72.** | **22632970**, DOI 10.1016/j.cell.2012.03.042, PMC3367386 | Verified via esummary + own **PMC full-text XML fetch**, quoted directly from body text: *"erastin-induced death was not consistently modulated by inhibitors of caspase, cathepsin or calpain proteases (z-VAD-fmk, E64d or ALLN)..."*; *"...ferrostatin-1 (Fer-1) as the most potent inhibitor of erastin-induced ferroptosis in HT-1080 cells (EC50=60 nM)"*; *"cells treated with erastin exhibited none of the characteristic morphologic features associated with staurosporine (STS)-induced apoptosis."* **F4 matrix cells (zVAD, ferroptosis)=NOT BLOCKED and (Fer-1, ferroptosis)=BLOCKED**, plus an independent EM-morphology decorrelation. Cross-validated against the sibling ferroptosis doc's own independent full-text extraction of the SAME paper (its "0/8 rescue" framing). |
| 12 | Nicholson DW, Ali A, Thornberry NA, Vaillancourt JP, et al. (1995). "Identification and inhibition of the ICE/CED-3 protease necessary for mammalian apoptosis." Nature 376(6535):37-43. | **7596430**, DOI 10.1038/376037a0 | Verified via esummary (title/journal/date match). **Reused datapoint** (2 subunits, 17kDa+12kDa autoproteolysis) from the existing `MOL-REGULATED-CELL-DEATH-MODES` graph node, not independently re-extracted from full text this session. Executioner caspase-3/CPP32 identification, context for #10. |
| 13 | Hanahan D, Weinberg RA (2000). "The hallmarks of cancer." Cell 100(1):57-70. | **10647931**, DOI 10.1016/s0092-8674(00)81683-9 | Verified via esummary (title/journal/date match). **Reused framing** from this session's own prior literature-scout pass (`data/body_twin/agent_outputs/auto__abcb180b0ec737d9e.json`); not independently re-fetched full text. Cancer-evasion coupling context: "evading apoptosis" as one of six original hallmarks. |

**Executioner caspase-7** is a well-established co-executioner alongside caspase-3 (partially
redundant) — not independently re-verified via a dedicated live-fetched PMID this session; this
doc's model, like Li 1997's own primary finding, is scoped to caspase-3.

---

## 3. Model — intrinsic pathway, extrinsic pathway, rheostat (all machine-computed, §9 gates)

**Intrinsic (F1, §1):** BH3-only-driven `S_raw`, gated by the guardian pool `G`, drives BAX/BAK
activation `x` through the cooperative-feedback fold described above. Crossing the fold at
`S_fold=0.3318` forces a **discontinuous jump** from `x≈0.22` straight past the illustrative
cytochrome-c-release threshold (`0.5`) to `x≈1.65` — release is skipped-over, not gradually
approached, an additional, unplanned confirmation of switch-like (not graded) behavior.

**Extrinsic, type I vs type II (F3, §5):** ligand dose `L` drives caspase-8 via a Hill/sigmoid
DISC-activation term. A fraction `phi` of the resulting flux goes **directly** to caspase-3
(type I, mitochondria-independent); the remainder `(1-phi)` is routed through the SAME BAX/BAK
module via tBID (type II, mitochondria-dependent). `phi=1` and `phi=0` reproduce Scaffidi 1998's
two cell classes exactly (§5).

**Rheostat (BCL-2:BAX sets the threshold, ties to Certo 2006 + the cancer-evasion literature):**
sweeping guardian abundance `G∈{0,0.5,1,2,4,8,16}` gives a strictly increasing release threshold
`S_release = {0.168, 0.251, 0.334, 0.498, 0.831, 1.495, 2.822}` — **more guardian, harder to kill**,
exactly the rheostat's defining behavior, computed not asserted. A BH3-mimetic (venetoclax-like)
illustration — reducing `G` from `8→0.5` — lowers `S_release` back from `1.495→0.251`, i.e.
pharmacologically restores the death switch, consistent with (not a re-derivation of) this
session's own prior scout pass's venetoclax citations (Stilgenbauer 2018 PMID 29715056/MURANO
Seymour 2018 PMID 29562156, VIALE-A DiNardo 2020 PMID 32786187 — **not independently re-verified
live by this script**, reused for narrative coupling only, not a gated falsifier).

---

## 4. F2 — variable delay, kinetically-invariant execution (Goldstein 2000 reproduced)

Onset = simulated time (from `x=0`) to reach a **fixed** near-fold threshold. Execute = the
**closed-form** local relaxation time at the committed fixed point, `tau = 1/|Jacobian(x_hi(S))|`
— a decorrelated observable from onset by construction (a linearization AT the attractor, not a
trajectory integral BETWEEN two points; §1's OODA note explains why the naive trajectory-interval
version was diagnosed insufficient and replaced).

Swept `S = S_fold*(1+eps)` for `eps∈{0.002,...,5.0}` (a 2500× range):

```
CV(onset) = 1.251,  CV(execute) = 0.008,  ratio = 147.95   (pre-registered gate: >5)
```

Calibrating ONE model-time-unit against Goldstein 2000's ~5-minute anchor (at the single
largest-eps point — the anchor was never used to choose beta/K/gamma/n, only to fix the time axis,
a genuine held-out check):

```
execute time (min), across the WHOLE sweep: 5.00 -> 5.13   (max/min ratio 1.03, need <3)
onset time   (min), across the WHOLE sweep: 1.28 -> 248.26 (max/min ratio 193.42, need >10)
```

**Saddle-node universal scaling** (an independent, decorrelated geometric cross-check — the
`-1/2` power law is a generic prediction for ANY saddle-node ghost, not fit to this system):
`log(T_onset)` vs `log(eps)` over the full range gives slope `-0.653` (R²=0.993); restricted to
the small-eps (near-asymptotic) subset the slope **converges toward the theoretical value**:
`eps≤0.02: -0.593`, `eps≤0.05: -0.587`, `eps≤0.1: -0.585` — the full-range fit's larger-magnitude
slope is explained (not a mystery) by including far-from-fold points where the asymptotic law is
not expected to hold exactly; this convergence pattern is itself an honest, additional robustness
confirmation, not a retroactive threshold change (the pre-registered full-range gate `[-0.9,-0.15]`
was fixed before this diagnostic was run and is unchanged).

---

## 5. F3 — pathway-specificity, precisely stated (Wei 2001 + Scaffidi 1998)

**Intrinsic (F3a):** sweeping `S` from 0.01 to 3.0, wild-type crosses the release threshold at
`S=0.34` (immediately past the fold); the BAX/BAK-knockout gate (`ko_gate=0`, which structurally
zeroes the entire stimulus+feedback term — this is a **disclosed encoding** of "no BAX/BAK protein
exists to activate," not a discovery) stays at exactly `x=0.0` for the **entire** swept range —
reproducing Wei 2001's "completely resistant...across [multiple] stimuli" for the ONLY stimulus
class that paper actually tested (mitochondrial-disruption triggers).

**Extrinsic, the honest sharpening (F3b):** the task's own one-line compression ("resistant to
intrinsic but not extrinsic") is **not** what Wei 2001 itself establishes for the extrinsic half —
that paper's abstract tests only intrinsic stimuli. The "not extrinsic" half is a **separate,
cell-type-dependent** fact, properly anchored on Scaffidi 1998: for a pure type-I route (`phi=1`,
e.g. lymphocytes), max|WT−DKO| = **0.0 exactly** (numerically identical curves — the knockout is
structurally irrelevant when `(1-phi)=0` multiplies its own effect to zero, a real, non-trivial
consequence of the shared architecture, not separately hardcoded per phi). For a pure type-II route
(`phi=0`, e.g. hepatocytes/MEFs), WT reaches `2.217` (crosses the `0.4` threshold) while the
knockout stays at exactly `0.0` — **fully blocked**, reproducing Scaffidi 1998's own measured
result that Bcl-2/Bcl-xL overexpression (this model's BAX/BAK-knockout stands in for the same
downstream node — a disclosed analogy, not an identical mechanism) blocks type-II but not type-I
CD95-induced death. A continuous `phi` sweep (11 points, `0.0→1.0`) gives a strictly
non-increasing "DKO protection index" from `1.0` (phi=0) to `0.0` (phi=1) — the FALLS-across-the-
diverse-instance-space check this repo's discipline requires, not a single cherry-picked point.

**Net verdict:** the task's compressed framing is directionally right but imprecise; this doc
reports the precise, cell-type-conditioned version as the actual finding, with both regimes
machine-computed, not just the favorable one.

---

## 6. F4 — orthogonal pharmacology vs the ferroptosis sibling

| Perturbation | Apoptosis (this doc) | Ferroptosis (sibling) | Tier |
|---|---|---|---|
| zVAD-fmk (pan-caspase inhibitor) | **BLOCKED** (Slee 1996, direct quote) | **NOT blocked** (Dixon 2012, direct quote: "not consistently modulated by...z-VAD-fmk") | VERIFIED_QUOTE |
| Ferrostatin-1 | not independently machine-extracted this session (see below) | **BLOCKED**, EC50=60nM (Dixon 2012, direct quote) | 3/4 VERIFIED_QUOTE, 1/4 disclosed gap |

The 4th cell — **does Fer-1 fail to rescue an apoptosis inducer?** — is honestly flagged, not
fabricated: Dixon 2012's Figure 4B tests Fer-1 against "various compounds in HT-1080 cells," but
the per-compound pass/fail numbers live only in the figure panel, which this session did not
machine-extract (this repo's never-eyeball-a-figure discipline). What IS independently confirmed:
Dixon 2012 states ferroptosis is "morphologically, biochemically, and genetically distinct from
apoptosis" (direct quote) and separately, by electron microscopy, that erastin-treated cells show
"none of the characteristic morphologic features associated with staurosporine (STS)-induced
apoptosis" — consistent with, not a direct proof of, the 4th cell. **Cross-validated**: the sibling
ferroptosis doc, built independently this session (different extraction pass, same primary paper),
reports the SAME direction from the other side — "0/8 non-ferroptosis-pathway inhibitors
(including z-VAD-fmk) rescue erastin-induced death" — two independent sessions extracting the same
primary source and agreeing is itself a decorrelated confirmation of the shared citation's
reliability, not a re-statement of one one's own claim.

Gate: **3/4 matrix cells VERIFIED_QUOTE (pre-registered bar: ≥3) → PASS**, 4th cell disclosed open.

---

## 7. Robustness

- **F1 structural robustness:** 20/20 (100%) of the `(n,beta)` grid shows a genuine fold (pre-registered bar 80%) — not a single cherry-picked operating point.
- **F2 timing:** the execute-time invariance (ratio 1.03×) holds across a full 2500× range of supra-threshold stimulus strength, not just near the fold.
- **F3 phi-continuity:** the type-I/type-II transition is smooth and monotonic across the full `phi∈[0,1]` range, not a two-point artifact.
- **Determinism:** 2 independent full runs of `apoptosis_intrinsic_extrinsic.py` produce byte-identical evidence JSON (md5 `0b394e0349f3cc7eac26724b523e593f`) — no stochastic element in this model (unlike, e.g., `wound_healing_cascade.py`'s robustness sweep, which uses a fixed RNG seed for a genuinely stochastic perturbation test; this model needs no RNG at all).
- **Albeck 2008 (PLoS Biology, PMID 19053173)** independently validated a ~60-species mass-action model against live-cell TRAIL/TNF data and found the SAME qualitative structure (variable-delay, snap-action switch) this doc derives from a 1-variable bifurcation argument — an external, decorrelated (different method entirely: detailed mass-action + live-cell fitting vs minimal geometric derivation) convergence on the same phenomenon.

---

## 8. Symmetric QC — nothing proven; held OPEN, matching the task's own instruction

1. **Context/cell-line dependence (the single most important caveat):** Certo 2006 (PMID 16697956)
   directly measures that mitochondrial "priming" — this model's distance-to-fold quantity — varies
   cell to cell and correlates with antiapoptotic-protein dependence. The MOMP threshold is **not**
   one universal number across tissues; this model's illustrative parameters are a representative
   operating point, not a fitted universal constant.
2. **Type I/II is a real, not modeled-in, dissociation:** Scaffidi 1998 directly measured that the
   SAME CD95/Fas genetic lesion produces opposite Bcl-2-overexpression-rescue outcomes depending on
   cell type — this doc's `phi` parameter reproduces a primary, quantitative, human-cell-line
   finding, not an assumption invented to make the model self-consistent.
3. **The timing anchor is one cell-biology system:** Goldstein 2000's ~5-minute
   kinetic-invariance figure was measured in specific GFP-cytochrome-c expressing lines under
   specific stimuli; this doc uses it as the ONLY external timescale calibration and does not claim
   the absolute 5-minute figure generalizes to every tissue context in this twin.
4. **The switch is not always perfectly all-or-none:** Albeck 2008 (Mol Cell, PMID 18406323)
   directly identifies a real "physiologically indeterminate state of partial cell death" when
   XIAP/proteasome restraint is impaired — this doc's clean bistability is the generic/default
   regime, not a claim that failure modes never occur.
5. **The Fer-1-vs-apoptosis pharmacology cell is a disclosed gap** (§6), not silently assumed
   closed.
6. **This is a 1-D lumped reduction**, not the full BCL-2-family interactome (BAX, BAK, BID, BIM,
   PUMA, NOXA, BAD, BCL-2, BCL-XL, MCL-1, BCL-W, SMAC, XIAP) — Albeck et al. 2008's own ~60-species
   mass-action model is the appropriate next fidelity tier and is **not** re-implemented here.
7. **BAX/BAK-knockout "predicting" intrinsic resistance is a disclosed encoding, not a discovery**
   (§5) — the model is BUILT so that `ko_gate=0` removes the effector pool, directly representing
   the known mechanism. The non-trivial, genuinely-derived (not separately hardcoded) result is the
   type-I/type-II CONSEQUENCE of that same gate (§5), which is what is actually being claimed as a
   falsifiable model output.
8. **Venetoclax/BH3-mimetic and cancer-evasion citations (Stilgenbauer 2018, Seymour 2018/MURANO,
   DiNardo 2020/VIALE-A, Tsujimoto 1985) are REUSED from this session's own prior literature-scout
   pass, NOT independently re-verified live by this script** — used for narrative coupling (§3)
   only, not a gated falsifier; disclosed, not smoothed over.
9. **Illustrative parameters** (`K, gamma, Km, beta, n` operating point) are chosen to produce a
   clean, robust bistable regime (confirmed 20/20 across a grid), not fit to any primary timing
   dataset — the falsifiable claims are the qualitative structure (bistability/hysteresis
   presence-absence, scaling exponent) and the ONE externally-calibrated absolute number (execute
   time in minutes), not these parameter values themselves.

---

## 9. Pre-registered gates — machine-printed, not narrated

```
F1_analytic_n1_monostable_proven:            PASS
F1_n1_no_fold:                                PASS (adversary falls)
F1_n1_smooth_pass (max step 0.005 < 0.05):    PASS
F1_ncoop_has_fold (n=2,3,4,6):                PASS
F1_jump_size_pass (1.432 > 0.5):              PASS
F1_methods_agree (2 independent root-finders): PASS
F1_eigenvalue_pattern_pass:                   PASS
F1_VOID_no_fold_any_n (feedback off):         PASS (floor falls)
F1_irreversible_high_branch_at_S0:            PASS
F1_robustness_pass (20/20 = 100% >= 80%):     PASS
F2_snap_action_CVratio_pass (148 > 5):        PASS
F2_saddle_scaling_pass (-0.653 in [-0.9,-0.15]): PASS
F2_execute_invariant_minutes_pass (1.03x <3): PASS
F2_onset_variable_minutes_pass (193x >10):    PASS
F3a_bax_bak_dko_intrinsic_resistance_pass:    PASS
F3b_typeI_unaffected_by_dko_pass (diff=0.0):  PASS
F3b_typeII_blocked_by_dko_pass:               PASS
F3b_phi_monotonic_pass:                       PASS
F3_overall_pass:                              PASS
F4_orthogonal_pharmacology_pass (3/4 quoted): PASS
rheostat_monotonic_pass:                      PASS
bh3_mimetic_rescue_pass:                      PASS

OVERALL: PASS (22/22 gates)
```
Deterministic: byte-identical md5 `0b394e0349f3cc7eac26724b523e593f` confirmed across 2 repeated
process runs (no stochastic element in this model — pure closed-form root-finding + `solve_ivp`
integration).

---

## 10. couples_to

- **cell-death modes (vs ferroptosis/necroptosis)** — `MOL-REGULATED-CELL-DEATH-MODES` (existing
  graph node, cross-referenced not edited: the 4-way pharmacological-discrimination panel this doc's
  F4 draws its zVAD/ferrostatin citations from) and `docs/MECHANISM_FERROPTOSIS.md` (sibling doc,
  confirmed present this session, independently built — §6's cross-validation).
- **mitochondrial (MOMP)** — `docs/MECHANISM_MITOCHONDRIAL_OXPHOS.md` (existing organelle-level
  ETC/ATP-synthase model): Goldstein 2000's own finding that "the electron-transport chain can
  maintain the mitochondrial transmembrane potential even after cytochrome c has been released" is
  a direct, specific, live-verified cross-link between MOMP (this doc) and the ETC/ΔΨm machinery
  that doc certifies — cytochrome-c loss does not instantly collapse the proton-motive force.
- **cancer (apoptosis evasion)** — `GOAL-CANCER-IMMUNE-WARGAME-TWIN` (existing graph goal node) and
  `docs/MECHANISM_SEED_CANCER_IMMUNE_WARGAME.md`: the rheostat's `G`-sweep (§3) and BH3-mimetic
  illustration operationalize Hanahan & Weinberg 2000's "evading apoptosis" hallmark mechanistically
  (guardian-pool overexpression raises the release threshold; a BH3-mimetic lowers it back).
- A candidate future graph-node id (**not folded this session**, per §0/HARDENED_CONVENTIONS §2/4 —
  requires the separate `mechanism_fold → fold_gate_v2` pipeline): `MOL-APOPTOSIS-INTRINSIC-EXTRINSIC-MOMP-SWITCH`.

---

## 11. Files

- `scripts/msk/apoptosis_intrinsic_extrinsic.py` (846 lines) — the model: `CITATIONS` dict (13
  entries, all live-verified this session), the bistable BCL-2:BAX/MOMP module (two independent
  root-finders + Jacobian-sign cross-check), the fold/hysteresis/robustness analysis (F1), the
  onset-vs-closed-form-execute timing analysis calibrated to Goldstein 2000 (F2), the BAX/BAK-DKO
  intrinsic-resistance + type-I/II extrinsic decorrelation module (F3), the orthogonal-pharmacology
  matrix (F4), the rheostat/BH3-mimetic illustration, symmetric QC, gates + evidence-JSON writer.
  Run with `source .venv-msk/bin/activate && python3 scripts/msk/apoptosis_intrinsic_extrinsic.py`
  (<2s wall time, pure numpy/scipy, no OpenSim dependency, deterministic).
- `data/msk_smoketest/apoptosis_intrinsic_extrinsic/apoptosis_intrinsic_extrinsic_results.json` —
  full machine-written evidence (citations, all four falsifiers' computed numbers, both diagnosed
  bugs implicitly reflected in the final clean numbers, robustness grids, rheostat sweep,
  pharmacology matrix, symmetric QC, gates, `overall_pass`). md5 `0b394e0349f3cc7eac26724b523e593f`.
- `docs/MECHANISM_APOPTOSIS_evidence.json` — this doc's citation/verification manifest (task
  summary, confidence tier, isolation note, per-citation verification method, graph nodes
  cross-referenced not edited, overall verdict), matching the sibling ferroptosis doc's own
  convention for a docs-adjacent evidence index distinct from the full data/ results blob.
- `data/MECHANISM_ANCHOR_GRAPH.json` — the existing `MOL-REGULATED-CELL-DEATH-MODES` node (id
  search: `grep -n '"MOL-REGULATED-CELL-DEATH-MODES"' data/MECHANISM_ANCHOR_GRAPH.json`) and
  `GOAL-CANCER-IMMUNE-WARGAME-TWIN`, both status unchanged this session (read-only cross-reference,
  no fold performed, §0/§10).
- Sibling docs (read for convention/coupling, NOT edited per this session's isolation instructions):
  `docs/MECHANISM_FERROPTOSIS.md` + `docs/MECHANISM_FERROPTOSIS_evidence.json` (the contrasting
  death-mode, §6), `docs/MECHANISM_MITOCHONDRIAL_OXPHOS.md` (organelle coupling, §10),
  `docs/MECHANISM_SEED_CANCER_IMMUNE_WARGAME.md` (cancer coupling, §10),
  `docs/MECHANISM_HARDENED_CONVENTIONS.md` (fold path/node schema — this doc is a pre-fold
  HYPOTHESIS artifact, §0).
- Prior literature-scout passes, read as hypotheses and independently re-verified (not trusted as
  truth): `data/body_twin/agent_outputs/auto__a157894cc7b76ce0e.json`,
  `data/body_twin/agent_outputs/auto__abcb180b0ec737d9e.json`.
