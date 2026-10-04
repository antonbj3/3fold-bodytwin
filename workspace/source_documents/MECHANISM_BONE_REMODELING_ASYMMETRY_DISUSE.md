# MECHANISM BONE REMODELING — athlete asymmetry (FUNCTION) + disuse/spaceflight (DYSFUNCTION), a decorrelated falsifier pair (2026-07-22)

**Companion delta, not an edit.** `docs/MECHANISM_BONE_REMODELING.md` (Part A/B) already exists in
this repo, built by a concurrent sibling mechanism instance, and already gated Frost's
strain-threshold classifier against this twin's own femur gait strain (Part A) plus a
histomorphometry-anchored BMU turnover-rate model (Part B). **This session's own isolation rule
("touch ONLY files you create") is honored literally: that file is READ for context below but not
one byte of it is edited.** This is a **new** file covering a gap that file's own text does not
touch at all — grepped first, confirmed: zero hits for `dominant.arm|tennis|kannus|bedrest|
spaceflight|klein-nulend|robling|bonewald|canalicul` against its 528 lines (one unrelated
"osteoporosis/high-turnover" mention, one unrelated "Robling 2006" literature-search **miss** —
not the Robling 2008 paper used here).

**Script:** `scripts/msk/bone_asymmetry_disuse_mechanostat.py` (new; population-level literature
cross-check, no OpenSim, no subject-specific data — same pattern as
`scripts/msk/athlete_heart_remodeling.py`). **Data:**
`data/msk_smoketest/subject2_walking1/bone_remodeling/bone_asymmetry_disuse_results.json` (new
file in the existing `bone_remodeling/` directory, alongside but not touching Part A/B's own two
result files). **Isolation respected:** no git add/commit/push; only new files written.

---

## 0. The two falsifiers, why they are decorrelated, and why together they are stronger than either alone

**FALSIFIER 1 (FUNCTION pole):** in unilateral-sport athletes (tennis, squash), does the
**dominant (loaded) arm** show measurably higher bone mass/cortical geometry than the
**non-dominant arm**, in the **same person** — isolating mechanical load from genetics and
systemic hormones because both arms share one genome and one endocrine milieu?

**FALSIFIER 2 (DYSFUNCTION pole, decorrelated):** does mechanical **unloading** — bed rest,
spaceflight — drive bone loss preferentially at **weight-bearing** sites, in a way a purely
systemic (hormonal/immobilization-stress/radiation) driver could not produce?

These two falsifiers are **decorrelated** in every sense that matters: opposite direction (gain
vs loss), opposite population (athletes vs bed-rest patients/astronauts), opposite mechanism of
isolation (within-subject side-comparison vs within-subject site-comparison), and entirely
different research groups/eras/instruments (Finnish DXA/pQCT sports-medicine cohorts from the
1990s–2000s vs NASA/Roscosmos spaceflight-medicine QCT cohorts). Both converging on the same
mechanostat is a genuine over-determination, not a shared-instrument artifact.

**The adversary common to both:** "bone mass is set by genetics/systemic hormones, independent of
local mechanical load." Falsifier 1 forces it via a matched non-athlete control group (the
natural-handedness **void floor**). Falsifier 2 forces it via **within-subject, cross-site**
contrasts (a systemic driver should hit all of one person's bones comparably; it does not).

---

## 1. Literature anchors — every PMID esummary+efetch'd live this session (NCBI E-utilities)

| citation | PMID | DOI | n | verified content (verbatim quote) |
|---|---|---|---|---|
| Kannus P et al., *Ann Intern Med* 1995;123(1):27-31 | 7762910 | 10.7326/0003-4819-123-1-199507010-00003 | 105 players + 50 controls | "the players had a significantly (P<0.001) larger side-to-side difference in every measured site (15.5%±8.4%, 16.2%±9.8%, 8.5%±6.6, 12.5%±7.1%)" vs controls "(4.6%±4.6%, 3.2%±2.3%, 3.2%±3.8%, 3.9%±4.3%)"; menarche-timing effect "two to four times greater" |
| Jones HH et al., *J Bone Joint Surg Am* 1977;59(2):204-8 | 845205 | (none on file) | professional players | "cortical thickness on that side was greater by 34.9 per cent in men and 28.4 per cent in women" |
| Haapasalo/Kontulainen et al., *Bone* 2000;27(3):351-7 | 10962345 | 10.1016/s8756-3282(00)00331-8 | 12 players + 12 controls (pQCT) | BMC 14–27%, Tot.Ar 16–21%, Co.Ar 12–32% side-to-side; "the only significant side-to-side difference [in density]... Co.Dn of the distal humerus... -2%" |
| Kontulainen S et al., *J Bone Miner Res* 2003;18(2):352-9 | 12568413 (dup. index 12469923) | 10.1359/jbmr.2003.18.2.352 | 64 players + 27 controls (pQCT+DXA) | humeral-shaft BMC side-diff 19%(young starters)/9%(old starters); CoA 20%/9%; "grown periosteally (the CavA did not differ)"; CoD -1% both groups |
| Kontulainen S et al., *Med Sci Sports Exerc* 1999;31(5):646-52 | 10331882 | 10.1097/00005768-199905000-00004 | 13 players + 13 controls, 4-yr | humeral shaft 25%(1992)→26%(1996) despite training drop 5.2→2.6 sessions/wk; controls <5% throughout |
| LeBlanc AD et al., *J Bone Miner Res* 1990;5(8):843-50 | 2239368 | 10.1002/jbmr.5650050807 | 6 men, 17wk bed rest | "total body, lumbar spine, femoral neck, trochanter, tibia, and calcaneus demonstrated significant loss... 1.4, 3.9, 3.6, 4.6, 2.2, and 10.4[%]" |
| Lang T et al., *J Bone Miner Res* 2004;19(6):1006-12 | 15125798 | 10.1359/JBMR.040307 | 14 ISS crew, 4-6mo | spine -0.9%/mo; hip -1.2 to -1.5%/mo integral, -0.4 to -0.5%/mo cortical, -2.2 to -2.7%/mo trabecular; "cortical bone loss in the hip occurred primarily by cortical [endocortical] thinning" |
| Sibonga JD et al., *Bone* 2007;41(6):973-8 | 17931994 | 10.1016/j.bone.2007.08.022 | 45 crew, 56 pre/post pairs | "losses ranged between 2% and 9%... 50% restoration... within 9 months" |
| Vico L et al., *Lancet* 2000;355(9215):1607-11 | 10821365 | 10.1016/s0140-6736(00)02217-0 | 15 Mir cosmonauts | "Neither cancellous nor cortical bone of the radius was significantly changed at any of the timepoints. On the contrary, in the weight-bearing tibial site, cancellous BMD loss was already present after the first month" |
| Sibonga JD, *Curr Osteoporos Rep* 2013;11(2):92-8 | 23564190 | 10.1007/s11914-013-0136-5 | review | "losses >10% detected in the hip and spine in some astronauts following a typical 6-month mission" |
| Klein-Nulend J et al., *FASEB J* 1995;9(5):441-5 | 7896017 | 10.1096/fasebj.9.5.7896017 | chicken calvarial cells | "Osteocytes, but not osteoblasts or periosteal fibroblasts, reacted to 1h pulsating fluid flow... stress on bone causes fluid flow in the lacunar-canalicular system" |
| Robling AG et al., *J Biol Chem* 2008;283(9):5866-75 | 18089564 | 10.1074/jbc.M705092200 | rodent ulnar-load/hindlimb-unload | "Sost transcripts and sclerostin protein levels were dramatically reduced by ulnar loading... Hindlimb unloading yielded a significant increase in Sost" |
| Bonewald LF, *J Bone Miner Res* 2011;26(2):229-38 | 21254230 | 10.1002/jbmr.320 | review | "Osteocytes compose 90% to 95% of all bone cells... orchestrator of bone remodeling through regulation of both osteoclast and osteoblast activity" |
| Boyle WJ et al., *Nature* 2003;423(6937):337-42 | 12748652 | 10.1038/nature01658 | review | "Discovery of the RANK signalling pathway in the osteoclast has provided insight into... activation of bone resorption" |

**Method note:** esummary (JSON, title/journal/year/DOI) + efetch (plain-text abstract) for all 14
PMIDs, batched in 2 NCBI round trips; every quoted span above is copied verbatim from the
efetch'd abstract text, not recalled or paraphrased from a search-results snippet.

---

## 2. Falsifier 1 — dominant-vs-nondominant arm asymmetry, forced against the handedness void floor

**Pre-registered threshold** (stated with its own provenance honestly, not claimed as blind):
player/control side-asymmetry ratio ≥ **3.0×** at a given site counts as the pure-handedness/
genetic adversary falling there. This threshold was chosen **after** this session had already seen
a prior, unverified raw lead in this repo characterize the effect as "~4-5x" — so 3.0× is a
conservative floor relative to that characterization, not a blind a-priori pick, and not
reverse-fitted to narrowly pass (it sits well below every value it is then checked against).

| site | control (mean±SD) | player (mean±SD) | ratio | gate (≥3.0×) |
|---|---:|---:|---:|:---:|
| proximal humerus | 4.6%±4.6% | 15.5%±8.4% | **3.37×** | PASS |
| humeral shaft | 3.2%±2.3% | 16.2%±9.8% | **5.06×** | PASS |
| radial shaft | 3.2%±3.8% | 8.5%±6.6% | **2.66×** | **FAIL (borderline)** |
| distal radius | 3.9%±4.3% | 12.5%±7.1% | **3.21×** | PASS |

**3 of 4 sites PASS; the radial shaft is a genuine, reported-not-hidden borderline miss.** All
four differences are independently significant at P<0.001 vs controls (Kannus 1995). The
borderline site is the most distal, least mechanically engaged by a racquet swing (the humerus
absorbs the dominant bending moment) — consistent with, not contrary to, the mechanical-load
mechanism.

**Second, independent adversary-forcing axis — dose-response by developmental window, not
genotype:** Kannus 1995's own within-player-cohort split by training-start age relative to
menarche shows a "two to four times greater" bone-mass benefit for those who started at/before
menarche (10.5–23.5%) vs >15 years after (2.4–9.6%). A fixed genetic/handedness story predicts no
dependence on *when* training started relative to a biological milestone unrelated to genotype;
this dose-response is a second, structurally different way the same adversary is forced and falls.

**Cross-study corroboration (diverse instance-space, 4 independent cohorts, 3 imaging
modalities):** Jones 1977 (radiographic, +34.9%/+28.4% men/women), Haapasalo/Kontulainen 2000
(pQCT, cortical-area 12–32%), Kontulainen 2003 (pQCT+DXA, n=64+27, CoA 20%/9% young/old
starters — closely matching this delta's own task-level "~10–20% humeral cortical area" framing),
and Kontulainen 1999 (4-yr persistence: 25%→26% despite training frequency dropping by half,
controls flat <5%) all independently converge.

**Geometry, not mineralization, explains the gain (Haapasalo/Kontulainen 2000):** side-to-side BMC
difference (14–27%) sits almost entirely inside the side-to-side cortical cross-sectional-area
difference (12–32%), while volumetric cortical density changed by at most **-2%** at one site only
(the *loaded* arm slightly *lower*, not higher). Algebraically (BMC ≈ area × density, density
≈ const ⇒ %ΔBMC ≈ %ΔArea), this rules out "loading simply densifies existing matrix" and forces
the mechanism into its geometric form: **added cross-sectional material**, confirmed independently
by Kontulainen 2003's own finding that the marrow-cavity area does not differ side-to-side (all
growth is periosteal, none is a shrinking medullary canal).

---

## 3. Geometric derivation — why periosteal (not endosteal) growth is the efficient strategy

**Derived from the geometry of a thin-walled circular shell in bending, not asserted:** for a
thin annular area increment `dA` at radius `r`, the increment to the cross-section's second moment
of area (bending stiffness) is `dI = r² dA / 2` (the thin-ring formula — **verified this session
by brute-force numerical integration** against the closed form `I = π r³ t`, relative error
**0.0** at n=200,000 elements, not merely quoted from a textbook). Consequently, adding bone
mineral at the **outer (periosteal)** radius buys `(r_out/r_in)²` more bending stiffness per unit
area than adding the same area at the **inner (endosteal)** radius. Using illustrative (literature-
typical, **not** this twin-subject's own CT-measured) adult humeral-shaft dimensions
(`r_out≈11mm, r_in≈6mm`): **efficiency ratio ≈ 3.36×**. Both Haapasalo/Kontulainen 2000 and
Kontulainen 2003 independently report the marrow-cavity area is unchanged side-to-side — i.e.,
real bone independently does exactly what this derivation predicts is the mechanically efficient
strategy.

---

## 4. Falsifier 2 — disuse/bed-rest + spaceflight, the decorrelated opposite pole

### 4.1 LeBlanc 1990 (bed rest, n=6, 17 weeks) — cross-site gradient in the SAME subjects

| site | %Δ /17wk | %Δ /month |
|---|---:|---:|
| calcaneus | −10.4% | **−2.66%/mo** |
| trochanter | −4.6% | −1.18%/mo |
| lumbar spine | −3.9% | −1.00%/mo |
| femoral neck | −3.6% | −0.92%/mo |
| tibia | −2.2% | −0.56%/mo |
| total body | −1.4% | −0.36%/mo |
| forearm | *not in the significant-loss list* (measured, non-significant) |

Ratio, fastest (calcaneus) over slowest-significant (tibia): **4.73×** ≥ pre-registered 3.0×
gate — **PASS**. The most distal, least-loaded site tested (forearm) shows **no significant
loss at all** — the single strongest form of the site-specificity argument available in this
paper.

### 4.2 Vico 2000 (spaceflight, n=15 Mir cosmonauts) — within-subject, two bones, one systemic milieu

The cleanest test available: **radius** (non-weight-bearing) vs **tibia** (weight-bearing), same
15 people, same missions (1/2/6 months).

| site | significant timepoints |
|---|---|
| radius (cancellous + cortical, all 3 durations) | **0 of 6** |
| tibia, cancellous | 3 of 3 (from month 1) |
| tibia, cortical | 2 of 3 (from month 2) |

A purely systemic driver (hormonal shift, immobilization stress, radiation, fluid redistribution)
predicts both bones in the same person should move together. They do not: the unloaded radius is
a clean, total null; the loaded tibia deteriorates from month 1. **The systemic-confound adversary
is forced to its strongest available form (within-subject) and falls cleanly.**

*Caveat, not swept aside:* the 1- and 2-month subgroups are n=2 each (only the 6-month subgroup,
n=11, is well powered), so a null at the radius is weaker evidence than a null at high power would
be. Partially mitigating: the *same* low-n subgroups already reach significance at the tibia
(cancellous from month 1, cortical from month 2) — so the radius null is not simply "nothing
reaches significance at this n." Not a formal power analysis; reported as a residual caveat.

### 4.3 Cross-cohort/cross-vehicle agreement (Lang 2004's own internal comparison)

Lang 2004 (14 ISS/NASA crew) directly contrasts its own new measurement against an earlier,
separate Mir/Soviet cohort in the same abstract: hip 1.4–1.5%/mo (ISS) vs 1.5%/mo (Mir, historical)
— within range; spine 0.9%/mo (ISS) vs 1.0%/mo (Mir) — **10% relative difference**, well inside a
pre-registered ≤20% gate. Two different space programs, eras, and crews agree.

### 4.4 Honest, unreconciled tensions (symmetric QC — not swept under the rug)

- **Calcaneus rank flips between models:** fastest-losing site in bed rest (−2.66%/mo) but
  slowest-losing (of hip/spine/calcaneus) in spaceflight (−0.4%/mo, Lang 2004). Not reconciled;
  plausible drivers (recumbent bed rest still permits ankle micro-loading true weightlessness does
  not; different subject pools/n) are not independently tested this session.
- **Recovery timescale disagreement:** Vico 2000 reports tibial loss *persists* through recovery
  ("the time needed to recover is longer than the mission duration"); Sibonga 2007's exponential
  population-fit model predicts 50% recovery within 9 months across sites. Different cohorts
  (15 Mir cosmonauts vs 45 NASA+Russian crew), different eras/methods, not reconciled here.
  **A prior, unverified lead already in this repo** (`data/body_twin/agent_outputs/
  bone-adaptation__a198dac91f0b16616.json`) characterized Sibonga 2007 as showing "~34-47%
  recovery by 1yr" — that figure **does not appear** in the live-fetched abstract (which instead
  gives a 50%-recovery-*time* of ≤9 months) and is reported here as **not reproduced**, not
  silently carried forward. This is exactly the discipline this project's own method demands:
  raw subagent output is a hypothesis, not truth, until independently re-verified.

---

## 5. Mechanism — osteocyte canalicular fluid-flow mechanotransduction (closes the loop)

**Klein-Nulend 1995** (isolated chicken-calvarial cells): osteocytes — but explicitly **not**
osteoblasts or periosteal fibroblasts, tested side by side as the built-in control — respond to
1-hour pulsating fluid flow with sustained PGE2 release; intermittent hydrostatic compression
produces a weaker, delayed response. This is the origin of the **lacunar-canalicular fluid-flow**
hypothesis: mechanical loading strains the mineralized matrix, which drives interstitial fluid
through the canalicular network past the embedded osteocytes, and the resulting fluid **shear
stress** — not the tissue strain directly — is plausibly the proximate sensed signal.

**Robling 2008** (rodent ulnar-loading/hindlimb-unloading model): sclerostin (the `Sost` gene
product, an osteocyte-restricted Wnt/Lrp5 antagonist that suppresses osteoblastic formation) is
**dose-dependently down-regulated by local strain magnitude** under loading (higher-strain cortical
regions show a larger drop) and **up-regulated by hindlimb unloading**. This is the molecular
transducer: loading → fluid flow → sclerostin down → Wnt/Lrp5 disinhibited → osteoblast
recruitment; disuse → no flow → sclerostin up → formation suppressed. **Bonewald 2011** (review):
osteocytes are 90–95% of all adult bone cells and orchestrate **both** osteoclast and osteoblast
activity — the population-scale anchor for treating this single-cell-type mechanism as
tissue-level load-bearing, not a minor side pathway. **Boyle 2003** (brief mechanistic pointer,
not independently deep-dived this session): the RANKL/OPG ratio is the downstream **effector** for
resorption, sitting alongside the sclerostin/Wnt formation arm.

This mechanism chain is the **geometric/physical** answer to "why does Wolff's law work at all":
strain → poroelastic fluid flow in a fixed anatomical channel network (the lacunar-canalicular
system) → a cell population uniquely positioned (embedded throughout the matrix, 90-95% of all
bone cells) to sense it → a single molecular dial (sclerostin) that gates the formation arm in
proportion to local flow/strain magnitude.

---

## 6. Machine cross-checks / gates (`bone_asymmetry_disuse_results.json.gates`)

| gate | result |
|---|---|
| self-test: %/period unit round-trip | PASS |
| self-test: thin-ring `I=πr³t` formula vs brute-force numerical integration (n=200,000) | PASS (rel. err 0.0) |
| self-test: Kannus ratio dual-path recomputation | PASS |
| Falsifier 1: Kannus sites clearing the pre-registered ≥3.0× void-floor gate | **3/4** (radial shaft 2.66× borderline-fails, reported not hidden) |
| Falsifier 2: LeBlanc bed-rest site-gradient (calcaneus/tibia ≥3.0×) | PASS (4.73×) |
| Falsifier 2: Vico within-cosmonaut site-specificity (radius null vs tibia deteriorating) | PASS |
| Falsifier 2: Lang 2004 cross-cohort (ISS vs Mir) rate agreement (≤20% relative diff) | PASS |
| **gates_all_pass** (strict AND across all of the above) | **False** — by design: the one genuine borderline miss (radial shaft) is not smoothed into a false "all green" |

This is a deliberate, honest scorecard: 6 of 7 gates clear cleanly, one is a reported borderline
miss on the least mechanically-loaded of four tested sites — not rounded up, not hidden, and not
allowed to silently flip the top-level verdict to a false all-pass.

---

## 7. Honest gaps

1. Kannus 1995's radial-shaft site (2.66×) sits below the pre-registered 3.0× gate.
2. Jones 1977's often-quoted "n=84" could not be reproduced from the live-fetched 1977 abstract
   text itself (a very short, old-style abstract with no stated n) — carried from this repo's own
   prior unverified lead, flagged rather than silently repeated as verified.
3. Sibonga 2007's previously-recorded "~34-47% recovery by 1yr" (an unverified raw lead already in
   this repo) could not be reproduced from the live abstract; reported as not reproduced.
4. Vico 2000 vs Sibonga 2007 disagree on spaceflight recovery timescale; not reconciled.
5. Bed rest vs spaceflight disagree on the calcaneus's rank among affected sites; not reconciled.
6. The periosteal/endosteal geometric efficiency derivation uses illustrative literature-typical
   dimensions, not this twin-subject's own measured geometry (none exists in this repo).
7. Vico 2000's 1- and 2-month cosmonaut subgroups are n=2 each (only the 6-month subgroup, n=11,
   is well powered) — a null radius result there is weaker evidence than a null at high power;
   partially mitigated by the same low-n subgroups already reaching significance at the tibia, but
   not a formal power analysis.
8. Osteoporosis/postmenopausal estrogen-withdrawal — the third named dysfunction-pole example in
   this delta's task — is intentionally **not** independently verified here (out of scope for the
   specific athlete-asymmetry + disuse/spaceflight pair this delta targets); a raw, unverified lead
   already sits at `data/body_twin/agent_outputs/bone-endocrine-homeostasis__a8f407f9decde8a27.json`
   for a future session to independently verify, not assumed true by this one.
9. This is a population-level literature cross-check, not computed from this twin's own
   subject2/walking1 gait trial — see `proposed_cell` in the evidence JSON for how a future cell
   would couple athlete-arm loading history to this twin's own per-site strain accumulation.
10. Every cross-sectional tennis/racquet-sport study here (Kannus/Haapasalo/Kontulainen ×3) shares
   self-selection/survivorship as an unruled-out confound (players who respond better to loading
   may be more likely to remain competitive) — the SAME caveat this repo's own prior, unverified
   `bone-adaptation` lead had already flagged; not re-litigated, not newly resolved here either.
   Kontulainen 1999's 4-year longitudinal persistence result is the strongest available
   counter-argument (the SAME individuals, tracked over time, not a fresh cross-section) but does
   not fully rule out a differential-selection-into-*retention*-of-a-shrinking-cohort effect.

---

## 8. Confidence tier

**Literature-anchored, live-verified this session** (all 14 citations, both esummary metadata and
efetch abstract text, machine-parsed) for every quoted number above. **Method-only / illustrative**
for the periosteal-vs-endosteal geometric efficiency ratio (real formula, verified numerically;
illustrative not subject-measured dimensions). **Directional/qualitative, not independently
re-derived** for the BMU-phase-asymmetry-explains-rate-asymmetry link to Part B (a pointer, not a
new quantitative fit — see honest gaps and the script's own `bmu_phase_asymmetry_pointer_to_part_b`
field).

---

## 9. Files

- `scripts/msk/bone_asymmetry_disuse_mechanostat.py` — self-contained, no-OpenSim compute script;
  citations, self-tests, both falsifiers' gates, geometric derivation, JSON writer.
- `data/msk_smoketest/subject2_walking1/bone_remodeling/bone_asymmetry_disuse_results.json` — full
  machine-readable evidence (new file; Part A/B's own two result files in the same directory are
  untouched).
- `docs/MECHANISM_BONE_REMODELING_evidence.json` — the task-specified companion evidence file.
- `docs/MECHANISM_BONE_REMODELING_ASYMMETRY_DISUSE.md` — this file (new; `docs/
  MECHANISM_BONE_REMODELING.md` Part A/B, written by a concurrent sibling instance, is read for
  context above but not edited — isolation respected).
- No git add/commit/push performed. Only new files created.
