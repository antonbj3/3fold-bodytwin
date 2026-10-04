# MECHANISM MENISCUS LOAD DISTRIBUTION — hoop-stress mechanism + meniscectomy→OA dysfunction pole (2026-07-22)

Couples to the KNEE-CELL flagship, `docs/MECHANISM_KNEE_CARTILAGE_MATERIAL.md` (cartilage geometry/material
this mechanism loads), `docs/MECHANISM_CARTILAGE_CONTACT.md` / `docs/MECHANISM_JAM_CONTACT_DECORR.md`
(contact force/pressure OUTPUT), `docs/MECHANISM_KNEE_LIGAMENTS.md` (same knee, same citation-verification
discipline — see §7 for an equivalent citation-framing catch), and the function↔dysfunction organizing
axis (meniscectomy-driven OA as the dysfunction pole, alongside ACL/cartilage failure modes). No existing
meniscus cert was found in this repo (`grep -ril meniscus docs/` hits are optics-meniscus-lens and
tear-film-meniscus documents, unrelated; one MSK-ORTHOPEDIC entry in `MECHANISM_CELL_INVENTORY.md` notes a
dryad-meniscus-CADAVER-dataset was dropped for a *different* hip/shoulder-extension cell — this is a fresh
build, not a rebuild). Script: `scripts/msk/meniscus_hoop_stress_model.py`. Evidence:
`docs/MECHANISM_MENISCUS_LOAD_DISTRIBUTION_evidence.json`. Pure literature + first-principles statics; no
OpenSim model touched; no git operations.

---

## Headline

**The hoop-stress mechanism is confirmed and the "passive spacer" adversary is forced to fail on 3
independent, diverse tests — not by narration, by direct measurement.** A wedge-shaped meniscus under
axial load redirects a fraction of that load radially outward (elementary wedge statics,
`F_radial ∝ tan(θ)`); a closed/anchored circumferential fiber loop reacts that radial thrust as hoop
tension (`T = F_radial_per_length × r`, the same identity as pipe hoop stress). This is not a hypothesis
I am asserting from geometry alone — **Jones et al. 1996 (PMID 11415635) directly measured the mechanism
itself**: intact human medial menisci carry real circumferential strain (1.54–2.86% across
posterior/middle/anterior thirds), a **complete radial tear "completely defunctioned" that strain (→~0)**,
while an equal-mass **longitudinal tear only *redistributes* it** (still functional). The strongest fair
form of the adversary — "only tissue *mass* removed matters, not fiber *continuity*, so any in-situ cut
that excises nothing should act like the intact spacer" — is forced against three independent,
diverse-instance-space tests (different species, different menisci, different methods) and **falls on
all three** (§3). The downstream consequence — contact area/stress with vs. without an intact meniscus —
replicates across 3 independent measurements at **2.21×–4.00× more contact area and 2.00×–3.35× less peak
stress** (Fukubayashi & Kurosawa 1980, Baratz 1986, Rivarola et al. 2026), bracketing the task brief's own
"triples area / halves stress" framing exactly (§2). The dysfunction pole (meniscectomy → OA) is
confirmed with an effect size that **exceeds** the brief's own "~2–3×" guess in the classic single cohort
(Roos 1998: RR=14.0) and **matches it closely** in a modern pooled meta-analysis (OPTIKNEE 2022: OR=3.14)
(§6). One quantitative thread is left honestly open, not papered over: three studies that varied radial
tear WIDTH disagree on exactly where the effect *begins* (§5) — a real, disclosed, unresolved
cross-study discrepancy, not a forced-fit precision claim.

---

## 0. Pre-registration (stated before the adversary was forced)

- **Claim C**: hoop-tension conversion (wedge → radial thrust → circumferential fiber tension) governs
  meniscal load-spreading; contact-area increase / peak-stress decrease with an intact meniscus lands in
  a pre-registered **2×–4× band** (the task brief's own "triples/halves" framing, tested against the
  measured range, not assumed).
- **Adversary ¬C (forced to its strongest fair form)**: the meniscus is a passive, inert wedge/spacer —
  contact mechanics depend on meniscal tissue **mass/volume present**, not fiber **continuity**. Under
  this adversary, an in-situ tear that excises *zero* tissue should behave like the intact meniscus
  *regardless of orientation* (radial or longitudinal) — only actual excision should matter.
- **Falsifier**: (a) real cadaveric contact-area/peak-stress WITH vs. WITHOUT an intact meniscus;
  (b) tear-ORIENTATION dependence — radial (severs the hoop path) must raise stress/cause gapping far
  more than equal-mass longitudinal (parallel to the hoop fibers); (c) a zero-tissue-removed root/anchor
  tear must approach total-meniscectomy-level dysfunction if continuity (not mass) governs; (d)
  meniscectomy must show elevated long-term OA relative risk, replicated across ≥2 independent designs.
- **Symmetric-QC commitment**: every ratio below is computed by `scripts/msk/meniscus_hoop_stress_model.py`
  from raw numbers taken verbatim from PubMed abstracts — never hand-arithmetic in prose, never eyeballed
  from a figure. Every PMID was live-verified via NCBI eutils (title/journal/author cross-check) this
  session — the task brief itself contained one citation-framing error, caught and reported (§7), the
  same discipline `docs/MECHANISM_KNEE_LIGAMENTS.md` applied to two wrong given-PMIDs in an earlier build.

---

## 1. Geometric derivation — wedge statics → ring/hoop statics (first principles, dimensionally checked)

The meniscus cross-section is a wedge (thick peripherally at the capsular attachment, thin at the free
inner edge). Two elementary-statics steps convert axial compression into circumferential tension —
**derived here, not asserted**, and a real units bug was caught and fixed while deriving it (an earlier
draft mixed a lumped force with a line load; §"caveat" below documents the fix):

**Step A — wedge force-redirection (per unit circumferential length).** A total axial load
`F_axial_total` distributed around a ring of circumference `2πr` gives a local axial line-load
`q_axial = F_axial_total / (2πr)` [N/m]. A wedge element with its loaded face inclined at angle θ to the
horizontal, in equilibrium under a purely vertical line load, develops a reaction normal to that face of
magnitude `q_axial/cos(θ)`; its horizontal (radially outward) component is
`q_radial = q_axial · tan(θ)` [N/m] — the same reason a chock or an arch voussoir under vertical load
develops outward thrust. Pure Newtonian statics, no citation needed.

**Step B — ring/hoop equilibrium (Barlow's-law form).** A thin ring of radius `r` carrying a uniformly
distributed outward radial line-load `q_radial` develops circumferential tension `T` found by cutting the
ring across a diameter and balancing forces on one half:
`∫₀^π q_radial·r·sin(φ) dφ = 2T  ⟹  T = q_radial · r` [N]. This is the identical relation behind
thin-walled pressure-vessel hoop stress (`σ_hoop = P·r/t`) and is exactly the mechanism
**Shrive, O'Connor & Goodfellow 1978 (PMID 657636)** describe as *"a simple system ... to show how the
menisci can bear load."*

**Combining A+B**: `T = F_axial_total · tan(θ) / (2π)` — **the radius cancels**. For a uniformly
axisymmetric ring, total hoop tension depends only on total axial load and wedge angle, not ring size (a
bigger ring spreads the same load over more circumference, exactly compensating). This radius-independence
is a genuine, falsifiable, machine-checked property of the idealization — swept over `r = 10–30mm` at
fixed θ=20°, `T` stayed constant to **1.2×10⁻¹⁶ relative spread** (floating-point-exact; see
`geometric_model.radius_independence_check` in the evidence JSON). At a representative total load of
1000N, illustrative tensions are **T ≈ 43N (θ=15°) to 92N (θ=30°)** — reported as order-of-magnitude only:
no live-verified circumferential tensile modulus or failure-load number was found this session to convert
this into an independent stress cross-check (Tissakht & Ahmed 1995, PMID 7738050, confirmed real/on-topic
but its abstract does not itself state the modulus value — an honest, disclosed gap, not a recalled
number dressed as verified).

**Why this predicts the orientation-dependence (§3) geometrically, not just empirically**: a ring is a
closed (or firmly anchored open-arc) tension-carrying loop — like a chain, cutting it *anywhere* fully
across breaks the load path entirely, regardless of *where* around the circumference the cut is (this is
why a root/anchor tear, §4, is mechanically equivalent to a mid-substance radial tear — both eliminate
the closed-loop condition). A **radial** cut runs across the ring — a full-thickness radial tear severs
every circumferential fiber crossing that meridian, exactly like cutting a hose-clamp. A **longitudinal**
cut runs parallel to the hoop fibers (splits the wedge into inner/outer strips along its own length) —
each strip can still independently carry hoop tension around its own smaller loop, just redistributed,
not abolished. This is the first-principles reason radial and longitudinal tears are NOT interchangeable
even though both remove zero tissue mass — and it is exactly what was measured (§3).

---

## 2. Falsifier A — contact area/stress WITH vs. WITHOUT an intact meniscus (measured, machine cross-checked)

| Study (PMID) | species/n | condition | contact area | peak stress |
|---|---|---|---|---|
| Fukubayashi & Kurosawa 1980 (6894212) | human cadaver, n=7, 1000N, 0° | **intact** | 1150 mm² | 3 MPa |
| | | **total meniscectomy** | 520 mm² | 6 MPa |
| Baratz et al. 1986 (3755296) | human cadaver, n=4, ~1780N | **partial meniscectomy** | −10% | +65% |
| | | **total meniscectomy** | −75% | +235% |
| Rivarola et al. 2026 FE (41539441) | MRI-derived FE, validated r=0.91 vs cadaveric/computational benchmarks | **intact** | 110±8 mm² | 1.2±0.2 MPa |
| | | **complete (100%) radial tear** | 35±6 mm² | 3.3±0.4 MPa |
| | | **anatomic repair** | 95±7 mm² (86% restored) | 1.4±0.3 MPa (ns vs. intact) |

**Machine-computed ratios** (`compute_ratios()` in the script — every number below is code output, not
hand arithmetic):

- Fukubayashi & Kurosawa: area ratio (intact/meniscectomy) = **2.21×**; peak-stress ratio
  (meniscectomy/intact) = **2.00×** (exactly halved, matching the task brief's "roughly halves" almost
  perfectly).
- Baratz total meniscectomy: area ratio = **4.00×**; stress ratio = **3.35×**. Baratz partial
  meniscectomy: area ratio = **1.11×**; stress ratio = **1.65×** — a genuine **dose-response** (partial
  < total on both axes, confirmed monotonic by code: `baratz_dose_response_monotonic = True`).
- Rivarola FE: area ratio (intact/100%-tear) = **3.14×**; stress ratio (100%-tear/intact) = **2.75×**.
- **Cross-study range** (over-determination, not a cherry-picked single number): area ratio spans
  **2.21×–4.00×**, stress ratio spans **2.00×–3.35×**. The task brief's "triples area (≈3×) / halves
  stress (2× ratio)" framing falls **inside both measured ranges** (`task_brief_triples_area_halves_stress_within_measured_range`
  = `{True, True}`).
- **Internal consistency check** (independent arithmetic on Fukubayashi & Kurosawa's own raw numbers, not
  quoting the paper's own ratio back at itself): mean pressure = load/area gives 0.87 MPa (intact) vs.
  1.92 MPa (meniscectomy) — both far below the reported *peak* pressures (3/6 MPa), confirming the
  pressure distribution is sharply non-uniform in both conditions. The **peak/mean concentration factor**
  is nearly the same with or without the meniscus (3.45× vs. 3.12×, <15% relative difference — code flag
  `fk_stress_concentration_factor_roughly_conserved = True`). This means the mechanism is primarily about
  the meniscus **enlarging the load-bearing footprint**, not about locally reshaping the stress-concentration
  *pattern* — a genuinely distinct, falsifiable sub-claim that the raw numbers happen to support.

---

## 3. Falsifier B — the adversary forced to its strongest form, and where it falls

**Strongest fair adversary**: *"the meniscus is a passive spacer; only tissue mass removed matters, not
fiber continuity — an in-situ tear that excises nothing should behave like the intact meniscus regardless
of orientation."* Three independent, diverse-instance-space tests (different species, different menisci,
different measurement modalities) were forced against this — `adversary_test()` in the script, verdicts
computed, not narrated:

| Test | Adversary predicts | Measured | Verdict |
|---|---|---|---|
| **Orientation, equal retained mass** — Sukopp et al. 2024 (PMID 37986646), porcine n=12, RSA gapping | No difference (same mass either way) | Longitudinal: **no** gapping/pressure change vs. native. Radial: **significant** gapping + pressure rise, normalized by suture repair | **ADVERSARY FALLS** |
| **Zero-mass-removed root/anchor tear** — Allaire et al. 2008 (PMID 18762653), human n=9 | Behaves like intact (full spacer mass retained) | Root tear: **+25% peak pressure (p<0.001)**, statistically **indistinguishable** from total meniscectomy | **ADVERSARY FALLS** |
| **Direct mechanism measurement** — Jones et al. 1996 (PMID 11415635), human n=19, strain gauges | Hoop strain is irrelevant to spacer function | Complete radial tear: **"completely defunctioned"** the meniscus (hoop strain → ~0) | **ADVERSARY FALLS (at the mechanism itself, not a downstream proxy)** |

`adversary_test()['summary']`: **3/3 tests, adversary falls on all 3**, across a genuinely diverse
instance space (human ×2 tests + porcine ×1; medial ×2 + lateral cross-check in §4; RSA-gapping,
pressure-film, and strain-gauge as three independent measurement modalities; two independent labs —
Pittsburgh/Harner group for Allaire, Ulm University for Sukopp).

Cross-corroboration on the **lateral** side (a different meniscus, same logic): **LaPrade/Jansson et al.
2014 (PMID 24647503)**, human cadaver n=8 — lateral meniscus posterior-root avulsion or a radial tear
adjacent to the root (again, zero tissue excised) significantly decreased contact area and increased
mean/peak contact pressure at nearly every flexion angle tested, restored by in-situ repair. Same
qualitative result, opposite compartment — the mechanism is not a medial-only artifact.

---

## 4. Why root/anchor tears ≈ total meniscectomy (the geometric point, machine-confirmed)

§1 derived that a ring's hoop-tension load path breaks if cut *anywhere* fully across — root avulsion (a
tear at the peripheral anchor) and a mid-substance complete radial tear both eliminate the same
closed-loop condition. This is not just an analogy: **Allaire 2008 measured no statistically significant
difference** between root-tear and total-meniscectomy peak contact pressure, and **LaPrade/Jansson 2014**
found the equivalent on the lateral side. The geometric model predicts this convergence; the data confirm
it independently on both compartments.

---

## 5. Radial-tear-width dose-response — confirmed direction, DISPUTED threshold location (honest, not papered over)

Every dataset that varied tear severity shows a **monotonic, non-decreasing** relationship between how
much of the hoop path is cut and how much contact mechanics degrade (`fe_dose_response_monotonic = True`,
`baratz_dose_response_monotonic = True`, `freutel_stress_area_tradeoff_check = True`). But the exact
**width at which the effect becomes detectable** is NOT consistent across the three studies that measured
it as a swept parameter:

| Study | method | 50% width | 60% width | 90% width | 100% width |
|---|---|---|---|---|---|
| Bedi 2010 (20516315) | human cadaver, dynamic gait simulator | — | **no significant change** | some rise (+1.3 MPa posterocentral) | (not tested; partial meniscectomy tested instead, +1.4 MPa further) |
| Rivarola 2026 FE (41539441) | FE, static 1000N axial | **+75% peak stress** (2.1 vs. 1.2 MPa) | — | — | +175% (3.3 vs. 1.2 MPa) |
| Jones 1996 (11415635) | human cadaver, strain gauge | strain **already reduced** anteriorly | — | — | strain → ~0 ("completely defunctioned") |

A **fourth** dose-response confirmation, using a different measured quantity entirely (attachment force,
not contact pressure): **Freutel et al. 2015 (PMID 24671386)**, porcine n=6, graded partial meniscectomy
(intact → 50% → 75% → 75%+anterior-horn PM) — posterior **attachment force** fell by up to **55%**
(p=0.003) while contact area fell by up to **23%** (p=0.01) and max contact pressure rose by up to **40%**
(p=0.02), all monotonic with resection extent (`freutel_stress_area_tradeoff_check = True`). Critically,
circumferential strain in the *remaining* peripheral rim was **not** significantly affected by partial
resection — consistent with §1's picture that the outer peripheral rim is the critical hoop-bearing belt:
as long as some of it survives intact, it keeps carrying normal strain; what degrades is the *amount* of
tissue engaging that mechanism. The authors' own conclusion: partial meniscectomy "reflects the impaired
ability of the meniscus to transform axial joint load into meniscal hoop stress" — the identical mechanism
this document derives geometrically in §1. **Disclosed, not concealed**: Freutel 2015 shares its
institute AND two named authors (Seitz AM, Ignatius A) with Sukopp 2024 (§3's orientation test) — both
from the Institute of Orthopaedic Research and Biomechanics, University of Ulm, confirmed live via PubMed
author-affiliation fields. This is **one Ulm research lineage across two papers (2015+2024), not two
independent porcine labs** — counted as such throughout this document (see §8 gap 6), not silently
inflated into extra "independent" replication.

A toy geometric model (`fiber_continuity_fraction()` + `predicted_stress_ratio()` in the script) was
**calibrated on Bedi's own pattern** (flat through 60%, some rise by 90% → breakpoint fit at width-fraction
0.75) and then **tested held-out** against Rivarola's independent 50%-width point (never used in
calibration): the calibrated model predicts **ratio=1.0 (no change)** at 50% width; Rivarola **measured
1.75**. **MISMATCH** — the model, calibrated on one real dataset, fails a genuine held-out prediction
against a second independent dataset. This is reported as an open disagreement, not resolved: 2 of 3
dose-response sources (Rivarola, Jones) favor an earlier onset (~50%); 1 of 3 (Bedi) favors a later one
(~90%). Plausible reconciling factors — dynamic multidirectional gait loading (Bedi, with real tissue
viscoelasticity/friction and capsular/ligamentous secondary restraints) vs. static near-frictionless axial
FE loading (Rivarola); differing definitions of "% width"; ordinary small-N cadaveric variability — are
disclosed as open, not adjudicated here. What IS robust regardless of the threshold-location dispute: the
strong-form passive-spacer adversary (§3) already falls from the *zero-tissue-removed* tests alone,
independent of exactly where the width-threshold sits.

---

## 6. Dysfunction pole — meniscectomy → osteoarthritis

- **Fairbank 1948 (PMID 18894618)** — the originating, eponymous qualitative description (ridge
  formation, femoral-condyle flattening/squaring, joint-space narrowing after meniscectomy). This is a
  pre-abstract-era PubMed record: **no abstract text is indexed**, so no numeric figure is live-extractable
  from it — reported honestly as qualitative-only, with quantitative confirmation deferred to the two
  entries below.
- **Roos et al. 1998 (PMID 9550478)** — 21-year follow-up, isolated meniscectomy (n=107) vs. matched
  controls (n=68): any radiographic change 71% vs. 18% (simple prevalence ratio, code-computed:
  **3.94×**); Kellgren-Lawrence ≥2 ("definite" OA) 48% vs. 7% (simple ratio: **6.86×**); the paper's own
  matched-pairs statistic: **RR=14.0 (95% CI 3.5–121.2)**. This **exceeds** the task brief's own "~2–3×"
  framing substantially — reported as a correction, not silently reconciled to fit the brief.
- **OPTIKNEE 2022 meta-analysis (PMID 36455966)** — modern, multi-cohort, random-effects pooled estimate:
  total medial meniscectomy (in an ACL-reconstruction context) **OR=3.14 (95% CI 2.20–4.48)**; partial
  meniscectomy **OR=1.87 (95% CI 1.45–2.42)**. This pooled, tighter-CI estimate **matches the task brief's
  "~2–3×" framing closely** — triangulating against Roos's higher, wide-CI single-cohort figure (small
  sample, isolated-injury context vs. OPTIKNEE's ACL-associated, multi-cohort context — different
  populations, same qualitative and roughly-overlapping-magnitude conclusion).
- **Partial > total preservation**, confirmed in **two independent outcome types** (over-determination):
  acute contact mechanics (Baratz 1986: partial −10%area/+65%stress vs. total −75%area/+235%stress) AND
  10+-year OA risk (OPTIKNEE: partial OR=1.87 vs. total-medial OR=3.14) — same ordering, two decades and
  two entirely different outcome measures apart.

---

## 7. The medial-~50%/lateral-~70% load-share figure — verification status (honest)

The task brief's opening claim ("menisci transmit ~50% medial to ~70% lateral of tibiofemoral load") is
**directionally and approximately confirmed**, but not from one clean single-number primary source:

- **Jones et al. 1996's own abstract background statement** (independent of that paper's new strain-gauge
  result): *"the menisci transmit approximately 50% of the load through the knee."*
- **Shrive, O'Connor & Goodfellow 1978 (PMID 657636)**: menisci bear **≥45%** of total joint load in
  partially-degenerate human knees and **~75%** in healthy pig knees — brackets the "50–70%" range from a
  different angle (species/degeneration axis, not medial/lateral compartment axis).
- **Ahmed & Burke 1983 (PMID 6688842)** is the paper most often cited for the specific medial/lateral
  split — confirmed live as real and on-topic ("a significant fraction of the joint compressive load is
  transmitted through the menisci... total meniscectomy causes a drastic alteration") but **its own
  abstract does not state the specific 50%/70% numeric split** (would require full-text, not fetched this
  session — an honest, disclosed gap, not a fabricated number).

Net: the qualitative claim and rough magnitude are corroborated by 2 independent live-verified sources;
the precise medial-vs-lateral percentage split is a widely-repeated figure in the review literature that
this session could not trace to a specific primary-abstract number.

**Separately, a citation-framing catch** (same discipline as `docs/MECHANISM_KNEE_LIGAMENTS.md`'s two
wrong given-PMIDs): the task brief cites "Fukubayashi 1980 / Kurosawa 1980" as though these were two
separate papers. Live verification (PMID 6894212) confirms this is **one single paper co-authored by
both** (Fukubayashi T, Kurosawa H, *Acta Orthop Scand* 1980) — reported prominently, not silently
corrected.

---

## 8. Symmetric QC — what this does NOT prove (honest gaps, full list)

1. **Threshold-location disagreement (§5)** — unresolved across 3 independent dose-response studies;
   the toy calibrated geometric model fails its own held-out test against Rivarola's 50%-width point.
2. **No live-verified circumferential tensile modulus** — Tissakht & Ahmed 1995 (PMID 7738050) confirmed
   real/on-topic but its abstract doesn't state the modulus; the geometric model's illustrative tension
   values (43–92N) are order-of-magnitude only, not stress- or failure-load-anchored.
3. **Ahmed & Burke's specific 50%/70% split not abstract-confirmed** (§7) — corroborated in direction and
   magnitude by two OTHER independent sources, not traced to this specific paper's primary numbers.
4. **Fairbank 1948 is qualitative-only** here — no abstract text indexed in PubMed (pre-abstract era);
   quantitative confirmation of the same dysfunction pole rests on Roos 1998 + OPTIKNEE 2022 instead.
5. **The task brief's own "~2–3×" OA-risk framing undersold the classic single-cohort figure** (Roos:
   RR=14.0) while matching the modern pooled meta-analysis (OPTIKNEE OR=3.14) — both reported, not merged.
6. **Small-N primary studies** (n=3–19 cadaveric specimens each, standard for this literature) — mitigated
   by replication across 10 studies spanning human cadaver (×7 distinct groups: Fukubayashi&Kurosawa,
   Baratz, Ahmed&Burke, Allaire, Bedi, LaPrade/Jansson, Jones), porcine (×2 papers, but **one Ulm research
   lineage, NOT independent** — Sukopp 2024 and Freutel 2015 share 2 named authors, confirmed live), and
   FE (×1, Rivarola). This is weaker over-determination than an initial count of "porcine ×2 independent
   labs" would suggest — corrected here after live-checking author affiliations (§3, §5); the CORE
   adversary-falls verdict (§3) does not depend on Freutel at all, so this correction does not weaken it.
7. **Axisymmetric-ring idealization** — real menisci are open C-shaped arcs under measurably non-uniform
   load (Jones 1996: 2.86%/2.65%/1.54% ant/mid/post strain); the radius-independence result is a property
   of the idealized model, not independently checked against a non-axisymmetric analytic/FE solution here.
8. **Abstract-level extraction only** — all numbers sourced from PubMed abstracts via NCBI efetch, not
   independently re-extracted from each paper's full-text Methods/Results tables (adequate for every claim
   used here except the two explicitly disclosed gaps above).
9. **A units bug was found and fixed in this session's own first-draft derivation** (a lumped force fed
   directly into a line-load ring formula) — disclosed as a demonstration that the dimensional-consistency
   check in §1 was actually exercised, not merely asserted.

---

## 9. Couples to

- `docs/MECHANISM_KNEE_CARTILAGE_MATERIAL.md` — the cartilage geometry/material this hoop-spread load
  ultimately loads (femoral/tibial thickness maps, elastic-foundation modulus).
- `docs/MECHANISM_CARTILAGE_CONTACT.md` / `docs/MECHANISM_JAM_CONTACT_DECORR.md` — the downstream
  contact-force/pressure OUTPUT this mechanism modulates.
- `docs/MECHANISM_KNEE_LIGAMENTS.md` — same knee model family; same live-PMID-verification discipline
  (that doc caught two wrong given-PMIDs; this one catches a two-papers-treated-as-one framing error and
  an undersold OA relative-risk figure).
- `docs/MECHANISM_CUTTING_ACL_MECHANISM.md` — the athlete knee-load axis (landing/cutting/squat loads)
  this hoop mechanism must survive without radial-tear-scale hoop-continuity loss.
- The function↔dysfunction organizing axis — meniscectomy-driven OA (§6) as a dysfunction pole
  structurally analogous to other MECHANISM cells' injury→degeneration mechanisms.
- MSK-ORTHOPEDIC cell entry in `docs/MECHANISM_CELL_INVENTORY.md` — that cell dropped a dryad-meniscus-
  cadaver dataset as off-regime for a *hip/shoulder* extension; this document is knee-native and does not
  depend on that dataset.

---

## 10. Files

- `scripts/msk/meniscus_hoop_stress_model.py` — raw literature data table (14 PMIDs), machine-computed
  ratios (`compute_ratios()`), dimensionally-checked wedge→ring geometric derivation
  (`wedge_hoop_model()`, with a live radius-independence sweep), the calibrated/held-out threshold model
  (`fiber_continuity_fraction()`, `predicted_stress_ratio()`), and the 3-test forced-adversary check
  (`adversary_test()`). Run: `.venv-msk/bin/python3 scripts/msk/meniscus_hoop_stress_model.py`.
- `docs/MECHANISM_MENISCUS_LOAD_DISTRIBUTION_evidence.json` — full machine-readable evidence: pre-
  registration, all 15 PMID live-verification records (title/journal/author cross-checks), the complete
  computed-evidence block (raw data + ratios + geometric model + adversary test verdicts), headline
  verdicts, and this document's honest-gaps list in machine-readable form.

No git commit, no git push performed (isolation respected). No existing file modified.
