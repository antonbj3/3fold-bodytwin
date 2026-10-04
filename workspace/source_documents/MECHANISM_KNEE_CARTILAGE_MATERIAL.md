# MECHANISM KNEE CARTILAGE — GEOMETRY (thickness map) + MATERIAL (elastic foundation) characterization (2026-07-22)

Extends `docs/MECHANISM_KNEE_LIGAMENTS.md` (ligaments) and `docs/MECHANISM_JAM_CONTACT_DECORR.md` /
`docs/MECHANISM_CARTILAGE_CONTACT.md` (contact force/pressure OUTPUT) one layer further down: the
cartilage GEOMETRY and MATERIAL substrate those legs' contact model depends on. Extracts the
`lenhart2015.osim` model's actual `Smith2018ContactMesh` cartilage thickness map (tibial/femoral,
+patella as a bonus 3rd limb) and elastic-foundation modulus, and validates EACH against its OWN
external anchor, live-verified — not against each other, and not against a tautology gate.

Script: `scripts/msk/cartilage_thickness_material.py`. Evidence:
`data/msk_smoketest/cartilage_thickness_material/cartilage_thickness_material_results.json`.
Inputs read in place (unchanged, mtime-reverified): `scratchpad/knee_lig/lenhart2015.osim` + its 6
STL meshes, and the already-completed `data/msk_smoketest/jam_contact_decorr/results/joint-mechanics/
gait_driven.h5`. No git operations; no existing file edited.

---

## Headline

**GEOMETRY and MATERIAL get OPPOSITE verdicts, and the model's own citation trail explains why.**

| substrate | model value | anchor (pre-registered band) | verdict |
|---|---|---|---|
| Femoral cartilage thickness (weight-bearing, gait peak) | 2.42mm median / 2.46mm mean | 2.0–3.0mm (in-vivo MRI) | **PASS** (in band, all 4 statistics) |
| Tibial cartilage thickness (weight-bearing, gait peak) | 1.89mm median / 2.10mm mean | 2.0–4.0mm (in-vivo MRI) | **BORDERLINE** (mean in band, median 5.6% below floor) |
| Cartilage elastic modulus (all 3 surfaces, uniform) | 5.0 MPa | 0.5–1.5 MPa (cadaveric indentation) | **FAIL** (3.3×–10× above band) |

The GEOMETRY is this donor's own dynamic-MRI-derived bone/cartilage surface reconstruction
(subject-specific, real) and lands inside or at the edge of the population MRI-morphometry band.
The MATERIAL is a single scalar (5 MPa) that Lenhart2015's own Methods cite to **two computational
finite-element papers (Blankevoort & Huiskes 1991; Caruntu & Hefzy 2004), not an empirical cartilage
mechanical-testing study** — a chain of modeling-convenience reuse, not a tissue measurement — and it
sits 3.3×–10× above the direct-indentation/aggregate-modulus literature. This is exactly the
geometry-vs-material asymmetry the task asked to keep separate: **the geometry is well-grounded, the
material is a borrowed, twice-removed, order-of-magnitude-adjacent stiff assumption.**

---

## 1. What was extracted from the model (raw XML values, live-grepped)

`scratchpad/knee_lig/lenhart2015.osim`, `ContactGeometrySet` (lines ~13805–13936), 3
`Smith2018ContactMesh` objects (femur_cartilage, tibia_cartilage, patella_cartilage), **identical**
material properties on all three:

| property | value | note |
|---|---|---|
| `elastic_modulus` | 5,000,000 Pa = **5.0 MPa** | uniform per triangle, all 3 meshes |
| `poissons_ratio` | 0.45 | uniform |
| `thickness` (flat default) | 0.003 m = 3.0 mm | **NOT actually used** (see below) |
| `use_variable_thickness` | **true** | overrides the flat default |
| `mesh_back_file` | the matching `*-bone.stl` | ray-cast target for thickness |
| `min_thickness` / `max_thickness` | 0.001 m / 0.01 m = **1–10 mm** | clip bounds on the ray-cast |
| contact law | Bei & Fregly 2003 lumped-parameter elastic foundation | combines both meshes' thickness+modulus in a pair |

**Because `use_variable_thickness=true`, the model does NOT use a flat 3mm cartilage layer** — it
computes a per-triangle thickness live, "by calculating the distance along a normal ray cast from
the center of each triangle in `mesh_file` to intersection with `mesh_back_file`" (the property's
own doc-comment, quoted verbatim), clipped to [1,10]mm. That per-triangle field is the actual
GEOMETRY substrate and is not stored anywhere in the `.osim` or in the already-completed JAM h5
output (checked live: `gait_driven.h5`'s `Smith2018ArticularContactForce` groups store
`triangle_proximity`/`triangle_pressure` per contacting-mesh-PAIR during simulation — the
instantaneous gap between the two OPPOSING cartilage surfaces — not the static cartilage-to-bone
layer thickness). So it had to be re-derived directly from the raw STL geometry — §2.

---

## 2. Method — reproducing the model's own documented ray-cast, not a proxy

`scripts/msk/cartilage_thickness_material.py`:

1. **STL parsing, hand-rolled** (no `trimesh`/`numpy-stl`/network available in `.venv-msk` or system
   Python — confirmed live). 5 of 6 meshes are ASCII (MeshLab); **`lenhart2015-R-tibia-bone.stl` is
   BINARY** (MeshMixer format — caught because a naive `grep -c "facet normal"` silently returned 0;
   `file` confirmed `data`/binary). Binary parser cross-checked arithmetically: file size 1,844,584
   bytes → `(1844584-84)/50 = 36890.0` triangles exactly (no remainder) — matches the parser's own
   count, an independent, from-first-principles confirmation the binary reader is correct.
2. **Vectorized Möller–Trumbore ray-triangle intersection**, one ray per cartilage-triangle center,
   candidates pruned via a `scipy.spatial.cKDTree` over bone-triangle centroids (K=40, retried at
   K=150 for any triangle still missing a hit). Direction (does the cartilage mesh's own face normal
   point toward or away from the bone) was **decided empirically per mesh pair, not assumed**: the
   correct direction won by a factor of ~17× in hit-rate for all 3 independent mesh pairs (e.g.
   femoral: 94.4% hit rate one way vs 5.6% the other; tibial 94.4% vs 5.6%; patellar 98.6% vs 1.4%) —
   a clean, non-ambiguous separation, reported alongside as the forced-adversary comparison (the
   "losing" direction's stats are in the evidence JSON, not discarded).
3. **Cross-check machinery** (all reported in the JSON, not narrated):
   - Facet counts independently parsed in Python **match the already-solved JAM h5 output's own
     per-triangle array width** exactly for all 3 cartilage meshes (29782/8602/6052) — confirms the
     STL files read here are the same meshes the actual prior JAM run used.
   - Zero degenerate (zero-area) triangles in any of the 3 cartilage meshes.
   - An independent, algorithm-free lower-bound cross-check (nearest-centroid Euclidean distance,
     free from the same KDTree query) agrees closely with the ray-cast thickness for all 3 pairs
     (e.g. femoral nearest-centroid median 2.00mm vs ray-cast median 2.07mm; tibial 1.29mm vs
     1.40mm) — two independent geometric measures, same answer.
   - Residual "no local hit" rate (~5.6% femoral/tibial, 1.4% patellar) is excluded from statistics
     (not zero-filled), concentrated at the expected location (peripheral/rim triangles where the
     ray misses within the 5cm cutoff) — disclosed as a bounded limitation, not silently absorbed.
4. **"Loaded/weight-bearing subset"** — defined from the model's OWN solved gait contact mechanics
   (triangles with pressure > 1kPa at the peak-force frame), not an arbitrary hand-picked anatomical
   region. **OODA catch, forced not skipped**: a naive `argmax` of contact-force magnitude over the
   full time series picked up a documented ~0-0.05s startup-transient artifact (same gotcha already
   flagged in `scripts/msk/cartilage_contact.py`'s `SETTLE_TIME_S`) — caught because femur_cartilage's
   naive peak (t=0.04s) disagreed with tibia_cartilage's naive peak (t=0.54s) for the SAME `tf_contact`
   force pair, which by Newton's third law must peak at the same instant on both sides. Fixed by
   excluding t<0.10s before taking `argmax`; **after the fix, all of tf_contact's peak times agree at
   t=0.54s**, exactly matching the already pre-registered/verified peak time in
   `docs/MECHANISM_JAM_CONTACT_DECORR.md` — an external cross-check PASS.

---

## 3. GEOMETRY result vs external anchor (in-vivo MRI morphometry)

Pre-registered band (supplied by the task before any number here was computed): **femoral 2.0–3.0mm,
tibial 2.0–4.0mm**.

| | femoral | tibial | patellar (bonus) |
|---|---:|---:|---:|
| whole-mesh median (clipped mm) | 2.07 | **1.40** | 2.17 |
| whole-mesh mean (clipped mm) | 2.02 | **1.67** | 2.16 |
| loaded-subset median (gait peak, clipped mm) | 2.42 | 1.89 | 1.94 |
| loaded-subset mean (gait peak, clipped mm) | 2.46 | 2.10 | 1.93 |
| fraction of whole mesh at the model's own 1mm clip floor | 18.3% | **37.1%** | 16.3% |
| in pre-registered band? | **YES, all 4 stats** | mean-in / median-just-below | (no pre-registered band; bonus) |

**Femoral: clean PASS.** All four statistics (whole-mesh and loaded-subset, median and mean) land
inside 2.0–3.0mm.

**Tibial: borderline, not a clean pass or fail.** The loaded-subset MEAN (2.10mm) is inside the
2.0–4.0mm band; the loaded-subset MEDIAN (1.89mm, a more outlier-robust statistic) sits 5.6% below
the lower edge; the whole-mesh statistics are both below band. 37.1% of the tibial cartilage surface
sits at the model's own 1mm clip floor (vs 18.3% femoral) — this donor's tibial cartilage genuinely
reads thin, not just at the periphery.

**Forced adversary for the femoral PASS** (a leaning-positive result — is it spurious?): (a) clipping
to [1,10]mm cannot itself manufacture a narrow 2–2.5mm concentration — it only bounds extremes, so
landing mid-range is a real signal; (b) an independent, algorithm-free distance metric (nearest-
centroid) agrees closely, ruling out an artifact specific to the ray-cast implementation; (c) a
wholesale model unit-scale bug (e.g. cm-vs-m) is independently ruled out — `docs/
MECHANISM_KNEE_LIGAMENTS.md`'s own segment-length check (femur 0.374m, tibia 0.403m, anatomically sane
adult-in-meters) and this model's own BW=608N (~62kg) both cross-confirm meters-scale correctness.
**The adversary falls — the femoral pass is not spurious.**

**Forced adversary for the tibial borderline result** (a leaning-negative-ish result — is the model
really thin, or is the comparison unfair?): the same hit-rate/cross-check machinery that validated
femoral applies identically to tibial (94.4% hit rate, close nearest-centroid agreement, no
degenerate triangles) — **no evidence of a tibia-specific algorithm bug**. The more likely
explanation is genuine single-donor anatomical variation (the task's own reminder: "lenhart2015 is a
single generic donor," n=1, being compared to a population band) compounded by this analysis's
"loaded subset" being a coarse >1kPa-anywhere-in-the-gait-contact-patch definition that does not
separate medial from lateral tibial plateau (the JAM h5 does carry a finer 6-region breakdown,
`regional_contact_area`/`regional_mean_pressure`, that could refine this further — **not done here**,
flagged as a cheap follow-up, not chased for scope/lean reasons).

---

## 4. MATERIAL result vs external anchor (cadaveric indentation)

Pre-registered band: **0.5–1.5 MPa** (aggregate/equilibrium modulus).

Model value: **5.0 MPa**, uniform across femoral/tibial/patellar cartilage.

`ratio_to_upper_bound = 3.33×`, `ratio_to_lower_bound = 10.0×`. A tertiary cross-check (Wikipedia's
"Cartilage" article, live-fetched, citing Mansour 2013 / Patel et al. 2019 / Korhonen et al. 2002)
independently states aggregate modulus "typically 0.5 to 0.9 MPa" — landing inside the pre-registered
band and putting the model's value at **5.56× above that upper figure**. **FAIL — not borderline.**

**Critical finding from the model's OWN citation lineage** (checked live, PMC full text of Lenhart
2015 itself, PMID 25917122): the paper's Methods state cartilage was "assumed to have an elastic
modulus of 5 MPa" citing references **10 and 12** — which are:

- Ref 10: Blankevoort L, Huiskes R. *Ligament-bone interaction in a three-dimensional model of the
  knee.* J Biomech Eng. 1991;113:263-269. — **a computational FE knee model**, not a cartilage
  material-testing study (the same paper already used for this repo's own ligament work).
- Ref 12: Caruntu DI, Hefzy MS. *3-D anatomically based dynamic modeling of the human knee...* J
  Biomech Eng. 2004;126:44-53. — **also computational**, not empirical.

So the 5 MPa traces through a **chain of computational-modeling reuse** (Blankevoort & Huiskes 1991
→ Caruntu & Hefzy 2004 → Lenhart 2015), not a direct citation to a cartilage indentation/mechanical-
testing study. This is a materially different, more critical provenance than "measured on tissue and
found to be 5 MPa" — and it is the mechanistic reason the material and geometry claims get opposite
verdicts: **the geometry is this donor's own real MRI-derived surface reconstruction; the material is
a twice-removed modeling assumption.**

**Forced steelman for the model (before accepting the FAIL at face value)**: two legitimate
mitigating arguments were checked, not just asserted:

1. **Regime mismatch, not fully resolved.** The live-verified Shepherd & Seedhom 1999 paper (PMID
   10342624) is titled *"the 'instantaneous' compressive modulus"* — a stiffer, fast/undrained
   poroelastic regime, distinct from the equilibrium/aggregate modulus (~0.5-1.5 MPa) most often
   quoted. Gait-rate loading is fast, so an "effective" gait-relevant modulus could legitimately sit
   above the equilibrium value. **This session could not extract Shepherd & Seedhom's own topline
   instantaneous-modulus number live** (abstract fetches returned qualitative comparisons only,
   full text paywalled) — so this fork is disclosed as genuinely open, not resolved in the model's
   favor by assumption.
2. **Lumped elastic-foundation model form.** The contact law (Bei & Fregly 2003) explicitly combines
   thickness and modulus of BOTH contacting meshes into a lumped effective stiffness — a known
   simplification (independent-spring foundation, no inter-triangle shear) that can require an
   adjusted effective modulus to reproduce realistic pressure/area, distinct from a literal claim
   about tissue material modulus. This is a real methodological point, but does **not** change the
   fact that 5 MPa, taken as a literal material property and checked against decorrelated tissue-
   testing literature, is a genuine, disclosed mismatch — it only explains a plausible motive.

**What this does NOT do**: it does not re-litigate `docs/MECHANISM_CARTILAGE_CONTACT.md`'s own finding
that this model's PRESSURE output (6.2/2.8 MPa gait) matches the same paper's own published pressure
figures — that is a same-donor, same-paper internal-consistency check (a good no-bug confirmation),
**not** an independent external validation of the 5 MPa modulus itself, and the two are kept
deliberately decorrelated here (an aggregation-cell trap this repo's own method notes warn about:
internal self-consistency is not external validation).

---

## 5. Confidence tier and falsifier

**Confidence tier, split as the task requested:**
- **Geometry (thickness map): in-vivo-anchored.** Anchors are MRI morphometry studies (Cohen/
  Ateshian 1999, PMID 10367018, in-vitro-calibrated in-vivo MRI; Frobell/Eckstein 2010, PMID
  20496431, n=1003 in-vivo MRI, Eckstein F senior author — the task's named "Eckstein" anchor)
  plus a cadaveric cross-check (Shepherd & Seedhom 1999, PMID 10343537). All 3 PMIDs' identity/
  topic/authors were confirmed LIVE this session (NCBI eutils esummary + PubMed abstract fetch).
- **Material (elastic modulus): cadaveric.** Anchors are Shepherd & Seedhom 1999 (PMID 10342624,
  instantaneous modulus) and Athanasiou et al. 1991 (PMID 2010837, aggregate modulus, distal
  femoral cartilage, cross-species) — both live-identity-verified.
- **Disclosed verification gap (honest, not hidden)**: exact topline mm/MPa figures from the primary
  papers' own Results sections were **not** extracted verbatim live this session — several abstracts
  state only comparisons ("knee > ankle," "ankle > knee") rather than a single topline absolute
  number, the Frobell/Eckstein abstract reports Z-scores/percent rather than mm, and full texts were
  paywalled (Wiley HTTP 402). The numeric bands used as the decision threshold are the ones the
  **task itself pre-registered** before any measurement here was computed, cross-checked once
  against a tertiary source (Wikipedia, citing named modern reviews) that landed inside the same
  band. This is reported explicitly so the reader can weigh it, not smoothed over.

**Falsifier for the GEOMETRY verdict**: a finer regional split (using the JAM h5's own existing
6-region breakdown to isolate medial vs. lateral tibial plateau) that shows the medial/central
weight-bearing region alone sits materially higher than the 1.89-2.10mm found here would upgrade
tibial from "borderline" to "pass"; conversely, access to the Frobell/Eckstein 2010 full text
showing this donor's implied population percentile is itself outside 2 SD would downgrade it.

**Falsifier for the MATERIAL verdict**: live extraction of Shepherd & Seedhom's own instantaneous-
modulus topline figure showing it is ≥ ~3-5 MPa at physiological gait strain rates would soften (not
eliminate — the citation-lineage finding stands regardless) the FAIL to a regime-mismatch explanation
rather than an unexplained stiff assumption.

---

## 6. Honest gaps (full list)

1. **Single generic donor** (task's own explicit reminder) — lenhart2015 is one MRI-imaged healthy
   subject; geometry and material are both n=1, compared against population-level (geometry) or
   small-cadaveric-series (material) anchors. Not resolved by this leg; a structural limitation of
   the source model, already flagged elsewhere in this repo's knee docs.
2. **Exact primary-source topline mm/MPa figures not extracted live** (§5) — paper identity/topic
   fully verified live; the specific numbers rest on the task's own pre-registered bands plus one
   tertiary cross-check, not a verbatim live re-derivation from the named papers' Results sections.
3. **Tibial "loaded subset" is a coarse >1kPa-anywhere definition**, not a medial/lateral split — the
   h5's own 6-region breakdown could refine this; not done here (lean/scope).
4. **~5.6% (femoral/tibial) to 1.4% (patellar) of triangles have no valid ray-cast hit** within the
   5cm/K=150 search — excluded from statistics (not zero-filled), concentrated at mesh
   rim/periphery as expected; a bounded, disclosed limitation, not a silent gap-fill.
5. **Regime mismatch for the material anchor not fully resolved** (§4.1) — instantaneous vs
   equilibrium modulus fork disclosed, not closed, for lack of a live-extracted topline number.
6. **Poisson's ratio (0.45) not independently anchor-checked** — out of the task's explicit scope
   (thickness + elastic modulus only), reported from the XML for completeness but not validated
   against a decorrelated external source here.
7. **Static geometry only** — this characterizes the thickness map and modulus as fixed model
   inputs; it does not re-run or re-validate the dynamic contact-pressure OUTPUT (already covered,
   decorrelated, in `docs/MECHANISM_CARTILAGE_CONTACT.md`/`docs/MECHANISM_JAM_CONTACT_DECORR.md`).
8. **"Li" (named by the task alongside Eckstein for the geometry anchor) was NOT independently
   verified this session** — WebSearch's session budget was exhausted (2000/2000) before a
   Li-specific query could be run, and no PMID for a "Li knee cartilage thickness" paper is
   reported anywhere above. This is disclosed explicitly rather than silently substituted:
   Cohen/Ateshian 1999 + Frobell/Eckstein 2010 + Shepherd & Seedhom 1999 were used as the live-
   verified geometry anchors instead, all independently on-topic and PMID-confirmed, but they are
   not "Li." Similarly, the "Eckstein" match (Frobell et al. 2010) is Eckstein as **senior/last
   author** of a multi-author OAI consortium paper, not a solo-authored Eckstein paper — reported
   precisely, not rounded up to a cleaner-sounding match.

---

## 7. Files

- `scripts/msk/cartilage_thickness_material.py` — hand-rolled ASCII+binary STL parsers, vectorized
  Möller–Trumbore ray-cast (with empirical direction decision + K-neighbor retry), settle-window-
  corrected gait-peak detection, literature-anchor tables, verdict computation.
- `data/msk_smoketest/cartilage_thickness_material/cartilage_thickness_material_results.json` — full
  machine-readable evidence: per-mesh-pair thickness stats (whole-mesh + loaded-subset, raw + clipped,
  winning AND losing ray direction), cross-checks (facet-count-vs-h5, degenerate-triangle count,
  nearest-centroid agreement), literature anchor table with per-citation live-verification notes, and
  the geometry/material verdict blocks.
- Inputs read in place, unchanged (mtime-reverified in the JSON's `provenance` block):
  `scratchpad/knee_lig/lenhart2015.osim`, `scratchpad/knee_lig/Geometry/*.stl` (6 meshes),
  `data/msk_smoketest/jam_contact_decorr/results/joint-mechanics/gait_driven.h5`.

No git commit, no git push performed (isolation respected).
