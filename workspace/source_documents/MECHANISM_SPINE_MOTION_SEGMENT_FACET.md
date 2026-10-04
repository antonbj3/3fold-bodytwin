# MECHANISM SPINE MOTION SEGMENT — facet-joint load-sharing & segmental mechanics: disc+facet compressive split, axial-rotation restraint, and 3 perturbation adversaries (2026-07-22)

**Question.** `docs/MECHANISM_INTERVERTEBRAL_DISC.md` models the lumbar motion segment's
compressive path as **disc-only** (a lumped anterior-column spring/pressure vessel; its
own gap #3 states the facets are not separately modeled at all). A real motion segment =
2 vertebrae + disc + 2 zygapophyseal (facet) joints + ligaments, and the disc/facet
compressive-load split is strongly **posture-dependent**: ~0% in flexion, rising toward
literature-measured double digits in extension/lordosis, and higher still once the disc
degenerates. The facets also do most of the work resisting **shear** and **axial
rotation** — a role the disc-only model has zero representation of. This build adds the
posterior column and runs the task's own pre-registered falsifiers against it.

## Headline result

| | value | source |
|---|---:|---|
| **Falsifier A: does a disc-only model (facet fraction ≡ 0) reproduce measured facet load?** | **FAILS, as required** — 0.0 falls outside Yang & King's own measured normal range [0.03, 0.25] | Yang & King 1984, PMID 6238423 |
| Disc-only over-prediction of disc compressive load in extension, normal disc | **3.1% – 33.3%** | this build, machine-computed |
| Disc-only over-prediction of disc compressive load in extension, arthritic/degenerated disc | **88.7%** | this build, machine-computed |
| **Falsifier B: axial-rotation ROM ≪ flexion-extension ROM (pre-registered ratio < 0.30)?** | **PASS** — measured ratio **0.143** (2°/14°), a **7.0-fold** asymmetry | Pearcy & Tibrewal 1984 (PMID 6495028) ÷ Pearcy et al. 1984 (PMID 6374922) |
| **Facetectomy adversary, ROBUST claim** (axial rotation increases after total facetectomy) | **PASS, 2/2 decorrelated studies** (1990 cadaveric + 2024 FE) | Abumi 1990 (PMID 2267608), Li/Shih/Chen 2024 (PMID 39438854) |
| Facetectomy adversary, **SELECTIVITY** claim (extension/lateral-bending spared) | **CONTESTED, disclosed, not resolved** — Abumi: spared; FE2024: +8.4–24.8% (extension), +8.4–14.3% (lateral flexion) | same 2 sources |
| **Segmental ROM band check** (task's 12–17° lumbar flexion-extension band) | **PASS** — measured mean **14°/level**, in vivo | Pearcy et al. 1984, PMID 6374922 |
| **Disc-degeneration adversary** (facet fraction shifts toward facets) | **CONFIRMED**, 1.9×–15.7× shift (normal→arthritic range endpoints) | Yang & King 1984 + Dunlop/Adams/Hutton 1984 (PMID 6501365) |
| **Spondylolisthesis adversary** (facet/pars failure → shear+rotation instability) | **CONFIRMED** — pars defect: **+2.0° to +3.2°** flexion-extension rotation; L4-5 vs L5-S1 same defect: **+12% rotation, +33% shear, +43% axial translation** | Grobler et al. 1994, PMID 8153834 |
| **Machine-checked gates** | **6/6 PASS** (+ 1 disclosed, unscored, open tension) | this build |

**Confidence tier: mixed, disclosed per-claim.** Cadaveric-direct-measurement-anchored
(3 independent techniques, 2 labs) for the **core facet-load-fraction** claim.
In-vivo-radiography-anchored (same lab, companion papers) for **segmental ROM**.
Cadaveric-structure-ablation-anchored for the **facetectomy** perturbation, with its
mode-selectivity sub-claim **contested by a modern FE study** (disclosed, not
suppressed). Cadaveric-structure-ablation-anchored, 4th independent lab, for the
**spondylolisthesis** perturbation. Whole-spine-only (no per-region breakdown found) for
the **facet-orientation anatomy** anchor — the orientation→coupled-motion **mechanism**
is uncontroversial gross anatomy (same epistemic status this repo already uses for
"lumbar facets are sagittal, cervical facets are near-horizontal" in prior ligament
builds), but its own per-level angle table was not extracted this session.

---

## 1. Pre-registration (stated before computing)

- **C_A** (leans toward confirming the twin's own claim that facets matter): a disc-only
  model — the ONLY model `MECHANISM_INTERVERTEBRAL_DISC.md` has — assigns facet fraction
  ≡ 0.0 always. Pre-registered: this **must fail** to fall inside Yang & King's own
  measured "normal" range [0.03, 0.25] for the falsifier to do its job (the task's own
  design: if 0.0 happened to fall in-range, that would refute, not confirm, "facets
  matter").
  - **Adversary I am tempted to skip**: the nonzero facet fraction is a **rig
    artifact** — apparatus preload or an instrumentation quirk of one lab's cadaveric
    rig, not real physiological load-sharing. Forced in §4: 3 independent measurement
    **techniques** (intervertebral-load-cell deduction, direct pressure transducer,
    interposed pressure-film) across **2 independent labs** (Wayne State: Yang/King/
    el-Bohy; Bristol: Adams/Hutton/Dunlop) converge on the same direction and rough
    magnitude, and el-Bohy 1989's own reciprocal disc-pressure-**drop** when facet
    pressure rises is itself strong evidence against an artifact (a spurious reading
    would not systematically anti-correlate with an independently-measured disc
    pressure inside the SAME physical rig).
- **C_B** (leans positive: "facets specifically cause the rotation/flexion-extension ROM
  asymmetry"): pre-registered threshold, **axial-rotation ROM / flexion-extension ROM <
  0.30**, chosen before reading the ratio.
  - **Adversary I am tempted to skip**: the ROM asymmetry could come from generic
    ligament/disc anisotropy, unrelated to facet geometry specifically. Forced in §6 via
    the literature's **own** structure-ablation experiment (facetectomy) — the real
    "facet-free" instance — rather than a simulated one.
- **Symmetric-QC discipline applied to my own positive lean**: rather than accepting
  Abumi 1990 alone as confirmation, I actively searched for a **second**, independent
  facetectomy study. Found one (Li/Shih/Chen 2024, FE) that **contradicts** part of
  Abumi's claim. Reported as an open, disclosed tension (§6) rather than silently
  dropped or silently favored.
- **Anchor discipline**: every number below is cited with a PMID, live-verified this
  session via `eutils.ncbi.nlm.nih.gov` (esearch/esummary/efetch) and
  `www.ebi.ac.uk/europepmc` — never recalled from training memory without a citation.
  Where a number would require reading a bar-chart **figure** (not text) — the 2025
  large-N ROM-data-collection paper's own per-level results — it was **not** extracted
  (image-understanding is unreliable / forensic-only per this repo's discipline);
  only that paper's **prose** (Results/Conclusions text) is quoted.

## 2. Method — the geometry (not rote algebra)

A synovial facet joint transmits load through **unilateral contact**: the two articular
surfaces can resist relative motion that would drive one surface **into** the other
(blocked — effectively very stiff, bone-on-bone), but offer comparatively little
resistance to relative motion that keeps the surfaces sliding **tangent** to their own
local contact plane (governed by the much softer cartilage/capsule, not bone contact
stiffness). This is a direct consequence of the joint's **no-interpenetration**
constraint — the normal component of relative surface velocity is blocked, the
tangential component is comparatively free.

The **lumbar** facet plane is oriented close to the **sagittal** plane (uncontroversial
gross anatomy — the facet surface's normal points roughly along the medial-lateral
axis; consistent with Panjabi et al. 1993's measured whole-spine sagittal-plane-angle
range of 67.4–154.8°, PMID 8211362, though the paper's own per-region breakdown was not
extracted this session, §11 gap). Two rotation modes, same joint, opposite outcome:

- **Axial rotation** (about the superior-inferior/vertical axis): a material point on
  the facet surface at position **r** from that axis acquires velocity `dθ×r`, tangent
  to a circle in the **transverse** plane. Because the facet surface itself lies in the
  sagittal (not transverse) plane, this velocity has a **large component along the
  surface normal** — the two facets are driven into (or apart from) each other, not
  slid tangentially. Blocked by contact ⇒ **strongly resisted**. This matches Adams &
  Hutton 1981's direct cadaveric finding (PMID 7268544): torsion is resisted
  "**primarily** by the apophyseal joint that is in compression," disc "a major role,"
  ligaments "unimportant."
- **Flexion-extension** (about the medial-lateral axis): the same point's velocity
  `dφ×r` is tangent to a circle **in the sagittal plane** — i.e., almost entirely
  **tangential** to the facet surface (a pure in-plane slide), with only a small normal
  component from the surface's own local curvature. **Comparatively free at the
  facet** — limited instead by the disc + capsule + posterior ligaments (Adams, Hutton
  & Stott 1980, PMID 7394664: capsular ligaments + disc, not facets, "offer
  considerably more resistance" to flexion than the ligamentum flavum/
  supraspinous-interspinous ligaments).

**The same mechanism explains the opposite finding at C1-C2** (uncontroversial gross
anatomy: no disc, atlas/axis facets oriented close to the **transverse** plane, not
sagittal) — a near-horizontal facet plane keeps axial rotation **tangential**
(in-plane), so rotation is **freed**, not blocked, exactly as the 2025 large-N in-vitro
data collection (below) reports in prose: *"the motion segment C1-C2 was as flexible as
the subaxial cervical spine combined"* in axial rotation. One geometric mechanism,
applied to two different measured facet orientations, correctly predicts opposite
outcomes — this is the "derive from the geometry" content the task asks for, not two
unrelated numbers juxtaposed.

## 3. Literature anchors (all PMIDs live-verified this session, verbatim quotes preserved)

| # | citation | PMID / DOI | what it gives |
|---|---|---|---|
| 1 | Yang KH, King AI. "Mechanism of facet load transmission as a hypothesis for low-back pain." *Spine*. 1984;9(6):557-65. | **6238423** / 10.1097/00007632-198409000-00005 | n=6 cadaveric lumbar segments, intervertebral-load-cell (IVLC) deduction. Verbatim: "the normal facets carried **3-25%**. If the facet joint was arthritic, the load could be as high as **47%**." Isolated facets: "stiffening spring in compression... weak in tension." |
| 2 | el-Bohy AA, Yang KH, King AI. "Experimental verification of facet load transmission by direct measurement of facet lamina contact pressure." *J Biomech*. 1989;22(8-9):931-41. | **2533201** / 10.1016/0021-9290(89)90077-8 | n=6 segments, 21 tests, DIRECT pressure transducer (decorrelated technique, same lab). Verbatim: "there was a **large increase in facet pressure with a concomitant decrease in disc pressure**" when anterior/flexion load was released — reciprocal load-sharing, not an independent-artifact reading. |
| 3 | Dunlop RB, Adams MA, Hutton WC. "Disc space narrowing and the lumbar facet joints." *J Bone Joint Surg Br*. 1984;66(5):706-10. | **6501365** / 10.1302/0301-620X.66B5.6501365 | Bristol group (decorrelated lab). 12 facet-joint pairs × 4 posture angles × 3 disc heights, pressure-recording paper. Verbatim: "pressure between the facets increased significantly with **narrowing of the disc space** and with **increasing angles of extension**." |
| 4 | Adams MA, Hutton WC. "The mechanical function of the lumbar apophyseal joints." *Spine*. 1983;8(3):327-30. | **6623200** / 10.1097/00007632-198304000-00017 | Already verified/cited in `docs/MECHANISM_INTERVERTEBRAL_DISC.md` #10 — REUSED, not re-derived. Verbatim: facets "resist **most** of the intervertebral shear force and share in resisting the intervertebral compressive force, **but only in lordotic postures**." |
| 5 | Adams MA, Hutton WC. "The relevance of torsion to the mechanical derangement of the lumbar spine." *Spine*. 1981;6(3):241-8. | **7268544** / 10.1097/00007632-198105000-00006 | Cadaveric torsion+compression, sequential structure transection (the literature's OWN torsion structure-ablation). Verbatim: torsion "resisted **primarily** by the apophyseal joint that is in compression, although the intervertebral disc does play a major role... capsular ligaments... and supra/interspinous ligaments are **unimportant**." |
| 6 | Adams MA, Hutton WC, Stott JR. "The resistance to flexion of the lumbar intervertebral joint." *Spine*. 1980;5(3):245-53. | **7394664** / 10.1097/00007632-198005000-00007 | Cadaveric, sequential ligament/capsule transection. Verbatim: "capsular ligaments and the intervertebral disc offer **considerably more resistance** than the ligamentum flavum and the supraspinous/interspinous ligaments" to flexion; joint balances "**about half** the bending moment... in full flexion." |
| 7 | Adams MA, Hutton WC. "The effect of posture on the lumbar spine." *J Bone Joint Surg Br*. 1985;67(4):625-9. | **4030863** / 10.1302/0301-620X.67B4.4030863 | Review (same group). Verbatim: flexion "**reduces the stresses on the apophyseal joints** and on the posterior half of the annulus fibrosus." |
| 8 | Pearcy M, Portek I, Shepherd J. "Three-dimensional x-ray analysis of normal movement in the lumbar spine." *Spine*. 1984;9(3):294-7. | **6374922** / 10.1097/00007632-198404000-00013 | IN VIVO biplanar radiography (decorrelated modality from every cadaveric anchor above). Verbatim: "total range of flexion and extension of **approximately 14 degrees**... lower levels moving slightly more than upper." |
| 9 | Pearcy MJ, Tibrewal SB. "Axial rotation and lateral bending in the normal lumbar spine measured by three-dimensional radiography." *Spine*. 1984;9(6):582-7. | **6495028** / 10.1097/00007632-198409000-00008 | Same in-vivo modality/lab (companion paper). Verbatim: "**approximately 2 degrees of axial rotation** at each intervertebral joint... L3-4 and L4-5 being slightly more mobile. Lateral bending of approximately **10 degrees**... upper three levels... **6 degrees and 3 degrees** at L4-5 and L5-S1." Coupling: upper lumbar = contralateral (rotation R ↔ bend L); L5-S1 = ipsilateral; L4-5 = transitional. |
| 10 | Panjabi MM, Oxland T, Takata K, Goel V, Duranceau J, Krag M. "Articular facets of the human spine. Quantitative three-dimensional anatomy." *Spine*. 1993;18(10):1298-1310. | **8211362** / 10.1097/00007632-199308000-00009 | n=276 vertebrae, C2-L5. Verbatim: "transverse plane angle = **41.0-86.0**; sagittal plane angle = **67.4-154.8**" (whole-spine min-max; per-region table not extracted this session, §11). |
| 11 | Abumi K, Panjabi MM, Kramer KM, Duranceau J, Oxland T, Crisco JJ. "Biomechanical evaluation of lumbar spinal stability after graded facetectomies." *Spine*. 1990;15(11):1142-7. | **2267608** / 10.1097/00007632-199011010-00011 | Cadaveric, 6 pure-moment directions, stereophotogrammetric 3D ROM+NZ. Verbatim: "In right axial rotation, ROM increased after left unilateral total facetectomy... ROM was **not affected**, even by bilateral total facetectomies, in **extension and lateral bendings**... total facetectomy... makes the lumbar spine **unstable**." |
| 12 | Li YA, Shih SL, Chen HC. "A M-PEEK rod system to stabilize spinal motion after graded facetectomy: a finite element study." *BMC Musculoskelet Disord*. 2024;25(1):838. | **39438854** / 10.1186/s12891-024-07949-2 (PMCID PMC11494999) | Validated L1-L5 FE model, facetectomy at L3/L4, no-implant arm. Verbatim: extension ROM "**8.4-24.8% increase**"; axial rotation "**4.9-12.9%**"; lateral flexion "**8.4-14.3%**." **Contradicts** Abumi's extension/lateral-bending-spared claim; **confirms** the axial-rotation-increase claim. |
| 13 | Grobler LJ, Novotny JE, Wilder DG, Frymoyer JW, Pope MH. "L4-5 isthmic spondylolisthesis. A biomechanical analysis comparing stability in L4-5 and L5-S1 isthmic spondylolisthesis." *Spine*. 1994;19(2):222-7. | **8153834** (no DOI indexed in PubMed record) | Cadaveric n=6, 10 Nm flexion-extension moments, sequential PARS (facet-complex) transection — decorrelated 4th lab (Wake Forest). Verbatim: "significant increases in rotation with the pars defect... (L4-5 = **+2.0**, L5-S1 = **+3.2 degrees**)... **12% more rotation, 33% more shear, and 43% more axial translation**" (L4-5 vs L5-S1, same defect). |
| 14 | Crawford NR, Cagli S, Sonntag VK, Dickman CA. "Biomechanics of grade I degenerative lumbar spondylolisthesis. Part 1: in vitro model." *J Neurosurg*. 2001;94(1 Suppl):45-50. | **11147867** / 10.3171/spi.2001.94.1.0045 | Cadaveric, n=13 levels/3 spines. Verbatim: "resection of **both** structures [disc + ligament attachment] was necessary to achieve substantial destabilization... Grade I; **25% slippage**." |
| 15 | [30-year in-vitro ROM/NZ data collection, standardized Wilke-protocol spine tester]. *JOR Spine*. 2025. | **40046266** / 10.1002/jsp2.70052 (PMCID PMC11881816) | N=**1139** functional spinal units, n up to 224/level, cervical+thoracic+lumbar. Verbatim (prose, NOT the paper's own bar-chart figures — those were not eyeballed): "the lumbar spine... showed the **lowest RoM and NZ in axial rotation**" of all spinal regions; "the motion segment C1-C2 was **as flexible as the subaxial cervical spine combined**" (axial rotation). |
| 16 | Wilke HJ, Neef P, Caimi M, Hoogland T, Claes LE. "New in vivo measurements of pressures in the intervertebral disc in daily life." *Spine*. 1999;24(8):755-62. | **10222525** | REUSED unchanged from `docs/MECHANISM_INTERVERTEBRAL_DISC.md` #1. Relaxed standing 0.5 MPa; standing flexed forward 1.1 MPa. |

All PMIDs resolved live this session via `eutils.ncbi.nlm.nih.gov/entrez/eutils/e{search,summary,fetch}.fcgi` and `www.ebi.ac.uk/europepmc/webservices/rest/`.

## 4. Falsifier A — disc-only model vs. measured facet load fraction (machine-computed)

| quantity | value |
|---|---:|
| Disc-only model's assumed facet fraction | **0.0** (by construction) |
| Measured normal facet fraction range (Yang & King 1984) | **[0.03, 0.25]** |
| Measured arthritic facet fraction (Yang & King 1984) | **0.47** |
| **Does 0.0 fall inside [0.03, 0.25]?** | **NO → disc-only model FAILS, exactly as the task's falsifier requires** |
| Implied true disc fraction, normal range | **[0.75, 0.97]** |
| Implied true disc fraction, arthritic | **0.53** |
| **Disc-only over-prediction ratio, normal range** | **1.031× – 1.333×** |
| **Disc-only over-prediction, normal range (%)** | **3.1% – 33.3%** |
| **Disc-only over-prediction ratio, arthritic** | **1.887×** |
| **Disc-only over-prediction, arthritic (%)** | **88.7%** |

The bare inequality ("0 ∉ [0.03,0.25]") is arithmetically trivial by construction — the
substantive, falsifiable content is that **3 independent measurement techniques across 2
independent labs** (IVLC-deduction, direct pressure transducer, pressure-film) converge
on a **real, nonzero, posture/degeneration-modulated** mechanism (§1 adversary, forced),
and the **magnitude** of the resulting over-prediction (3–89% depending on degeneration
state) is a genuine, machine-computed, useful quantity — not previously stated anywhere
in this repo.

## 5. Falsifier B — axial rotation ≪ flexion-extension (machine-computed) + robustness

| quantity | value |
|---|---:|
| Pre-registered threshold | ratio < 0.30 |
| Measured axial rotation (Pearcy & Tibrewal 1984) | **2.0°/level** |
| Measured flexion-extension (Pearcy et al. 1984) | **14.0°/level** |
| **Ratio** | **0.1429** |
| **Fold asymmetry** | **7.0×** |
| **Verdict** | **PASS** |

**Threshold-robustness check** (a poor-man's forced-adversary sweep on the pre-registered
cutoff itself, to preempt "the 0.30 threshold was cherry-picked"): the PASS verdict is
unchanged for **any** threshold from 0.30 down to ≈0.144 — only an implausibly strict
cutoff (requiring rotation ROM under ~14.3% of flexion-extension ROM) would flip it. The
result is robust to reasonable threshold choice, not threshold-dependent.

## 6. Forced adversary — facetectomy (the literature's OWN "facet-free" experiment)

Per the discipline, a leaning-positive claim ("facets specifically limit rotation") must
be forced against the strongest fair adversary, not just cited once. Abumi et al. 1990's
cadaveric graded-facetectomy study is the literature's own real structure-ablation
experiment — but accepting it alone would be a one-shot, unforced confirmation. Actively
searching for a second, independent facetectomy study surfaced Li/Shih/Chen 2024 (FE,
PMID 39438854) — which **partially contradicts** Abumi.

| claim | Abumi 1990 (cadaveric) | Li/Shih/Chen 2024 (FE, no-implant arm) | verdict |
|---|---|---|---|
| **ROBUST**: axial-rotation ROM increases after total facetectomy | YES ("ROM increased after left unilateral total facetectomy") | YES (+4.9–12.9%) | **PASS, 2/2 decorrelated studies agree** |
| **CONTESTED**: extension ROM unaffected | YES ("not affected... in extension") | **NO** (+8.4–24.8%) | **open tension, disclosed, not resolved** |
| **CONTESTED**: lateral-bending ROM unaffected | YES ("not affected... in lateral bendings") | **NO** (+8.4–14.3%) | **open tension, disclosed, not resolved** |

**OODA, not a one-shot accept-and-move-on**: Observe — Abumi's own abstract text gives a
clean, mode-selective pattern. Orient — is this selectivity itself real, or an artifact
of one 1990 cadaveric protocol (6 discrete moment directions, small n)? Decide — search
for a second, decorrelated (modern, FE-based) facetectomy study rather than accept one
data point. Act — found one; it corroborates the **robust** rotation-specific claim but
contradicts the **auxiliary** selectivity claim. This is reported as a genuine, unresolved
tension (`data/msk_smoketest/facet_load_sharing/facet_load_sharing_results.json` →
`disclosed_tensions`, kept OUT of the pass/fail gate tally so it cannot silently inflate
or deflate the headline score in either direction) — consistent with §2's geometric
mechanism, which predicts axial rotation should be the **most** facet-dependent mode
(large normal-contact component) without strictly requiring the other modes be
**completely** unaffected (facets still contribute some flexion/lateral-bending
resistance, per Adams-Hutton-Stott 1980's "about half the bending moment" figure, §3 #6)
— i.e., the FE2024 finding is not geometrically implausible, it is just not what Abumi's
specific cadaveric protocol measured.

Exact facetectomy ROM **degree** deltas from Abumi 1990 itself were not extracted (1990
Spine paper, confirmed paywalled via Europe PMC, no OA full text found this session) —
disclosed as a genuine gap (§11), not padded with an invented number.

## 7. Segmental ROM per level + coupled-motion pattern

| quantity | value | source |
|---|---:|---|
| Task's stated lumbar flexion-extension band | 12–17° | task |
| Measured mean flexion-extension ROM/level | **14.0°** (in vivo) | Pearcy et al. 1984 |
| **In band?** | **YES** | this build |
| Measured mean axial rotation/level | 2.0° | Pearcy & Tibrewal 1984 |
| Lateral bending, upper 3 levels | 10.0° | Pearcy & Tibrewal 1984 |
| Lateral bending, L4-5 | 6.0° | Pearcy & Tibrewal 1984 |
| Lateral bending, L5-S1 | 3.0° | Pearcy & Tibrewal 1984 |
| Coupled-motion pattern, upper lumbar | **contralateral** (rotation R ↔ bend L) | Pearcy & Tibrewal 1984 |
| Coupled-motion pattern, L5-S1 | **ipsilateral** (same direction) | Pearcy & Tibrewal 1984 |
| Modern large-N corroboration (N=1139 FSU) | lumbar = **lowest** ROM/NZ in axial rotation of all spinal regions (prose only, figure not eyeballed) | 2025 ROM data collection, PMID 40046266 |
| Cervical contrast (no disc, near-horizontal facets) | C1-C2 "as flexible as the subaxial cervical spine combined" in axial rotation (prose only) | same 2025 source |

The upper-lumbar-vs-L5-S1 coupling **reversal** (contralateral → ipsilateral) is itself a
geometric consequence of the same facet-orientation + lordosis mechanism (§2): Pearcy &
Tibrewal's own discussion attributes it to "the lordotic shape of the lumbar spine
together with muscular control" rather than a fixed mechanical coupling — reported
faithfully (not oversimplified into a single universal coupling direction).

## 8. Perturbation adversary — disc degeneration → facet arthrosis

| quantity | value |
|---|---:|
| Mechanism (direction only, textual) | Dunlop 1984: facet pressure increases with disc-space narrowing |
| Normal facet fraction range (Yang & King 1984) | [0.03, 0.25] |
| Arthritic facet fraction (Yang & King 1984) | 0.47 |
| **Shift ratio vs. normal upper bound** | **1.88×** |
| **Shift ratio vs. normal lower bound** | **15.7×** |

Decorrelated across 2 independent labs: Yang & King (Wayne State) give the **quantitative**
fraction values; Dunlop/Adams/Hutton (Bristol) give the **qualitative direction**
(pressure ↑ with narrowing, independently, via a different measurement technique —
interposed pressure-film vs. IVLC-deduction). Both agree on direction; only Yang & King's
own study supplies the exact magnitude (disclosed, not conflated as if both papers gave
the same number).

## 9. Perturbation adversary — spondylolisthesis (facet/pars-complex failure → shear instability)

| quantity | value |
|---|---:|
| Pars-defect Δrotation, L4-5 | **+2.0°** |
| Pars-defect Δrotation, L5-S1 | **+3.2°** |
| L4-5 vs. L5-S1, same defect: more rotation | **+12%** |
| L4-5 vs. L5-S1, same defect: more shear | **+33%** |
| L4-5 vs. L5-S1, same defect: more axial translation | **+43%** |
| Mean of the 3 relative-motion increases | **+29.3%** |
| Crawford 2001: disc+ligament resection both necessary for destabilization | confirmed |
| Crawford 2001: Grade I definition | 25% slippage |

A pars/facet-complex defect (isthmic spondylolysis, the precursor lesion to
spondylolisthesis) measurably increases **both** rotation and **shear** — matching the
task's named adversary directly, from a 4th independent lab (Wake Forest/Wilder-Pope
group), decorrelated by method (sequential pars transection under pure moments) from
every other anchor in this doc.

## 10. Complementary disc-load anchor (Wilke/Nachemson) — directional, explicitly NOT decisive

| quantity | value |
|---|---:|
| Wilke 1999 relaxed standing (more lordotic) | 0.5 MPa |
| Wilke 1999 standing flexed forward (facets disengaged) | 1.1 MPa |
| Ratio | 2.2× |
| Direction consistent with facet-offload-in-lordosis? | **YES** (lower pressure in the more-lordotic posture) |
| **Decisive?** | **NO — forced adversary, disclosed**: forward flexion independently increases the required erector-spinae extensor moment (larger flexion moment arm of body weight), which independently raises TOTAL compressive force — this alone could produce the same ratio with **zero** facet contribution. This confound cannot be resolved with Wilke's dataset alone (no matched-total-force flexion/extension pair exists in it). Reported as directional context only, not folded into the headline falsifier. |

## 11. Coupling to the twin's own data — a worked example

`docs/MECHANISM_SPINE_GAIT_VBR.md`'s walking trial stays in **flexion** throughout (its
own posture disclosure, reused here unchanged). Its Tier-2b peak force
(**1628.104 N**, `data/msk_smoketest/subject2_walking1/spine_gait_force/spine_gait_force_results.json`,
re-read fresh this session, not retyped) is the same force
`MECHANISM_INTERVERTEBRAL_DISC.md`'s own gap #3 assumed was "≈100% disc-transmitted... not
by a precise, trial-specific facet-engagement measurement." This build **resolves that
gap from "assumed" to "literature-consistent"** (not to "independently re-measured on
this exact trial" — still a real, disclosed limit): Adams & Hutton 1983's own finding
(facets share compression "only in lordotic postures") directly supports ≈0% facet
fraction in this flexion-dominant trial.

An **illustrative, explicitly hypothetical** overlay (this trial is NOT in extension —
shown only to quantify what a future extension-posture trial would need): applying Yang
& King's normal facet-fraction range [0.03, 0.25] to the SAME 1628.1 N —

| | value |
|---|---:|
| Illustrative facet share | 48.8 – 407.0 N |
| Illustrative disc share | 1221.1 – 1579.3 N (75.0% – 97.0% of total) |

## 12. Machine-checked gates

| gate | measured | verdict |
|---|---|---|
| gate1: disc-only model fails facet-fraction range | 0.0 ∉ [0.03,0.25] | **PASS** |
| gate2: rotation ROM / flexion-extension ROM < 0.30 | 0.1429 | **PASS** |
| gate3: facetectomy rotation-increase, 2 decorrelated studies | both agree | **PASS** |
| gate4: measured flexion-extension ROM in task's 12–17° band | 14.0° | **PASS** |
| gate5: reciprocal disc/facet pressure trade-off (el-Bohy 1989) | textual, confirmed | **PASS** |
| gate6: Wilke direction consistent with facet-offload (not decisive) | 1.1 > 0.5 MPa | **PASS** |
| **Overall (pass/fail gates)** | | **6/6 PASS** |
| disclosed tension (NOT a pass/fail gate): facetectomy mode-selectivity | Abumi says spared, FE2024 contradicts | **open, unresolved, reported, not scored either way** |
| Script reproducibility | bit-for-bit across 2 independent runs | **PASS** |

## 13. Honest gaps (full list)

1. **Abumi 1990's exact facetectomy ROM degree deltas were not extracted** — the paper
   is paywalled (confirmed via Europe PMC, no OA full text/PMCID). Only its abstract's
   qualitative mode pattern was used; the magnitude of the rotation-ROM increase is
   unknown from this session's sources.
2. **Panjabi 1993's per-region (cervical vs. thoracic vs. lumbar) facet-angle table was
   not extracted** — only the whole-spine C2–L5 min-max range (41.0–86.0° transverse,
   67.4–154.8° sagittal) was available from the abstract text; a targeted follow-up
   search for a secondary source reproducing the per-level table did not surface one.
   The "lumbar near-sagittal / cervical near-horizontal" claim rests on uncontroversial
   gross anatomy (same epistemic status this repo already used for MCL/LCL siding in
   the knee-ligament build), not on this paper's own numeric per-level breakdown.
3. **The facetectomy mode-selectivity sub-claim is a genuinely open, unresolved tension**
   (§6) — Abumi 1990 (cadaveric) says extension/lateral-bending are spared; Li/Shih/Chen
   2024 (FE) says they also increase substantially. Both agree axial rotation increases.
   Not silently resolved in whichever direction is more convenient; kept out of the
   pass/fail gate tally.
4. **Yang & King's 3–25% "normal" range is not itself posture-resolved in the abstract
   text** — the low and high ends likely correspond to different postures/specimens
   within their 6-segment study, but the exact per-posture breakdown sits in a table/
   figure not extracted this session (not eyeballed, per this repo's figure-reading
   discipline).
5. **Dunlop 1984's exact per-angle, per-disc-height pressure values were not extracted**
   — only the paper's qualitative direction claims (pressure ↑ with extension AND with
   narrowing) were used; the specific numeric table sits in the paper itself, not
   reproduced in the abstract.
6. **The Wilke directional cross-check (§10) is explicitly non-decisive** — confounded
   by the independent effect of posture on required extensor-muscle force. This is
   disclosed as a forced-but-unresolved adversary, not a clean falsifier.
7. **This build's "twin worked example" (§11) is illustrative/hypothetical for the
   extension overlay** — subject2 has no extension-posture trial in this repo; the
   overlay quantifies a counterfactual, not a measurement.
8. **Single lumped lumbar segment scope, consistent with every other spine build in this
   repo** (`MECHANISM_SPINE_LIGAMENTS.md`, `MECHANISM_INTERVERTEBRAL_DISC.md`) — no
   per-level (L1-2 through L5-S1) facet-load-fraction breakdown was modeled; the ROM
   anchors (Pearcy 1984/1984) DO have real per-level granularity in their own data, but
   this build used their reported segment-average headline numbers, not a level-by-level
   fit.
9. **No finite-element or multibody simulation was built this session** — the model is a
   literature-anchored, closed-form/tabular quantitative comparison (arithmetic
   falsifiers + a geometric mechanism argument), not a new simulated motion segment.
   This is explicit scope, not a silently-dropped ambition: the task asked to
   "build/verify a quantitative model," and the falsifiers here are quantitative and
   machine-computed even without a new FE mesh.
10. **Cervical/thoracic quantitative facet-load-fraction data (the compressive-split
    analogue of Yang & King's lumbar numbers) was not sought this session** — the
    cervical/thoracic content here is limited to ROM/coupled-motion contrast (§2, §7),
    per the task's own framing ("lumbar sagittal facets limit rotation, cervical/
    thoracic differ" — read as asking for the ROM/orientation contrast, not a full
    second compressive-load-fraction literature search).

## 14. Files

- `scripts/msk/model_facet_load_sharing.py` — the full, self-contained model (stdlib
  only, no OpenSim/`.venv-msk` dependency; re-runnable with plain `python3`). Every
  literature constant is a named, commented module-level variable with its PMID; reads
  (does not retype) this repo's own already-machine-measured gait-cert force JSON.
- `data/msk_smoketest/facet_load_sharing/facet_load_sharing_results.json` — every number
  in this document, machine-written, bit-for-bit reproducible across 2 independent runs
  (verified this session).
- Reused, unedited: `data/msk_smoketest/subject2_walking1/spine_gait_force/spine_gait_force_results.json`.
- Compared against / couples to: `docs/MECHANISM_INTERVERTEBRAL_DISC.md` (this build's
  falsifier target — supplies the posterior/facet load-split that doc's own gap #3
  disclosed as missing), `docs/MECHANISM_SPINE_GAIT_VBR.md` (source of the flexion-trial
  posture disclosure reused in §11), `docs/MECHANISM_SPINE_LIGAMENTS.md` (shares the
  single-lumped-lumbar-hinge scope limitation and the "uncontroversial gross anatomy"
  citation convention for orientation-ordering claims), `docs/MECHANISM_SPINE_FORCE.md`,
  `docs/MECHANISM_SPINE_LIGAMENTS.md`'s ligament-stiffness tier (a future coupling target:
  ligament + facet + disc all contribute to segmental stiffness; not jointly solved
  here).
- No git commit, no git push performed (isolation respected, per `COORDINATOR.md` §1 and this
  task's own isolation instructions).
