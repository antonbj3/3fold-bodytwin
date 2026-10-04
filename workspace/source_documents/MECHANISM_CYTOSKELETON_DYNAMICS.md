# MECHANISM CYTOSKELETON DYNAMICS — actin treadmilling + microtubule dynamic instability,
# the same asymmetric-end kinetics mechanism, forced across two polymers (2026-07-22)

**Status: HYPOTHESIS awaiting independent QC. Confidence tier: IN-VITRO-ANCHORED** (purified-protein
kinetics / video and dark-field light microscopy — every rate constant and frequency below comes
from reconstituted, purified-protein systems; §6 discloses the in-vitro-vs-in-vivo gap explicitly,
not swept under). This is a **NEW cell**, checked live this session
against `data/MECHANISM_ANCHOR_GRAPH.json` before building (same dedup discipline as
`docs/MECHANISM_DNA_REPAIR_KINETICS.md`): a **SEED-DESIGN stub already exists**,
`AUTO-VERIFY-CYTOSKELETON-NUMERIC-TABLES-FULL` (`legs=[]`, no datapoints, `status:OPEN`), and its
claim text explicitly **names** "Pollard1986 rate-constant table" and a dynamic-instability review
(PMID 25823928) as **unexecuted targets**. This cell IS that execution — not a duplicate. A
second, already-populated node, `MOL-CYTOSKELETON-DYNAMICS`, covers a genuinely **different axis**
(cell-mechanical stiffness/AFM Young's-modulus vs malignancy — cofilin knockdown, microfluidic
transit velocity) and its own `regime_note` states verbatim: *"did not chase Pollard1986 or
Tinevez2009 specifically this pass (STILL OPEN, left for that stub)"* — this cell is that disclosed
gap.

Script: `scripts/msk/cytoskeleton_dynamics.py`. Evidence (machine-written):
`data/msk_smoketest/cytoskeleton_dynamics/cytoskeleton_dynamics_results.json` (md5
`f920a9b984e9d9aefa788619c94e64b7`, **byte-identical across 2 independent runs** — determinism
confirmed; pure numpy, no RNG, closed-form algebra only; all floats independently re-verified
finite by walking the parsed JSON in a fresh process, not by trusting the script's own internal
check). <1s wall time, no OpenSim dependency.

**Every citation below was verified LIVE this session** via direct `curl` against Europe PMC REST
(search + `fullTextXML` for open-access primary sources) and cross-checked title/journal/volume/
pages/year. **Two of my own initially-recalled identifiers were wrong** — a live, in-session
demonstration of the ~62%-recall-drift risk this task was pre-warned about: Hyman 1992's PMID
(recalled 1421573 → that ID actually resolves to an unrelated osteopontin paper; corrected to
**1421572** via a live title search), and Pollard 1986's full numeric Table I could **not** be
independently extracted from the primary source this session (disclosed in full below, §6) despite
two separate live fetch routes.

---

## 0. Falsifiers (pre-registered) + verdict up front

| # | Falsifier | Verdict |
|---|---|---|
| F1a | Actin barbed-vs-pointed Cc **asymmetry** (Pollard 1986, textbook-tier rate constants) clears a pre-registered ≥2× bar | **PASS**, ratio **5.10×** (Cc_barbed=0.1207 µM, Cc_pointed=0.6154 µM) |
| F1b | **Forced adversary**: a symmetric-rate-constant (no end-asymmetry) version of the identical model must predict exactly zero net flux | **FALLS, PASS** — flux = **0.0 exactly** (closed-form algebraic necessity — see §6 epistemic-tier caveat) |
| F1c | Flux-balance self-consistency at the derived treadmilling steady state | **PASS** (residual < 1e-9 — code-correctness check, see §6) |
| F1d | Independent-method (single-filament fluorescence, Fujiwara 2007) cross-check of nucleotide-state-dependent Cc, same order of magnitude | **PASS** |
| F2a | MT growth velocity (Walker 1988's own measured rate constants, converted via MT helical geometry) lands in the task's "several µm/min" regime across its own tested tubulin range | **PASS**, 3/4 points in [1,10] µm/min (0.68, **1.66, 2.32, 3.47** µm/min at 7/10/12/15.5 µM tubulin) |
| F2b | Over-determination: Cc derived from Walker's own kon/koff matches Walker's own independently-stated "~5 µM" | **PASS**, 1.1%/7.0% relative error (plus/minus ends) — internal-consistency check, see §6 |
| F2c | Catastrophe frequency ↓ / rescue frequency ↑ with tubulin concentration (Walker 1988's own data) | **PASS** (live-quoted directly, qualitative) |
| F2d | Cross-polymer distinctiveness: the SAME formula, on MT's own near-symmetric rate constants, predicts a treadmilling flux ≥10× smaller than the catastrophic-shrinkage rate scale (explaining why dynamic instability, not treadmilling, dominates pure-tubulin MT) | **PASS**, ratio **701.7×** |
| F3 | Decorrelated hydrolysis-necessity: 2 empirical legs (GMPCPP/MT; ATP-hydrolysis-dependence/actin) + 1 theoretical leg (Hill 1980, general thermodynamic proof) all sign-agree | **PASS**, 3/3 |

`overall_pass = True` (10/10 gates). **But see §6 before treating this as "10 independent
confirmations"** — 3 of the 10 gates (F1b, F1c, and half of F2b) are closed-form
model-correctness checks, not independently falsifiable empirical data points; this is disclosed
explicitly, not smoothed over.

---

## 1. Geometric mechanism — ONE two-end birth-death model, both polymers

Per this repo's geometric-thinking discipline: treadmilling is not a separate law for actin and a
separate law for microtubules — it is **one** structural fact about a polymer with two chemically
distinct ends, each independently obeying second-order association / first-order dissociation
kinetics (confirmed as the right model FORM by Pollard 1986's own data: elongation rate is linear
in monomer concentration at both ends, live-verified, §2). For an end with on-rate `kon`
(µM⁻¹s⁻¹) and off-rate `koff` (s⁻¹):

- **Critical concentration**: `Cc = koff/kon` — the monomer concentration at which net addition at
  *that end alone* is zero.
- **Net flux at any concentration C**: `J(C) = kon*C − koff`.
- **Treadmilling steady state** `C*`: solve `J_barbed(C*) + J_pointed(C*) = 0` (barbed-end growth
  exactly balances pointed-end loss, filament length constant, subunits flow barbed→pointed) —
  closed form: **`C* = (koff_b + koff_p) / (kon_b + kon_p)`**. This is *not* the average of the two
  Cc's; it is the point where the two ends' fluxes cancel, first predicted theoretically by
  `wegner_1976` a decade before Pollard measured the actual rate constants.

**The forced adversary, made explicit (§6 gives the honest epistemic weight)**: if the two ends
were kinetically identical (`kon_b=kon_p`, `koff_b=koff_p`), `C*` collapses algebraically to the
single shared Cc, and `J(C*) = kon*(koff/kon) − koff ≡ 0` — **by definition**, not by simulation.
The real, *falsifiable* content is not this identity — it is that **Pollard's measured rate
constants are NOT equal at the two ends** (an empirical fact that could have come out symmetric
and didn't), which is what pries `C*` open into a genuine window strictly between two *different*
Cc's (`0.1207 < 0.1705 < 0.6154` µM) where one end grows while the other shrinks. The algebra is
the reason *why* asymmetric rate constants are what generate a treadmilling window; it is not
itself new data.

---

## 2. Actin — barbed/pointed asymmetry and the treadmilling flux

**Pollard TD (1986). J Cell Biol 103(6 Pt 2):2747-2754. PMID 3793756, DOI
10.1083/jcb.103.6.2747, PMCID PMC2114620.** Full abstract verified verbatim, live, via two
independent fetch routes (Europe PMC `fullTextXML`, PMC HTML article page). EM elongation assay
using *Limulus* acrosomal-process nuclei, growth rates 0.05–280 s⁻¹; **elongation rate is linear in
ATP-actin concentration from the critical concentration to 20 µM at both ends** (confirms the
2-parameter birth-death model form is not an assumption — it is what this paper's own data
requires). Quoted directly: compared with ATP-actin, "ADP-actin associates slower at both ends,
dissociates faster from the barbed end, but dissociates slower from the pointed end," and at steady
state "there will be slow subunit flux from the barbed end to the pointed end" — treadmilling,
stated directly, with "the pointed end...relatively stable."

**Honest gap, stated up front**: the paper's own abstract text defers the actual numeric rate
constants to *"(table; see text)"* — and **both** live fetch routes this session returned only the
abstract/metadata section (~7–12 Kb), not the body Results section containing Table I (an
OCR/tagging gap in this 1986 scan, not a paywall — the article is fully open access). A further
live search for a secondary source explicitly re-quoting the table did not succeed. The values used
below (**kon_barbed=11.6, koff_barbed=1.4 µM⁻¹s⁻¹/s⁻¹; kon_pointed=1.3, koff_pointed=0.8**) are the
**textbook-tier** numbers universally attributed to this paper across the field, matching this
task's own pre-registered number — disclosed as textbook-tier, not primary-live-extracted, this
session.

| quantity | barbed end | pointed end | ratio (pointed/barbed) |
|---|---:|---:|---:|
| kon (µM⁻¹s⁻¹) | 11.6 | 1.3 | — |
| koff (s⁻¹) | 1.4 | 0.8 | — |
| **Cc = koff/kon (µM)** | **0.1207** | **0.6154** | **5.10×** |

Asymmetry confirmed against a pre-registered ≥2× bar (**PASS, 5.10×**, more than double the
threshold). **Task's own pre-registered framing, checked precisely, not smoothed to fit**: the
task states "the ~0.1-0.2 µM critical-concentration difference driving treadmilling." Checked both
readings — Cc_barbed **alone** (0.1207 µM) sits inside that range; the Cc_pointed−Cc_barbed
**difference** (0.4947 µM) does not. Both reported exactly (§0/JSON); the task's framing most
plausibly intended the barbed-end value itself.

**Treadmilling steady state**: `C* = (1.4+0.8)/(11.6+1.3) = 0.1705 µM` (strictly between the two
Cc's, as required). Net flux `J = 0.578 subunits/s` (barbed↔pointed balance residual < 1e-9,
confirming the solver). Converted via the actin helical-lattice geometry (2.7 nm axial rise/monomer,
**Holmes KC, Popp D, Gebhard W, Kabsch W (1990). Nature 347:44-49, PMID 2395461, DOI
10.1038/347044a0**, live-verified; textbook-tier constant, not independently re-extracted as an
explicit number from the live abstract): **treadmill rate = 0.0937 µm/min = 5.62 µm/h.**

**Independent-method over-determination (F1d)**: **Fujiwara I, Vavylonis D, Pollard TD (2007).
PNAS 104(26):8827-8832. PMID 17517656, DOI 10.1073/pnas.0702510104, PMCID PMC1885587.** Modern
single-filament fluorescence microscopy (not bulk EM), same senior author 21 years later,
explicitly investigating "the origins of the difference in the critical concentration at the two
ends of actin filaments in the presence of ATP" (direct quote). Measured: saturating phosphate
reduces Mg-ADP-actin Cc from **1.8 → 0.06 µM**, almost entirely via reduced off-rates at both ends;
barbed-end on-rate increases only 15% with phosphate, still 3-fold below ATP-actin's. Different
nucleotide states than Pollard's ATP-actin (disclosed, not a strict replication), but the same
0.06–1.8 µM order of magnitude as the ATP-actin Cc's above — same-order-of-magnitude check PASSES.

**A third independent single-filament method directly observes the flux-balance structure
itself**: **Fujiwara I, Takahashi S, Tadakuma H, Funatsu T, Ishiwata S (2002). Nat Cell Biol
4(9):666-673. PMID 12198494, DOI 10.1038/ncb841.** Quoted directly: "during the steady-state phase,
a treadmilling process of elongation at the barbed end and shortening at the pointed end occurs, in
which both components of the process proceed at **approximately the same rate**." This is a
directly-observed confirmation of the `J_barbed = −J_pointed` structure §1's `C*` derivation
predicts algebraically — validating the model's *shape*, not merely a fitted number.

---

## 3. Microtubule — dynamic instability, rate constants live-extracted from primary full text

**Mitchison T, Kirschner M (1984). Nature 312(5991):237-242. PMID 6504138, DOI 10.1038/312237a0.**
Full abstract verified verbatim (discovery paper): "microtubules in vitro coexist in growing and
shrinking populations which interconvert rather infrequently. This dynamic instability is a
general property of microtubules." Qualitative, as expected for a short Nature letter — the
rate-constant table is Walker 1988's contribution.

**Walker RA, O'Brien ET, Pryer NK, Soboeiro MF, Voter WA, Erickson HP, Salmon ED (1988). J Cell
Biol 107(4):1437-1448. PMID 3170635, DOI 10.1083/jcb.107.4.1437, PMCID PMC2115242.** Unlike
Pollard 1986, **the full body text WAS successfully fetched live this session** (open access,
Europe PMC `fullTextXML`) — every number below is a **direct quote from the primary Results
section**, not a textbook-tier value:

| quantity | plus end | minus end |
|---|---:|---:|
| kon (µM⁻¹s⁻¹) | 8.9 | 4.3 |
| koff (s⁻¹) | 44 | 23 |
| Cc = koff/kon (µM), derived here | **4.944** | **5.349** |
| Cc, paper's OWN stated value | "~5 µM" (both ends) | "~5 µM" (both ends) |
| rapid-shortening rate (dimers/s, conc.-independent) | 733 | 915 |

Porcine brain tubulin, MAP-free, 37°C, tubulin concentrations 7–15.5 µM tested. Quoted directly:
"as the tubulin concentration was increased, catastrophe frequency **decreased** at both ends, and
rescue frequency **increased dramatically** at the minus end...the frequency of catastrophe was
slightly greater at the plus end, and the frequency of rescue was greater at the minus end" — the
classic GTP-cap-size qualitative signature, confirmed directly from primary text (F2c). The paper
also quotes an earlier population-assay estimate (Gard and Kirschner: kon+=1.4 µM⁻¹s⁻¹, ~6× lower)
and **explains the discrepancy itself**: population assays underestimate rates because net growth
is diluted by nucleation lag and in-incubation catastrophes — a disclosed, mechanistically-explained
cross-method difference, not a contradiction.

**F2a — growth velocity, converted via MT helical-lattice geometry** (13 protofilaments × 8 nm
axial rise/dimer = 1625 dimers/µm, standard textbook constant):

| tubulin (µM) | flux (dimers/s) | growth velocity (µm/min) | in task's "several µm/min" band [1,10]? |
|---:|---:|---:|---|
| 7.0 (Walker's own lower bound) | 18.3 | 0.676 | no (below) |
| 10.0 | 45.0 | **1.662** | yes |
| 12.0 | 62.8 | **2.319** | yes |
| 15.5 (Walker's own upper bound) | 93.95 | **3.469** | yes |

3/4 points clear the pre-registered ≥50%-in-range bar (**PASS**) — reported exactly, including the
one point (7 µM, Walker's own lowest tested concentration) that falls just under the disclosed
[1,10] µm/min interpretation of "several," not hidden.

**F2b — over-determination, disclosed epistemic weight**: Cc derived here purely from Walker's own
kon/koff (1.1% / 7.0% relative error vs the paper's own stated "~5 µM") — confirms no
transcription/unit error entered the rate constants used, but both numbers originate in the *same*
paper; this is an internal-consistency check, not an independent replication (§6).

---

## 4. Decorrelated hydrolysis-necessity check (F3) — 2 empirical legs + 1 theoretical proof

| leg | system | citation | manipulation | effect (direct quote) |
|---|---|---|---|---|
| **GMPCPP** | microtubules | **Hyman AA, Salser S, Drechsel DN, Unger E, Mitchison TJ (1992). Mol Biol Cell 3(10):1155-1167. PMID 1421572** (recall-corrected, see banner), DOI 10.1091/mbc.3.10.1155, PMCID PMC275680 | slowly-hydrolyzable GTP analog replaces GTP | depolymerization rate collapses **5000×** (500→0.1 s⁻¹); "GMPCPP also **completely suppresses** dynamic instability"; GMPCPP's own hydrolysis rate is 4×10⁻⁷s⁻¹ (negligible); "GTP hydrolysis...is not required for normal polymerization but is essential for depolymerization and thus for dynamic instability" |
| **ATP-hydrolysis-dependent fluctuations** | actin | **Duttagupta M, Riley AT, Brieher WM (2025). PNAS 122. PMID 41264255, DOI 10.1073/pnas.2505077122, PMCID PMC12663992** | direct hydrolysis/Pi-release dependence test (disclosed substitute, §6) | "barbed-end fluctuations depended on adenosine triphosphate (ATP) hydrolysis and release of inorganic phosphate (Pi) and were blocked by phalloidin and barbed-end capping agents"; frequent ±2-6-subunit and rare ±50-150-subunit excursions; ~1/3 of steady-state ATP consumption goes to these fluctuations, 2/3 to treadmilling; interpreted as "a dampened form of dynamic instability in pure actin" |
| **Thermodynamic theory** | actin OR microtubules (explicit in the title) | **Hill TL (1980). PNAS 77(8):4803-4807. PMID 6933529, DOI 10.1073/pnas.77.8.4803** | closed-form limit: free energy of hydrolysis X→0 | "the monomer flux and the ATP flux can both be expressed in terms of X and rate constants...**both fluxes approach zero as X leads to 0**...this limit corresponds to ATP equilibrium" — a thermodynamic **proof**, not an empirical correlation |

**3/3 legs sign-consistent**: hydrolysis (not mere nucleotide *binding*) is mechanistically required
for net polymer dynamics, across 2 different polymers, 2 different empirical methods/eras, and 1
general theoretical argument that covers both by name. This is the same 3-leg-decorrelated-check
structure `docs/MECHANISM_DNA_REPAIR_KINETICS.md` §4 used, applied here across mechanism *types*
(2 empirical + 1 theoretical) rather than 3 empirical gene systems.

**Disclosed substitution**: the task's own pre-registration named "AMP-PNP" as the actin-side
non-hydrolyzable analog. Despite multiple targeted live searches this session (exact phrase,
alternate spellings, general nonhydrolyzable-analog terms), the specific historical AMP-PNP-on-actin
experiment was **not independently located**. `duttagupta_2025` is used as a mechanistically
equivalent, more modern substitute (direct test of hydrolysis-dependence via a different reagent
class) — disclosed here and in the evidence JSON, not a silent swap.

---

## 5. Cross-polymer distinctiveness (F2d) — why actin treadmills and pure-tubulin MT does not

The **same** `C*`/`J` machinery, applied to Walker 1988's own **near-symmetric** MT rate constants
(Cc_plus/Cc_minus ratio = **1.08×**, vs actin's **5.10×**), predicts a theoretical pure-tubulin
treadmilling flux of **1.174 dimers/s (2.60 µm/h)** — compare this to the **same paper's own**
catastrophic-shrinkage rate scale (**824 dimers/s**, average of 733/915): the treadmilling flux is
**701.7× smaller** than the shrinkage-event velocity scale, clearing a pre-registered ≥10×
"dynamic-instability-dominated" bar by a 70-fold margin. This is the **geometric** reason (not a
separate empirical fact) that pure-tubulin microtubules are dynamic-instability-dominated rather
than treadmilling-dominated, while actin — with its much larger Cc asymmetry and no comparably fast
catastrophic-shrinkage mode — treadmills smoothly instead. It is the *same* asymmetric-end-kinetics
mechanism in both cases; what differs is the **size of the asymmetry relative to the polymer's own
shortening-velocity scale**.

**Plausibility cross-check, disclosed weight**: the computed 2.60 µm/h pure-tubulin theoretical flux
is same-order-of-magnitude as **Horio T, Hotani H (1988). Cell Motil Cytoskeleton 11(4):229-236.
PMID 2972399, DOI 10.1002/cm.970100127** — a **directly measured** MAP-stabilized treadmilling flux
of **0.9 µm/h** (real-time darkfield observation of a dynein-decorated marker block). This is a
different (MAP-containing) system, not a same-system replication — disclosed as a plausibility
check, not a strict validation. Critically, the SAME paper found: "microtubules undergo treadmilling
and do NOT exhibit any dynamic instability" once MAP-stabilized — "treadmilling can take place in
the steady state **only after** microtubules have been stabilized by MAPs." Independent
confirmation of the dynamic-instability phenomenon itself, by a wholly different visualization
instrument (dark-field microscopy vs Walker's video-enhanced light microscopy), comes from **Horio
T, Hotani H (1986). Nature 321:605-607. PMID 3713844, DOI 10.1038/321605a0**: both ends "alternate
quite frequently between the two phases in a stochastic manner," with "no correlation...either among
individual microtubules or between the two ends of a single microtubule," and the two ends have
"remarkably different characteristics."

---

## 6. Symmetric QC — honest gaps and honest epistemic-tier accounting, held OPEN, not swept

**Which gates are genuinely empirical vs. closed-form model-correctness checks (stated plainly, not
buried)**:
- **Empirical/externally-anchored** (real, data-contingent findings that could have failed):
  F1a (Pollard's measured asymmetry ratio), F1d (Fujiwara 2007 independent method), F2a (Walker's
  measured growth velocities), F2c (Walker's measured catastrophe/rescue direction), all 3 F3 legs
  (Hyman 1992, Duttagupta 2025, Hill 1980).
- **Closed-form model-correctness / definitional checks** (necessary infrastructure, not
  independent data points): **F1b** (symmetric adversary flux ≡ 0 by the algebraic definition of
  Cc — see §1's honest framing), **F1c** (the C* solver satisfies its own defining equation to
  float precision — a code-correctness check), and half of **F2b** (Cc derived from Walker's own
  kon/koff vs Walker's own stated value — both numbers originate in one paper, an internal-
  consistency check, not an independent replication). These are legitimate and necessary (they
  confirm the model correctly encodes the claimed mechanism, catching a class of bug the DNA-repair
  precedent's md5-determinism check catches for a different failure mode) — but **10/10 PASS should
  not be read as "10 independent empirical confirmations."** Only 6 gates carry that weight.

**Other honest gaps**:
1. **Pollard 1986's own numeric Table I was NOT independently extracted from primary full text this
   session** — 2 separate live fetch routes (Europe PMC `fullTextXML`, PMC HTML article page) both
   returned only the abstract/metadata section; the abstract's own printed text defers to "(table;
   see text)"; a further search for a secondary source re-quoting the table did not succeed. The
   textbook-tier values used match the task's own pre-registered number and are universally
   attributed to this paper, but are disclosed as textbook-tier, not primary-live-extracted.
2. **IN-VITRO ≠ IN-VIVO (the task's own pre-registered caveat, held OPEN)**: every rate constant
   here is purified-protein kinetics in one specific buffer/temperature. In vivo, formins and the
   Arp2/3 complex nucleate/processively cap/elongate actin, profilin/cofilin gate monomer
   availability and severing, and +TIPs (EB1/3, XMAP215, stathmin, katanin, spastin) heavily
   modulate MT catastrophe/rescue/severing — none are in this model. This repo's own
   `MOL-CYTOSKELETON-DYNAMICS` node already shows cofilin knockdown alone nearly **doubles** max
   cell length (52.16→99.33 µm, P<0.0001, live-verified there) — regulators, not bare polymer
   kinetics, dominate cell-scale behavior in vivo.
3. **Dynamic-instability parameters are highly assay/tubulin-source dependent (task's own
   caveat)**: Walker 1988 used porcine brain tubulin, MAP-free, 37°C, one buffer; Horio & Hotani
   used calf brain tubulin under darkfield conditions — a partial cross-source check, not a
   systematic survey. No second independently-parametrized rate-constant table was live-verified.
4. **The 2-state GTP-cap model is a simplification** — Brouhard GJ (2015). Mol Biol Cell
   26(7):1207-1210. PMID 25823928, DOI 10.1091/mbc.e13-10-0594, PMCID PMC4454169 (the exact PMID
   the anchor-graph stub named, independently re-confirmed live this session) states directly: the
   canonical explanation "has been recently subverted, particularly...how GTP-tubulin forms
   polymers and why GTP hydrolysis disrupts them" — a primary source actively flagging the
   simplification, not merely this cell's own caveat. Its own contrast, quoted directly: "polymers
   such as F-actin will grow continuously as long as the subunit concentration is high enough,"
   unlike MT's sudden catastrophic shrinkage.
5. **Treadmilling and dynamic instability are not a universal, absolute dichotomy** — Horio & Hotani
   1988 found near-mutual-exclusivity specifically in their MAP-content-gated system; Duttagupta
   2025 shows pure actin *also* has a bounded (±50-150-subunit), dampened instability-like mode
   coexisting with treadmilling. §5's clean actin-vs-MT distinction is a first-order,
   literature-supported picture, not an absolute partition.
6. **The AMP-PNP actin experiment named in this task's own pre-registration was not independently
   located live this session** despite multiple targeted searches (§4) — `duttagupta_2025` is a
   disclosed substitute, not a silent swap.
7. **Wegner A (1976). J Mol Biol 108(1):139-150. PMID 1003481, DOI
   10.1016/s0022-2836(76)80100-3** — PMID/DOI/journal/pages verified live; abstract text itself
   unavailable (pre-abstract-indexing era) — historical/conceptual citation only, no quantitative
   content from it is used in any gate.
8. **Geometric conversion constants** (2.7 nm/actin-monomer; 13-protofilament/8 nm-dimer MT
   lattice) are standard structural-biology textbook facts (cited to Holmes 1990 for actin) — not
   independently re-extracted as explicit numbers from the live-fetched abstract this session.
9. **`TASK_GROWTH_RANGE_UM_MIN=(1,10)`** is this session's own disclosed interpretation of "several
   µm/min" — not a literature-quoted numeric range.
10. **Single buffer/temperature/species condition per paper** — no claim of universality across
    ionic strength, temperature, or tubulin/actin source.
11. **11 of 13 citations carry a full abstract or full-text live confirmation this session**;
    `wegner_1976`'s abstract was not retrievable (bibliographic metadata only); `holmes_1990`'s
    specific 2.7 nm number was not independently extracted from its own abstract (textbook-tier).

None of the above are swept into `overall_pass`; every gate that composes it is disclosed above
with its true epistemic weight, matching this repo's symmetric-QC convention.

---

## 7. Pre-registered gates — machine-printed, not narrated

```
F1.F1a_asymmetry_confirmed                               PASS  (5.10x >= 2.0x gate)
F1.F1b_forced_symmetric_adversary_falls                  PASS  (flux=0.0 exactly; model-correctness, Sec 6)
F1.F1c_flux_balance_selfconsistent                       PASS  (residual < 1e-9; solver-correctness, Sec 6)
F1.F1d_fujiwara_independent_method_same_order            PASS
F2.F2a_growth_velocity_in_task_range                     PASS  (3/4 points, 75% >= 50% gate)
F2.F2b_cc_overdetermination_matches                      PASS  (1.1%/7.0% relerr <= 15% tol; Sec 6)
F2.F2c_catastrophe_rescue_direction_literature_confirmed PASS  (live-quoted, qualitative)
F2.F2d_cross_polymer_di_dominance_confirmed              PASS  (701.7x >= 10.0x gate)
F3.F3_all_legs_sign_consistent                           PASS  (3/3)
F3.F3_all_legs_fall_confirming_necessity                 PASS  (3/3; GMPCPP 5000x suppression)

OVERALL: PASS (10/10) -- 6 empirical/externally-anchored, 3 model-correctness/definitional,
1 (F2b) split between the two -- see Sec 6 for the honest breakdown.
```
Deterministic: byte-identical md5 `f920a9b984e9d9aefa788619c94e64b7` confirmed across 2 independent
process runs (no RNG — closed-form algebra + a small literature parameter grid only). All floats
in the evidence JSON independently re-verified finite by parsing the written file in a fresh
process (not by trusting the script's own internal check).

---

## 8. couples_to — the disclosed relationship to the rest of the graph

Per this repo's convention (`docs/MECHANISM_HARDENED_CONVENTIONS.md` §4b), this cell **couples_to**:

- **`AUTO-VERIFY-CYTOSKELETON-NUMERIC-TABLES-FULL`** — this cell **executes** that stub's
  previously-unexecuted targets (Pollard 1986 rate-constant table + PMID:25823928). Not itself
  folded into the graph this session (§ below).
- **`MOL-CYTOSKELETON-DYNAMICS`** — the cell-mechanical/malignancy axis this cell's polymer-kinetics
  axis sits underneath (cofilin/Rho/Rac-driven cell shape and stiffness are downstream consequences
  of the same actin dynamics modeled here, at a coarser, cell-scale readout). Reused one datapoint
  by reference (cofilin knockdown cell-length effect, §6 item 2) as disclosed, not re-derived.
- **`cell-biology (cell division/motility/transport)`** (task's own requested coupling): microtubule
  dynamic instability is the literal mechanistic basis for mitotic-spindle search-and-capture
  (Mitchison & Kirschner's own original motivating application, not built as a separate cell this
  session — a prospective pointer, matching this repo's established disclosure convention for
  not-yet-instantiated coupling targets); intracellular transport (kinesin/dynein walking along
  microtubule tracks) depends on the *stable* (non-catastrophizing) subset of the same MT population
  this cell models as dynamically unstable — the two behaviors coexist on different MT
  subpopulations in the same cell (stabilized vs dynamic), not a contradiction. Cell migration/
  motility couples through `MOL-CYTOSKELETON-DYNAMICS`'s own already-populated legs (Latrunculin-B/
  microfluidic-transit datapoints, live-verified there).
- **`muscle (actin overlaps the contractile apparatus, distinct pool)`** (task's own requested
  coupling, checked live this session): `grep`-checked `docs/MECHANISM_CROSSBRIDGE_MODEL.md` and
  `docs/MECHANISM_CONTRACTION_DYNAMICS.md` for tropomyosin/troponin/CapZ/nebulin/tropomodulin —
  **no hits**; neither doc currently addresses thin-filament (sarcomeric actin) STABILITY
  mechanisms, so this is a genuinely open disambiguation, not a contradiction to reconcile. The
  actin this cell models (cytoplasmic/cortical, ATP-driven barbed/pointed treadmilling, Cc~0.1-0.6
  µM) is mechanistically DISTINCT from sarcomeric thin-filament actin, which is capped
  (CapZ/tropomodulin), side-bound (tropomyosin/troponin), and length-set (nebulin) specifically to
  **suppress** treadmilling and hold a fixed length for the crossbridge cycle those docs model —
  same monomer chemistry, opposite regulatory outcome (dynamic vs static), a real and useful
  contrast, not yet built as an executed cell.

Per `docs/MECHANISM_HARDENED_CONVENTIONS.md` §2/§4, promoting any of this into canonical
`data/MECHANISM_ANCHOR_GRAPH.json` nodes/edges requires the separate `mechanism_fold → fold_gate_v2`
pipeline — **not performed this session**, matching the established precedent of every sibling
`MECHANISM_*` doc (DNA-repair kinetics, wound-healing, thyroid axis) stating its results in
prose/evidence-JSON only.

---

## 9. Files

- `source_repository/scripts/msk/cytoskeleton_dynamics.py` — the model: `CITATIONS`
  dict (13 entries, tiered live-verified, each with PMID+DOI+exact quoted finding + disclosed
  tier), the shared 2-end birth-death kinetic machinery (`critical_concentration`, `end_flux`,
  `treadmill_c_star`, `symmetric_adversary`), F1 (actin, §2), F2 (microtubule, §3), F3 (decorrelated
  hydrolysis-necessity, §4), the cross-polymer distinctiveness computation (§5), the symmetric-QC
  list (§6, including the honest empirical-vs-definitional gate accounting), the gates block, and
  the evidence-JSON writer with a finite-value check. Run with `source .venv-msk/bin/activate &&
  python3 scripts/msk/cytoskeleton_dynamics.py` (<1s wall time, pure numpy, no OpenSim dependency,
  deterministic — no RNG, confirmed byte-identical across 2 runs).
- `source_repository/data/msk_smoketest/cytoskeleton_dynamics/cytoskeleton_dynamics_results.json`
  — full machine-written evidence (all 13 citations, F1/F2/F3 computed numbers, geometric
  constants, pre-registered thresholds, symmetric QC, gates, `overall_pass`). md5
  `f920a9b984e9d9aefa788619c94e64b7`, confirmed byte-identical across 2 independent runs; all
  floats independently re-verified finite in a fresh process (0 non-finite values found).
- Read but NOT modified (isolation: touch only files created this session):
  `source_repository/data/MECHANISM_ANCHOR_GRAPH.json` (the
  `AUTO-VERIFY-CYTOSKELETON-NUMERIC-TABLES-FULL` stub this cell executes, and the
  `MOL-CYTOSKELETON-DYNAMICS` sibling node this cell complements — both confirmed distinct,
  §0/§8), `source_documents/MECHANISM_HARDENED_CONVENTIONS.md` (fold path/node
  schema/doc conventions — this doc is a pre-fold HYPOTHESIS artifact, §8),
  `source_documents/MECHANISM_DNA_REPAIR_KINETICS.md` (read for doc-structure/
  epistemic-tier-disclosure convention precedent, not re-litigated),
  `source_documents/MECHANISM_CROSSBRIDGE_MODEL.md`,
  `source_documents/MECHANISM_CONTRACTION_DYNAMICS.md` (grep-checked for the
  muscle/distinct-pool coupling, §8; not modified), `source_repository/COORDINATOR.md`
  (isolation/startup conventions).
