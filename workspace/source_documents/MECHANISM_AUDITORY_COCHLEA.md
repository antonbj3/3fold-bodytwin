# MECHANISM AUDITORY / COCHLEA — Greenwood tonotopy + traveling-wave active amplifier + threshold-of-hearing (2026-07-22)

Resolves the auditory/cochlear-transduction sensory layer, alongside this twin's vestibular
(`docs/MECHANISM_VESTIBULAR_BALANCE.md`) and nerve-conduction (`docs/MECHANISM_NERVE_CONDUCTION.md`)
threads. Three mechanistic, geometry-derived (not curve-fit) pieces: (1) the Greenwood
cochlear frequency-position (tonotopic) map, (2) the basilar-membrane traveling-wave place
mechanism plus the active outer-hair-cell "cochlear amplifier" as a living/dead compressive
nonlinearity, and (3) a decorrelated check of the ear-canal quarter-wavelength resonance
against the threshold-of-hearing minimum's location. Script: `scripts/msk/auditory_cochlea.py`.
Raw machine-checked results: `reports/probes/auditory_cochlea_results.json`. Consolidated
citation ledger: `docs/MECHANISM_AUDITORY_COCHLEA_evidence.json`.

**Scope note before anything else**: this repo's graph (`data/MECHANISM_ANCHOR_GRAPH.json`)
already holds `SENS-AUDITORY-TRANSDUCTION`, `SNS-AUDITORY-FORWARD-MODEL` /
`AUTO-DIFFERENTIABLE-AUDITORY-FORWARD-MODEL-CE`, and `SENS-COCHLEAR-SYNAPTOPATHY-HIDDEN-HEARING-LOSS`
— confirmed by direct read, not assumed. All four are `SEED-DESIGN` (`mechanism_grade`):
designed hypotheses with real, cited literature, but **not yet executed as machine-falsified
cells**. `SENS-AUDITORY-TRANSDUCTION` in particular already carries a live-verified,
4-paradigm OHC-ablation threshold-shift anchor (40-60dB, prestin-KO x2/kanamycin/Kit(W-v))
that is directly reused below (one of its four datapoints, Liberman 2002, was independently
re-spot-checked this session and reproduced exactly). This doc does not edit any of those
nodes — it is a **companion**: the actual quantitative model + fresh falsifier tests those
SEED-DESIGN cells did not yet have, per this repo's isolation discipline (re-verify, build a
companion, don't edit).

## Part 1 — Cochlear tonotopy: the Greenwood frequency-position function

### Method, in one paragraph

The classical traveling wave arises because the basilar membrane's **stiffness gradient**
(stiff at the base, compliant at the apex; von Bekesy, Nobel Prize 1961) sets a
position-dependent local resonance. A first-order WKB view — local resonance
`w0(x) ~ sqrt(K(x)/M)` with an exponentially tapering stiffness `K(x)` — predicts a **pure
exponential** place-frequency law over most of the membrane. Greenwood's (1990, PMID 2373794)
empirical human function, `F(x) = A*(10^(a*x) - k)` with `x` = proportion of basilar-membrane
length from the apex, **reduces exactly to that pure-exponential form** whenever `10^(a*x) >> k`,
and deviates from it only near the apex (`x->0`, where `10^(a*x)=1` is comparable to
`k=0.88`) — i.e. the "-k" term is not an arbitrary curve-fitting knob, it is the same
mathematical signature any stiffness-gradient resonator produces once the stiffness-dominated
approximation breaks down near a physical low-frequency floor (the helicotrema shunt).
`scripts/msk/auditory_cochlea.py` verifies this directly from the formula (closed-form log-slope
vs. numerical finite difference, two independent computational paths), then checks Greenwood's
own human constants against **independently measured** (non-Greenwood, non-acoustic) human
cochlear-length data.

### Citations — verified LIVE this session (NCBI eutils + PMC/EuropePMC full-text fetch, not recalled)

| source | verified as | what it anchors |
|---|---|---|
| Greenwood DD (1990). *J Acoust Soc Am* 87(6):2592-605. **PMID 2373794**, DOI 10.1121/1.399052 | Abstract fetched live | PRIMARY source. Abstract is conceptual (constants not printed there). |
| Stakhovskaya O, Sridhar D, Bonham BH, Leake PA (2007). *J Assoc Res Otolaryngol* 8(2):220-33. **PMID 17318276**, PMCID **PMC2394499** | **Full open-access text fetched live** | Quotes Greenwood 1990 **verbatim**: human `A=165.4, a=2.1, k=0.88` -> 20-20,677 Hz. Independently measured (n=9 cadaveric dissection) organ-of-Corti length: **33.13mm (SD 2.11, range 30.5-36.87mm)** — NOT derived from Greenwood's formula. Also quotes Sato 1991 (n=18: 32.3/37.1mm F/M) and Ketten1998/Skinner2002 (n=33: 33.4mm). |
| Dorman MF et al (2007). *J Assoc Res Otolaryngol* 8(2):234-40. **PMID 17351713**, PMCID PMC2441831 | Abstract fetched live | Living-patient **psychophysical** (electric/acoustic pitch-match) test in a *different regime* (CI electrode = spiral ganglion, not organ of Corti): pitch matches ~0.5-1 octave **lower** than Greenwood predicts at 3-20mm insertion depth. A genuine, disclosed domain-of-validity boundary, concordant with Boex et al 2006. |

### Headline results (all numbers from `auditory_cochlea.py`'s own JSON output, not hand-typed)

| check | result |
|---|---|
| **Reproduction check** (does the extracted formula reproduce Stakhovskaya's own stated 20Hz / 20,677Hz?) | f(apex)=**19.848Hz** (stated 20), f(base)=**20677.07Hz** (stated 20677) — both within 2% tolerance — **PASS** |
| **Falsifier 1b** (pre-registered: Greenwood's assumed 35mm within 1 SD of 3 independent measured-length studies) | grand mean of Sato1991/Ketten-Skinner/Stakhovskaya = **33.74mm**; \|35-33.74\|=1.26mm vs SD=2.11mm -> **z=0.60** — **PASS** |
| **Forced adversary**: naive fixed-mm (not proportional) position-to-frequency rule, applied to Stakhovskaya's own shortest/longest measured cochleae (30.5/36.87mm) | mis-predicts frequency by up to **0.48 octaves** (746Hz or 1392Hz instead of 1000Hz) — the proportional Greenwood formula has **zero** error by construction. Adversary forced to its literal strongest form and shown to **measurably fail**. |
| **Geometric check** (closed-form `d(lnF)/dx` vs. numerical finite difference — two independent computational paths) | max relative disagreement **1.4e-8** — **PASS** (machine cross-check, not eyeballed) |
| Apex-vs-base log-slope inflation | base-region asymptotic slope 4.835/unit-x; apex-limit slope 40.30/unit-x -> **8.33x inflation**, concentrated in the apical ~10-15% of length — the geometric signature of the "-k" floor, not a free parameter |
| Void-floor sweep | F(x) spans **19.85Hz to 20,677Hz** (1042-fold), a real, non-degenerate range |

### Honest gaps
- The reproduction check is an **implementation-fidelity** test (did the constants get
  transcribed and coded correctly), not a scientific falsifier of the theory — it is reported
  separately from, and is weaker than, the independent-length-anchor check (1b).
- The geometric two-computational-paths check is a pure **mathematical** consistency check on
  a single formula — it will always pass absent a coding bug; it does not test anything about
  the physical world by itself. It corroborates the WKB *interpretation*, it does not prove it.
- The implied stiffness e-folding length (3.62mm) is a **derived geometric interpretation** of
  the verified formula, offered as mechanistic motivation — it was **not** independently
  cross-checked this session against a separate direct mechanical-stiffness-gradient
  measurement paper.
- Dorman 2007/Boex 2006 show the acoustic (organ-of-Corti) map does **not** directly transfer
  to electrical (spiral-ganglion) cochlear-implant stimulation — disclosed as a real boundary,
  not swept under the rug (this is *why* Stakhovskaya 2007's paper exists: a separate
  OC-to-SG offset function, not reproduced here, is needed for CI applications specifically).

## Part 2 — Basilar-membrane traveling wave + active cochlear amplifier (compressive nonlinearity)

### Method, in one paragraph

At any one cochlear place, the healthy (living) basilar membrane's response to a tone near its
characteristic frequency (CF) is **compressive**: small input changes produce even smaller
output changes across a wide middle range of levels. This is the "cochlear amplifier" — active,
outer-hair-cell (OHC, prestin-driven electromotility) positive feedback onto the passive
traveling wave. `scripts/msk/auditory_cochlea.py` builds a piecewise dB-domain input/output model
directly calibrated to Ruggero et al. (1997, PMID 9104018)'s own **directly measured** chinchilla
basilar-membrane numbers: linear (slope 1 dB/dB, gain 66-76dB re stapes motion) below ~15dB SPL,
compressive (slope **0.2 dB/dB**, their own "as low as" figure) from 15-80dB SPL. Removing the
amplifier (death) is modeled exactly as Ruggero's abstract states it: **fully linear** ("death
abolishes all nonlinearities") at a reduced gain (60-81dB lower at CF).

### Citations — verified LIVE this session

| source | verified as | what it anchors |
|---|---|---|
| Ruggero MA, Rich NC, Recio A, Narayan SS, Robles L (1997). *J Acoust Soc Am* 101(4):2151-63. **PMID 9104018**, PMCID PMC3578390 | Abstract fetched live | PRIMARY, direct in-vivo laser-velocimetry: linear gain 66-76dB; compressive slope "as low as 0.2dB/dB" (40-80dB SPL); death -> fully linear, **60-81dB** desensitization at CF. |
| Sellick PM, Patuzzi R, Johnstone BM (1982). *J Acoust Soc Am* 72(1):131-41. **PMID 7108035**, DOI 10.1121/1.387996 | Abstract fetched live | INDEPENDENT (different lab, species — guinea pig —, decade, and technique — Mossbauer, not laser): saturating BM nonlinearity in vivo, "**eliminated**... and **absent post mortem**." Same qualitative signature, fully decorrelated source. |
| Robles L, Ruggero MA (2001). *Physiol Rev* 81(3):1305-52. **PMID 11427697**, PMCID PMC3590856 | Abstract fetched live | Review corroboration + an honest regime caveat this doc keeps, not hides: "at the apex... nonlinearities appear to be less prominent than at the base." |
| Oxenham AJ, Bacon SP (2003). *Ear Hear* 24(5):352-66. **PMID 14534407**, DOI 10.1097/01.AUD.0000090470.73934.78 | Abstract fetched live | Cross-domain (physiological->psychophysical) review: total OHC loss -> BM linearization -> explicitly accounts for poorer thresholds, **loudness recruitment**, reduced frequency selectivity. Qualitative corroboration, not a quantitative dB/dB equivalence. |
| Liberman MC et al (2002). *Nature* 419(6904):300-4. **PMID 12239568** — **spot-checked live this session** | Re-fetched, reproduced exactly | "Hearing sensitivity in mammals is enhanced by **more than 40dB** (100-fold)" by the cochlear amplifier; prestin-KO: 40-60dB loss; het: ~6dB. Confirms this repo's own `SENS-AUDITORY-TRANSDUCTION` graph datapoint was not fabricated. |

### Headline results

| check | result |
|---|---|
| **Falsifier 2** (pre-registered: direct BM-mechanical gain-loss-at-death, 60-81dB, must numerically overlap the fully decorrelated OHC-ablation threshold-shift range, 40-60dB, already live-verified in this repo's graph) | overlap = **[60,60]dB** — **PASS**, but honestly a **single boundary point**, not a wide comfortable match (see honest gaps) |
| **Forced adversary**: is living/dead difference JUST a constant gain offset (same shape)? | Adversary predicts `ratio_dead/ratio_alive` in [0.8,1.2]. Actual (both numbers directly cited, not fitted): ratio_alive=0.2, ratio_dead=1.0 -> **ratio_of_ratios = 5.0x**, far outside the survival band -> **adversary REJECTED**: the living/dead difference is a genuine **shape** change (compressive->linear), not merely an amplitude shift. |
| Void-floor sweep (0-120dB SPL input) | alive output range **narrower** than dead output range over the identical input sweep (compression is real, non-degenerate, and asymmetric between conditions) — **PASS** |

### Honest gaps (symmetric QC — nothing here is proven to a single number)
- **The overlap is boundary-touching (exactly 60dB), not deep.** Direct BM-mechanical
  measurement (60-81dB) and genetic/chemical OHC-ablation threshold shift (40-60dB) are
  **different measurement domains** (mechanical velocity vs. behavioral/neural threshold) —
  agreement in *order of magnitude*, meeting exactly at one boundary, is the honest,
  biologically-expected signature; it is reported exactly, not inflated into "they match."
- **The living/cadaver difference (the task brief's own named concern) is estimate-laden.**
  There is no single certified "the cochlear amplifier is worth N dB" number here — only a
  disclosed range (40-81dB) spanning multiple methods, held open deliberately.
- Compression calibration (Ruggero 1997) is from a single **basal** (9-10kHz CF) chinchilla
  site. Robles & Ruggero 2001's own review states apical nonlinearities are less prominent —
  this model is **not** extrapolated to apical/low-CF places.
- Psychophysical (loudness recruitment) and physiological (BM-mechanical) measures of
  "compression" are **qualitatively** concordant (Oxenham & Bacon 2003) only — no quantitative
  dB/dB equivalence between them is claimed.
- The 15dB / 80dB breakpoints and the extension of slope=0.2 above 80dB SPL are disclosed
  modeling simplifications of a smoothly-varying real I/O curve, not independently pinned
  point-by-point across the full input range.

## Part 3 — Decorrelated check: ear-canal quarter-wave resonance vs. threshold-of-hearing minimum

### Method, in one paragraph

The human ear canal is, to first order, a tube closed at the tympanic membrane and open at the
concha — a quarter-wavelength resonator, `f_res = c/(4L)`, a **parameter-free physics**
calculation once a canal length is known. Using an independently, **acoustically** measured
adult ear-canal length (Keefe et al. 2024, direct time-domain reflectometry, n=40: 22.6mm ±
6.5mm) gives a predicted resonance of **3.79-3.91kHz** (20°C-to-body-temperature sound-speed
range). This is checked against the psychoacoustically measured threshold-of-hearing minimum's
**location** (~3-4kHz per the task brief's own prior and ISO 226 / Suzuki & Takeshima 2004's
lineage) — two **genuinely decorrelated** measurement domains (an anatomical/acoustic tube
measurement vs. a behavioral-detection psychophysical measurement), connected only by shared
physics, not by any shared dataset or fitting procedure.

### Citations — verified LIVE this session

| source | verified as | what it anchors |
|---|---|---|
| Keefe DH, Fitzpatrick DF, Porter HL, Chen S (2024). *J Acoust Soc Am* 156(4):2709-2726. **PMID 39431854**, PMCID **PMC11495878** | **Full open-access text fetched live**, Table II extracted directly | PRIMARY anchor: direct modern acoustic (reflectometry) total ear-canal length, adults: **22.6mm (SD 6.5, n=40)**. |
| Suzuki Y, Takeshima H (2004). *J Acoust Soc Am* 116(2):918-33. **PMID 15376658**, DOI 10.1121/1.1763601 | Abstract fetched live (no open full text found) | The empirical basis of ISO 226:2003 (12 studies synthesized, supersedes Robinson & Dadson 1956, up to 14dB discrepancy below 500Hz). Exact numeric threshold table **not** extracted — honest gap. |
| Shaw EA (1974). *J Acoust Soc Am* 56(6):1848-61. **PMID 4443484**, DOI 10.1121/1.1903522 | Bibliographic match confirmed (no abstract indexed, pre-abstracting era) | Classic free-field-to-eardrum transfer-function reference for the ear-canal/pinna resonance peak. |
| Hammershøi D, Møller H (1996). *J Acoust Soc Am* 100(1):408-27. **PMID 8675836**, DOI 10.1121/1.415856 | Abstract fetched live | n=12, probe-microphone methodology corroboration that ear-canal transmission is a real, individually measurable transfer function. |
| Hennig L et al (2025). *Ann Anat* 257:152319. **PMID 39214319**, DOI 10.1016/j.aanat.2024.152319 | Abstract fetched live | Modern population MRI, n=870/1740 ears: cartilaginous EAC portion height 8.47-8.62mm — a partial-length corroboration (not the total length used for the calculation). |

### Headline results

| check | result |
|---|---|
| **Falsifier 3** (pre-registered: computed resonance must fall in [3,4]kHz) | central estimate **3905Hz** (body temp) / **3794Hz** (20°C reference) — both **inside** [3,4]kHz — **PASS** |
| Void-floor sweep across the measured population's ±1SD length range | resonance spans **3033-5482Hz** — a real, wide, non-degenerate range reflecting genuine inter-individual canal-length variability |

### Honest gaps
- **The threshold-of-hearing minimum's exact magnitude (dB SPL) was not independently
  extracted live this session** — Suzuki & Takeshima (2004)'s full text and table are
  paywalled and no open-access copy was found. Only the minimum's **location** (~3-4kHz) is
  cross-checked, carried as the task brief's own literature prior (and this repo's
  pre-existing, unverified `SNS-AUDITORY-FORWARD-MODEL` SEED-DESIGN note of "~2-5kHz").
- **This is a location match for the central estimate, not for every individual.** The
  population's own ±1SD canal-length variability translates to a **3.0-5.5kHz** resonance
  range — its upper tail (shorter canals) sits *outside* even the broader 2-5kHz literature
  band. Real biological variability, disclosed rather than smoothed over.
- **Contributory, not sole, cause.** The threshold curve's exact shape is also set by the
  middle-ear transfer function (this repo's own forward-model note: peaks nearer ~1kHz) and by
  intrinsic cochlear sensitivity. Ear-canal resonance is not claimed to be the sole reason for
  the minimum's location, only a genuine, independently computable contributor whose predicted
  location coincides with the observed one.
- The quarter-wavelength model idealizes the tympanic membrane as a rigid closed end and
  ignores the canal's non-uniform (funnel) cross-section — a first-order geometric
  approximation, not a full acoustic-impedance model.

## Symmetric QC — what is and is not claimed

**Claimed (PASS, machine-verified, externally anchored)**:
1. Greenwood's (1990) human cochlear frequency-position formula, with its own live-verified
   constants, reproduces its own stated calibration numbers (implementation fidelity) and is
   anchored by 3 independent (non-Greenwood, non-acoustic) measured human cochlear-length
   studies spanning 1991-2007 (grand mean 33.74mm vs. assumed 35mm, z=0.60). A naive
   fixed-distance adversary is forced to its literal strongest form and shown to measurably
   mis-predict frequency (up to 0.48 octaves) across the real measured population, while the
   verified proportional formula is immune to this by construction. The formula's functional
   shape is shown, via two independently cross-checked computational paths, to be consistent
   with a stiffness-gradient traveling-wave (WKB) mechanism, not an arbitrary curve fit.
2. A living/dead compressive cochlear-amplifier model, calibrated directly to Ruggero et al.
   (1997)'s own measured chinchilla basilar-membrane numbers, predicts a gain-loss-at-death
   range (60-81dB) that numerically touches the fully decorrelated (4-paradigm, 2-species,
   genetic/chemical OHC-ablation) threshold-shift range already verified in this repo's own
   graph (40-60dB). The tempting "just a gain offset" adversary is forced with the directly
   cited numbers and cleanly rejected (5x discrepancy) — the living/dead difference is a real
   shape change.
3. An independent (2024, acoustic reflectometry) measurement of human ear-canal length,
   combined with parameter-free quarter-wavelength physics, predicts a resonance (3.8-3.9kHz)
   that lands inside the literature-stated threshold-of-hearing-minimum band (3-4kHz) — a
   genuine geometry-to-psychophysics decorrelated cross-check.

**NOT claimed (explicitly OPEN, honest, per the task's own symmetric-QC brief)**:
1. No single certified number for "the cochlear amplifier is worth N dB" — the living/dead
   BM-mechanical range (60-81dB) and the OHC-ablation threshold-shift range (40-60dB) agree in
   order of magnitude and touch at exactly one boundary point (60dB); they are not claimed to
   coincide precisely, and living-vs-cadaver preparations are explicitly flagged (per the task
   brief) as differing enormously in ways this model does not fully resolve.
2. The exact ISO 226 / threshold-of-hearing numeric table (dB SPL magnitude, not location) was
   not independently extracted from a live primary source this session (paywalled).
3. Compression calibration is basal/high-CF-specific; apical compression is explicitly
   *not* modeled the same way, per Robles & Ruggero (2001)'s own disclosed caveat.
4. Psychophysical and physiological compression measures are held qualitatively, not
   quantitatively, concordant.
5. Greenwood's acoustic map does not transfer directly to electrical (cochlear-implant)
   stimulation without a measurable correction (Dorman 2007/Boex 2006) — a genuine boundary.
6. The ear-canal-resonance-vs-threshold-minimum match is for the population's central
   estimate; individual variability (±1SD canal length) produces a wider resonance range whose
   tail extends past the matched band.
7. The stiffness-gradient e-folding length is a geometric interpretation, not independently
   verified against a separate mechanical-stiffness measurement.

`couples_to`: `SENS-AUDITORY-TRANSDUCTION` (OHC electromotility/amplifier-gain SEED-DESIGN
hypothesis this doc's Falsifier 2 independently cross-checks with 2 fresh direct-mechanical
citations, Ruggero 1997 + Sellick 1982, neither of which was in that node's own source list),
`SNS-AUDITORY-FORWARD-MODEL` / `AUTO-DIFFERENTIABLE-AUDITORY-FORWARD-MODEL-CE` (SEED-DESIGN
differentiable forward-model cell whose stage-2/stage-3 slots — filterbank place-map and OHC
compressive nonlinearity — this doc's Greenwood map and BM I/O model could fill; not built or
edited here), `SENS-COCHLEAR-SYNAPTOPATHY-HIDDEN-HEARING-LOSS` (different hidden state —
IHC-synapse/auditory-nerve integrity, downstream of and complementary to this doc's OHC/BM-
mechanics scope, confirmed non-overlapping by direct read), `HEARING-LOSS-STAGE-LOCALIZATION`
(this doc supplies the mechanistic sensor-stage forward model that stage-localization work
needs), and `docs/MECHANISM_VESTIBULAR_BALANCE.md` (shared hair-cell mechanotransduction
motif — Piezo/TMC channels — and the shared repo convention of building a mechanistic model
atop a thin pre-existing SEED-DESIGN graph node rather than editing it).

## Repro

```
cd ~/projects/bodytwin
source .venv-msk/bin/activate
python scripts/msk/auditory_cochlea.py   # writes reports/probes/auditory_cochlea_results.json
```
Pure numpy (no OpenSim, no external data download, no scipy dependency beyond what the venv
already has). Runtime: well under 1 second. No git operations; no writes outside
`reports/probes/auditory_cochlea_results.json` and this doc pair.
