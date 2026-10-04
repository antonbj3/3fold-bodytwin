# MECHANISM MUSCLE ENERGETICS — the ATP/PCr/oxidative-flux layer (2026-07-22)

Builds the twin's first **CERTIFIED high-energy-phosphate metabolism model**: resting
[PCr]/[ATP]/[Pi], the creatine-kinase (CK) near-equilibrium constraint that buffers [ATP] while
[PCr] does the work, and — the key measurable — the post-exercise PCr recovery time constant
**tau** as a 31P-MRS in-vivo readout of mitochondrial oxidative capacity. Couples the ATP-demand
side to this twin's own already-measured muscle-force/metabolic work (`metabolic_cost.py`, the
Static-Optimization thread) and the O2-supply side to the already-built VO2/Fick chain
(`respiratory.py`, `cardiac_output.py`, `muscle_perfusion.py` — **no `MECHANISM_BLOOD_OXYGEN_
TRANSPORT.md` exists in this repo**; that role is currently played by these three docs together,
checked directly by `find`/`grep` before writing this, not assumed). Script:
`scripts/msk/muscle_energetics.py`. Evidence:
`data/msk_smoketest/subject2_walking1/muscle_energetics/muscle_energetics_results.json`.

**No re-solve, no OpenSim.** Section C reads `metabolic_cost_results.json` and
`respiratory_results.json` read-only and does pure arithmetic; Sections A/B are self-contained
closed-form algebra (numpy not even required) against live-verified literature constants.

## Pre-registered falsifier (stated in the task before this script was written)

Does modeled PCr-recovery kinetics reproduce the **measured** 31P-MRS tau (approx. 30-50s
healthy muscle; longer in deconditioned/diseased), and does the implied oxidative ATP flux,
converted to O2 via the P/O ratio, connect to the whole-body VO2/blood-O2 Fick chain? Symmetric
adversary, forced not assumed: a **deliberately wrong-muscle-group** Qmax (forearm, much less
oxidative than leg muscle) must produce a **wrong, out-of-band** tau — proving the model is
falsifiable, not a tautology that "works" for any input.

**Result: PASS, 13/13 machine gates, non-knife-edge on every one** (margins reported below, not
just pass/fail). One real bug was caught and fixed by the model's own self-check gate during
this build (a micromolar-to-molar unit error, `1e-9` where `1e-6` was needed) — reported here
rather than silently corrected, per this project's own machine-cross-check discipline.

## Citations — every PMID verified LIVE this session via NCBI eutils `efetch`/`esummary`

`WebSearch` was session-quota-exhausted (2000/2000 calls used) before this task started — the
same documented constraint `MECHANISM_CROSSBRIDGE_MODEL.md`/`MECHANISM_RESPIRATORY.md`/
`MECHANISM_CARDIAC.md` already disclose hitting. Fallback: direct `curl` to
`eutils.ncbi.nlm.nih.gov` (`esearch` -> `esummary`/`efetch`, `rettype=abstract`), the same
repo-precedented route. Every number below is quoted from a live-fetched abstract this session,
not recalled.

| # | Citation | PMID / DOI | Role |
|---|---|---|---|
| 1 | Lawson JW, Veech RL (1979). Effects of pH and free Mg2+ on the Keq of the creatine kinase reaction and other phosphate hydrolyses and phosphate transfer reactions. *J Biol Chem* 254(14):6528-37. | **36398** (no DOI on record — pre-DOI-era JBC, confirmed via esummary's own empty `articleids` doi field) | **Primary CK equilibrium constant**, quoted directly: K_CK = [ATP][Cr]/([ADP][PCr][H+]) = 3.78e8 M^-1 (free [Mg2+]=0), **1.66e9 M^-1 (free [Mg2+]=1mM, physiological)**, pH 7.0, 38C. |
| 2 | Veech RL, Lawson JW, Cornell NW, Krebs HA (1979). Cytosolic phosphorylation potential. *J Biol Chem* 254(14):6538-47. | **36399** (no DOI on record, same pre-DOI issue) | Companion paper — confirms CK (+adenylate kinase) near-equilibrium in **intact** brain/muscle tissue; states free cytosolic [ADP] is **"probably 20-fold lower than measured cell ADP content"** — the primary-source statement that free ADP sits below the direct-measurement floor. |
| 3 | Kemp GJ, Meyerspeer M, Moser E (2007). Absolute quantification of phosphorus metabolite concentrations in human muscle in vivo by 31P MRS: a quantitative review. *NMR Biomed* 20(6):555-65. | **17628042**, `10.1002/nbm.1192` | **Primary resting-concentration anchor** (their own new data, n=4, calibrated 31P-MRS, human calf): **[PCr]=33+/-2 mM, [Pi]=4.5+/-0.2 mM, [ATP]=8.2+/-0.4 mM**. Also: freeze-clamp biopsy PCr reads ~20% *lower* than 31P-MRS (CK-mediated breakdown during freezing — a disclosed cross-modality artifact). |
| 4 | Meyer RA (1988). A linear model of muscle respiration explains monoexponential phosphocreatine changes. *Am J Physiol* 254(4 Pt 1):C548-53. | **3354652**, `10.1152/ajpcell.1988.254.4.C548` | **The foundational model.** Rat gastrocnemius (species flagged), tau intensity-*independent* across stimulation rates (mean 1.44 min). Their own words: *"a simple first-order electrical analog model... assumes equilibrium of the creatine kinase reaction, which is modeled as a **chemical capacitor**, with capacitance proportional to the total creatine level, and PCr level proportional to the cytosolic free energy of ATP hydrolysis."* |
| 5 | Arnold DL, Matthews PM, Radda GK (1984). Metabolic recovery after exercise and the assessment of mitochondrial function in vivo in human skeletal muscle by means of 31P NMR. *Magn Reson Med* 1(3):307-15. | **6571561**, `10.1002/mrm.1910010303` | Founding human 31P-NMR paper. **The primary source for this task's own pre-registered open caveat**: severe-vs-mild exercise gave slower *apparent* PCr recovery, traced to intracellular **acidosis** (pH shifts the CK equilibrium itself), not necessarily lower mitochondrial capacity — recommends pH-corrected free-[ADP] recovery as the more specific index. |
| 6 | Ryan TE, Southern WM, Reynolds MA, McCully KK (2013). A cross-validation of near-infrared spectroscopy measurements of skeletal muscle oxidative capacity with phosphorus magnetic resonance spectroscopy. *J Appl Physiol* 115(12):1757-66. | **24136110**, `10.1152/japplphysiol.00835.2013`, PMC3882936 | **Primary tau anchor.** n=16 healthy young adults (22.5+/-3.0y), calf, ~10s plantar-flexion, [PCr] 34.25+/-2.86 -> 25.89+/-2.96 mM (24% depletion). **tau(31P-MRS)=31.5+/-8.5s; tau(NIRS muscle-VO2-recovery, a decorrelated optical modality, same subjects)=31.5+/-8.9s** (r=0.88-0.95) — two independent physical observables agreeing on one hidden state. Qmax=1.16+/-0.28 mM ATP/s via their own Qmax=[PCr]_rest·k_PCr. |
| 7 | Conley KE, Jubrias SA, Esselman PC (2000). Oxidative capacity and ageing in human muscle. *J Physiol* 526(Pt 1):203-10. | **10878112**, `10.1111/j.1469-7793.2000.t01-1-00203.x`, PMC2269983 | **Primary deconditioning/aging anchor** (this task's own explicit ask). Vastus lateralis, n=9 adults(38.8y)/40 elderly(68.8y). Explicitly names **"the linear model of oxidative phosphorylation described by Meyer"**, uses "Oxidative capacity = k_PCr·[PCr]_rest" — the *same* formula this script uses, found independently. [PCr]_rest: adult 29.24+/-1.60, elderly 27.11+/-0.88 mM. **Qmax: adult 1.16+/-0.147, elderly 0.61+/-0.04 mM ATP/s (47.4% lower)**; biopsy-confirmed mitochondrial volume density adult 3.6+/-0.11% vs elderly 2.9+/-0.15%. |
| 8 | Hinkle PC (2005). P/O ratios of mitochondrial oxidative phosphorylation. *Biochim Biophys Acta* 1706(1-2):1-11. | **15620362**, `10.1016/j.bbabio.2004.09.004` | P/O ratio: classic consensus **2.5 (NADH-linked substrates), 1.5 (succinate)**; revised (H+/ATP=10/3 stoichiometry) 2.3 (NADH), 1.4 (succinate). **Task-brief correction**: the task's own stated range (2.3-2.7) is only partly supported — Hinkle's verified abstract tops out at 2.5, not 2.7; this script uses the verified 2.3-2.5 band. |
| 9 | Jeneson JAL, Wiseman RW, Kushmerick MJ (1997). Non-invasive quantitative 31P MRS assay of mitochondrial function in skeletal muscle in situ. *Mol Cell Biochem* 174(1-2):17-22. | **9309660** (no DOI on record) | **Forced adversary** — forearm-flexor Qmax=0.24+/-0.06 mM ATP/s, via a *different* protocol (steady-state flux-vs-driving-force fit, not recovery kinetics). Deliberately the WRONG muscle group. Also: cytosolic ATP-hydrolysis midpoint free energy dGp,0.5 = **58.1+/-1.2 kJ/mol**. |
| 10 | Wallimann T, Wyss M, Brdiczka D, Nicolay K, Eppenberger HM (1992). Intracellular compartmentation, structure and function of creatine kinase isoenzymes... the 'phosphocreatine circuit' for cellular energy homeostasis. *Biochem J* 281(Pt 1):21-40. | **1731757**, `10.1042/bj2810021`, PMC1130636 | Landmark review, cited for the **spatial** PCr-shuttle concept (context only — no numeric quote extracted; abstract text not returned by `efetch` this session, disclosed, not fabricated). |
| 11 | Harris RC, Soderlund K, Hultman E (1992). Elevation of creatine in resting and exercised muscle of normal subjects by creatine supplementation. *Clin Sci (Lond)* 83(3):367-74. | **1327657**, `10.1042/cs0830367` | Corroborating quote: *"No changes were apparent in the muscle ATP content"* even under a real creatine-supplementation perturbation — supporting evidence for ATP stability, not a source of any gated number. |

**Self-correction, disclosed rather than hidden** (same discipline this repo's other docs already
apply): the task brief's own framing implied a 2015 Kemp/Meyerspeer/Moser *NMR Biomed* review; a
live search found no such record. The real paper is the **2007** one above (same three authors,
same journal) — found and substituted before use, not silently guessed.

**Second discrepancy vs the task brief, caught in this session's own final review pass, not
omitted**: the task states resting [PCr] ~= 25-30 mM. Checking all four verified values directly:
**calf muscle reads systematically above that band** (Kemp 2007: 33.0mM; Ryan 2013: 34.25mM — both
independently verified, both above 30mM), while **vastus lateralis/quadriceps sits inside it**
(Conley 2000: adult 29.24mM, elderly 27.11mM — both within [25,30]). This looks like a genuine,
muscle-group-dependent effect (consistent with calf/plantar-flexor muscle's known different
fiber-type composition vs quadriceps), not a citation error — both calf papers agree closely with
each other (33.0 vs 34.25, 3.7% apart) and both quadriceps/vastus-lateralis numbers agree closely
with each other (29.24 vs 27.11) — but it means the task's own stated 25-30mM band is a better
match to thigh muscle than to the calf muscle this doc's primary tau anchor (Ryan 2013) actually
used. Flagged, not smoothed over.

## Section A — resting state + the creatine-kinase near-equilibrium constraint

**Geometry, not heuristics**: the system lives on two conservation-law manifolds — total creatine
(Cr_tot = [Cr]+[PCr]) and total adenine nucleotide (TAN = [ATP]+[ADP], ignoring the much smaller
AMP/myokinase pool, a disclosed simplification) — intersected with the CK equilibrium surface
K_CK = [ATP][Cr]/([ADP][PCr][H+]). Solving these three equations together gives a **closed form**
(no root-finding needed):

```
[ADP] = [ATP][Cr] / (K_CK · [PCr] · [H+])                         (equilibrium surface)
[ATP] = TAN / (1 + [Cr]/(K_CK·[PCr]·[H+]))                         (closed form after substitution)
```

This is the **"recover the unobserved from an over-determined constraint"** move: free [ADP] is
*not* 31P-MRS-visible (per Veech et al. 1979, ~20-fold below measured total cellular ADP) but is
algebraically recoverable from the MRS-**observable** manifold {[PCr], [Cr], [ATP], pH} via the
equilibrium constraint — structurally the same shape (not asserted numerically identical) as this
project's own sigma_min/over-determination framing elsewhere: the cheap, directly-measurable
quantities pin down the one the instrument cannot see.

**Machine self-check** (a pure algebra check, not physiology): plugging [PCr]=[PCr]_rest back into
the closed form must return [ATP]=[ATP]_rest exactly. **First run: 0.173% residual — a real bug**
(resting free-ADP, computed in micromolar, was converted to molar with `1e-9` instead of `1e-6`).
Fixed at source; **re-run: 0.00e+00 relative error** (floating-point exact). Caught by the gate
built for exactly this purpose, not by inspection — the payoff of building a self-check before
trusting downstream numbers.

**Resting free [ADP], derived (not measured — no citation claims a direct measurement of this
number, it is a consequence of the equilibrium constraint applied to Kemp 2007's measured
[PCr]/[ATP])**: **14.22 uM**, using an illustrative total-creatine pool Cr_tot=42.5mM (a standard
literature figure, **not independently re-verified via a live primary-source fetch this
session** — flagged, used only in this section, never in the tau falsifier below). This resting
ADP estimate is itself sensitive to the assumed Cr_tot (a +/-20% Cr_tot sweep moves it from 1.50
to 26.94 uM) — reported honestly as an order-of-magnitude ("low micromolar") result, not a precise
point estimate.

### The buffering demonstration — machine-computed, not asserted

Sweeping [PCr] from rest down to 90% depletion (a severe, near-exhaustive bout), holding Cr_tot
and TAN conserved throughout:

| PCr depletion | [PCr] (mM) | [ATP] with buffer (mM) | ATP change | [ADP] (uM) | [ATP] with **NO** buffer (mM) | feasible? |
|---:|---:|---:|---:|---:|---:|---|
| 0% | 33.00 | 8.2000 | 0.00% | 14.2 | +8.20 | yes |
| 10% | 29.70 | 8.1929 | -0.09% | 21.3 | +4.90 | yes |
| 25% | 24.75 | 8.1789 | -0.26% | 35.3 | -0.05 | **no** |
| 50% | 16.50 | 8.1370 | -0.77% | 77.2 | -8.30 | **no** |
| 75% | 8.25 | 8.0138 | -2.27% | 200.4 | -16.55 | **no** |
| 90% | 3.30 | 7.6657 | **-6.52%** | 548.6 | -21.50 | **no** |

**Forced adversary (not a hand-wave)**: the "NO buffer" column asks what would happen to [ATP] if
the *same* amount of high-energy phosphate (X = PCr_rest - PCr_new) were drawn directly from ATP
hydrolysis with no CK/PCr system at all. It goes **infeasible (negative) by ~25% PCr-equivalent
depletion** — i.e. without the buffer, the entire resting ATP pool (8.2mM) would be exhausted by
an amount of contractile work that, WITH the buffer, barely dents [ATP] (-0.26% at that same X).
**The PCr/CK system extends the usable high-energy-phosphate reserve by a factor of ~4.0x**
(PCr_rest/ATP_rest = 33.0/8.2) while keeping [ATP] within single-digit percent of rest even at
90% pool depletion. Robustness check (forced variant, not assumed): repeating the 90%-depletion
computation at Cr_tot +/-20% (34.0 / 42.5 / 51.0 mM) gives ATP changes of -5.29% / -6.52% / -7.71%
— all comfortably inside the same conclusion, confirming this is not an artifact of the one
illustrative Cr_tot value chosen.

Corroborated qualitatively by Harris, Soderlund & Hultman (1992, PMID 1327657): *"No changes were
apparent in the muscle ATP content"* even under a real creatine-supplementation perturbation of
the total creatine pool.

## Section B — the tau falsifier (the task's central question)

**Model** (Meyer 1988's linear/circuit model, the same formula both Ryan 2013 and Conley 2000
independently use): monoexponential recovery `PCr(t) = PCr_end + delta_PCr·(1-e^(-t/tau))`, with
`tau = [PCr]_rest / Qmax` (the initial-rate relation — at t=0, ADP is near its most-elevated,
maximally mitochondria-stimulating level).

**Pre-registered discriminating test**: combine a Qmax measured in **one** independent study (a
*different* muscle, lab, subject pool, and protocol than the tau measurement) with a resting
[PCr], and check whether the **predicted** tau reproduces Ryan/McCully (2013)'s **directly
measured** tau — without ever touching that paper's own internally-computed Qmax.

| source of Qmax | muscle / study | Qmax (mM ATP/s) | PCr_rest used | **tau predicted (s)** | ratio to measured | within 1 SD of measured (31.5+/-8.5s)? |
|---|---|---:|---|---:|---:|---|
| Conley et al. 2000, **adult** | vastus lateralis (independent study) | 1.16 | Kemp 2007 calf, 33.0mM | **28.4** | 0.903 | **YES** |
| Conley et al. 2000, **adult** | vastus lateralis (independent study) | 1.16 | Ryan 2013's own calf, 34.25mM | **29.5** | 0.937 | **YES** |
| Conley et al. 2000, **elderly** | vastus lateralis, deconditioning direction | 0.61 | Ryan 2013's own calf, 34.25mM | **56.1** | 1.782 | prolonged, correct direction |
| Jeneson et al. 1997, **FORCED ADVERSARY** | forearm flexors, wrong muscle group | 0.24 | Ryan 2013's own calf, 34.25mM | **142.7** | 4.531 | **falls outside 2 SD** — correctly rejected |

**Reading this result**: a Qmax measured in a *different leg muscle*, *different lab*, *different
protocol*, combined with an independently-measured resting PCr, predicts a tau within **~6-10% of
the directly-measured value** — both variants land inside 1 measurement SD. This is genuinely
non-circular: Ryan/McCully's own 1.16 mM ATP/s (same-paper Qmax, computed *from* their own tau) is
never used for this comparison; Conley 2000's adult Qmax (1.16 mM ATP/s, coincidentally
numerically close but from an entirely separate study/lab/muscle/protocol) is. **The forced
adversary correctly fails**: substituting a wrong-muscle-group Qmax (forearm — genuinely, and by a
well-known, cited physiological reason, far less oxidative than habitually-loaded leg muscle) gives
a tau of 142.7s, more than **4.5x** the measured value and outside 2 measurement SDs — proving the
model is sensitive to a real confound, not a tautology that accepts any input. **The elderly
comparison correctly reproduces the task's own qualitative deconditioning prediction**
(56.1s > 28-30s, i.e. prolonged) using a real, independently-measured 47.4%-lower Qmax — not an
assumed slowdown.

**What this does NOT resolve** (pre-registered, held open per the task): a genuinely
**disease**-specific (not merely aging) tau-prolongation number was searched for directly this
session — mitochondrial myopathy, chronic renal failure/uremic myopathy (PMID 8621800, Marrades
et al. 1996, which DOES show a real, quantitative disease-related 31P-MRS bioenergetic
abnormality: resting [PCr]/[Pi] 8.9+/-1.2 in chronic renal failure vs 10.9+/-1.5 in healthy
controls, p=0.01, though not a directly-quotable tau-in-seconds figure from the live-fetchable
abstract), and Kemp/Taylor/Styles/Radda (1993, PMID 8457430, notes proton efflux during recovery
"is increased in some cases of mitochondrial myopathy") — all confirm disease-related bioenergetic
abnormalities are real and documented, but **no clean, single disease-cohort tau-in-seconds number
was extracted from a live-fetchable abstract this session**. The **aging** comparison (Conley
2000) stands as the quantitative "deconditioned" exemplar; the **disease** side is corroborated
qualitatively (real, cited, disclosed) but not quantitatively gated. Disclosed, not papered over.

**Held OPEN per this task's own pre-registration** (Arnold, Matthews & Radda 1984's own finding):
tau's dependence on fiber type, training status, and **pH-correction method** is real and
unresolved by this script — intracellular acidosis during severe exercise slows the *apparent*
PCr recovery independently of true mitochondrial capacity; this script does not model pH drift
during recovery at all (constant pH=7.0 assumed throughout).

## Section C — P/O ratio -> O2 -> couples to this twin's own VO2/Fick chain

Converts Qmax (mM ATP/s) to an O2 consumption ceiling via `ATP_per_O2 = 2 x P/O` (O2 has 2 O
atoms). **DOC FIX 2026-07-28d (stale magnitude)**: this section previously named Hinkle's classic
P/O=2.5 (NADH-linked) as the operative band. The mechanism graph's current best point estimate is
now **P/O(NADH)=2.78+/-0.04** (`docs/MECHANISM_MITOCHONDRIAL_OXPHOS.md`, self-corrected 2026-07-26,
cross-checked to 0.2% against Mookerjee 2017's decorrelated whole-chain P/O_max=2.79) — Hinkle's
2.5 remains a correctly-cited historical consensus figure but is no longer this graph's operative
value. At P/O=2.78, ATP_per_O2=5.56 and the ceiling recomputes to 1.16/5.56 = **0.209 mM O2/s per
liter of active muscle** (a ~10% FALL from the old 0.232, not a rise — a higher P/O means fewer O2
molecules are needed per ATP, so the O2 *flux* ceiling for a fixed ATP-flux Qmax falls, even though
the ATP obtainable per unit O2 rises). `scripts/msk/muscle_energetics.py`'s own `PO_RATIO_NADH=2.5`
constant (L359) was NOT changed this session (out of scope for this doc-only fix — recomputing it
live would require re-running the script and re-verifying all 13 gates); the 0.209 figure above is
a hand-recomputation for this doc's narrative only, disclosed as such, not yet the script's live
output.

**Coupling to the twin's own already-measured work** (subject2/walking1, "combined-corrected"
config — `metabolic_cost.py`'s own literature-anchored tendon+mass-corrected variant, **not**
re-solved here): walking metabolic power = **400.4 W** (5.1197 W/kg x 78.2kg — reusing
`respiratory.py`'s own already gross-ified number, verified consistent: 3.9197 W/kg net +
1.2 W/kg basal = 5.1197 W/kg, matching exactly, avoiding the net-vs-gross definitional mismatch
`MECHANISM_METABOLIC_CALORIMETRY.md` already warns about). Converting via ATP free energy
(dGp,0.5=58.1 kJ/mol, Jeneson 1997) gives an implied whole-body ATP turnover of **6.891 mmol/s**.

Normalizing by this repo's **own already-established, literature-anchored active muscle mass**
(`metabolic_cost.py`'s `LITERATURE_ANCHORED_MUSCLE_MASS_KG = 17.2` kg — reused verbatim, **not**
the raw Fmax-derived 45.6kg, which that same script's own sanity gate already flags as
implausible at 58.3% of body mass vs its 5-30% band) and a standard muscle density (1.06 kg/L):
**0.425 mM ATP/s per liter of active muscle — 36.6% of the adult-leg Qmax ceiling.** Comfortably
submaximal (gate: <50%), as physiologically required for a self-selected-speed walk, and does not
exceed the ceiling (gate: <100%) — a real, non-vacuous constraint the twin's own number could in
principle have violated, and did not. **Recomputation note (DOC FIX 2026-07-28d)**: this 36.6%
figure is 0.425 (mM ATP/s, from the Jeneson dGp thermodynamic conversion) / 1.16 (mM ATP/s, Qmax
measured directly from 31P-MRS PCr-recovery kinetics) — an ATP-flux-to-ATP-flux ratio that does
**not** involve the P/O ratio anywhere in its derivation. The P/O correction above (2.5->2.78)
therefore leaves this gate's 36.6% margin **exactly unchanged**, not "more comfortable" — that
would require Qmax itself to be re-derived from an O2 supply constraint via P/O, which it is not.

**Internal consistency cross-check (disclosed as internal, not a second external validation)**:
the same power, routed through P/O instead of the caloric-equivalent-of-O2 method
`respiratory.py` already uses (E_O2=20.9 kJ/L), gave an O2 flux of 1.853 L/min at the old P/O=2.5
vs `respiratory.py`'s own already-computed **1.149 L/min** — a ratio of 1.613 (implied efficiency
1/1.613=62.0%). **DOC FIX 2026-07-28d (hand-recomputed at P/O=2.78, script constant unchanged —
same disclosure as above)**: O2 flux falls to **1.666 L/min** (higher P/O = fewer O2 molecules
needed for the same ATP turnover), ratio vs `respiratory.py`'s 1.149 L/min falls to **1.450**, and
the implied thermodynamic efficiency **rises to 1/1.450 = 69.0%** (still physically required to be
<100%; still a real, non-tautological direction check — captures *less* than full combustion
energy, just less-less than the old estimate implied). This is not an
error: it is the **implied thermodynamic efficiency of trapping substrate-combustion free energy
into ATP's phosphoanhydride bonds**, computed by comparing two independent, already-cited
stoichiometric conventions (P/O ratio + ATP free energy vs the caloric equivalent of O2) that
should differ by exactly this factor if oxidative phosphorylation is not 100% efficient.
**Gate**: falls in (0,1) at either P/O value — a weak sanity gate mainly ruling out a units/sign
bug, not an independently-cited validation of the specific 62.0%/69.0% figure (no dedicated
efficiency-measurement primary source was verified this session for this number — disclosed, not
asserted as externally anchored). Both sides of this comparison share the SAME upstream
`metabolic_cost.py` output, so this is internal-consistency evidence, not an independent
replication (this repo's own common-mode discipline, applied to itself).

## Section D — timescale separation (couples to the muscle-force / SO thread)

Gait cycle duration (subject2/walking1, already established elsewhere in this repo) = **1.57s**.
tau (measured) = 31.5s. Ratio = **0.0498** — the gait cycle is **~20x shorter** than the oxidative
recovery time constant. Structurally analogous (not numerically identical) to `muscle_fatigue.py`'s
own singular-perturbation argument (`L=10 >> F,R~1e-3`, a fast recruitment servo nested inside slow
fatigue/recovery kinetics): here, CK near-equilibrium re-establishes on a sub-second enzyme-kinetics
timescale (fast) while mitochondrial PCr resynthesis unfolds over tens of seconds (slow) — nested
inside a stride cycle shorter than either. **Consequence for the twin**: during steady-state
submaximal walking, PCr should settle to a mild, roughly-steady *partial*-depletion operating
point rather than fully depleting and recovering each stride — consistent with Section C's finding
that walking's implied ATP flux sits at 36.6% of Qmax (a sustainable, not stride-limited, load).

## Confidence tier

Per **this repo's own established ledger convention** (`docs/MECHANISM_TRUST_LEDGER.md`), **not**
the task brief's own generic phrasing ("in-vivo-anchored"): that legend explicitly reserves
`in-vivo-anchored` as a term-of-art for the OrthoLoad implant-telemetry joint-force family
specifically — `MECHANISM_METABOLIC_CALORIMETRY.md` already resolved this identical tension for a
different (indirect-calorimetry) real-human in-vivo dataset by downgrading to the next tier, and
this doc follows that same established precedent rather than re-litigating it. **Tier:
cadaveric-or-published-plausibility** — every anchor here (Kemp 2007, Ryan 2013, Conley 2000,
Hinkle 2005) is a genuinely external, decorrelated, **published, real in-vivo-human-measurement**
dataset (31P-MRS and NIRS are directly non-invasive in-vivo techniques, arguably a strong member of
this tier), but is compared against this twin's own quantities only in Section C/D (the primary
Section B falsifier is a pure cross-literature check with no twin-side model output at all).

## Honest gaps

1. **pH-correction held OPEN** (pre-registered): Arnold/Matthews/Radda 1984's own finding that
   apparent tau confounds true mitochondrial capacity with intracellular-acidosis recovery. Not
   modeled here (constant pH=7.0 throughout).
2. **Fiber-type and training-status dependence held OPEN** (pre-registered): not modeled; Conley
   2000's aging comparison is the only "varies with condition" axis quantified here.
3. **Disease-specific (not merely aging) tau-prolongation**: searched for directly, not found as a
   clean quotable number this session (Section B) — a disclosed gap, not a fabricated substitute.
4. **Cr_tot=42.5mM is illustrative**, not independently re-verified via a live primary-source fetch
   this session — used only in Section A's qualitative buffering demonstration (shown robust to
   +/-20%), never in the Section B tau falsifier (which needs only PCr_rest and Qmax).
5. **P/O ratio band — UPDATED 2026-07-28d (was stale)**: this script's live constant
   (`PO_RATIO_NADH=2.5`, `muscle_energetics.py` L359) still uses Hinkle 2005's classic 2.3-2.5
   verified band and was NOT changed this session (out of scope for a doc-only fix). The mechanism
   graph's own current best estimate has since moved to **P/O(NADH)=2.78+/-0.04**
   (`docs/MECHANISM_MITOCHONDRIAL_OXPHOS.md`, self-corrected 2026-07-26, cross-checked to 0.2%
   against Mookerjee 2017's decorrelated whole-chain P/O_max=2.79) — the task brief's originally-
   queried upper bound of 2.7 is now essentially confirmed (from a different, independent route
   than Hinkle's own abstract: a structural/flux cross-check, not a literature re-read). Section C
   above hand-recomputes what changes if the script's constant is updated to 2.78: the O2-ceiling
   figure and the internal-efficiency cross-check both move; the 36.6%-of-Qmax gate does NOT (it
   is P/O-independent by construction). Re-pointing the live script constant to 2.78 and
   re-verifying all 13 gates is left as an OPEN follow-up, not done here.
6. **Cross-muscle-group comparison, disclosed**: the tau falsifier's independent leg combines
   vastus lateralis Qmax (Conley) with calf PCr_rest (Ryan/Kemp) — the same class of caveat this
   repo's own `muscle_perfusion.py`/`MECHANISM_VASCULATURE.md` already flags for applying
   quadriceps flow data to calf.
7. **Section C/D coupling is single-subject, single-trial** (subject2/walking1) and downstream of
   ONE shared `metabolic_cost.py` instrument chain — internal consistency, not independent
   replication, exactly as this repo's own common-mode discipline requires disclosing.
8. **No `docs/MECHANISM_BLOOD_OXYGEN_TRANSPORT.md` exists** in this repo (checked via `find` before
   writing this doc) — the Fick-chain role is currently played jointly by `MECHANISM_RESPIRATORY.md`
   (VO2), `MECHANISM_CARDIAC.md` (Q=VO2/(a-vO2diff), Q=HR*SV), and `MECHANISM_VASCULATURE.md`
   (per-muscle flow) — coupled to here, not fabricated as one doc that doesn't exist.
9. **Lumped, well-mixed-compartment model only** — no spatial PCr-shuttle/microcompartmentation
   (Wallimann et al. 1992's own subject), no explicit glycolytic/anaerobic ATP contribution
   (assumed negligible at this submaximal steady-state intensity, not separately measured).
10. A real **unit-conversion bug** (mu M->M, `1e-9` vs `1e-6`) was caught by this script's own
    self-check gate during this build and fixed at source — reported here per this project's
    machine-cross-check discipline, not silently corrected.

## Files

- `scripts/msk/muscle_energetics.py` — the full model: CK equilibrium + closed-form ATP
  buffering algebra, the tau falsifier (3 independent-Qmax variants + 1 forced adversary), the
  P/O-to-VO2 coupling (reads `metabolic_cost_results.json`/`respiratory_results.json` read-only),
  timescale separation, all 13 gates, self-contained (no OpenSim dependency).
- `data/msk_smoketest/subject2_walking1/muscle_energetics/muscle_energetics_results.json` — full
  machine-readable evidence: every citation, every constant, the full buffering sweep, the Cr_tot
  robustness sweep, the tau falsifier table, the P/O/VO2 coupling numbers, all 13 gates
  (regenerated fresh each run).
- Inputs read read-only, never modified: `data/msk_smoketest/subject2_walking1/metabolic_cost/
  metabolic_cost_results.json`, `data/msk_smoketest/subject2_walking1/respiratory/
  respiratory_results.json`.
- Reused context, not modified: `docs/MECHANISM_METABOLIC_COST.md`, `docs/MECHANISM_RESPIRATORY.md`,
  `docs/MECHANISM_CARDIAC.md`, `docs/MECHANISM_VASCULATURE.md`, `docs/MECHANISM_MUSCLE_FATIGUE.md`
  (the singular-perturbation analogy), `docs/MECHANISM_CROSSBRIDGE_MODEL.md` (sibling ATP-consuming
  mechanistic layer), `docs/MECHANISM_TRUST_LEDGER.md` (tier vocabulary).

No git commit, no git push performed (isolation respected, per this repo's own COORDINATOR.md /
`MECHANISM_HARDENED_CONVENTIONS.md`).
