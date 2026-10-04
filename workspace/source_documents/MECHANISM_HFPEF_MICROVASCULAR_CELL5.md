# HFpEF cell 5 — does the microvascular chain discriminate, as the paradigm claims?

Pre-registration: `data/hfpef_diastolic_v1/PREREG_CELL5.md`, committed **before** any cell-5 anchor
value was in hand. Code: `scripts/tissuetwin/hfpef_cell5_microvascular.py`.
Evidence: `data/hfpef_diastolic_v1/cell5_evidence.json`. Anchors: `anchors/cell5_*.json`.
Upstream: cell 1 (`MECHANISM_HFPEF_DIASTOLIC.md`), cell 2 (`MECHANISM_HFPEF_MATERIAL_CELL2.md`).


## ★ The microvascular impairment is measurable; the DISCRIMINATION is not

Cell 5's upstream half rests on a ratio of ratios — HFpEF's coronary flow reserve relative to its
own control, against HFrEF's relative to its own — and that ratio, 0.932, has never carried a
detectability number. The dispersions were in the anchors the whole time, stated as **text** in an
`sd_or_iqr` field where no instrument could use them. Structured out of the records' own strings:

| study | comparison | values | n | \|t\| |
|---|---|---|---|---|
| Dryer | CFR, control vs HFpEF | 3.84 ± 1.89 vs 2.55 ± 1.60 | 14 / 30 | **2.21** |
| Dryer | IMR, control vs HFpEF | 19.7 vs 26.7 | 14 / 30 | **2.19** |
| Takafuji | CFR, control vs HFrEF | 4.03 ± 1.47 vs 2.87 ± 0.86 | 26 / 26 | **3.47** |

**Both diseases have a measurable microvascular deficit inside their own studies.** That is a real
positive and it had never been computed here.

But the quantity cell 5's conclusion actually turns on is the *difference between* those deficits:
0.664 in HFpEF against 0.712 in HFrEF, ratio **0.933** — reproduced exactly from the raw numbers.
Compared directly, HFpEF's 2.55 against HFrEF's 2.87 gives **|t| = 0.95**, and it is a **cross-study**
pairing besides, so it is inadmissible for two separate reasons.

Three more comparisons became computable when the parser learned to classify a dispersion *statistic*
rather than just grab its number — `SD 1.89`, `±0.62 (SEM, per journal convention)` and
`911–1188 (25th–75th percentile)` all needed different conversions, and reading the SEM as an SD
would have understated spread by √n:

| study | quantity | values | \|t\| |
|---|---|---|---|
| Mohammed | capillary density | 1044 vs 788 | **8.66** |
| Mohammed | microvessel density | 1316 vs 967 | **7.85** |
| van Heerebeek | **PKG activity** | 5.11 vs 9.18 | **6.51** |

The third matters most. Cell 5's downstream half rests entirely on van Heerebeek's PKG measurement —
the paper this cell reports as never independently replicated. Its internal evidence is now known to
be **strong**: |t| 6.51 within its own study. What is missing from that leg is *replication*, not
*significance*, and those are different failures. The rarefaction numbers behind G8 are likewise
strongly resolvable.

That is what `CHAIN_FAILS_TO_DISCRIMINATE` means, now as a number rather than a band: the deficits
are each real, each strongly measured, and the data still cannot tell HFpEF from HFrEF apart.

**And the null has a bound, which it did not before.** "Does not discriminate" reports that a
difference was not *detected*; what size of difference would have had to be there is a separate
number:

    observed HFpEF-vs-HFrEF CFR gap        0.32 units      |t| 0.95
    smallest gap resolvable at 80% power   0.94 units      3.0x the observed
                                                           33% of the HFrEF mean

So the verdict is **silent about a modest discrimination, not evidence against one** — the same
correction G4's titin null needed. A chain that discriminated by, say, half a CFR unit would look
exactly like this in these data.

Two independent reasons the contrast is inadmissible anyway, both kept on the record rather than
only the more convenient one: it is **cross-study** (Dryer against Takafuji, different labs and
protocols, where a sibling falsifier measured control groups for one quantity differing **2.563×**
between labs), and **each arm is itself underpowered** — Dryer's own control-vs-HFpEF positive
reaches |t| 2.21 while its 80 % floor is 43 % of the control mean against an observed 33.6 %. Dryer's own power adds the third leg — the
smallest CFR reduction it could have detected at 80 % power is **43 %** of the control mean, against
an observed 33.6 %, so even its positive is underpowered for the effect it found.


## G2 — the stiffness budget cannot be formed: PASS by abstention
The pre-registered decisive number was S_micro = f_density × f_PKG × f_titin. It cannot be built:
family D — the change in passive tension per unit PKG manipulation — has **zero quantitative record
for reduced-EF anywhere**. The chain has no comparator arm at its final step, which the
pre-registration fixed in advance as a PASS by abstention rather than a failure.

That is itself worth stating plainly: **the step that would make this chain a quantitative
explanation of stiffness has never been measured in the disease the chain is supposed to be
distinguished from.**

## G3 — the discrimination test, which the pre-registration called the real point
Two adjacent links *do* have both arms, and together they make a checkable prediction: if
microvascular dysfunction **causes** the PKG deficit, the disease with the worse flow reserve must
have the worse PKG.

Each arm is compared only against a control measured in the **same study with the same quantity** —
CFR, MFR, MBFR and MPRi are not one quantity, so cross-study absolute values are reported but never
differenced.

| link | HFpEF (ratio to own control) | HFrEF (ratio to own control) | HFpEF/HFrEF |
|---|---|---|---|
| **upstream** — coronary flow reserve | 0.664 | 0.712 | **0.933** |
| **downstream** — PKG activity | 0.695 | 1.133 | **0.614** |

Against cell 1's own resolution band [0.70, 1.43]:

- the **upstream** step does **not** discriminate — microvascular impairment is essentially equal in
  the two diseases (0.933, inside the band);
- the **downstream** step discriminates **strongly** — the PKG deficit is HFpEF-specific (0.614,
  outside the band), and reduced-EF PKG activity is not merely spared but sits *above* its non-HF
  comparator (1.133).

### The downstream step has TWO quantities, and they are not one measurement
The first version of this cell pooled them, and that pooling was deciding the verdict. Separated:

| downstream quantity | HFpEF | HFrEF | ratio | discriminates? |
|---|---|---|---|---|
| **PKG activity** — homogenate assay at saturating cGMP, i.e. PKG **capacity** | 0.557 | 1.254 | **0.444** | yes |
| **PKG activity index** — pVASP/VASP, an **in vivo** phosphorylation readout | 0.833 | 1.012 | **0.824** | no |

Capacity is what PKG *could* do if fully stimulated; the index is what it *is* doing in the intact
heart. A myocardium can carry high capacity and low signalling — that is exactly what elevated PDE5
would produce — so these are different quantities and the pooled median of 0.6135 averaged across a
distinction that decides the answer. Same class as the CFR/IMR polarity error earlier in this cell,
one level subtler: not opposite directions, but different measurements.

**So the conclusion has to be stated over both readings, and it survives either way — but differently:**

- on the **in vivo** readout, *neither* link discriminates (upstream 0.93, downstream 0.82). The
  chain is then internally **consistent** — and fails to deliver the HFpEF-specificity it was
  invoked to explain;
- on the **capacity** assay, the downstream discriminates (0.44) while the upstream does not.
  Internally **inconsistent**.

Under both readings the chain does not supply a discriminating mechanism. What changes is *which
way* it fails. `CHAIN_FAILS_TO_DISCRIMINATE_UNDER_EITHER_PKG_READOUT`.

**Superseded:** the earlier headline `CHAIN_IS_INTERNALLY_INCONSISTENT` came from the pooled
downstream median and is not defensible once the two quantities are separated. It is kept in the
evidence file under `pooled_verdict_superseded` rather than deleted.

Cell 1's measured material-stiffening ratio, 0.83, sits inside the same band — i.e. it agrees with
the upstream step and disagrees with the downstream one. That is hypothesis **H1** from the
pre-registration: the route is real but is not what differentiates the two diseases.

## The pairing-rule falsifier — how much of that verdict is my choice of pairings?
The upstream comparison is n = 1 study per arm, which is thin for a claim this strong, and the
paired-within-study rule was the strictest defensible one — but strictness chosen after seeing the
answer is not a control. So the upstream ratio was recomputed under three admissibility rules:

| rule | admissible pairs | upstream ratio span | reading |
|---|---|---|---|
| **R1** same study *and* same quantity | 1 × 1 | **0.932** | does not discriminate |
| **R2** same quantity, control from any study | 4 × 2 | **0.847 – 1.024** | does not discriminate |
| **R3** any reserve quantity, any study | 6 × 9 | **0.681 – 3.121** | discriminates, mostly toward *reduced-EF* being worse |

**Under the two defensible rules the upstream step does not discriminate, and the cell's conclusion
stands.** Under R3 — which mixes CFR with MBFR and MPRi, related but not identical measurements —
the span straddles: most of it points to reduced-EF having *worse* flow reserve (up to 3.12×), which
does not rescue the chain but contradicts it harder, since that disease would then have worse
perfusion *and* better PKG.

But the bottom edge, **0.681**, grazes below the band, so a pairing exists that supports the chain,
and that has to be said. It is the weakest one available: it divides a preserved-EF ratio built from
two different studies' CFR by a reduced-EF ratio built from one study's CFR over *another* study's
MFR — cross-study and cross-quantity at both ends. The verdict is recorded as
`A_PAIRING_RULE_EXISTS_THAT_SUPPORTS_THE_CHAIN` rather than as an unqualified negative.

So the honest form of the headline: **the two adjacent links order the diseases differently under
every same-quantity comparison; only by mixing measurement types can the upstream step be made to
support the chain, and the same mixing more often contradicts it.**

## What this rests on, stated as weakly as it deserves
- **Upstream: n = 1 study per arm.** Only studies carrying their own control can enter a paired
  comparison, and only one per arm does.
- **Downstream: n = 2 records but ONE paper** (van Heerebeek 2012). The discrimination claim's
  strongest quantitative support is one lab's one cohort.
- **That paper is doing triple duty.** It spans cell 5's families C and D *and* already served
  cell 2's family B. Kolijn likewise spans C and D; Zile also crosses from cell 2. One cohort
  carrying anchor roles across two cells is the strongest form of the tautology the
  pre-registration named — recorded as n_eff=1, not laundered.
- **cGMP concentration itself has no recoverable number for any group in any arm**, so the link
  between the microvascular step and the PKG step is not measured at all — it is inferred.

## An error of mine, caught in the first run
The first pairing pooled Dryer's CFR (higher = better perfusion) with Dryer's **IMR** — index of
microcirculatory **resistance**, where higher is worse — as if both were ratios-to-control in the
same direction. The median of 0.664 and 1.355 came out at **1.0097**, a number that reads as "HFpEF
perfusion is normal" and is pure sign error. Polarity is now declared by name: only reserve-type
quantities may enter a ratio, and resistance indices are excluded explicitly rather than by hoping
they look different.

## G8 — the orphaned anchor finally does work, and it decorrelates the upstream link
The brief's third anchor, myocardial capillary density, was harvested and consumed by nothing —
found by the pooling audit rather than noticed. It cannot extend the discrimination test, because
family A has no reduced-EF arm with a number. But it can do something the cell needed: give an
**independent measurement family** for the same upstream impairment, which until now rested on
coronary flow reserve alone at n = 1 study per arm. Biopsy morphometry and coronary physiology share
no method and no cohort.

| quantity (kept separate — they count different things) | healthy | HFpEF | ratio | implied diffusion distance |
|---|---|---|---|---|
| capillary density | 1044 | 788 | **0.755** | ×**1.151** |
| microvessel density (capillaries + arterioles) | 1316 | 967 | **0.735** | ×**1.167** |
| *functional*, coronary flow reserve | — | — | **0.664** | — |

### The agreement criterion was borrowed; it is now derived
The first version declared the families to agree because their ratio-of-ratios fell inside **cell
1's** resolution band — a band defined for *material stiffening* ratios and reused because it was to
hand. Nothing physical says a structural count ratio and a functional reserve ratio must be
numerically equal, so "inside a borrowed band" is not a test.

There is a derivable relation. Treat the microvascular bed as parallel conductances: if vessel
number falls by a factor f, maximal conductance falls ~f while resting flow is autoregulated to
unchanged demand, and flow reserve is max/rest — so the prediction is **CFR_ratio = density_ratio**.

| structural quantity | predicted CFR ratio | observed | unexplained factor | share of the log loss explained by rarefaction |
|---|---|---|---|---|
| capillary density | 0.755 | 0.664 | **0.880** | **68.7 %** |
| microvessel density | 0.735 | 0.664 | **0.904** | **75.3 %** |

**Vessel rarefaction accounts for roughly 70–75 % of the measured flow-reserve loss**; the remaining
factor of 0.88–0.90 is what vessel count does not explain — vasodilator capacity, endothelial
function, whatever else lives there. That is a measured residual, not a story.

**The assumption that most weakens it, stated rather than buried:** coronary flow reserve is
dominated by **arteriolar** resistance, not by capillary number. Capillary rarefaction and
arteriolar reserve are different vascular beds, so a close match may be partly coincidental. The
residual is reported as measured; no causal reading is claimed.

The diffusion-distance factor is pure geometry with no free parameter — mean intercapillary spacing
goes as 1/√(areal density) — so the measured rarefaction implies a **15–17 % longer diffusion path**
in HFpEF myocardium.

The two density quantities are kept **separate** rather than pooled. They were nearly averaged
together, and the pooling audit that caught it exists because this cell had already made that
mistake twice.

## The replication hunt — a searched negative, and a direction conflict
Cell 5 named two measurements that would decide it. Both were hunted (PubMed, Europe PMC, direct
fetches, plus two review articles checked specifically for citations beyond the two known sources).

**Neither exists.** **Zero** labs independent of Amsterdam/Paulus and the linked Bochum network have
measured PKG activity, a PKG activity index, or cGMP concentration **with a number** in reduced-EF
human myocardium. Family D's reduced-EF cell is likewise still empty, so G2 stays an abstention.

That changes the standing of the cell's weakness rather than removing it: *unreplicated* was an
inference before, and is now a **searched** negative. **The field's discrimination claim has never
been independently replicated at its decisive step.**

Two independent labs did measure *adjacent* quantities in reduced-EF human tissue, and one of them
matters:

- **Pokreisz 2009** (Leuven, independent, and with **true non-failing controls** — the best
  comparator design found anywhere in this harvest) found **PDE5, the cGMP-*degrading* enzyme,
  elevated** in explanted dilated and ischaemic human left ventricles. Elevated PDE5 predicts
  **lower** cGMP/PKG in reduced-EF — the **opposite** direction to the single paper this cell
  depends on, where reduced-EF PKG activity is the *highest* of all three groups.
- **Persoon 2018** (Regensburg, independent) found cGMP/PKG-I down-regulated with unloading in
  dilated but not ischaemic cardiomyopathy. Wrong comparator for this purpose.

Neither carries a recoverable number, so neither settles anything. But the only independent evidence
bearing on the *direction* of the decisive step points against the paper that supplies it.

## What would decide it
- **A PKG-activity measurement in reduced-EF myocardium from a second, independent lab.** The whole
  discrimination rests on one cohort; a replication either way settles it.
- **A PKG dose-response on reduced-EF cardiomyocytes** — family D's empty cell. That single
  measurement would also convert G2 from an abstention into a real stiffness budget.
