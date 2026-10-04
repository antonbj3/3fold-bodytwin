# HFpEF cell 3 — the brief's knife-edge: is HFpEF several diseases, mechanically?

Code: `scripts/tissuetwin/hfpef_cell3_dispersion.py`. Evidence:
`data/hfpef_diastolic_v1/cell3_dispersion_evidence.json`. Threshold pre-registered in
`data/hfpef_diastolic_v1/PREREG.md`, section "THE KNIFE-EDGE THE BRIEF NAMES, MADE MEASURABLE".

The brief: *"HFpEF is probably several diseases with one name. If the clusters cannot be distinguished
mechanically, THAT is the finding — report it, do not force a model."*

## ★ The question this cell cannot reach, now measured elsewhere

Everything above is the spread **between cohort means**. This doc has always said cohort means are
by construction blind to multimodality *within* a cohort; harvest pass 4 landed a measurement of
exactly that, on the most direct stiffness parameter in the project.

Lam CS (PMID 17404159), one study, three groups, EDPVR exponent β:

| group | n | β | dispersion | CV |
|---|---|---|---|---|
| community controls | 617 | 5.96 | ± 0.06 | 1.0 % |
| hypertension, no HF | 719 | 6.05 | ± 0.41 | 6.8 % |
| HF, normal EF | 244 | 7.09 | ± **3.55** | **50.1 %** |

**The means barely separate. The dispersions separate by tens.** CV ratio HFpEF/control = **49.74**.

A 1 % CV for a fitted biological exponent across 617 subjects is tight enough that the dispersion
could have been a standard *error* mislabelled as a deviation, which would have made the ratio
31.28 instead. Rather than hedge, I fetched the source: **PMC2001291 Table 2's footnote reads
verbatim "Data are mean ± SD; Comparisons adjusted for age and sex, as well as body surface area
(BSA) where appropriate"**. The paper rules the SEM reading out itself. Both readings crossed cell
3's pre-registered 2× line by more than an order of magnitude anyway, but the number that stands is
49.74 and it stands for a stated reason.

The same fetch confirmed α is **not tabulated** here either — the fifth consecutive record to omit
the one parameter that would lift every stiffness value in this project off its floor.

This does **not** contradict cell 3. Cell 3 measured between-cohort spread of geometry and found
none; this is within-cohort spread of a stiffness exponent, and they are different questions. What
it does mean is that the knife-edge result should never be quoted as "HFpEF is not heterogeneous" —
only as "HFpEF cohort *means* are not more spread than HFrEF's".

Four things it does not establish, all recorded rather than argued away: the β is **noninvasive**
(echo E/e′-derived LVEDP), a different evidence tier from every invasive record; it is a **power-law**
β and must not be ratioed against the exponential family cell 1 uses; dispersion in a sicker group
can be measurement noise as much as biology, and nothing here separates them; and `n_eff = 1`, one
study sharing an author lineage with the other power-law study on file.

## The threshold was pre-registered; the quantity is substituted, and that is declared
The pre-registration fixed a dispersion ratio of **2x** as the line, on **invasive** kappa per
cohort. Three literature passes established that the invasive family cannot bear per-cohort kappa
values at all (one citation-only record per arm, floors only), so that exact test cannot run.
Substituted: **model-predicted kappa per cohort** — still a kappa, still per cohort, predicted
rather than measured — with the 2x threshold carried over unchanged. Dispersion on the two raw
geometric inputs is reported beside it, because a dispersion result living in a single derived
quantity is not a finding.

## Measured

| quantity | HFpEF spread | HFrEF spread | healthy spread | ratio HFpEF/HFrEF | verdict (2x) |
|---|---|---|---|---|---|
| wall fraction phi, all modalities | **3.761x** | 1.646x | 1.504x | **2.284** | ★ crosses 2x |
| wall fraction phi, **CMR only** | 1.270x | 1.646x | 1.504x | **0.771** | not over-dispersed |
| predicted kappa, all | 3.077x | 2.365x | 2.113x | 1.301 | not over-dispersed |
| predicted kappa, **CMR only** | 1.178x | 2.365x | 2.113x | **0.498** | not over-dispersed |
| absolute cavity volume, all | 2.048x | 1.197x | 1.419x | 1.712 | not over-dispersed |
| absolute cavity volume, **CMR only** | 1.155x | 1.197x | 1.419x | 0.965 | not over-dispersed |

★ **This result changed, and the change is the honest headline.** With cell 1's cohort coverage
completed — two records had been kept out by the spelling of a JSON key — HFpEF's wall-fraction
spread across **all modalities** rose to **3.761x**, ratio **2.284**, which **crosses the
pre-registered 2x line for the first time**. The verdict is now
`OVER_DISPERSION_IS_A_MODALITY_ARTEFACT`, and it rests entirely on the modality control: restricted
to CMR, HFpEF's spread is **1.270x**, *below* HFrEF's 1.646x, ratio **0.771**.

So the knife-edge conclusion survives — HFpEF is not mechanically over-dispersed — but it no longer
survives without the control. The excess is carried by echocardiographic cohorts, now two of them:
the Bangladesh study (phi 1.596) and Venkateshvaran (phi **2.422**, LV mass 199.42 g against an EDV
of 78.41 mL, both absolute and mutually consistent — an extreme concentric cohort, not a units
error, and it was checked for one because Liu's 1.0604 turned out to be exactly that).

Echo and CMR do not measure the
same volumes, and that single modality difference was the whole apparent signal.

On absolute cavity volume HFpEF is if anything slightly *less* dispersed than HFrEF (ratio 0.965).

## What this instrument structurally cannot see — the load-bearing caveat
These are **cohort means, not patients**. Sub-phenotypes of a syndrome appear as multimodality
*within* a cohort, and the between-cohort spread of means is blind to it by construction: five
cohorts each containing two distinct phenotypes in the same proportion would show **zero** excess
dispersion here. So this is a negative on the **between-cohort mechanical axis only**.

It does not say HFpEF is one disease. It says chamber mechanics does not separate these cohorts —
which is precisely the finding the brief asked to be reported rather than modelled around. Two
further limits belong with it:

- HFpEF sub-phenotypes are usually proposed on **comorbidity / metabolic** axes (obesity, diabetes,
  hypertension), not on chamber geometry. A null on the mechanical axis is not a null on those.
- n = 5–7 cohorts per group. The spread statistic is max/min, chosen because an IQR is barely
  defined at that n, and it is the *most* outlier-sensitive choice — so it favours finding
  over-dispersion, and still found none.

## A second, independent instrument: WITHIN-cohort dispersion
Cohort means are blind to within-syndrome spread by construction. Within-cohort standard
deviations were re-surfaced from fields already present in the geometry anchors, giving a
coefficient-of-variation comparison at a level cell 3's first instrument cannot reach.

**The indexing bias had to be controlled first, and it runs toward a false positive.** Per-patient
indexing to body surface area removes body-size variance, so an indexed CV is systematically
smaller than an absolute one — and the reduced-EF arm on disk is 4/4 indexed while the healthy arm
is mostly absolute. Comparing across that would deflate the reduced-EF CV and make HFpEF look more
dispersed than it is: exactly the direction that would manufacture the answer being watched for.

| quantity | HFpEF median CV | HFrEF median CV | healthy | ratio HFpEF/HFrEF | verdict (2x) |
|---|---|---|---|---|---|
| EDV, **indexed vs indexed** | 0.2425 (n=3) | 0.2881 (n=4) | 0.1467 (n=1) | **0.842** | not over-dispersed |
| LV mass, **indexed vs indexed** | 0.1935 (n=2) | 0.3144 (n=4) | 0.2405 (n=1) | **0.615** | not over-dispersed |
| LV mass, mixed indexing (biased) | 0.3000 (n=4) | 0.3144 (n=4) | 0.2405 (n=3) | 0.954 | — |
| wall thickness | 0.1345 (n=2) | **no data (n=0)** | 0.1652 (n=2) | — | **untestable** |

**HFpEF is not more heterogeneous within cohorts either.** The bias control was worth running: on
LV mass, mixed indexing gives 0.954 and the clean comparison gives 0.615, so the artefact was
inflating HFpEF's apparent dispersion by a factor of 1.55.

**But leave-one-out separates what survives from what does not, and I had overclaimed.** At n = 2–4
cohorts per arm a median is one or two numbers, so the direction was tested by deletion:

| quantity | full ratio | range under any single deletion | pre-registered verdict robust? | direction robust? |
|---|---|---|---|---|
| EDV | 0.842 | 0.550 – **1.060** | yes (max 1.060 ≪ 2) | **no** |
| LV mass | 0.615 | 0.236 – **0.995** | yes | only just — one deletion reaches parity |

Dropping a single **17-patient** cohort (Shao) takes EDV to 1.060 and mass to 0.995. So the
*pre-registered* verdict — not over-dispersed, ratio far below 2 — survives every deletion, and
that is what cell 3 is entitled to claim. The *directional* reading "HFpEF is less dispersed than
HFrEF" does not survive and is withdrawn: without that one small cohort both quantities sit at
parity.

And the gap is now exact: **wall thickness has no reduced-EF arm at all (n=0)**, which is precisely
the axis the residual live possibility lives on.

## Published phenogroups, checked for circularity — and the strongest positive is circular
A separate harvest went after the data cell 3 is blind to: patient-level or within-cohort
distributions, and **published HFpEF phenogroups reporting chamber geometry per cluster**. Three
independent phenomapping papers were found with per-cluster geometry
(`anchors/geometry_patient_level.json`); the check is in
`anchors/phenogroup_circularity_check.json`.

| paper | signal | checked | verdict |
|---|---|---|---|
| Kyodo 2023 (n=365 + 273 validation) | wall thickness across phenogroups, **p<0.0001, replicated** | its Supporting Table S1 lists the 24 clustering variables: **"Echocardiography at admission IVST, PWT, LVDd, LVDs"** | **CIRCULAR** |
| Fazzini 2024 (n=591) | LV EDV index p=0.52, LV mass p=0.093 — **null on both** | paper names which variables were *excluded* from clustering; geometry is not on that list | null is the informative direction |
| Tamaki 2024 (n=211) | LVEDD p=0.009; LV mass index p=0.098 | input list not retrieved | unresolved — treated as not independent |

**The one twice-replicated positive is a restatement of its own construction.** Clusters built
partly *on* IVST and PWT differ in IVST and PWT; the validation cohort was phenomapped "using the
same 24 variables", so it reproduces the construction rather than testing it. The finding may be
clinically useful; it cannot serve as independent evidence of within-syndrome mechanical
multimodality.

This is the **second** time in this cell that the strongest-looking positive turned out circular —
Sinning 2011 split its DHF/no-DHF cohorts by the measured EDPVR itself, confirmed twice on raw text.

**Net effect on cell 3: the verdict is unchanged and better supported.** The nulls are on
end-diastolic volume and mass — exactly the two quantities cell 3 tested — and they arrive from a
third, independent direction: within-cohort clustering, after between-cohort spread and
within-cohort CV. All three levels agree. A null
cannot be manufactured by the construction the way a positive can. The residual live possibility is
that HFpEF's mechanical heterogeneity lives specifically in **wall thickness / concentricity**,
which no non-circular published phenogrouping has yet tested.

## The residual axis, now tested non-circularly — a split decision
A targeted harvest closed the gap the section above named. Two phenogroupings were found whose
clustering inputs contain **no geometry at all**, so a geometric difference between their clusters
cannot be manufactured by the construction. The pivotal one is a POSITIVE, so its input list and
its numbers were re-verified here against primary text — the same check that caught Kyodo
(`anchors/noncircular_phenogroup_verification.json`).

**Mitic 2023 (n=75, echo) — verified non-circular.** Verbatim: *"we applied hierarchical clustering
to the concentrations of six biomarkers (sST2, Gal-3, GDF-15, Syn-1, BNP, and cystatin C)"*. Its
clusters differ in wall and mass but **not** in chamber size:

| | remodeled (n=24) | non-remodeled (n=51) | p |
|---|---|---|---|
| LV mass index, g/m² | 142.88 ± 19.48 | 101.71 ± 19.94 | **<0.001** (1.40×) |
| IV septum, mm | 13.19 ± 1.85 | 11.44 ± 1.20 | **0.027** |
| Posterior wall, mm | 11.38 ± 0.52 | 10.09 ± 1.00 | **0.003** |
| chamber size | 54.13 ± 3.23 | 51.47 ± 5.22 | 0.238 (null) |

Two things the verification adds that the p-values do not. First, a **unit error in the source
table**: its "End-diastolic volume" is given in *millimetres* — 54 mm is a diameter, not a volume —
so that null carries no volume information. Second, four of the six clustering biomarkers are
cardiac **remodeling/fibrosis** markers and the authors name the clusters "Remodeled" and
"Non-remodeled": finding more LV mass in the remodeling cluster is not circular, but it is close to
what the input choice predicts. It is the weakest kind of surprise.

**The vote on the residual axis:**

| quantity | positive | null | detail |
|---|---|---|---|
| LV mass | 1 | 2 | Mitic n=75 p<0.001 vs Tamaki n=211 p=0.098, Fazzini n=591 p=0.093 |
| wall thickness | 1 | 0 | Mitic only; no other non-circular source reports it |
| **relative wall thickness** (concentricity proper) | 0 | 1 | Tamaki **p=0.608** |
| chamber size | 1 | 1 | Tamaki LVEDD p=0.009; Mitic diameter p=0.238 |

**Split decision, weighted toward the nulls by n.** The single positive is the smallest study
(n=75, echo, self-described preliminary) clustering on fibrosis biomarkers; the two larger studies
are null on mass; and the one measure of concentricity proper is the flattest null in the whole
evidence base. Cell 3's negative is not overturned — and the axis is no longer *untestable*, which
is the real change: it is now tested, and the answer is genuinely mixed.

## Cell 4 — what the model PREDICTS for that contrast, and it is not what the name suggests
`scripts/tissuetwin/hfpef_cell4_phenogroup_prediction.py`, evidence
`cell4_phenogroup_prediction_evidence.json`. Mitic's two clusters are the first within-syndrome
geometric contrast in the whole anchor base that is not a restatement of its own clustering, so
cell 1's machinery was run on it unchanged.

The one assumption — the source reports chamber size in millimetres, i.e. a diameter — was refused
rather than converted, and swept instead across its whole defensible range: from *the chambers are
the same size* (taking the p=0.238 null at face value) to *shape-preserving* V ~ d³ = 1.163.

**Two** assumptions of mine are swept, not one — the second was caught on re-reading: the absolute
wall-fraction anchor. The first version placed both clusters around the HFpEF cohort mean
phi = 0.756, which is not even consistent with Mitic's own masses, and g(phi) is nonlinear. It is
now swept across the full measured HFpEF range (phi is index-invariant, so every cohort qualifies —
0.644 to 1.596, not just the four carrying absolute volumes).

| assumption corner | K = kappa(remodeled) / kappa(non-remodeled) |
|---|---|
| chambers equal | 0.951 – **0.972** – **0.997** |
| shape-preserving (V ~ d³) | 0.836 – **0.846** – 0.858 |
| **all 324 combinations** (6 phi anchors × 9 volume ratios × 6 materials) | **0.836 – 0.997** |

**A 1.405× LV mass difference between two non-circular HFpEF phenogroups predicts NO increase in
passive chamber stiffness.** The wall term is weak and points the compliant way, and whatever
chamber-size difference exists points the same way, so the answer never exceeds 1 anywhere in the
swept space. The honest strength of it: *never stiffer* is robust; *meaningfully more compliant* is
not — at the tightest corner (chambers equal, high phi anchor) the difference is only 0.3 %.

A structural check came along for free: K is **exactly** independent of the absolute chamber size
(invariance ratio 1.000000000001), which is what scale-freedom requires and what G4c measured in
cell 1 — the same premise showing up again in a cell that did not assume it. The "remodeled" phenotype, named for its
fibrosis-biomarker profile and carrying 40 % more myocardium, is not predicted to be a stiffer
chamber.

**And a closure test that also demotes this cell, which is the honest part.** At the chambers-equal
corner the contrast is a pure wall-fraction change, so cell 1's separately measured wall exponent
predicts cell 4's ln K analytically. Cell 4 never used that exponent — it ran the full thick-wall
solver over cohort geometries — so agreement is a genuine closure between two independently built
cells. Measured over 36 combinations: **max |residual| 0.00015, which is 0.32 % of the ln K span
itself** (mean 0.00004). It closes.

Which means cell 4 is **not independent corroboration** that the wall lever is weak. It *is* that
lever, applied to a real, non-circular, measured contrast. Its contribution is that the contrast
exists — not that it confirms anything cell 1 did not already say.

**The falsifier is concrete:** an invasive EDPVR measured in a biomarker-defined *remodeled* HFpEF
subgroup that came out stiffer than its non-remodeled counterpart would break this.

## What would decide the stronger question
**Per-patient** end-diastolic volumes and masses from a single HFpEF cohort, rather than cohort
tables — which converts the test from between-cohort spread to within-cohort multimodality, where a
several-diseases signal would actually live. The harvest established that no open dataset currently
delivers it: ACDC, Sunnybrook (its CSV was downloaded and inspected — `PatientID`, `Gender`, `Age`,
`Pathology` only), M&Ms, PhysioNet MIMIC-IV echo (per-case volumes but no disease label and no
mass), Zenodo/Figshare/Dryad/OSF all fail on at least one of open-download, disease-label, or mass.
UK Biobank and BioLINCC/TOPCAT hold the fields behind an application. **This needs an access
request, not more searching.**

**A second decider, now partly delivered:** non-circular phenogroupings exist and give a split
answer (above). What would break the tie is one more of them at n in the hundreds, clustering on
axes that are NOT remodeling biomarkers, reporting **relative wall thickness** — the concentricity
measure — rather than absolute septal thickness.

**And a gap the harvest could not close after a 16-angle search:** reduced-EF wall thickness now has
5 records where it had 0, but the one CMR record reports an IQR and the one SD-bearing record is
echo, so a modality-matched CMR-vs-CMR dispersion comparison is still impossible from published
literature. That is a real property of the literature, not a shortfall of searching.
