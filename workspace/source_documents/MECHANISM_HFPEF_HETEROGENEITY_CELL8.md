# HFpEF cell 8 — where the heterogeneity lives

Cell: `data/hfpef_diastolic_v1/`. Pre-registration: `PREREG_CELL8.md`, committed before the run.
Code: `scripts/tissuetwin/hfpef_cell8_heterogeneity_location.py`. Evidence: `cell8_evidence.json`.
Upstream: cell 3 (`MECHANISM_HFPEF_KNIFE_EDGE_CELL3.md`), the within-cohort dispersion falsifier.

## The question cell 3 could not reach

Cell 3 asked whether HFpEF is mechanically over-dispersed and answered **no** — from the spread
*between* cohort means, in geometry. Harvest pass 4 then supplied a *within*-cohort dispersion for
the EDPVR stiffness exponent. Cell 3 already computed within-cohort CVs for geometry. Both halves of
a comparison were sitting on disk and had never been put against each other.

    R = CV_within(stiffness) / CV_within(geometry),  computed per arm

## The result

| stiffness source | geometry set | R disease | R control | ratio of R |
|---|---|---|---|---|
| Lam, noninvasive | indexed-only | 2.297 | 0.052 | **44.17** |
| Lam, noninvasive | mixed indexing | 2.056 | 0.053 | 39.05 |
| Popovic, hybrid invasive | indexed-only | 0.315 | 0.058 | 5.40 |
| Popovic, hybrid invasive | mixed indexing | 0.282 | 0.059 | **4.77** ← conservative |

**All four readings give the same verdict** — `HETEROGENEITY_IS_CONCENTRATED_IN_THE_MATERIAL_PARAMETER`
— against a pre-registered threshold of 2.0. The cell quotes the **conservative 4.77**, chosen by
rule (the smallest ratio across stiffness sources and geometry sets) rather than by taste, because
the noninvasive method is the one the replication falsifier showed diverging on the disease side.

Underneath, the raw dispersions:

| arm | CV stiffness (noninvasive / invasive) | CV geometry (median) |
|---|---|---|
| controls | 0.0101 / 0.0113 | 0.1936 |
| HFpEF | 0.5007 / 0.0687 | 0.2180 |

Read plainly: **in health the stiffness exponent is roughly twenty times more tightly determined
than the chamber's geometry; in HFpEF that inversion collapses.** The quantity that goes loose in
disease is the material one, by between 4.8× and 44× depending on the method used to measure the
stiffness.

This is the project's central finding — that the stiffening is material rather than geometric —
restated in a second and independent dimension. Cells 1, 2, 5, 6 and 7 reached it through the
*magnitude* of the residual. This reaches it through the *variability*, and the two were derived
from different anchors.

## Why the comparison is a ratio of ratios

The obvious objection is that a sicker cohort measured with a noisier instrument shows larger
dispersion in everything. That objection is pre-registered as G3 and it is answered by the form of
the statistic: a common multiplicative noise factor `k` in the diseased arm appears in **both**
CV_stiffness and CV_geometry and **cancels exactly in R**. It cannot produce a ratio-of-R different
from 1. What would survive the null is a factor acting on one quantity and not the other — an EDPVR
fit that degrades in sicker patients while their imaging does not. That remains possible and is not
excluded by anything here.

## ★ Replication: the direction holds, the magnitude does not

Cell 8's G3 named the one thing that would survive its noise null — *a factor acting on one quantity
and not the other* — and did not exclude it. The sharpest form of that objection is that the control
arm's 1.0 % CV is not biology but an artefact of one fitting method, which would make R_control
meaningless and the 44× with it. A second study tests it, and it was already on disk: **Popovic**
(PMID 36994635) fits the same power-law exponent by a **different method** — hybrid invasive (PCWP by
micromanometer catheter plus echo volume) against Lam's fully noninvasive echo-derived LVEDP. Its
IQRs had never been extracted; they are now, from PMC10330082 Table 3.

| | Lam (noninvasive, SD) | Popovic (hybrid invasive, IQR/1.349) |
|---|---|---|
| controls | 1.01 % | **1.13 %** |
| HFpEF | 50.07 % | **6.87 % / 7.19 %** |
| disease/control CV ratio | **49.74** | **6.08 / 6.36** |

**The leg that could have killed cell 8 failed.** Two independent methods reproduce the ~1 % control
dispersion, agreeing within **1.12×** — a fitted exponent really is that tightly determined in a
healthy chamber, and it is not one lab's artefact.

**But the magnitude does not replicate: 49.74 against 6.08, a disagreement of 8.2×.** Both still
cross cell 8's pre-registered 2× line, so the verdict
`HETEROGENEITY_IS_CONCENTRATED_IN_THE_MATERIAL_PARAMETER` stands in a second study and a second
method. The *size* of the effect does not.

And the divergence sits exactly where the surviving hypothesis predicts it: the **control** side
agrees across methods while the **disease** side does not, which is the shape a noninvasive LVEDP
surrogate degrading in sicker patients would produce. That is not proof it is happening — but cell 8
named this hole, and the replication attempt found the data leaning into it rather than away.

**The figure that should be quoted is the conservative one.** In this falsifier's own quantity —
the disease/control CV ratio — that is 6.08 against Lam's 49.74. Carried through into cell 8's
`ratio_of_R` it is **4.77 against 44.17**. Cell 8 now computes all four readings itself and selects
the smallest by rule, so no prose has to reconcile two numbers by hand.

### The conservative arm is internally coherent

Promoting Popovic's reading to the headline puts weight on three numbers from one table. They can be
checked against nothing but themselves. Under `P = α·V^β`, doubling the pressure fixes the volume
ratio exactly:

    30/15 = (V30/V15)^β    ⟹    V30/V15 = 2^(1/β)

β, V15 and V30 are reported separately and were not fitted to this equation jointly. α never enters.

| group | β | V30/V15 observed | 2^(1/β) predicted | error |
|---|---|---|---|---|
| controls | 5.90 | 1.1273 | 1.1247 | **+0.23 %** |
| HFpEF 50–64 % | 6.08 | 1.1327 | 1.1208 | +1.07 % |
| HFpEF ≥65 % | 6.15 | 1.1275 | 1.1193 | +0.73 % |

Worst error 1.07 %, against a declared 2 % tolerance — which is what rounding alone produces, since
the volumes are reported to whole millilitres and 1 mL in 102 is already ~1 %. The table is coherent
with the functional form it declares. This does **not** validate the single-beat Klotz assumption
behind those volumes; it validates that the numbers belong to one fit.

## What this does not establish

- **n_eff = 2** for the stiffness half after the replication above — Lam (PMID 17404159) and
  Popovic (PMID 36994635), which share a Borlaug author lineage, so the two are not fully
  independent. It was n_eff = 1 when this cell was written and the verdict text said so.
- The two halves come from **different studies**. No paper on file reports both a stiffness and a
  geometry dispersion, so the arms are matched across studies — the cell's central weakness, stated
  in G1 rather than smoothed over.
- The stiffness dispersion is **noninvasive** (echo-derived LVEDP) and a **power-law** exponent. No
  ratio is formed between it and the exponential-family β cell 1 uses for K_meas; the model-form
  scan exists because those two were once pooled.
- The dispersion statistic was verified at source: PMC2001291 Table 2's footnote states
  *"Data are mean ± SD"*, which closes the standard-error reading that would have given 31.28.

## What would decide it

One study reporting, for the same patients, a dispersion in the EDPVR stiffness parameter **and** a
dispersion in LV volume and mass. That single design collapses G1's weakness and takes n_eff off 1.
