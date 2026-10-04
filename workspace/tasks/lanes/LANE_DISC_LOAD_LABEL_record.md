# LANE_DISC_LOAD_LABEL — is the disc's 2,09× a model error or a load label error?

Results directory `results/LANE_DISC_LOAD_LABEL/`. Decision **outside the eye**.

## The measured state, calculated by the coordinator 2026-10-04
`tasks/assembly/disc_load_decision.py` records that our pressure multiplier `k = 1,3–1,5` is
**1,809 to 2,560×** above the multiplier implied by the in-vivo reference (`0,585959` to `0,718629`),
and the source `MECHANISM_INTERVERTEBRAL_DISC.md` itself calls it FAIL.

The numbers behind it, each recalculated:

- disc area `A = 1800 mm²`, nominal stress `σ = 0,9045 MPa`
- the model's **reference load** is therefore `F = σ·A =` **1628,1 N**
- reference: Wilke, in-vivo nucleus pressure during walking, `0,53–0,65 MPa`
- `k_implied = P/σ` gives 0,585959 and 0,718629

## The obstacle, exactly
The comparison assumes the reference is measured at the same load as the model's reference. That is not shown.
I calculated which load would make our own `k` compatible with the measured pressure:

| measured pressure | k | load required | corresponding σ |
|---|---|---|---|
| 0,53 MPa | 1,5 | **636,0 N** | 0,3533 MPa |
| 0,53 MPa | 1,3 | 733,8 N | 0,4077 MPa |
| 0,65 MPa | 1,5 | 780,0 N | 0,4333 MPa |
| 0,65 MPa | 1,3 | **900,0 N** | 0,5000 MPa |

The range **636–900 N** is not unreasonable for walking, while 1628 N is more associated with lifting. If Wilke's
walking pressure was measured under a load in that range **there is no model error at all** — the entire 2,09× is then a
comparison made at the wrong load, the same disease as the meniscus load labels (see
`LANE_MENISCUS_LOAD_LABEL`, where the products span a factor 26).

## The operation
Find out which load Wilke's walking measurement actually corresponds to, from the source itself — body weight,
the subject's mass, segment level and how the load across L4/L5 was calculated. Give the load in newtons with
the source's line. Then decide which of the two statements stands: *the model overestimates pressure* or
*the comparison was made at the wrong load*. Update `disc_load_decision.py` so that the load is an **input**
and never a constant hidden in a nominal stress.

## Strongest control
Same geometry, same reference, but the comparison made at 1628,1 N in both terms — that means looking up
a published nucleus pressure at a load near 1628 N instead of during walking. If no such measurement
exists, that is **the result**, and the edge must then say that reference and model are not comparable, not that
the model is wrong.

## Falsifier
If Wilke's walking load is above 1200 N, the FAIL verdict holds and the multiplier must be recalibrated.
If it is below 950 N, FAIL is invalid as worded and must be withdrawn.

## Regler
Status `PENDING_INDEPENDENT_REVIEW`, ingen utsaga om biologisk validering, inget excluded_category material,
Nothing that names the collaborative. Each published source with PMID eller DOI, volume and pages.
