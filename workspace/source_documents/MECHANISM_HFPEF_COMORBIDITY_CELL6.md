# HFpEF cell 6 — the comorbidity drivers as quantified exposures

Cell: `data/hfpef_diastolic_v1/`. Pre-registration: `PREREG_CELL6.md`, committed before any cell-6
anchor value was in hand. Code: `scripts/tissuetwin/hfpef_cell6_exposures.py`.
Evidence: `cell6_evidence.json`. Anchors: `anchors/cell6_exposure_slopes.json` (21 slope records).
Upstream: cell 1 (`MECHANISM_HFPEF_DIASTOLIC.md`).

This cell had no document until now; its result lived in the claims ledger alone.

## The question

Cell 1 measured that passive chamber stiffness has two geometric levers and only two — cavity size at
an exponent of exactly −1, and wall fraction at a weak, compliant-direction exponent. The dominant
comorbidities act on exactly those levers: obesity through body size and therefore cavity volume,
hypertension through concentric hypertrophy and therefore wall fraction. So the pre-registered
question was not whether these exposures associate with HFpEF — they plainly do — but whether they
can reach the chamber **through the levers that exist**.

    d ln kappa / d(exposure) = -1 x d ln V/d(exposure) + (wall exponent) x d ln phi/d(exposure)

Nothing is fitted: cell 1's exponents are inputs, the exposure-to-geometry slopes are anchored, and
the product is a prediction.

## Verdicts

| gate | verdict |
|---|---|
| G0 EF-blind | PASS |
| G1 exposures, not labels | PASS |
| G2 direction test | MIXED_OR_STIFFENING |
| G2b lever ratio from exposure slopes | **SIZE_DOMINATES** |
| G3 magnitude vs the residual | **GEOMETRIC_ROUTE_IS_INACTIVE** |
| G3b slope-source falsifier | INACTIVE_UNDER_EVERY_SLOPE_SOURCE |

**The geometric route delivers 6.5 % of the measured log-stiffening, and points the compliant way.**
Against a pre-registered 20 % activity line, that is inactive with the sign wrong as well as the
magnitude.

## ★ The share now carries its slopes' uncertainty

The 6.5 % was computed from **one pair of slopes** as point estimates — Petersen's BMI-to-LV-mass and
BMI-to-EDV — and both state a 95 % CI in the anchor that had never been propagated. The
interval-propagation falsifier made this correction for cell 1; nobody had made it here.

    LV mass  +8.3 % per SD BMI   95% CI 7.6-8.9
    EDV      +4.8 % per SD BMI   95% CI 4.2-5.4

Walked over all eight CI corners — the two slopes move `d ln phi = d ln Vw − d ln V` in **opposite**
directions, so the corners must be walked rather than the midpoints combined:

| | share of the measured log-stiffening |
|---|---|
| point estimate | 6.54 % |
| across the CI corners | **5.87 – 7.77 %** |
| activity line | 20 % |

**2.57× of headroom at the worst corner.** The conclusion is robust to the uncertainty its own
sources state.

One correction to what the older prose implied: the sign is **not** opposite at every corner
(`all_corners_same_sign` is False, meaning some corners do share the residual's sign). The magnitude
result stands on its own; the sign is a corner-dependent detail and is recorded as such.

## What this does not establish

- The share rests on **one pair of slopes from one study**. Propagating Petersen's stated CI
  propagates Petersen's uncertainty and nothing wider — G3b's slope-source falsifier is the check
  that the verdict survives other sources, and it does.
- The exposure **delta** behind the share is a cross-cohort splice: an unweighted mean BMI of
  **30.89** across three HFpEF arms (Dattani 33.84, Liu 25.00) against **24.42** across two healthy
  arms, a delta of **6.47 kg/m²**. Cell 6 declares that splice separately rather than folding it
  into the interval.
- 18 of 21 slope records carry a usable CI or SE. Only the two decisive ones are propagated here;
  the rest remain a worklist.
