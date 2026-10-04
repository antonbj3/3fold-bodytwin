# PROOF_LANE_THREE_ROUTES — can three incompatible values of the same coefficient all be right?

## The measured state, computed by the coordinator
The nucleus pressure in the disc is written `P = k·σ` with `σ = F/A`, `A = 1800 mm²`. Three independent routes give `k`:

| route | k | from where |
|---|---|---|
| reference implies | **0,585959 – 0,718629** | `k = P/σ` from Wilke's in-vivo 0,53–0,65 MPa against σ = 0,9045 MPa |
| our model | **1,3 – 1,5** | empirical multiplier in `disc_load_decision.py` |
| published fiber equation | **4,877250096040** | `k = 1/g` with `g = 0,205033570210` at η = 0,16, γ² = 2, φ = 30° |

Spread across the three: **8,324×**. I checked the fiber identity myself: `1/0,205033570210 =
4,87725009604904`, that is, `k = 1/g` exactly. And it is brittle — the coefficient crosses **zero** at
`φ = 34,58743251190°`, only **4,5874°** from the paper's assumed 30°, and at 25° k is 2,686523862844.

Two further measured facts from the same week, both verified by me:

- **Force magnitude does not determine pressure.** Of 1440 pairs in the **same specimen** with identical axial
  preload and identical moment magnitude, 1420 differ by more than 0,002 MPa; largest gap
  **1,270237 MPa**, donor 12, T3–T4, extension versus flexion at ±7,5 N·m. It is 10,59× the width of
  the entire reference band 0,53–0,65 MPa.
- **Charge density is the missing quantity.** At the same force, different fixed charge densities give a
  pressure gap of **0,32164781492307004 MPa**, that is, 2,6804× the reference band's width. And the conditional
  pressure interval at σ·A = **1628,1 N** is `[0,2508975219612329; 0,6904589796073666] MPa`, which
  **contains** the reference band but is **3,663×** wider than it.

## The question
Is there a reading under which all three `k` values are **simultaneously true** — that is, where they measure
different quantities that just happen to bear the same name — or is at least one of them testably wrong?

I want this as a theorem, not as an argument. Specifically:

- what quantity each of the three `k` values is actually a ratio of, written out, so that the difference appears
  in the units or in the integrand and not in the prose;
- whether `k` can be a constant at all, given that `A` is a **total** area while pressure is measured
  **locally** — the same structural error that another edge in the net just failed on, where `κ = p₀A/W ≥ 1` always
  holds for a peak pressure and our four rows were at 0,1155–0,168;
- which of the three becomes identifiable **first** when charge density is included as an input, and
  which remains unidentifiable;
- and whether the fiber equation's zero crossing at 34,5874° makes it unusable as a calibration over
  the **entire** interval or only near the singularity — state the boundary in degrees.

## The strongest control
Setting `k = 1` directly, that is, pressure **is** the nominal stress. It is the
choice with the fewest assumptions and also lies between the reference and our model. If none of the three
routes beats `k = 1` against the reference, the entire multiplier is a fiction and that must be said.

## Falsifier
Construct a case yourself where two of the three routes provably measure the same quantity and still give
values that differ by more than the sum of their spreads. If you succeed, the "different quantities" reading is
ruled out and at least one value is wrong.

## Warning
A width-based rounding rule fails — your own `PROOF_LANE_PRECISION_PROP` let one fail on
all 1 200 constructed boundary cases. Use the endpoint rule when stating precision.

## Regler
`PENDING_INDEPENDENT_REVIEW`, ingen utsaga om biologisk validering, inget excluded_category material, inget som
names the collaboration. Each published source with PMID eller DOI, volume and pages.
