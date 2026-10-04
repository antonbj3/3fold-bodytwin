# LANE_COARSEST_SUFFICIENT — what is the COARSEST resolution that preserves the answer?

Resultatmapp `results/LANE_COARSEST_SUFFICIENT/`.

## The finding that makes the lane necessary, measured three times independently last night
Three lanes in three different areas of physics measured the same thing without knowing about each other:

| var | grov | fin | faktor |
|---|---|---|---|
| heat dose in tissue | 48,87 % error | 0,19 % at 10 µm | **254×** |
| transport across the blood brain barrier | 22,90 % misclumped | 0,126 % axially resolved, 256 segments | **181×** |
| structural sag | 0,0305 mm | 0,1346 mm over mesh refinement | **4,4×, and VERDIKTET TURN** |

The third is the most serious: at the same material density, the coarsest resolution passes one
deformation limit at which the two finer ones fall. A design exclusion with no specified resolution is
thus not an exclusion. And the fleet line then landed on a 90 percent cap **0,14 %
over** the limit, i.e. a margin fifty times smaller than the resolution spread — undecidable.

**Clumping is therefore not an approximation with a known cost with us. It has twenty percent error which
no one has declared.**

## What the lane is, and what it is NOT
This is operator adaptive resolution made concrete. He said: you should be able to choose resolution per
variable, and correct resolution of the body twin is at least as important as the robot side. The question is
not if resolution matters — it's measured three times. The question is **which resolution is the coarsest
sufficient for every greatness**, for to refine everything everywhere is both unnecessary and impossible.

And the machinery is already graded in the older project, so this is not a new build:
- `source_repository/scripts/physics_exp/decidability_abstention_atlas.py` with siblings —
  classifies a quantity as precision-limited, cross-quantity floored, or divergent. **Verified
  transferable to biology tonight:** it gave KLASS 2 for both bleeding rate and collagen ratio,
  the same answer as two hand counts, and it reproduced its own five programs exactly.
- `docs/inherited_memory/repr-router-compute-value-is-amr-sharpness-scaling-not-breakthrough.md` —
  the refinement value scales with the singularity strength: smooth ~1×, returning corner ~11×, crack
  ~1005× in nodes to the same accuracy. **And the lesson that applies directly to this lane: assess
  never an adaptive method on a happy test case.**
- `docs/inherited_memory/graded-abstention-paydownable-cliff-fundamental-action.md` — smoothness and
  repayment is TWO OBEROENDE axes, and gradual-but-fundamental exists.
- `docs/inherited_memory/dataset-mean-fidelity-is-area-diluted-certify-per-frame-worst-region.md` —
  an average value is diluted over area; certify the worst-case region and report the blindness rate.
Replicate lane-locally, don't write in the old project.

## Do like this
1. **Start with the three measured cases as calibration.** Your method must reproduce that heat dose
   requires 10 µm, that the transport requires axial resolution, and that the deflection vertices reverse. If
   if it doesn't, the method is wrong and you should discover it on the three and not on the others.
2. **For each quantity: find the coarsest resolution where the answer does not change.** The answer, not the number.
   A quantity whose VERDIKT is invariant under refinement needs no finer twin; one whose
   verdict vender needs it and cannot be certified without it. Report by magnitude: coarsest
   sufficient resolution, what determines sufficiency, and the order of convergence.
3. **The order of convergence is the diagnosis, not a detail.** A quantity that converges more slowly than
   its expected order carries a singularity, and then **diverges maximally under refinement while
   percentiles converge** — meaning the envelope must be defined on a percentile and not
   at a maximum. That's exactly the flotilla open question, so answer it while you're here.
4. **Separate the two axes** according to COMMON.md: smoothness (gradual to steepness) and amortization (finite
   random versus structural). Gradually ⇒ pay down with resolution. Dip ⇒ change reading, to
   refine a precipice never certifies.
5. **Measure the distribution before any system statement is made.** How many magnitudes has a coarse sufficient
   resolution that is coarser than the one we use, how many require finer, and how many are
   undecidable. Three numbers. It is the answer to the operator's question and it is not a number but a distribution.
6. **And connect to the cost.** A magnitude that requires 256 segment where we drive lumpy is a
   calculation debt with known size. Enter it as a factor against current cost, so the priority
   becomes possible.

## Strongest control and falsifier
- **Control:** the resolution each cell runs at today. The profit is the number of magnitudes where we drive for
  coarse for our own value, not the number of magnitudes improved by refinement — everything is improved
  of refinement and it is not interesting.
- **Falsifier:** if almost all quantities turn out to be invariant under refinement at current
  resolution, the three measured cases are exceptions and clumping is defensible. That would be reassuring
  results and should be reported outright, not downplayed.
- **Forbidden:** to evaluate an adaptive method on a happy test case; to report an average value that
  proof of sufficiency without the worst-case region; to write in `~/projects/bodytwin`.

## Delivery
`PORT.json`: the three distribution numbers, one row per quantity with coarsest sufficient resolution,
order of convergence, the classification of the two axes and the calculation cost factor. Plus the calibration line:
three out of three measured cases reproduced, or which missed.

Inga interna data. Allt PENDING_INDEPENDENT_REVIEW.
