# the proof lane: is there a test ON the summary, instead of case by case?

## The measured pattern, five chains, all with zero identity error

Five independent chains have shown the same thing: two states with **exactly identical summary** lead
to different downstream decisions. Each time it has been found by hand, one case at a time.

1. **Laser.** Weight pair q = (1/2, −1/2) summing to zero, summary S_N = Y = 0 exactly, and
   source emission **0 versus 4/3 Pa·s** with the declared threshold **2/3** exactly halfway between. The aggregate
   strictly underestimates in **14 412 of 15 413** instances, 93,51 %.
2. **External solver.** `exposure_kernel_identity_error = 0` with downstream gap **−0,0444081 strain**,
   i.e. **2,740e+05 × the reserve** 1,6204384e-07.
3. **Summary repair.** Registered task repaired to exactly 0, but
   `partial_future_count_gap_nmol = 39`, and calibration covers 158/415 = 38,07 % less than everything.
   Cooperative: **107 of 1350 rows** change verdict between two readings, 7,93 %, and the transport passes
   are a strict subset of the paired ones (229 − 122 = 107 exactly).
4. **Cells.** Everything quantized at one molecule: 0,1037837 µM per molecule in 0,016 µm³. The target gap is
   exactly **2,00000 molecules** and the smallest distinguishable contrast at one molecule’s error is **2**. Zero
   margin.
5. **The eye.** Identical mean thickness, index, hydration and surface allocation give **0,3056 D** in 4 of 4 pairs;
   rotationally symmetric control gives −1,36e-08 D.

## The question, and it is a theoretical question

Each time we have discovered insufficiency by finding a witness. That does not scale: 43
validation cells exist, 21 are accepted on a scalar test, and 35 are untested.

**Is there a test on the SUMMARY itself — a property of the mapping from state to
summary — that decides whether it can carry a given decision, without searching for witnesses?**

Four tasks:

1. **Formulate the condition exactly.** Given a state set X, a summary σ: X → S and a
   decision d: X → {0,1}, when does d factor through σ? The trivial answer is that d must be constant
   on σ’s fibers. Make it testable: which property of σ and d can be CHECKED without
   enumerating the fibers?

2. **Quantify gradual insufficiency.** Our cases are not binary. Laser gives 93,51 % strict
   underestimation, the cooperative case 7,93 % reversed rows, the cells zero margin exactly at
   the resolution limit. Is there a measure of HOW insufficient a summary is that is comparable
   between chains? Our five cases have different units and different decisions.

3. **Distinguish two causes.** A summary may be insufficient because it loses information
   (the fiber is too coarse) or because the decision lies at a threshold where the fiber happens to cross. Laser
   has the threshold exactly halfway between the two outcomes, the cells have the gap exactly equal to resolution.
   Is this the same phenomenon or two?

4. **State the cheapest proof.** For a given summary and decision: what is the minimum
   amount of work that decides sufficiency? One witness suffices to reject, but what is needed to
   SUPPORT it? We have never done the latter.

## Requirements

- Control: current practice, i.e. looking for witnesses case by case. Each proposed test is counted against
  how many of our five cases it would have caught in advance.
- Falsifier written out before each derivation.
- Say what cannot be decided. A limitation is a result, and an impossibility theorem is a
  stronger result than a weak test.
- Everything quantitative should carry unit and convention in the name.
- PENDING_INDEPENDENT_REVIEW. No tissue validation is claimed.
