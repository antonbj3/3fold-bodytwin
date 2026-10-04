# the proof lane: what does rank two mean, and what does rank three cost?

## What appeared in a correction

An the proof lane round corrected my claim about a calibration number. I had written that "158/415 = 38,07 %
of the calibration coverage is missing". The correct, according to the correction: **158/415 = 0,380723 is a
unsupported dimensionless data value gap during a rank two calibration**, and the corresponding
the stored quantity gap is **12640/83 = 152,2892 nmol**.

The same round reported for another chain: **exakt uN = 0** for the task at hand but
**AKN = (156000, −156000)ᵀ s⁻¹ ≠ 0**, and that the **dynamic rank grows from 2 to 3**.

And a third chain carries `conditional95_95_independent_run_count = 59` against 11 runs that exist,
as well as a cover from the triangle inequality covering 448 by 448 but whose center lies in 0 by 448 —
the envelope is 2,54× to 9,37× wider than the worst fault it is supposed to enclose.

## The question

Rank two returns. It is not a number but a structural property, and it determines what the calibration is
can carry.

1. **Why rank two?** What property of the calibration makes its rank exactly two, and what are the two
   the directions? State them in terms of the quantities being measured, not as abstract vectors.

2. **What does rank three cost?** The dynamic rank grows to 3 as soon as the task becomes dynamic
   (AKN ≠ 0 while uN = 0). How many independent measurements are required to support rank three, and which ones?
   This is the practical question: if rank three requires ONE new observable, we want to know which one.

3. **Are 0,380723 and 152,2892 nmol the same thing in two entities?** About the dimensionless gap and the
   stored quantity gap is the same deficiency expressed twice, only one of them should be reported. Decide it,
   and enter the conversion factor if they are the same.

4. **Connect the rank to the envelope width.** An envelope from the triangle inequality uses only the amplitudes
   and nothing about the angle. Is it a rank deficiency of the same kind? If the case is loose to
   the calibration is rank two, does rank three make it narrower — and by how much, in the same unit?

5. **59 runs against 11.** Is that requirement a consequence of the rank or of the dispersion? If the rank is
   the binding factor 59 may be the wrong number.

## Krav

- Control: today's reading, to treat the calibration as sufficient for the recorded task,
  which it demonstrably is (uN = 0 exactly).
- Falsifier before each derivation.
- Enity in each field name. State when two numbers are the same quantity in different units — that error has been made
  six times in the project today.
- Say what cannot be determined without the calibration matrix itself.
- PENDING_INDEPENDENT_REVIEW.
