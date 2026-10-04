# Steering r8 — resolve the 20 unresolved labels; the sample size is not what limits you

## What r7 settled, and what it did not
You ran the random sample I asked for, with an exact method and no missing-at-random assumption, and the
result is more careful than the headline:

| sample | fail | hold | unresolved | fail share of testable | 95 % finite-frame interval on the population |
|---|---|---|---|---|---|
| adversarial (unchanged) | 45 | 6 | 1 of 52 | **0.882** | — |
| random audit, n = 30 of M = 260 | 20 | 4 | **6** | **0.833** | **[0.485, 0.958]** |
| random source-native, n = 30 | 12 | 4 | **14** | **0.750** | **[0.238, 0.958]** |

So the point estimates under random sampling, 0.833 and 0.750, are close to the adversarial 0.882 — the
phenomenon is clearly not an artefact of choosing suspicious contracts. **But `majority_below_half`
is UNKNOWN for a real reason: the lower end of the interval is 0.485 and 0.238, so the data does not yet
exclude that fewer than half of all summaries are insufficient.** My falsifier was that a random share
under half would turn the routine requirement into a targeted one. That question is still open, and it
is the only thing standing between us and a certified system statement.

## This round: the unresolved labels, nothing else
**The interval width is driven by the unresolved union, not by n.** Six unresolved out of 30 stretch the
audit interval from a point range of [0.667, 0.867] to a population interval of [0.485, 0.958]. Another
thirty samples would add almost nothing; resolving the 6 and the 14 would collapse the interval.

1. **Resolve the 6 unresolved labels in the audit sample.** For each: is the summary sufficient or not,
   with the identity error and the downstream difference, or state exactly what makes it unresolvable.
2. **Then the 14 in the source-native sample**, same requirement. Twenty decisions in total.
3. **Recompute both intervals** with the same exact hypergeometric inversion and report whether
   `majority_below_half` is now resolved. That is the deliverable: one word, FALSE or still UNKNOWN,
   with the new interval.
4. **If a label genuinely cannot be resolved with available information**, say which information is
   missing and write it as an orderable item. An unresolvable label is not a failure; an unresolvable
   label reported as a fail or a hold would be.

## Control and falsifier
- **Control:** the adversarial distribution, unchanged at 0.882 over 52 contracts. The random sample is
  the weaker-informed comparison, which is the right direction for this claim class.
- **Falsifier:** if resolving the twenty leaves the interval still spanning 0.5, the population question
  cannot be answered from a frame of 260 and the honest move is to state the sample-level result only and
  stop claiming prevalence at all. Say that plainly rather than widening the frame.
- **Forbidden:** drawing a new sample before the twenty are resolved; assuming the unresolved are missing
  at random — your own method note rules that out and it is the reason these numbers are trustworthy;
  reporting 0.833 without its interval.
