# the proof lane: what precision does a number coming out of our operator chain carry?

## The measured problem

A swarm job took the number `2.220261437908497` from our own network and showed that it is exactly
**3397/1530**, hence the ratio of two four-digit readings. I verified:
`Fraction(x).limit_denominator(10**4)` gives 3397/1530 and `float(3397/1530) == x` is true. The number carries
four significant digits; the other eleven are floating-point representation.

The same error occurs in my own work three times today:
- a residual reported as **−0,000000** when the inputs carried 3e-5,
- a certificate interval quoted with ten digits that was **1,81 ULP in float32** wide,
- a clock bracket of 0,00576 s that was exactly **1509,7/2¹⁸**, hence the resolution of bisection.

And a Sol lane has tried to build a discriminator three times: the first filter caught 24 026 numbers that
were almost all 1/6, 4/3 and 8/3 in long format; the second caught 1/86400 and 1/60000, hence
unit conversions. Mechanically DISTINGUISHING the three cases has proved harder than expected.

## The question is propagation rather than labelling

1. **Formulate the propagation rule for our actual operator class.** The chains perform: sums of measured
   terms, ratios of measured terms, products with unit constants that are exact by definition,
   exponentiation with integer exponents, bisection to a tolerance, and linear solutions of small
   systems. What is the correct rule for how many significant digits the result carries, per operation? The
   classical digit rule is a rule of thumb and is wrong for ratios near one and sums with
   cancellation. Give the exact rule and state where the rule of thumb breaks.

2. **Distinguish three cases exactly.** (a) two short measurements divided, dangerous; (b) a simple fraction
   or a unit conversion, harmless; (c) an actual long computation or an irrational number,
   harmless. Is there a condition that determines which case a given double belongs to, without knowing
   the provenance? If the answer is no — which I suspect — state it as an impossibility theorem and specify
   the minimum provenance information that suffices.

3. **Cancellation and precision together.** Our toric margin rests on correlation −0,796 where
   cancellation removes 51,8 % of the error. A difference of nearly equal numbers loses significant digits. How
   many digits does a cancellation-dependent margin carry, and is that an argument AGAINST reporting
   it as a certificate?

4. **State the cheapest labelling.** We have tens of thousands of numbers. What minimum metadata per number makes
   precision computable afterwards — and which kind of number can never be rescued without
   rerunning the computation?

## Requirements

- Control: current practice, printing repr() of a double.
- Falsifier before each derivation.
- An impossibility theorem is a stronger result than a weak test. Say what cannot be done.
- PENDING_INDEPENDENT_REVIEW.
