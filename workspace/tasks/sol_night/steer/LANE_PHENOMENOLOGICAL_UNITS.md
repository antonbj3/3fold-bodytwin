# Styrning r17 — de 62 oklassade, inget annat

## Hindret
Styrningen till r16 var entydig: klassa de 62 UNKNOWN families, for the headline count 9,725 % is a lower
limit and span 9,7–13,9 % is the only honest answer until they're done.16 went instead to one
paket med 50 source-role and factor reference entries, and the span remains untouched.

## The operation this round
1. **Klassa de 62.** FYSISK eller FENOMENOLOGISK with the same criterion as for the others 1 429, eller
   `UNCLASSIFIABLE_WITHOUT` plus exactly what task is missing, three numbers, and the new span.
2. **Rapportera huvudtalet som intervall** i RESULTS.md and in the gate until they 62 One point is that:
   A lower limit applies is overreporting, and that is all I am asking of this round.
3. Om de 62 cannot rate with available information: say it straight out, put them as EN
   acquisition post with what's missing, and close the census section. It's a complete outcome.

## Control and falsifier
- **Kontroll:** the current situation; 9,725 % rapporterat som ett tal med 62 familjer oklassade.
- **Falsifierare:** om klassningen av de 62 moving the main number less than half a percentage point each
  uncertainty is insignificant and the census is complete — rapportera det rent ut.
- **Forbidden:** new packages, new families or new roles before the three numbers exist.

## ADDENDUM — the census missed a tuned constant that a swarm job just caught
A swarm job (Q115 + Q084, verified by me against its `results.json`) compared our epithelial renewal
front velocity to a held-out external measurement and reports `refutes_us: true`:

| | value |
|---|---|
| measured mouse duodenal front velocity | **9.00 ± 0.465 µm/h** control, **−0.798 ± 0.942 µm/h** during 0–10 h after Ara-C, 8.68 ± 0.655 beyond |
| our Q115 value | **30.0 µm/h**, i.e. `q115_over_anchor_ratio = 3.3333` |
| how our value arose | the job quotes the source: eta is *"chosen to give 0.03 mm/h"* |
| the test | parameter-free; nothing fitted, every anchor numeral held out |

Two things follow for this lane, and they matter more than the remaining 62 unknowns.

1. **A constant tuned to hit a chosen output is the purest case of a phenomenological debt, and the
   census did not flag it.** Find out why. Does the classifier only look at declared units and
   provenance fields, and not at whether a value was back-solved from a target? If so, that is a
   detectable pattern — a constant whose comment or docstring names the output it was chosen to
   produce — and it should become a census class of its own: `TUNED_TO_TARGET`.
2. **Sweep for the pattern across all 45 cells** and report the count. Phrases of the form "chosen to
   give", "set so that", "calibrated to match", "to reproduce" next to an assignment are the signature.
   Report how many constants carry it, in how many cells, and for each whether a reported output is
   sensitive to it. A tuned constant that nothing depends on is a label; one that carries an output is a
   debt that makes the output circular.

That second number is the one I want: **how many of our reported outputs depend on a constant that was
chosen to produce a desired value.** It bounds how much of the model is self-referential, which no other
measurement tonight does. Still deliver the 62 classifications, but this comes first if you must choose.

Control: the current census, which reports 9.725 % phenomenological and did not contain this case.
Falsifier: if no cell outside Q115 carries the pattern, it is a single defect and not a class — say so
plainly and the sweep is cheap insurance.
