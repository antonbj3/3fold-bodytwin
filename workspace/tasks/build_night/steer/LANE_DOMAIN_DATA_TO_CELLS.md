# Steering r4 — same physics, stated as the clinical question it actually is

## What happened in r3
The round ran to 216 650 tokens and was then refused by the provider's biological-risk filter, twice.
No `r3.json`, no new cell. One log out of 355 tonight, so this is isolated and not a systemic block.

**The cause is my own wording and it was also scientifically sloppy.** My earlier steer asked which
parameter change makes "a host surface be attacked". That is not how the question is posed in the
literature and it reads as seeking harm. The real question is **regulatory failure**: complement
over-activation on self surfaces is the mechanism of specific named diseases — paroxysmal nocturnal
haemoglobinuria (loss of CD55 and CD59 anchoring), atypical haemolytic uraemic syndrome (factor H
dysfunction) and C3 glomerulopathy — and under-activation is the mechanism of recurrent encapsulated
bacterial infection in complement deficiency. That framing is accurate, it is the clinical reason the
operator's seed matters, and it is what the cell should compute.

## This round
1. **Restate the cell's readout as a regulatory margin.** How far is a self surface from losing
   protection, in the units of the regulator that protects it: factor H activity, CD55 and CD59 surface
   density. The existing `rho_critical = 0.69012` is that margin seen from the recycling side; express it
   from the regulator side too, because that is the side a measurement can reach.
2. **Name the two clinical directions with their mechanism**, not as an abstraction: reduced regulator
   function on self surfaces (the haemolytic and glomerulopathy direction) and reduced terminal pathway
   function (the recurrent-infection direction). One number each, in the regulator's own unit.
3. **Then the three level-1 cells**, as the previous steering said: 45 of the 98 uncovered domains have a
   literature record with a unit, so pick the three with the tightest validity-range match and report per
   cell the consumed record with its locator, the computed number and the deviation.
4. Keep `FOLD_CD55_ILLUSTRATIVE` marked illustrative; it is not measured and must not carry a
   conclusion.

## Control and falsifier
- **Control:** the domain data's own numbers without our cell. The gain is what the cell computes that
  the data does not state.
- **Falsifier:** if the regulatory margin cannot be expressed in a regulator unit that someone measures
  clinically, the threshold stays uncertified and the cell's output is an information link, not a
  capability. Say that rather than converting it anyway.
- **Forbidden:** quoting `rho_critical` as certified while the source's own uncertainty is unmeasured; a
  fourth cell before the three are checked; reading the filtered records in `LIT_REFS_FILTER.md`.

## Extension 2026-10-04 11:40 — your searches are blocked, work from disk
The provider's content filter responds to your web searches with "flagged for possible biological
risk" and leaves no report. Six rounds have been lost like that over three lanes; one of them burned
129 754 tokens with no outcome. The driver booked them as ready because they had driven long enough, and filled up
on into the same wall. It is now fixed: a round whose log contains the string is counted as aborted.

**Change the operation, not the wording.** Do not reformulate the search to slip through the filter.
Instead, work from the material that is already on the table, and say in the report what you are looking for
NOT could do and what it would have given. A named unreachable source is a result; a
rewritten search that happens to pass is not verifiable.

You have published numbers in hand from previous rounds. Count on them. Do you need a source you can't
pick up, write it as acquisition target with quantity, unit and locator so that someone else can pick it up.
