# LANE_IMMUNE_OVERACTIVATION — the second failure mode: too much, not just too little

Results directory `results/LANE_IMMUNE_OVERACTIVATION/`. Decisions outside the eye.

## Measured state
Three decisions are built and executable, and I ran all three today:

- `complement_amplification_decision.py`: convertase lifetime **129,84 s** surface-bound and **111,09 s** in
  fluid phase (half-lives 90 and 77 s), C3 cleaved per convertase **231,1** at saturating C3 and
  **115,6** at Km.
- `complement_window_decision.py`: the host cell reverses at recycling fraction
  **0,6901150139456149**; the activator **never** reverses across the entire interval. Growth just below
  the turning point: host −0,00018768 per minute, activator +0,21710422 per minute. One sign change for
  the host over 2001 points, zero for the activator.
- `immune_window_decision.py`: the host is spared for ρ in **[0; 0,69]**, thus 691 of 1001 points; with
  two simultaneous requirements this shrinks to **510 of 1001**, the interval [0,181; 0,69].

The sources are locators in the code: PMID 6229580, doi `10.1016/j.molimm.2007.11.003`,
doi `10.1182/blood-2005-02-0782`, doi `10.1182/blood-2008-11-189944`.

## The obstacle, in two parts
1. **The edges are not anchored.** Six of the network's 163 edges concern the immune system, two are
   `TIGHT` — `T-E14` and `T-E15` — and **none** of the six carries an `external_reference` field. The numbers
   thus rest in code while the edges claim to be tight without pointing to a measurement.
2. **Only one failure mode is tested.** All three decisions ask whether the host is **spared**. Underactivation
   is the other half: a system that spares the host perfectly while allowing the activator to grow by
   0,217 per minute is not regulated, it is switched off. And overactivation — host damage — is the mode
   that actually causes clinical harm, and it is not scored at all.

## The operation
Make the window two-sided. A decision that returns **both** boundaries: the recycling fraction below
which the activator is no longer lysed (underactivation) and above which the host is damaged
(overactivation), in the same unit, against the same sources. Return the boundaries, not their description.

Then anchor `T-E14` and `T-E15`: set `external_reference` on each with whichever of the four
locators actually carries the number, and if none does, the status must drop from `TIGHT`.

## Strongest control
The current one-sided reading, thus only "is the host spared". If the two-sided window is not
narrower than [0; 0,69], two-sidedness has bought nothing — and if it is **empty**, that is the
most important result, because then no recycling fraction both spares the host and keeps
the activator in check.

## Falsifier
If none of the four locators measures the recycling fraction as a quantity — thus if 0,6901150139456149
is a model construction with no counterpart in any measurement — then `T-E14` and `T-E15` are not `TIGHT`
and the turning point must not support a decision. Say so with the locator you tested.

## Regler
`PENDING_INDEPENDENT_REVIEW`, ingen utsaga om biologisk validering, inget excluded_category material, inget som
names the collaborator. Each source with PMID eller DOI, volume and pages.
