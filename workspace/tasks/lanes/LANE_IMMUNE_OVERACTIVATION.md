# LANE_IMMUNE_OVERACTIVATION — the second failure mode: too much, not just too little

Results catalog `results/LANE_IMMUNE_OVERACTIVATION/`. Decisions outside the eye.

## Measured position
Three decisions are built and executable, and I ran all three today:

- `complement_amplification_decision.py`: convertase lifetime **129,84 s** surface bound and **111,09 s** in
  fluid phase (half-lives 90 and 77 s), C3 cleaved by convertase **231,1** at saturation C3 and
  **115,6** vid Km.
- `complement_window_decision.py`: host cell flips at recycle rate
  **0,6901150139456149**; the activator **never** reverses over the entire range. Growth just below
  tipping point: host −0,00018768 per minute, activator +0,21710422 per minute. A character change for
  values over 2001 points, zero for the activator.
- `immune_window_decision.py`: values are spared for ρ in **[0; 0,69]**, thus 691 of 1001 points; with
  two simultaneous demands shrink it to **510 of 1001**, the range [0,181; 0,69].

The sources are available as locators in the code: PMID 6229580, doi `10.1016/j.molimm.2007.11.003`,
doi `10.1182/blood-2005-02-0782`, doi `10.1182/blood-2008-11-189944`.

## The obstacle, in two parts
1. **The edges are not anchored.** Six of the grid's 163 edges touch the immune system, two of them stand
   `TIGHT` — `T-E14` and `T-E15` — and **none** of the six carries a `external_reference` field. The speech
   thus rests in code while the edges claim to be tight without pointing to a measurement.
2. **Only the one failure mode is tried.** All three decisions ask if the host is **spared**. Underactivation
   is the other half: a system that protects the host perfectly and at the same time allows the activator to grow with it
   0,217 per minute is not regulated, it is turned off. And overactivation — host damage — is that mod
   that actually does clinical damage, and it's not scored at all.

## The operation
Make the window two-sided. A decision that returns **both** limits: the recovery rate below
which the activator is no longer lysed (underactivation) and that above which the host is damaged
(overactivation), in the same unit, against the same sources. Return the boundaries, not their description.

Then anchor `T-E14` and `T-E15`: put `external_reference` on each with the one of the four
the locators that actually carry the number, and if none do, the status will drop from `TIGHT`.

## Strongest control
The one-sided reading that exists now, i.e. only "the values are spared". If the two-sided window does not
is narrower than [0; 0,69] the two-sidedness hasn't bought anything — and if it's **blank** it's
most important result, because then there is no recycling rate that both conserves values and lasts
aktivatorn i schack.

## Forgers
If none of the four locators measures the recycling rate as a magnitude — that is, if 0,6901150139456149
is a model construction without a counterpart in any measurement — so `T-E14` and `T-E15` are not `TIGHT`
and the turning point must not carry a decision. Say it with the locator you tried.

## Rules
`PENDING_INDEPENDENT_REVIEW`, no statement of biological validation, nothing like
names an individual. Any source with PMID or DOI, volume and pages.
