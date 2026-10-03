# LANE_CELLS_FROM_DOMAIN_DATA — 156 domains with data, 43 cells. Build from what already exists.

Resultatmapp `results/LANE_CELLS_FROM_DOMAIN_DATA/`.

## The finding the lane builds on
`~/projects/bodytwin/data/` has **156 physiology domains with data** and 1 309 completed working cells with
decisive numbers, of which 696 MEASURED and 431 REFUTED. The workspace has **43 cells** and cites none of
them. Nothing was lost; it was disconnected. The index is in `notes/OLD_CELL_INDEX.json` and
is queried with `python3 tasks/build_night/index_old_cells.py --quantity <text>`.

And old work lives in **four** places, not one: `data/` (1 620 cell records), `reports/` (821),
`scripts/physics_exp/` (374) and `docs/inherited_memory/` (417 distilled lessons). Swept 2/10:
**896 of 1 612 names are not called by anything in the workspace.** Map: `notes/OLD_WORK_UNCALLED.md`.

## Do this
1. **Choose domains by what makes an EDGE possible, not by count.** A new cell is only worth something
   if it shares a quantity with an existing cell at the same resolution level. Read `tasks/build_night/COMMON.md`
   about edges at the finest shared level and keying by quantity.
2. For each candidate domain: is there a decisive quantity with a unit and a number in the source? If yes, build a
   cell with the same contract as `tasks/free48/sources/Q036/Q036_model.py` — declared parameters with
   units, resolution level per quantity, dimensional check run every time, falsifier, and
   `external_referent`. Also look at `SURG_HEMOSTASIS_model.py` for how an identifiability certificate
   is written when a quantity cannot be determined.
3. **Deliver FEWER and COUPLED cells, not more.** Three cells that each share a quantity with an
   existing cell are worth more than twenty isolated ones. Say for each new cell which edge it enables, with
   quantity, unit and level.
4. **The 431 refuted cells are as useful as the measured ones** — a preserved negative prevents a
   later lane from rerunning a question that has already failed. If a candidate domain already has a REFUTED
   answer, do not build the cell; record the negative and move on.

## Control and falsifier
- **Control:** the 43 existing cells. The gain is the number of new cells that carry an edge to an existing
  quantity, not the number of new cells.
- **Falsifier:** if none of the 156 domains has a quantity shared with an existing cell at
  the finest shared level, the disconnection is not an accounting gap but a real
  incompatibility of quantities — report it directly, that would be a substantial finding.
- **Forbidden:** writing in `~/projects/bodytwin` (read view only); building a cell without
  a dimensional check; counting cells as progress.

## Delivery
`PORT.json`: per new cell, quantity, unit, resolution level, the edge it enables, and its
falsifier. Plus a list of candidate domains dismissed against a preserved REFUTED.

# Round 2 (koordinatorn, 3/10 00:15) — The oxygen roof breakage is the find, spread it
You pulled the same edge I pulled by hand and went further: the roof 5,4545 ATP per O2; and **12 av 12 nativa
fall above it, max 69,68 i.e. 12,8× over**. Den nativa formuleringen producerade ATP ur ingenting,
and the fault is only visible when the cell is connected to the molecular ceiling.
- **Changed operation:** determine how: BRETT Run the same roof control against all cells that carry one
  ATP- or oxygen quantity, not just those 12 the cases. Count how many cells and how many operating points
  Which breaks the roof, and how much. **Measure the distribution before there is a statement about the system.**
- **And say WHY shared distribution does not break.** 0 av 12 i den formuleringen mot 12 av 12 native
  a large difference; which term separates them? If shared propagation suppresses the violation without resolving it
  is the wrong wording, not the right one.
- **Connected find from my own hand tonight:** I swiped the same parameter over his source-backed tape
  [0,18333; 0,20370] genom Q088 and received **exakt ett teckenbyte av 18 driftpunkter**, vid 2,0 kPa
  syre, 0,6 mM/s work and 40 µm radie. Din 233,28 s for creatine phosphate at the same point is inside the
  mitt svep 225,9–235,1. The same operating point therefore carries both the character change and the roof breakage — Look at that.
- **Hinder:** `frozen_external_gate = FAIL` while the technical connection passes. Tell me which of them
  frozen criteria that are flawed and by which number.
- **Falsifierare:** If the ceiling break disappears when units are re-normalized, it was a unit trap and
  not physics — Try it first, that's the cheapest explanation.

# New round (coordinator, 3/10 00:55)
The distribution is measured and the statement is now permitted: 78 av 78 nativa punkter bryter syretaket, median
2,18× and max 12,8× over, while the split formulation breaks 0 av 78 med median 5,000 precis
under the ceiling. The unit invariance error is exactly 0,0So my single trap hypothesis is set by you.
the offence is limited to Q088, med 0 av 57 kinetiska punkter brytande.

The obstacle is now the interpretation, not the measurement. Two readings are open and they have different consequences:
either the native wording is wrong and is to be replaced, or the divided wording is one:
silence that hides the same error behind normalization. That the shared median lies on 5,000, that is
exactly at the roof rather than below it, is suspected — a formulation that lands just on its own
border is often built-in.

Modified operation: determines it. Contains the split formulation ceiling as a condition, explicitly
or implicitly? If yes, it suppresses the violation and must not be used as a correction. If no, it is a valid
The correction and the national wording shall be withdrawn.

And connect to my own measurement: I pulled last night a second edge to Q088, where the oxygen limit was set
at alveolar value 99,76 mmHg medan Q021 counting 83,998 at rest. With the calculated value goes
sensitivity[3] from 46,568 s till censurerad vid horisonten 600 s. So the same cell carries both
The roof breakage and an overly high oxygen limit. Try whether the two are connected: if the cell receives more oxygen than
Physiology allows, can a ATP-exchange over the roof is a consequence of that and not a separate bug.

Falsifier: if the ceiling breakage remains unchanged at Q021:s oxygen limit, are the two faults independent and
both shall be corrected separately.

# New Round (koordinatorn, 3/10 04:50)
Horisonten 10 ms med 100 ms and 1 s som FAIL is an honest delimitation and better than a claim of
Full validity. Keep it.

But one thing needs to be resolved now: the inherited frozen external criterion has fallen into FEM consecutive rounds on:
samma tal, 2,2727 % to a tolerance of: 2,0000. A glitch on 0,27 percentage points remaining
in five rounds is not a research result, it is a stuck gate. Determine which of three it is:
the tolerance is set too tight for the quantity, the inherited number is wrong, or our calculation is wrong.
Tell me which one and with what support, letting it stand a sixth round is bustywork.

Changed operation otherwise: what prevents 100 ms? Du klarar 10 ms without replaying history, and 100 ms
fall. Is it the length of the history, the number of carrier states, or that a slow process enters
mellan 10 and 100 ms? The latter is the interesting answer, because then there is a time scale in the system that
is not modeled, and it should be named.

Falsifierare: om 100 ms falls of pure calculation cost and not of physics, there is no scientific
limit without a resource limit, and then it shall be reported as a cost factor instead of as a
FAIL.

# OPERATOR DIRECTIV 3/10 09:25: EMOVE THE SOLUTION OF HELA BANAN
Operator's directives: Raise the resolution all over the track. You have moved the horizon 10 ms till 100 ms till
1 s and keeps the validation gate, which is the right kind of resolution increase in time.

Changed operation: do the same in RUMMET. Du har 42 spatial nodes; what happens at the double and the
quadruple? Report if they passed the gates holding, and especially if ANT-intervallet [−24 112,03;
−24 110,27] molekyler krymper eller flyttar. Ett intervall som FLYTTAR at refinement means that: 42
nodes were not enough and that previous rounds' numbers are resolution dependent.

Falsifier: if the range moves more than its own radius 0,877 molecules at doubled
Node density, previous results are not converged and should be relabelled.
