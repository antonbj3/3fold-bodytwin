# LANE_CELLS_FROM_DOMAIN_DATA — 156 domains with data, 43 cells. Build from what already exists.

Results folder `results/LANE_CELLS_FROM_DOMAIN_DATA/`.

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
