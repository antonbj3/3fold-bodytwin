# LANE_QUANTITY_KEY_BRIDGE — the source keys by DIMENSION, the cells by MODEL ROLE

Results folder `results/LANE_QUANTITY_KEY_BRIDGE/`.

## The root cause, now measured rather than named
Seed 2 was for the dense expert to be CONTENT that cells consume. After twelve rounds,
the source was well ordered — 188 records, replay-exact, N_eff per record, 940 forbidden questions declined,
82 of 155 older records identified as primary numerical observations — and **zero cells
consumed a single record**. The lane was closed with the root cause named but not addressed.

I measured it tonight and it is simpler than it looked:

**The source keys by PHYSICAL DIMENSION:** clearance_distance 46, pressure 28, time 22, force 20,
flow 16, proportion 8, temperature 7, dose 6, protocol_or_geometry 6, outcome_proportion 5,
concentration 4, tissue_damage 4, strain 2, current 2, signal_amplitude 2, density 1, volume 1.

**The cells key by MODEL ROLE:** `cut_depth`, `cut_length`, `plug_time`, `upstream_fraction`,
`vessel_density`, `vessel_radius`, `bleeding_rate`, `blood_loss` — and `R_thermal`, `R_recruit`,
`Tm`, `crosslink_density`, `gamma_fibril`, `gamma_tissue`, `hyp_fraction`, `recruitment`.

A record that says "this is a pressure" can never find a field named `upstream_fraction`.
It is a **level collision and not a vocabulary collision**, and it is solved with a mapping and not with
more records.

## What the bridge would give, calculated
I mapped the 16 cell variables to their dimension class and counted available records:

| cell variable | dimension class | records |
|---|---|---|
| `cut_depth`, `cut_length`, `vessel_radius` | clearance_distance | **46** |
| `plug_time` | time | **22** |
| `bleeding_rate` | flow | **16** |
| `upstream_fraction`, `hyp_fraction`, `recruitment`, `R_thermal`, `R_recruit` | proportion | **8** |
| `Tm` | temperature | **7** |
| `blood_loss` | volume | 1 |
| `vessel_density`, `crosslink_density` | density | 1 |
| `gamma_fibril`, `gamma_tissue` | energy_per_area | **0** |

**14 of 16 cell variables have at least one matching record. Two have zero — and those are exactly the cells' two
NULL quantities.** The source has **no class for energy per area at all**, which explains exactly why
the 188 records never touched the cells' real debts: the quantity missing from the model is the
quantity also missing from the source.

## Do this
1. **Build the mapping as code, not as a table in a report.** Every cell variable gets a
   dimension class and unit. The mapping belongs in `tasks/build_night/surgical_expert.py` as a
   `--for-cell <CELL>` mode that returns the records per variable. Read the script first; it already has
   `--quantity` and per-record `trust` and `N_eff`.
2. **Require three things of a record before it may be consumed**, and decline otherwise: (a) the dimension matches,
   (b) **the validity range overlaps the cell's** — a bursting limit measured on vessels with 4,5 to 7,5 mm
   diameter must not be fed into a cell whose vessels have at most 0,5 mm radius, and that error is already
   pointed out in the source's own material, (c) the record carries a measured value and not just threshold text.
3. **Deliver ONE actual consumption all the way through.** Choose the record that best matches `plug_time` or
   `cut_depth`, feed it into the cell, run, and report what moved. If nothing moves,
   the record is not load-bearing — say so and take the next one. Twelve rounds failed at precisely this, so it
   is the lane's only real test.
4. **And report N_eff with every consumed record.** A record with N_eff = 1 and high
   trust is not five independent records. The source has 33 records with N_eff = 1,0 and 155 older ones carrying
   `UNTRACED_NOT_ZERO`, thus unreviewed rather than zero — distinguish the two in the report.
5. **Energy per area must not be invented.** The two NULL quantities require a measurement that does not exist
   (acquisition item 8 in `notes/ACQUISITION_LIST.md`, two data laws). Add them as a class in the source's
   schema so that future harvesting can fill it, but do not invent a value.

## Strongest control and falsifier
- **Control:** the current situation, thus zero cells consuming a record after twelve rounds. The gain is
  the number of records actually fed into an executable cell that move something, not the number that match
  by dimension.
- **Falsifier:** if the 46 distance records all fail on validity range — wrong vessel calibre, wrong
  tissue, wrong scale — the dimension mapping is necessary but insufficient, and then
  the validity range is the real key. That would be a sharper finding than the bridge itself and must be
  reported as such.
- **Forbidden:** consuming a record whose validity range does not overlap the cell's; counting a
  dimension match as consumption; harvesting new records before a single one is consumed.

## Delivery
`PORT.json`: the mapping cell variable → dimension class → unit → number of records meeting all three
requirements, plus at least one actual consumption with what it moved and its N_eff.

No internal data. Everything PENDING_INDEPENDENT_REVIEW.
