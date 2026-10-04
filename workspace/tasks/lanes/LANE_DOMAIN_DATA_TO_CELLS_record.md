# LANE_DOMAIN_DATA_TO_CELLS — build cells from domain data that already exists and that no cell touches

Resultatmapp `results/LANE_DOMAIN_DATA_TO_CELLS/`.

## Why the lane exists, measured
The workspace has **43 cells**. The source project's `data/` has **505 directories**, and at least **19** names a
tissue or organ area. I checked seven of them and five bear a finished result with numbers:
`corneal_transparency` (20 kB), `complement_cascade` (28 kB), `bcell_affinity_maturation` (20 kB),
`capillary_starling` (24 kB), `tissuetwin_recovery` (96 kB). **None of the five is available as a cell with us.**
Two of them hit the operator's seeds straight, and they are the first two tasks of the lane.

And a warning to save you time: the **indexed cell archive** (1 620 entries in
`notes/OLD_CELL_INDEX.json`) is NOT the content source. Of the 1 274 with described decisive quantity is
**590, so 46,3 %, accounting formulations** (number of lines, artifacts, documents) and only **102,
thus 8,0 %,** carries a physical quantity without accounting words. The contents are in the `data/` directories, not
in the cell records. Therefore, extract `data/` and read the cell entries only when a directory is difficult to interpret.

The source project is readable but **read only**. Everything new is written in the workspace.

## Task 1 — cornea: give the cornea the conversion it lacks
`source_repository/data/corneal_transparency/corneal_transparency_results.json` carries exactly
the relation eye lane has as a fault, i.e. disturbed lattice order → scattering tail:

| magnitude | value |
|---|---|
| `C_PREF` | 16,3511 |
| lattice | a = 14,0 nm, φ_areal = 0,28, ρ2d = 4,5473e-4 /nm², d_hex = 50,3917 nm |
| index | n_fibril = 1,411, n_matrix = 1,365, Δn = 0,046 |
| thickness | L = 500 000 nm |
| **the cost of disorder** | `falsifier_ratio_poisson_over_physio` min **20,2466**, max **59,4856**, mean **38,1381** |
| **wavelength exponent** | **−4,8765** |
| validation | large disorder: median relative error 7,9271 %; perfect lattice should give 0 and doesn't — declared as sidelobe artifact from finite patch |

1. **Build the cell `CORNEA_SCATTER` in `tasks/free48/sources/`** with that chain: lattice order → scattering
   → transmission, with these parameters as declared input.
2. **And cross against the eyeline's own number.** `results/LANE_EYE_OPTICAL_TWIN/PORT_R7_V2.json` wears a
   meridional spreading tail at an arcminute of **1,6649e-4**, explicitly *not* converted from
   the axial phase model of the lane. The iris r7 standard values (d = 53 nm, r = 14,73 nm) are close
   of the domain data d_hex = 50,3917 nm and a = 14,0 nm, so the scales are comparable. **Say if the two numbers are
   compatible, and by what margin.** It is the lane's first real test.
3. **The exponent −4,8765 is steeper than the Rayleigh −4.** Determine if the difference is physical (correlated
   disorder) or an artifact of the finite patch, with the same argument as the validation-note
   uses. An answer here is worth more than a new parameter study.

## Task 2 — complement: operator immune seed has ready parameters
`source_repository/data/complement_cascade/complement_cascade_results.json` bears
the discrimination between host and activator, which IS the issue of under- and over-activation:

| magnitude | value |
|---|---|
| half-life unstabilized / activator / host | 3,0 / 30,0 / **1,05** min |
| properdin stabilization | **10,0×** |
| regulators | FOLD_H = 4,7619, CD55 = 6,0 (illustrative), MAC/CD59 = 2,3333 |
| C5 | Km = 1 400 nM, kcat = 0,288 /min, C5_0 = 400 nM |
| C3 | 6 684,49 nM (1,25 mg/mL, MW 187 000) |
| decay / inhibition / MAC | activator 0,0231 / 0,0200 / 0,3000, host 0,6601 / 0,0952 / 0,1286 |

1. **Build the cell `COMPLEMENT_DISCRIMINATION`** with the readout *how many times more MAC formation a
   activator area gets more than one host area*, and indicate which of the above numbers set the quota.
2. **And respond to the seed's two directions with the same cell:** which parameter change gives
   UNDERactivation, i.e. that an activator is treated as a host, and which gives OVERactivation, i.e
   that a host surface is attacked. Specify the margin to each turn in that parameter's own units. It is
   the form the operator requested: under- and over-activation as two limits around a work area.
3. **Mark `FOLD_CD55_ILLUSTRATIVE` as illustrativee throughout the report.** The name says that the number
   is not measured, and it must not carry a conclusion.

## Requirements that apply to both cells
- **A distribution gate, not a point-in-band test.** I measured last night that **21 of 43 cells carry a
  scalar point or band gate and only 7 carries any distribution statistics**, and a swarm job showed
  my own hemostasis cell that such a gate approves a model where **7,305 % of mass never closes**.
  Each new cell should report quantiles and tail mass, not just a band mean.
- **Declare the scope** in each input row: species, tissue, temperature, method. Of the night
  permeability theory was that published numbers can span 8–9 powers of ten depending on layer and
  driving force, and a number without range is not falsifiable.
- **The domain data is OUR EGEN previous calculation and thus no external reference.** It may be input data and
  point of comparison, never a validation.

## Strongest control and falsifier
- **Control:** the current state, i.e. that no cell touches any of the five directories. The win is the number
  executable cells that consume domain data and return a number, not the number of directories inventoried.
- **Falsifier:** if the cornea number and the eye lens 1,6649e-4 are incompatible by more than a power of ten is
  at least one of them wrong — and then it's the finding, not the cell. If so, report it as the main outcome.
- **Forbidden:** to write in `~/projects/bodytwin`; to inventory more directories before the two cells
  drive; to call the domain data a reference.

## Delivery
Two executable cells in `tasks/free48/sources/`, each its `PORT.json` with consumed numbers and their
areas of validity, the statement of compliance against the eyelash spread tail, and the two turning margins
for under- and over-activation.

No internal data. Everything PENDING_INDEPENDENT_REVIEW.
