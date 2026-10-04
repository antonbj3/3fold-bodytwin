# PROOF_LANE_NO_DISPERSION — can a scalar gate be decided without having the spread?

## The measured state, with numbers I verified myself
`results/LANE_DISTRIBUTION_GATE_43/ALL_43_GATES_R3.json` carries 43 cells. The scalar criteria
decide 41 of them (30 PASS, 11 FAIL, 2 UNKNOWN); after criterion containment 25 PASS, 11 FAIL,
7 UNKNOWN. With the spread carried, **5** are decided (2 PASS, 3 FAIL) and 38 remain UNKNOWN.

The 38 are not disagreement. **31 of them have the reason `"Matching predictive outcome measure absent or not
supplied to the native primary receiver"`** — the cell has no distribution counterpart to be tested
against. 6 are inherited design gates without biological admission, 1 is a one-sided risk entry.

Three cells actually flip: Q012, Q031, Q160, all `flip_direction = FALSELY_SAFE`. For Q160
I verified the SD/SEM rule against raw data: n = 13 published clusters, SD 102,80627271159125 µm²/s,
SD/√13 = 28,51332982315204 µm²/s exactly as recorded; the error mass 2/13 against the criterion ≥ 95 % per
specimen; Clopper-Pearson's one-sided 95 % interval [0,028; 0,410].

## The question
For the 31 cells we have the criterion, its frozen text and the scalar outcome — but **no
spread**. Is there a test **on the criterion and the operation themselves**, without the spread, that
decides whether any spread at all could flip the verdict?

Formally: for cell `c` with operation `F`, input `x̂` and criterion `G` (band, point or
direction), decide whether there is an admissible perturbation `δ` such that `G(F(x̂+δ))` differs from
`G(F(x̂))`. I suspect this is the same two-copy construction you used on
the insufficiency question: `d` factors through `σ` exactly when a formula is unsatisfiable. If that is
the case, the 31 can be divided into two classes without a single new measurement:

- **cannot flip** — the scalar verdict stands, and the cell needs no spread;
- **can flip** — the cell is undecided, and then I want to know **which** quantity must be measured and
  how narrow its spread must be for the verdict to stand.

## Strongest control
The trivial classification: say "can flip" for all 31. It is always correct and always useless.
Your result is only worth anything if it puts a **named number** of cells in "cannot flip"
with a certificate per cell.

## Falsifier
If you construct a cell where the test says "cannot flip" and an admissible spread still
flips the verdict, the construction is wrong and must be rejected, not patched. Build that
counterexample search yourself and run it against all cells you classify.

## What you have already corrected and must not redo
A width-based rounding rule fails: your own `PROOF_LANE_PRECISION_PROP` let one fail on
all 1 200 constructed boundary cases, and the real rule is that both enclosure endpoints must
round identically (`[1,38188861; 1,38218988]` carries `1,382`, three decimals, not four). Use
the endpoint rule here too when specifying how narrow a spread must be.

## Rules
Status `PENDING_INDEPENDENT_REVIEW`, no claim of biological validation, no excluded_category material,
nothing that names an individual person. Save the certificates as separate files and give the count in the final message.
