# LANE_ANCHOR_CELLS

Build **the cells between the data anchors** so that incomplete information gives a useful answer instead of UNKNOWN. Results directory `results/LANE_ANCHOR_CELLS/`. Form A.

## Why this lane, and it brings the whole evening together

Anton 2/10: remember the function with data anchors, building the cells between them so it can handle incomplete information.

It is the only thread that explains why five lanes closed tonight. Each time the outcome was that a measurement does not exist, and each time the system answered **nothing**:

- k = Γ_II/Γ_I at DEJ: the same bond inventory compatible with k from 1 to 10⁶. `native_k_upper = None`.
- Mode-resolved DEJ separation work: no published measurement exists, `negative_result: true`.
- Clinical update rate for the robot: unreported, `clinical_sufficiency_gate = UNKNOWN`.
- Existence of the enzyme residual: determined by raw degrees of freedom that publications do not print.
- Metformin capacity: unbounded above with probability 0,113 because the extraction ratio is too close to one.

**A twin that improves on does not answer UNKNOWN.** It answers with the set of worlds still compatible with what we know, and with which measurement shrinks the set most. That is the cell between the anchors, and it is what makes the anchors an instrument instead of a list.

We already have the parts, without them being an object: the FLOWX job's certificate saying finite versus unbounded per substance; the interval enclosures in the hidden-axis chain's round 4; the witness span [1 … 10⁶] in the bond inventory. Each is a cell in Anton's sense, built once and then discarded.

## Do this

1. **Define the cell as an object, not a report.** Minimum: which anchors bound it, which quantity it concerns, what is missing, **the set of still-compatible values** (interval, union of intervals, or "unbounded above" with the direction that is unbounded), and the measurement that would shrink the set most, with the expected shrinkage.
2. **Densify between anchors, not outside them.** A cell must never extrapolate outside the anchors that bound it without labeling that. If the set can only be specified outside the anchors' validity boxes, it is not a cell but a guess.
3. **Test against the five cases above.** They are tonight's own closures and are the right test set, since each gives UNKNOWN today. For each: can the cell give a useful set? Report per case what it gives, and be honest where it cannot.
4. **The direction of the unbounded part must be carried.** A lower bound and an upper bound have entirely different value: for k, a lower bound cannot falsify us but an upper bound can kill the mechanism for every rete-ridge wavelength we computed. The cell must therefore always say which direction is open, not just that it is open.
5. **It must be consumed.** A verdict nobody reads is worthless, measured four times today. Write at least one test that FAILS if a consumer treats a cell as UNKNOWN when it actually carries a finite set.

## Strongest control and falsifiers

Classify the claim according to COMMON.md: this is `capability`. There is no equally informed counterpart, because today's alternative is to answer UNKNOWN. The control is therefore our five actual closures and what they delivered, not a method competition.

- **Falsifier:** if the cell can only return UNKNOWN for all five cases, the construction adds nothing, and that must be said directly. If it can give a set for a single one of them, that is progress and must be reported as exactly that — one of five, not "solved".
- **Forbidden:** extrapolating outside the anchors' boxes without labeling; returning an interval without saying which direction is open; building the cell without a test that fails when it is ignored; changing the graph engine (pinned `352c6d3`, also at `/opt/agents/graph_engine`).

## Delivery

The cell object, its test, and a table of the five cases with what the cell gives for each. `PORT.json` so that the next lane can place a cell between two other anchors without rebuilding anything.

Everything PENDING_INDEPENDENT_REVIEW. No internal data.
