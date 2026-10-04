# PROOF_LANE_SCALE_BRIDGE — what must a tight edge across a scale jump declare?

## The measured state
The network has 168 edges and 163 variables. I built `tasks/assembly/resolution_audit.py` tonight and
measured: **112 of 163 variables are UNKNOWN** at the scale level, so only **29 edges have a known level at both
ends** and 139 cannot even in principle be tested for a scale jump. Of the 29, **11** cross a level.

One edge remained `TIGHT` across a scale jump without any declared bridge:
`T-E1-charge_inventory-vs-ionic_swelling_pressure`, **MOLECULE versus TISSUE, three steps apart**. (Two
others, `T-E14` and `T-E15` MOLECULE versus CELL, were lowered from TIGHT tonight when their axis proved
to be a model closure — so it was not the scale jump that defeated them.)

A TIGHT edge across a scale jump claims that a molecular quantity **determines** a tissue quantity. It
requires an averaging step: how many units, over what volume, with what assumption about their
arrangement. `T-E1` states none.

## The question
What must such an edge declare for the claim to be **well-posed** — not true, just
well-posed? I want it as a checklist supported by theorems, not as good practice. Specifically:

- the **minimum** set of declarations that makes a molecule-to-tissue edge well-posed, and proof that
  each is necessary through a counterexample that fails when it is omitted;
- when a scale jump may be called **TIGHT** at all: is there a condition under which the homogenisation error is
  zero, or is every scale jump inevitably approximate and hence at most `OPEN`?
- how the homogenisation error **propagates** to a decision: given a relative error `e` in
  the averaging step, what bound applies to the decision's outcome, and is the bound sharp?
- and the specific case: for charge inventory versus ionic swelling pressure, which quantities
  would need to be measured for the jump to become well-posed? We now have observed charge endpoints
  **0,026 and 0,261 mEq/g wet weight** that give the pressure intervals `[0,0233; 0,0386]` and
  `[0,6413; 0,6938] MPa`, separation **0,6026633024235674 MPa**, without overlap.

## The strongest control
Always setting `OPEN` on every edge that crosses a level. It is always correct and erases
the information that some jumps are actually well-determined. Your checklist is worth something only if it
places a **named count** of the eleven crossing edges on the right side.

## Falsifier
Construct an edge that passes your entire checklist and still gives a wrong decision because of
the scale jump. If you succeed, the checklist is insufficient and must be rejected.

## Rules
`PENDING_INDEPENDENT_REVIEW`, no biological validation statement, no excluded_category material, nothing like
names an individual. The endpoint rule applies to any precision.
