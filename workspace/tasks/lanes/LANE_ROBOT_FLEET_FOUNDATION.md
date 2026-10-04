# LANE_ROBOT_FLEET_FOUNDATION — what it takes to GENERATE a fleet of medical robots

Result folder `results/LANE_ROBOT_FLEET_FOUNDATION/`.

## The question
The operator: what is needed for us with generative CAD to be able to build a fleet of different medical
robots — and his own guess is that there will be a densification of a specific graph area, as well
that correct resolution of the body twin is at least as important. **Both parts are correct, and the
others are the binding.** The lane's task is to make it a numbered list of demands instead of a
conviction.

## What ALREADY EXISTS, verified last night — and it lives in TWO places
The operator corrected me here: I first looked in the CAD lane and he pointed to the field engine. **HHe has
right, and the division is the important one.** The field motor owns it GENERATIVE GEOMETRY, the CAD lane owns
ASSEMBLY and certification. For our question — derive geometry from a requirement — the field engine is it
primary.

### The field engine: generative geometry
`../3fold-field-engine/`:
- **`src/field_engine/cadbank_v1.py`** with `data/cadbank_v1` and `tests/test_cadbank.py` — a
  CAD bank with test.
- **`reports/topopt_first_assembly_observation.json`** plus `.npz` reference —
  **topology optimization on an assembly.** It is the generative primitive that actually responds to our
  query: give it loads and boundary conditions, get back a shape. The chain tissue → requirements → geometry goes
  through it, not through an assembly planner.
- `reports/recipe_extended_initial/*_cad_op_exec_v1.py.json` — CAD operations as executable steps, and
  `examples/factory_synth/cad`, plus a documented CAD pipeline in `docs/fig8/`.
**Read the topopt report first.** It determines if the requirement joint can be entered as loads and boundary conditions.

### the CAD lane: assembly and certification
`local_path`, branch `feat/cad-assembly-generative`, CadQuery as core:
- **`reports/certified_generative_assembly_loop.json`** — a generative loop that GROWS a detail into
  a geometric limit binds, and certifies the envelope with the binding cause NAMED: R_max
  10,0 mm, binding conditions "orthogonal-neighbour collision", predicted onset 10,0 mm, with a
  growth path where each step carries swept overlap volume in mm³. All gates pass.
- **`reports/certified_assembly_sequence.json`** — assembly order as PLANNER: it derives
  precedent out of swept overlap volume instead of having it given. 2 of 6 schemes feasible,
  pins before cap, cap-first always blocked, all gates green.
- **`reports/third_assembly_generalization.json`** — generalized to a real wheel assembly:
  43 solids, 12 classes, instantiation factor 3,58.
- **`data/ur10e_assembly_graph.json`** — a real robot arm imported as a graph: 13 links,
  12 leads, 6 revolute as expected.
Thus: topology optimization, a CAD-bank, certification of the enclosure, assembly planning and
class detection exists and is graded on a robot arm, a car component, and an assembly. **Read them
before you build anything, and start in the field engine topopt.**

## The three parts of the answer, and only one of them is built
**1. The generator: EXISTS, in two halves.** The field motor topology optimization takes loads and boundary conditions
and gives a shape; the CAD lanes loop takes geometry and provides a certified envelope with named binding
limit. The first half is the one we need for the claim stage, and it's the one I missed first.

**2. The claim joint: MISSING, and that's the gap.** Nothing drives the generator from tissue. Today lists
the robot branch components and asks what they can handle; a base would go the other way and derive the claim from
the tissue. **We have exactly one example of it working:** the peeling power at the skin's dermo-epidermal
limit, ≤0,1 mN, is a robotic requirement that comes from a tissue measurement and not from technical convenience —
and it turned out to be the binding mechanical requirement, in front of the tip force, as a pure
engineering analysis had pointed to. An example is not a method. The lane's core task is to make it happen
a method.

**3. Fleet: not tried.** A fleet is not many designs but one PARAMETERIZED FAMILY. The
the closest available is class and instantiation detection within ONE assembly (12 classes,
instantiation factor 3,58), not a family of interventions.

## Do this
1. **Take three interventions and derive the requirements from below the tissue.** Choose three that differ in physics, not in
   anatomy — for example, an incision in skin, a needle track in soft tissue, and a non-contact
   laser treatment. For each: which quantities set the requirement? Power, speed, accuracy,
   temperature, time.
2. **Classify each derived requirement into three boxes, and the third one is the important one:**
   - **CALCULATED** — we calculate the quantity today, with cell and unit specified.
   - **MEASURED BUT UNCOUPLED** — the number is in any of the four pools (see `COMMON.md`) but no cell
     counting it.
   - **UNMEASURED** — none of them. Then the demand is not derivable and the generator cannot be powered by it.
   **The numbers in the three boxes ARE the answer to the operator's question.** Report them before anything else.
3. **Test the resolution dependence, because it is his second point and it is decidable.** For each requirement in
   the box COUNTED: does the requirement change if the quantity is counted a level finer? Use the levels in `COMMON.md`
   (MOLECULE/CELL/FIBRIL/TISSUE/ORGAN/WHOLE_BODY). A requirement that is invariant under refinement needs
   no finer twin; one that moves does. **It divides the requirements list into the part that requires higher
   resolution and the one that doesn't ** — and that's a much sharper statement than "resolution is
   important".
4. **Put the claim against the generator for ONE cases.** Take a claim in the COUNTED box, feed it as binding
   condition in the certified generative loop, and see if the loop returns an envelope. If it does
   there is the chain tissue → requirements → geometry, and it is the lane's strongest possible delivery. About it
   don't do it, say which interface is missing.
5. **Fleet, first by step 4.** What varies between two medical robots for different interventions, and
   what is common? If the requirements in the COUNTED box are the same for several interventions and only the numbers differ,
   is the family parameterizable and the number of parameters is the dimension of the fleet. Count it.

## What the operator has already decided, and which should not be redone
- The threshold is a robot so good that the EU's regulations must be revised. The find there: **specifications
  doesn't move the regulatory framework, a damage measurement does.** All 15 mechanics nodes are regulatory inert.
- Binding boundary in the branch is registration, not sensory. **But the number is incorrectly delimited since last
  night** — 145 s is a kernel, the spatial part is unpublished, and a one-time overhead of
  866,75 s does not count. Do not build on 1 192 updates without reading `LANE_REGISTRATION_SCOPE`.
- Tissue physics belongs to the intelligence engine, the robot owns its precision and planning. But that one
  the division leaves **what a tool does with tissue** without an owner, and that's where surgery sits.
  That gap is lane's, and it should be named as such.

## Control and falsifier
- **Control:** the robot branch as it stands, with 0 of 202 ports carrying a computed number and 167 as
  requires recalculation. The gain is the number of claims in the box COUNTED, not the number of nodes.
- **Falsifier:** if almost all derived claims end up in the UNMEASURED box, is a theoretical basis for
  robotic construction not reachable from today's twin, and **then the list of requirements is instead one
  measurement order list** — a full and important outcome that must be reported plainly.
- **Prohibited:** to build a new generator; to write in `cad-to-simulation-J`, `3fold-field-engine`
  or `~/projects/bodytwin` — all three are read view; to count a component specification as one
  derived requirement.

## Delivery
`PORT.json`: the number of the three boxes, the list of requirements per procedure with quantity and unit,
the resolution dependency per requirement, and the outcome of the one generator attempt.

No internal data. Everything PENDING_INDEPENDENT_REVIEW.
