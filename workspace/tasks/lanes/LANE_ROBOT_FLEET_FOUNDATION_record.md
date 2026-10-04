# LANE_ROBOT_FLEET_FOUNDATION — what it takes to GENERATE a fleet of medical robots

Resultatmapp `results/LANE_ROBOT_FLEET_FOUNDATION/`.

## The question
The operator: what is needed for us to be able to build a fleet of different medical
robots with generative CAD — and his own guess is that there will be a densification
of a specific graph area, and that correct resolution of the body twin is at least
as important. **Both parts are true, and the second is the binding one.** Lanen's
job is to make it a numbered list of demands instead of a conviction.

## What REDAN EXISTS, verified last night — and it lives in TWO places
The operator corrected me here: I first looked in the CAD lane and he pointed to the
field motor. **He's right, and the sharing is the important one.** The field engine owns
the GENERATIVE GEOMETRY, the CAD lane owns ASSEMBLY and certification. For our question
— derive geometry from a requirement — the field engine is the primary one.

### The field engine: generative geometry
`../3fold-field-engine/`:
- **`src/field_engine/cadbank_v1.py`** with `data/cadbank_v1`
  and `tests/test_cadbank.py` — a CAD bank with tests.
- **`reports/topopt_first_assembly_observation.json`** plus `.npz` reference — **topology
  optimization on an assembly.** It's the generative primitive that actually answers
  our question: give it loads and boundary conditions, get back a shape. The tissue
  → requirements → geometry chain goes through it, not through an assembly planner.
- `reports/recipe_extended_initial/*_cad_op_exec_v1.py.json` — CAD operation as executable steps, and
  `examples/factory_synth/cad` , plus a documented CAD pipeline in `docs/fig8/` . **Read the topopt report
  first.** It determines if the requirement joint can be entered as loads and boundary conditions.

### CAD-lanen: assembly and certification
`local_path`, branch `feat/cad-assembly-generative`, CadQuery as core:
- **`reports/certified_generative_assembly_loop.json`** — a generative loop that GROWS a part
  to a geometric boundary bind, and certifies the enclosure with the binding cause NAMED: R_max
  10,0 mm, binding condition "orthogonal-neighbour collision", predicted onset 10,0 mm, with
  a growth path where each step carries swept overlap volume in mm³. All gates pass.
- **`reports/certified_assembly_sequence.json`** — assembly sequence as PLANNER:
  it derives precedence from swept overlap volume instead of having it given. 2 out
  of 6 orders doable, pegs before cap, cap-first always blocked, all gates green.
- **`reports/third_assembly_generalization.json`** — generaliserad to en verklig hjulmontering:
  43 solider, 12 klasser, instansieringsfaktor 3,58.
- **`data/ur10e_assembly_graph.json`** — a real robot arm imported as a graph:
  13 links, 12 joints, 6 revoluta as expected. Thus: topology optimization, a
  CAD bank, enclosure certification, assembly planning, and class detection exist
  and are graded on a robot arm, a car component, and an assembly. **Read
  them before you build anything, and start in the field engine topopt.**

## The three parts of the answer, and only one of them is built
**1. The generator: EXISTS, in two halves.** The field motor's topology optimization takes loads and boundary
conditions and gives a shape; The CAD lane's loop takes geometry and provides a certified envelope with named
binding boundary. The first half is the one we need for the claim stage, and it's the one I missed first.

**2. Requirements Link: MISSING, and that's the gap.** Nothing powers the generator from tissue. Today,
the robotics branch lists components and asks what they can do; a base would go the other way
and derive the demand from the tissue. **We have exactly one example of it working:** the peel
force at the dermo-epidermal boundary of the skin, ≤0,1 mN, is a robotic requirement that comes
from a tissue measurement and not from technical convenience — and it turned out to be the
binding mechanical requirement, rather than the tip force, that pure engineering analysis
had pointed to. An example is not a method. Lanen's core task is to make it a method.

**3. Fleet: not attempted.** A fleet is not many designs but a PARAMETERIZED FAMILY.
The closest thing is class and instantiation detection within ONE assembly
(12 classes, instantiation factor 3,58), not a family of interventions.

## Do this
1. **Take three interventions and derive the requirements from below the tissue.** Choose three that differ in physics, not in
   anatomy — for example, an incision in skin, a needle track in soft tissue, and a non-contact
   laser treatment. For each: which quantities set the requirement? Power, speed, accuracy,
   temperature, time.
2. **Classify each derived requirement into three boxes, and the third one is the important one:**
   - **COMPUTED** — we calculate the quantity today, with cell and unit specified.
- **MEASURED BUT UNCONNECTED** — the count is in any of the
four pools (see `COMMON.md` ) but no cell counts it.
   - **UNMEASURED** — none of them. Then the demand is not derivable and the generator cannot be powered by it.
   **AThe numbers in the three boxes ARE the answer to the operator's question.** Report them before anything else.
3. **Try the resolution dependence, because it is his second point and it is decidable.** For each requirement in
   the box CALCULATED: does the requirement change if the quantity is counted a level finer? Use the levels in `COMMON.md`
   (MOLECULE/CELL/FIBRIL/TISSUE/ORGAN/WHOLE_BODY). A requirement that is invariant under refinement needs
   no finer twin; one that moves does. **It divides the requirements list into the part that requires higher
   resolution and the one that doesn't ** — and that's a much sharper statement than "resolution is
   viktigt".
4. **Put the claim against the generator for ONE cases.** Take a claim in the COUNTED box, feed it as binding
   condition in the certified generative loop, and see if the loop returns an envelope. If it does
   there is the chain tissue → requirements → geometry, and it is the lane's strongest possible delivery. About it
   doesn't do it, tell me which interface is missing.
5. **The fleet, first by step 4.** What varies between two medical robots for different interventions, and
   what is common? If the requirements in the COUNTED box are the same for several interventions and only the numbers differ,
   is the family parameterizable and the number of parameters is the dimension of the fleet. Count it.

## What the operator has already decided, and which must not be redone
- The threshold is a robot so good that the EU regulations must be revised. The finding there: **specifications do
  not move the regulatory framework, a damage measurement does.** All 15 mechanical nodes are regulatory inert.
- Binding boundary in the branch is registration, not sensory. **But the number is incorrectly delimited since i
  night** — 145 s is a core, the spatial part is unpublished, and a one-time cost over
  866,75 s does not count. Do not build on 1 192 updates without reading `LANE_REGISTRATION_SCOPE`.
- Tissue physics belongs to the intelligence engine, the robot owns its precision and
  planning. But that division leaves **what a tool does with tissue** without an owner,
  and that's where surgery sits. That gap is lane's, and it should be named as such.

## Control and falsifier
- **Check:** the robot branch as it stands, with 0 of 202 ports carrying a calculated number and 167 requiring
  recalculation. The gain is the number of claims in the box COUNTED, not the number of nodes.
- **Falsifier:** if almost all derived requirements end up in the UNMEASURED box, a theoretical basis
  for robot construction is not reachable from today's twin, and **then the requirements list is instead
  a measurement order list** — a perfectly fine and important outcome to be reported plainly.
- **Prohibited:** to build a new generator; writing in `cad-to-simulation-J`
  , `3fold-field-engine` or `~/projects/bodytwin` — all three are read
  view; to count a component specification as a derived requirement.

## Deliverable
`PORT.json` : the number of the three boxes, the list of requirements per intervention with quantity
and unit, the resolution dependency per requirement, and the outcome of the one generator attempt.

No internal data. Everything PENDING_INDEPENDENT_REVIEW.
