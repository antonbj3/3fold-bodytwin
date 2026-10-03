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
`the public staging tree/3fold-field-engine/`:
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
`~/projects/cad-to-simulation-J`, gren `feat/cad-assembly-generative`, CadQuery as core:
- **`reports/certified_generative_assembly_loop.json`** — a generative loop that GROWS a part
  to a geometric boundary bind, and certifies the enclosure with the binding cause NAMED: R_max
  10,0 mm, binding condition "orthogonal-neighbour collision", predicted onset 10,0 mm, with
  a growth path where each step carries swept overlap volume in mm³. All gates pass.
- **`reports/certified_assembly_sequence.json`** — assembly sequence as PLANNER:
  the precedence out of swept overlap volume instead of getting it given. 2 av 6 arrangements that are feasible;
  sticks before lid, lid-first always blocked, all validation gates green.
- **`reports/third_assembly_generalization.json`** — generaliserad till en verklig hjulmontering:
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

Inga interna data. Allt PENDING_INDEPENDENT_REVIEW.

# New round (coordinator, 3/10 00:55)
You delivered the list of requirements counted: COUNTED 4, MEASURED BUT UNCONNECTED 3, UNMEASURED 8. That's
the answer to the question, and it means that a theoretical basis is not reachable from today's twin — the
list is mostly a measurement order list. Say it straight out of the gate, it's a perfectly good outcome.

And you nailed the resolution issue more strongly than I dared hope: the heat dose has 48,87%
error coarse vs. 0,19% at 10 µm, and at the same density, mesh refinement yields deflections
of 0,0305, 0,0804 and 0,1346 mm against the 0,0400 limit — the coarsest passes, the finer
two fall. THE DESIGN VERDICT REVERSES WITH THE RESOLUTION. It should be listed as the main
result, because it makes dissolution a certification condition and not an ambition.

Obstacles: complete_native_safe_requirements = 0 and certified_medical_envelopes = 0. None of
the four calculated requirements has yet yielded a certified envelope through the generator.

Changed operation: take ONE of the four CALCULATED requirements all the way through the field engine
topology optimization to a case, at the finest resolution you've shown is necessary. A case certified
at 10 µm is worth more than four certified at a resolution where the value reverses.

And the force proxy: 0,721 N against the hold 1,514 N is 52,4% wrong. Say whether the error is the shape of the
proxy or its resolution — you've already shown that resolution can carry a factor of 254 in another quantity.

Falsifier: if none of the four COUNTED requirements can be driven through the generator, the interface between requirements
and geometry is missing, and then it is the interface that must be specified — no more requirements.

# New round (the coordinator, 3/10 01:20)
You got a real ability last night: the generator can EXCLUDE a design family. The lower
limit for displacement is 0,06001 mm against the 0,0400 mm limit, so 1,50× over at best,
so the family falls on geometry and not on parameters. Being able to exclude is cheaper
and often more valuable than being able to choose, and that should be the main result.

Obstacles: certified_medical_envelopes = 0 and medical_fleet_dimension = None, and the remaining
four families are UNKNOWN. In addition, the bending law of the source failed in 4 out
of 4 instances, so the law we intended to inherit does not hold in any tested case.

Changed operation in that order: (1) Bending law first. A law that falls in all four tried instances
cannot carry four UNKNOWN families. Say which condition is broken — slender beam, small deformations,
linear material — and which shape holds instead. (2) Determine the four UNKNOWN using the same exclusion
method that determined the first one. Four excluded families is a better result than an unsafe
certified enclosure. (3) The raft dimension is still None: if all five families are settled,
the dimension is the number of parameters separating the survivors, and then it is counted.

And remember the resolution result from r1: the value flips between 0,0305 and 0,0804
mm deflection at the same density. Each exclusion must therefore be specified with the
resolution at which it applies, otherwise it is not an exclusion but an artifact.

Falsifiers: if no inflection within published assumptions reproduces the four instances of the source,
the source is wrong and not our application of it — report it straight, it shifts the blame right.

# New round (coordinator, 3/10 01:50) — a possible common cause of your bending team miss
You reported that the source bending law fell in 4 out of 4 instances. A swarm job last night pitted the Q052's
delivered mechanical coefficients against FOUR independently published samples under ISO 13314:2011, with
three replicates being: predicted modulus of elasticity 5,705 GPa vs. measured 3,000 for one sample and 11,048
vs. 2,900 for two others, thus 1,902 to 3,810 times too high, vs. that of the cell OWN tolerance of 25 percent.
No parameters were touched, so it's not a fit issue. I have recalculated the ratios from the raw data
and they are correct. The report is in results/SWARM48_GRAPH_REPORTS/BT-FW48-AUTO-c040386b4454fa/.

Try if it's the same four. If your flexural miss and the folded stiffness have the same source, a coefficient
correction is sufficient for both, and then your excluded design family is possibly ruled out on too stiff
a material rather than on geometry. That would be an important correction of your own main result.

Important: the job found that the FUNCTIONAL FORM survives a two-parameter calibration,
so the defect is in the coefficients and not in the structure. This means
that your generator can use the shape but not the supplied numbers.

Falsifiers: if your bending law falls on samples other than these four, the errors are independent
and both must be corrected separately — and then your geometric exclusion is unchanged.

# New round (coordinator, 3/10 02:05)
My crossing hypothesis was wrong and your round showed it: bending law passes 4 out
of 4 in the second question set with a max error of 1,1e-06, and source_law_refuted
is set to False. The two errors are independent, and the law holds where it is
correctly applied. Good thing you tried it instead of inheriting my guess.

The hurdle now is that the UP certification is not closed: your upper limits for the displacement
of the families span up to 2,3e8mm, and only one family has a meaningful upper limit of
0,86mm. A lower bound can exclude but never approve, so without closed upper bounds
no family can be certified — which is exactly why certified_medical_envelopes is 0.

Changed operation: close the upper limit of the three families that pass the lower one. What
makes the boundary unbounded — a singularity in the geometry, an unbound parameter, or a
load distribution without a maximum? Each response has its own action. A single certified
enclosure with a closed upper limit is worth more than multiple exclusions at night.

And list each exclusion with the resolution it applies to. You showed yourself
that the value flips between 0,030 and 0,080mm deflection at the same density,
so an exclusion without stated resolution is not an exclusion.

Falsifiers: if no upper limit can be deduced within published assumptions, the generator
can only exclude and never accept — and then it is that characteristic that should
be reported as the true limitation of the fleet, not the number of claims.

# New round (coordinator, 3/10 03:50) — margin is less than resolution spread
You tightened the upper limit 21,5 times, which is real work. But the result is a knife's edge
that cannot be judged: the best family's 90 percent upper limit is 0,040057 mm against the 0,0400,
limit, i.e. 0,14 percent over. And in Round 1 you showed yourself that the same quantity
moves from 0,0305 to 0,1346 mm over mesh refinement, a factor of 4,4. **A margin of 0,14 percent
cannot be resolved when the resolution alone adds a factor of four.** You also set new_physical_resolution_refinement
to False this round, so the scatter is not eliminated.

Changed operation, and it takes precedence over everything else: pin the resolution before you judge.
Run a convergence study on the BEST family alone — three or four levels of refinement, and report
if the upper bound converges and to what. If it converges below 0,0400 the family passes and you
have the project's first certified casing. If it converges over, it falls. If it does not converge,
the deformation limit is not determinable from this model, and THAT is the result.

A family resolved at converged resolution is worth more than five sharpened limits at unknown
resolution. And state the order of convergence: a quantity that converges as h in the first
order behaves differently from one that converges in the second. If the order is lower
than expected, there is a singularity in the geometry, and then the limit should be read
as a percentile value and not as a maximum—the graded finding of peak values ​​at singularities
says that maxima diverge under refinement while percentiles converge.

Falsifier: if the upper bound grows monotonically with refinement instead of converging, the maximum
is not a certifiable quantity for this geometry, and the envelope must be defined on a percentile.
Then say it plainly — it's a method statement that applies to all families and not just one.

# New round (coordinator, 3/10 04:35)
The convergence study reversed the outcome and gave the diagnosis: the upper limit is now
0,5 percent BELOW the deformation limit, a family passes at 40 micrometers, and the order
of convergence is 2,13 and 1,95 so second order as expected. No singularity, so maximum
IS in principle certifiable. Good work, and it was worth waiting for the verdict.

The obstacle is now a single thing: the circle coverage variant does not converge. It gives 0,013522, 0,015200,
0,016528, 0,016502 — rising and turning — while the frozen geometry converges nicely to 0,014166. As
long as it does not converge, maximum_upper_mm cannot be closed and no casing can be certified.

Changed operation: diagnose circle coverage. A quantity that rises and then reverses during refinement
has either a coverage that is not nested between levels, or a discretization that changes
which point is maximum. The first is an implementation error, the second is a reading problem
and should then be read as a percentile. Decide which, and state how you decided.

And if it's the readout problem: define the envelope on a percentile
instead of a maximum and THEN certify. A certified envelope at a declared
percentile is a true result; an uncertified maximum is nothing.

Falsifier: if the circle coverage does not converge even with nested levels, maximum is
not a certifiable quantity for that geometry, and the percentile definition is not a choice
but a requirement. Then say it as a method statement that applies to all families.

# New round (the coordinator, 3/10 05:05)
The diagnosis was right and the remedy right: non-nested support geometry, now rebuilt to nested levels
16, 32, 64, 128. But you read the four ascending numbers as if the maximum doesn't converge, and it does.
I calculated from your own numbers: the increments are 0,001678, 0,001302 and 0,000775 with quotients
of 0,776 and 0,595, so a geometric extrapolation gives the limit approximately 0,01842 mm — so 46 percent
of the 0,0400. deformation limit. It is slow convergence by a good margin, not divergence.

Changed operation: calculate limit value itself instead of leaving maximum_gate as UNKNOWN. Do
a Richardson extrapolation over the nested levels, specify the cutoff with its uncertainty, and
set the gate. If the margin of uncertainty is below 0,0400 you have the FIRST certified casing
of the project, and that is worth the night more than anything else the lane has done.

And use the order of convergence you already measured: 2,13 and 1,95 for the frozen geometry. If the outer
patch converges more slowly than the second order, there is a weaker singularity there, and then the
uncertainty in the extrapolation must be greater. Enter the order of the outer patch separately.

Falsifier: if the Richardson extrapolation does not stabilize over the four levels,
the ratios 0,776 and 0,595 are not a geometric pattern but noise, and then a fifth
level is required before the gate can be set. Then say it instead of guessing.
