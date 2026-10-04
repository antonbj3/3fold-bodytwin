# Launch announcement — dental, BodyTwin and manufacturing

You take over the research management for the dental project together with Anton. First, build an accurate overall picture of the project, the existing code, the graphs, the data and the researcher contact. Then work independently with concrete, testable steps. Communicate briefly and clearly in Swedish.

## The goal

Realize a coherent and verifiable chain from **real dental/anatomic geometry → field representation → materials and physics → simulation of intervention and function → generative design of construction and materials → fabrication → experimental validation**. Work iteratively: each step should be able to feed back requirements and results to the others.

**The whole chain should be developed as a coherent research graph.** Expand each step to concrete nodes and relationships: anatomy and tissues, materials and composition, physics and uncertainty, design methods, manufacturing techniques and verification. Link the factors affecting construction to mathematical models, design variables, goals and constraints. This is Anton's explicit direction. A first runnable case should provide feedback to the wider graph; the extent of the graph is not determined by what the first demonstration happens to support.

The project is being prepared for a possible collaboration with **Osayd Alawawda**, whose focus according to the project notes is prosthetics, digital dentistry, CAD/CAM, dental biomaterials and manufacturing. According to the existing brief, he first wants to see a real data set through the 3FOLD chain and then formulate the research question based on what the system can actually do. Check the original correspondence before attributing additional demands, resources or commitments to him.

The scientific ambition is to use the connection between the BodyTwin, the field engine and manufacturing for new useful relationships, representations or calculation methods. Search how geometry, material composition, biological and mechanical interfaces, load and manufacturing process influence each other. The result must be able to be compared with competent existing methods and provide a concrete basis for experiments.

**Material knowledge and simulation of the operations themselves are explicit innovation goals.** Seek new relationships between composition, structure, process and properties, better material models and testable material concepts. At the same time, develop the ability to represent how anatomy, tools and construction work together during a procedure. These tracks should get their own nodes, hypotheses, data needs and experiments and be linked to the generative design.

Exact prosthetic construction, materials, load cases, manufacturing process and research question are still open. The mapping is partly called "implants", but that title does not determine the final delimitation of the project. General BodyTwin and the collaborator's biomechanics questions have their own workspace.

## Keep the stress point

**Aim for a result that deserves to be tested by a scientist and built in a laboratory. Keep pressure on the load bearing restraint. When you reach a limit: make it measurable, open the mechanism behind it and try a way to move it. When an experiment fails: use it to formulate a sharper hypothesis. When it succeeds: subject it to a stronger counter test. Let each round leave us with better ability, clearer knowledge or a precisely defined gap.**

The stress point is what currently limits a credible end product or conclusion. It can be in geometric representation, error propagation, material data, interface model, load case, numerical method, manufacturing deviation or validation. Several constraints may be linked. Identify which investigation best separates their causes and what actually becomes possible if successful.

Drive your own hypotheses from first principles, our previous results and the primary literature. Break down supporting ideas recursively: **idea → mechanism → equation → operation → representation → assumptions and physical basis**. End each branch in a defined building block, verified knowledge, or explicit knowledge gap. Then try new combinations across domain boundaries. Search conceptually distant problems that share a real mechanism; show the equation mapping and what the transfer requires.

Work with high demands on utility, accuracy and full calculation cost. A positive result must withstand comparison for the same task, quality and coverage. Keep negative results and UNKNOWN. Do not change measures or thresholds afterwards to make an outcome positive. Breakthrough is the ambition; the evidence determines what we are allowed to say.

## Expandera grafen genom hela kedjan

Start from existing BodyTwin and manufacturing nodes and their underpinnings. Map overlaps, gaps and ambiguous identities before creating new nodes. Each area below should be broken down into usable building blocks with explicit connections to neighboring areas:

| Area | What can be represented in the graph |
|---|---|
| Anatomy and Geometry | anatomical regions, segmentation, shape variation, surfaces, volumes, thin layers, interfaces, coordinates and geometric uncertainty |
| Biological tissues | including gums and other relevant soft tissues, dental tissues, bone and supporting tissues; which are included is determined per construction and question |
| Material and composition | properties, composition, structure, spatial variation, material laws, parameters and calibration background of natural tissues and engineered materials |
| Material innovation | hypotheses about new compositions, microstructures, gradations, surfaces and interfaces; relationship between process, structure and properties; identification and testing of material models |
| Load and physics | load cases, motion, boundary conditions, contact, deformation and relevant candidates for time dependence, damage, transport and biological processes |
| Operations and interventions | condition before, during and after an intervention; tool–tissue–construction, movements, contact, changing geometry/topology and material state; the order of the elements and their verification basis |
| Generative design | form, topology, material distribution, composition and process choice as possible design variables; targets, constraints, algorithms and sensitivity |
| Manufacturing | relevant additive and subtractive methods, process parameters, orientation, accessibility, finishing, tolerances and the difference between intended and manufactured construction |
| Verification and experiment | measurement method, specimen, reference solution, calibration, separate validation, uncertainty, reproducibility and comparison with existing method |

The table is an expansion mission. It does not imply that these functions already exist, that all listed physics models are applicable, or that a single simulation must contain everything.

For each building block: specify identity, type, task, mechanism/model, input/output, entities, domain/region, scope, source, state of knowledge, uncertainty, and a way to test it. Missing data must be explicit. Distinguish concepts, hypotheses, variables, equations, implementations and experimental evidence. Relationships must state their meaning: for example, membership, model dependency, parameterization, interface, design constraint, or evidence. An assumed impact must be labeled as hypothesis until it has support.

The finished reading view has an existing type contract. New research nodes and new relationship types require a delineated proposal for additions in consultation with the graph session. Put the expansion under `notes/`, `tasks/` and `results/` with traceable connection to existing ID until loading and validation is found. Never type directly into generated graphs. This is a concrete next development step; the new domain nodes are not yet entered by this start message.

## Bring the factors into the generative design

Follow every relevant factor from underpinning to mathematical role. Determine if it is a design variable, a material/model parameter field, a boundary condition, a simulated condition, a target, a constraint, or an uncertainty. A factor can affect several links, which requires explicit relationships.

An organizing scheme is:

```text
x       = selectable geometry, material distribution/composition and manufacturing choices
theta   = anatomical, material and load-related conditions
u       = the states calculated by the selected simulation
F(u, x, theta) = 0                 the model's equations
g(u, x, theta) <= 0; h(u, x, theta) = 0   requirements and constraints
J(u, x, theta)                    several explicit design objectives
```

This is a structure to pose the problem, not a finished dental model. Contact and other conditions may require special differences or complementarity conditions. Discrete material/process choices and continuous shape parameters may need different search methods. Do not assume that all branches are differentiable or that goals can be weighted together without justification.

In particular, investigate the linkage **process selection → fabricated geometry and material properties → calculated response → design decision**. Also track geometry and model errors through the chain. When data is missing: formulate scenarios, intervals or a data collection task; do not set a probability distribution or clinical limit without support. Document which factors are included in each run, which have been delineated, and what could change the design decision.

For interventions, the schedule needs to be expanded with time, tool movement and changing geometry, contact and material condition. Represent the moments and their state transitions explicitly. Investigate how the simulated course of the intervention can provide goals or limitations for the construction's form, materials and accessibility. A final state of equilibrium is not sufficient to describe the entire process.

## Material knowledge as an innovation track

Starting from mechanisms, investigate what changes in composition, structure, grading, surface, interface or processing could produce a useful property. Formulate expected effect, possible goal conflicts, comparison material and an experiment that can disprove the hypothesis. Also seek innovation in how material properties are identified and modeled, including biological tissues.

Link material candidates to a possible manufacturing process and measurable properties. Maintain the distinction between the material intended and the material actually produced. An optimization outside the model's calibrated range provides a candidate to investigate; it does not establish that the material has the predicted properties. Build the feedback **material hypothesis → fabrication/sample → measurement → model update → new design** in the graph.

## Simulate the procedures themselves

Map relevant interventions and their sub-moments with connection to the selected anatomy and construction. Examples to assess are tissue handling, incision, drilling and insertion or placement of a construct. The examples are candidates for the domain map; no single intervention has yet been chosen as the first experiment.

For each step: define geometry and material before the step, the tool's representation, controlled and unknown variables, interactions, possible changes in geometry/topology and which state is left for the next step. Try which mechanical, thermal, transport-related or biological phenomena are needed for the question. Specify what can be measured during the course and what requires a separate model and validation afterwards.

Look for innovation both in the physics models and in the field engine's ability to represent and calculate the change. Link the results back to material selection, generative design and manufacturing. Interactivity or real-time requirements need an explicit use and measured trade-off against accuracy.

## Tissues and BodyTwin

**Yes, gums and other relevant tissues are included in the intended modeling work.** For each selected case, the geometry needs to be linked to region identity, properties/material layers, interfaces, loading and verification data. The model's level of detail is chosen according to what it should be able to answer. Distinguish geometric representation, mechanical response and possible biological time development.

Inventory what BodyTwin actually has in executable code, what only exists as a graph node and what remains. Follow its anatomy, regions, material coupling and uncertainty workflows where applicable. The fact that a component is in BodyTwin's graph does not mean that it already has a calibrated simulator. Check anatomical applicability for each construct; do not transfer the entire set of regions or boundary conditions between different dental cases of routine.

## Read in this order

Arbetsmapp: `local_path`.

1. `../AGENTS.md`, `START.md`, `../README.md`, `../FORMAT.md` and `../SETUP_STATUS.json`.
2. `external_research_path`. Read project intent and open questions carefully. Scientific figures and broad conclusions in older mapping need their own source checks before use.
3. `notes/ADDENDUM_DATA_COLLECTION_2026-09-22.md`, then the mapping directory's `_STATE_2026-09-17.md` and `_MASTERINDEX.md`. The addendum preserves anton-02's partially clipped submission and points out verification needs. Then read relevant sections in `01_imaging_geometry.md`, `02_biomechanics_fe.md`, `07_manufacturing_methods.md`, `08_materials_compositions.md`, `15_uncertainty_meshquality.md`, `17_fabricate_test_pairs.md` and `18_dental_ceramics_data.md` as needed. The master index does not automatically cover everything in the later mapping files.
4. `local_path`, then the workspace's `VALIDATION.json`. Use `./graph status`, `./graph show --id <ID>` and directed questions to orient yourself. Don't read tens of thousands of graph records or the entire corpus at random.
5. Read the real code that a candidate chain needs via `references/current_bodytwin` and `references/field_engine`. Differentiate between delivered implementation, bounded experiments and planned features. Then read the original correspondence with Osayd if available, without sending or changing anything.
6. Read `notes/ADDENDUM_COMPUTE_CELLS_2026-09-22.md` for published computational cells, private BodyTwin experiments, cad-to-simulation's cell registry, and later engine trials. Inventory relevant candidates and their actual contracts before new implementation is planned.

## What's available now

**The reading view of the joint graph is built and checked.** It is a basis for the research, not a certificate that all the included claims are verified.

- BodyTwin: 3 950 anchor nodes, mainly OPEN.
- Crafting: 97 anchor nodes from CS-I. It is the correct number in the dependency layer. The manufacturing's comprehensive basis is also in the relationships and measurement logs.
- Dental together: 4 047 anchor nodes, 624 dependencies, three constraint meshes with 48 variable names and 59 relationships.
- 13 log files: 3 799 line occurrences → 2 772 unique fold entries. Identical copies are merged with all origins preserved. A fold record does not have to be an independent physical measurement.
- Source index, explicit citation records, source files and stress map are connected as separate typed layers. A citation is not verified just because it is imported.
- `MERGED_GRAPH.json` shows the entire index; `GRAPH.json`, `RELATION_GRAPH.json`, `MEASUREMENTS.json`, `SOURCES_VIEW.json` and `EVIDENCE_BINDINGS.json` contain the details.
- `OPEN_INTEGRATION_ITEMS.json` shows unbound references, conflicts, deviating source values and data required for further numerical use.

The original files remain unchanged. Work through the readers and keep new tasks/results in `tasks/`, `results/`, `notes/`. Do not edit generated graphs or source graphs directly. The engine is version locked according to `../ENGINE.json`; the second graphing session develops the engine further. Coordinate API needs through concrete contracts and reproducible examples.

The grids contain 26 numerical marginal values but lack explicit σ. Units and margin definitions need to be reviewed prior to weighting. The same name does not automatically mean the same variable; two copies of a result do not provide two independent pieces of evidence. Common source provenance, area of validity and contradictions must accompany each connection. The engine's structural `rank` is an advisory heuristic.

`seed_graph.json` is an index from operator text to commits, and `substrate_graph.jsonl` is machine telemetry. They should not count as additional manufacturing skills. Older CS-D graphs are replaced generations.

## Data and actual preparedness

`data/corpus` leads to `external_mount`. The move is complete: **52 005 files, 113 330 792 472 bytes**, checksum verified. The old path works as a symlink. Separate FE geometry has also been mapped under `external_mount`; check actual structure and manifest before selecting files.

The mapping describes, among other things, Teeth3DS+, Open-Full-Jaw, OpenMandible, material and fatigue data as well as article and patent traces. The verification of the move applies to the existing files. It does not say that all planned downloads are complete, that the data belongs to the same individual, or that license and usability are reviewed for each new purpose. Check the identity, content and terms of use of the selected data.

Older mode files describe full disk and ongoing harvests. That data is historical: the storage move is complete and the harvesting processes are stopped. Read current catalogs and `external_research_path`. Plan additions based on actual gaps and available space; do not start old watchdogs out of routine.

The previous searches report gaps in coherent patient data and experimental reference. Treat them as results from that search scope. Verify the slot that will carry a new project before claiming that data or a method is missing in general.

## Working methods and coordination

- Field engine innovation work is left to a separate session; orientation file `local_path`. The second graph session owns the graph engine development. Dental owns the concrete application chain, its requirements and verification. Examine existing results before ordering the same experiment again.
- Start with a reproducible chain for an actual data set. Check incremental entities, geometry, region/material identity, boundary conditions, reference solution and export. Let the measured chain show where a research effort can make a difference.
- Separate source data, derivation, simulated result, measured result and hypothesis. Independently review positive results before posting them as supported.
- If you delegate: provide delimited tasks, clear file owners and a limited concurrency. Subagents should not start their own agent trees. Previously unchecked recursion caused a quota spike.
- Read relevant existing reports before new runs. Choose experiments that decide something. Calculate the full cost and compare to strong, simple baselines.
- Text, graphs and figures must be comprehensible and have correct units, axes and boundaries. Document run command, input, version and raw output.
- Keep private information private. Do not send emails or publish anything without Anton's express instructions. He writes Anton's correspondence himself.

## First re-reporting

After the orientation, give Anton a brief coherent picture:

1. What we are trying to enable and what Osayd has actually requested, with the source clear.
2. How the entire chain should be expanded in the graph, what is already in nodes or executable code and what additions to the graph contract are needed.
3. Which real data set is suitable for a first test and where the first load bearing stress point is located, with basis and residual uncertainties.
4. Two or three concrete continuations, one of which is recommended, with hypothesis, connection from factor to mathematical design role, strong baseline, decisive measurement and what a negative outcome would mean.

Save the orientation in `notes/ORIENTATION.md`. Tune this image with Anton before locking a large driving wave or a specific prosthetic research question. Continued work must be able to explain what obstacle it attacks and what the result would enable.
