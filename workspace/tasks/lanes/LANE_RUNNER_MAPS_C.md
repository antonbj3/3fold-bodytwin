# Task: maps/operators, geometry and unexplored ground (Sol C)

You're a lane runner. `lane-model`Agent with full mandate from the coordinator. (Anton)You own all the decisions. Work **oberoende** av andra agenter — skriv din egen tolkning.

## Goal
Build the **method level** that the boundary work rests on: define the **maps/operators** that connect subject spaces, their **geometric properties**, and map **the unexplored surface** (cross-subject ports that have no map yet). This is what allows diseases (Sol B) to be connected and solved, and the interface (Sol A) to become dense.

## Basic view (start from this, improve freely)
A "map" is an **OPERATOR/MAPPING**, not data — and a **geometric concept**: a structure-preserving mapping (morphism) between two representation spaces. Examples:
- Theta: x -> theta  (omics → boundary conditions/parameter set)
- M:     D_d -> y    (disease space → observable phenotype)
- graph <-> field <-> molecular (the T41 spectral unit is already such an operator)
- T43 contradiction engine, T44 method transfer, T45 operator composition, T48 emergent combination — treat as operators.

## Task
1. **Define maps explicitly as objects** with properties: domain/codomain (spaces), injectivity/surjectivity, which INVARIANTS are preserved, dimensional reduction, and how UNCERTAINTY/ERROR is transported through the mapping.
2. **Geometry:** the spaces (transcriptome, parameter, field, graph spaces), their metric/topology, and the operators' algebra (composition, inverse where it exists, isometries/near-isometries).
3. **Subject map:** list which subject spaces are connected (molecular signalling, systems biology/ODE, continuum physics/boundary conditions/transport, geometry/topology, statistics/information/identifiability, graph theory/representation).
4. **Unexplored ground (the explored surface):** identify ports between subjects that have no map yet. Candidates: metabolomics<->kinetics, epigenetics<->time constants, mechanotransduction<->transport, microbiome<->boundary conditions, pharmacokinetics<->BC. Each port is described as a method object with a discriminating test.
5. **Identifiability:** which maps are invertible/determinable from observable data, which are fundamentally underdetermined, and which smallest experiment decides.
6. **Prioritisation:** rank maps/ports by downstream capability and how many diseases/mechanisms they unlock.

## Inputs
- The graph and its working view: `./graph working rank`, `./graph working packet --id <ID>`, `notes/GRAPH_WORKFLOW.md`.
- T41–T50 in external_research_path and the v6 report BODYTWIN_MUTATIONS_IMPROVED_v6.md.
- First principles: external_research_path

## Outputs
1. `external_research_path` — insightful overall view: spaces, maps, invariants, error transport, unexplored ground.
2. Machine-readable list of maps/ports (JSON) with prioritisation.
3. `tasks/lanes/LANE_RUNNER_MAPS_RESULT.md` — decisions, what was created, what was NOT done.

## Constraints
- No fabricated data/citations; no automatic evidence admission. Follow AGENTS.md and START.md.
- No emails/pushes/credentials; no SEED results deleted.
