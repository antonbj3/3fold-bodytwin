# CX-PARAM — map what "parametric" can mean for the collaborator and what already exists (inventory, no new builds)

Background: the collaborator Rasmussen (the reference model) wrote that he wants to see how we model physiological geometry and "parameterization of shapes" (`references/john_context.md`, `notes/JOHN_REQUIREMENTS_AND_SOFTWARE.md` §4). The meaning is unclear. Anton is unsure how complete earlier work is, including scaling of muscle/bone to weight/load/force and the video pipeline. The goal is an honest map: interpretation → what exists → how complete → what is missing. No new models.

## Interpretations (add more if you find them)
P1 anatomical parameters (length, CCD, anteversion, torsion, tibia slope) → edit/instantiate
P2 statistical shape model (PCA modes) as parameters
P3 conditional instantiation from sparse measurements with uncertainty
P4 body parameters (height, mass, fat) → segment mass/inertia/muscle strength/fibre length (mass-scaling baseline laws; our alternatives)
P5 muscle parameters from the individual (PCSA, optimal fibre length, tendon slack length, muscle volume) from geometry/MR/volume
P6 tissue/material parameters (density → modulus, cartilage, layers)
P7 parameters as design variables for what-if analysis (osteotomy: CCD/AV → joint force), differentiability, emulator
P8 parameters from video/IMU (anthropometry, shape, landmarks)
P9 constructive/CAD-like parametric geometry (SDF primitives, field engine)

## Search, in this order (grep; `tool_find` does not index BodyTwin)
1. `~/projects/bodytwin/data/MECHANISM_ANCHOR_GRAPH.json` (claim/notes/scripts fields), `~/projects/bodytwin/scripts/msk/` (filenames, docstrings, content), `~/projects/bodytwin/docs/`, `~/projects/bodytwin/bt_memory/`.
3. Field engine and other workspaces under `~/projects/3fold-workspaces/` and `~/projects/3fold_staging/` (read only), for P8/P9.
Search words, in both Swedish and English: parametr, parameter, skalning, scaling, LengthMass, fat, fett, PCSA, fiber, tendon slack, muskelvolym, muscle volume, strength, styrka, SSM, PCA, mode, anteversion, CCD, torsion, osteotomi, what-if, emulator, gradient, differentiable, anthropometr, video, IMU, SDF, primitive.

## Deliverables: `results/CX-PARAM/` only
- `PARAM_MAP.md`, starting with `# CX-PARAM`. For each interpretation:
  - What exists: file path, result ID, node ID.
  - Status: working and tested / partial / script only / idea only / nothing found, each with evidence (test, number, audit ID).
  - The single most important missing piece, and how big it is.
- A final table: interpretation → our status → the reference model → relevance to the collaborator's publications (NATO 2024, femur population 2025, Sensors 2023) → recommended prio (1–3), with a one-line justification.
- `param_map.json`: the same data, machine-readable.

## Rules
- Only paths you have opened count as evidence. Never write "nothing exists" without listing which searches you ran. A negative claim needs a search log.
- Read-only everywhere except `results/CX-PARAM/`. No computation over 60 s. No emails, pushes or publishing. restricted model data and the collaborator's data stay internal.
