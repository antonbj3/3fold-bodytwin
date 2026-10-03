# CX-PARAM-DENTAL — parameters that generate design: the dental chain → transfer to BodyTwin and the collaborator (inventory, read-only)

Background: the collaborator wants to see "parameterization of shapes". The dental project (`local_path`, owned by another session, READ-ONLY for you) has parameter → design chains: crown (`results/DESIGN_crown`, `DS_CROWN_*`), implant/screw and LPBF fatigue (`results/K1_lpbf_implant_fatigue`, `BQ3_SCREW_MODULE_MECH`), lattice (`BQ3_LATTICE_CONST`), PDL (`BQ3_jobs/PDL_*`). Start in `START.md`, `GRAPH.json`, `CONSTRAINT_NETS.json`, `NEXT_ACTIONS.json`, and `tasks/AGENT_QUEUE.md`.

## Question
Which parts of the dental parameter → geometry → load → design loop exist, how complete are they, and which can be moved to BodyTwin as "parametric design from the individual's anatomy"?
Examples: an osteotomy wedge/plate from CCD/AV targets, an implant/stem shaped to the individual's femur, orthoses/insoles, a lattice with a stiffness gradient against stress shielding.

## Deliverables: `results/CX-PARAM-DENTAL/` only
- `DENTAL_PARAM_MAP.md`, starting with `# CX-PARAM-DENTAL`.
- For each dental design chain:
  - the parameters;
  - the generator (file path);
  - which load/strength model evaluates the design;
  - status (tested with a number / partial / script / idea), with evidence (file, test, result);
  - known errors (read `CORRECTED_CLAIMS.md`, `AUDIT.md`).
- A transfer table: dental component → BodyTwin counterpart (what exists there: `results/N7c/geomgr`, `results/CX-SOLVER2`, `results/BT-XF4`, `notes/RESULTS_INDEX.md`) → what is missing to close the loop → relevance to the collaborator (femur population/osteotomy 2025, NATO 2024 twin chain → injury risk) → prio 1–3.
- A proposal for the three most valuable joint experiments, each with a falsifiable criterion and null model. Proposals only, nothing built.
- `dental_param_map.json`: the same data, machine-readable.

## Rules
- Only paths you have opened count as evidence. Negative claims ("nothing exists") require a listed search log.
- Write nothing in dental or in `~/projects/bodytwin`. No computation over 60 s. No emails or pushes. No patient or personal data in the report.

## Entry points from the dental session (anton-69, 24/9) — start here
- Crown (most complete parametric shape generator): `cells/design/crown_design_geometry.py` (x = material, occlusal thickness, fissure rounding, cement gap, milling radius/compensation), `crown_design_problem.py` (x/θ/g/J + Pareto), `crown_design_pipeline.py`, `crown_case_sts3d.py` (parametric preparation on a real tooth: taper angle, reduction, chamfer). Reports: `results/DESIGN_crown/REPORT.md`, `results/CONTACT_damage/REPORT.md`, `results/LICENSE_gate/REPORT.md` (shape sensitivity F50 1.2–2.5× between teeth).
- Implant: `cells/physics/implant_thread_fe.py` (parametric thread, root radius), `iso14801_assembly_fe.py`; placement `cells/design/implant_site.py`, `design_implant_v2.py` (x = D, L from GUDID catalogue, pose, depth; set-valued robust rule), `design_margin.py` (uncertainty-aware margin, in progress).
- Overview: `notes/STATE_2026-09-24.md`, `notes/DENTAL_METHOD.md`.
- NOTE (dental): several results are REFUTED against anchors, see each REPORT. The parameterization works, but the predictions are not validated. Report this separately: generator status vs prediction status.
