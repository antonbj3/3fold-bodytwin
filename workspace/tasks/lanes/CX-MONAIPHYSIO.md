# CX-MONAIPHYSIO — how exactly can BodyTwin use NVIDIA's "monai-physio"?

Anton: "monai-physio — put a Sol agent on investigating exactly how we can use this from NVIDIA."

## Tasks
1. **Identify exactly what "monai-physio" is.** Search the web and GitHub: the NVIDIA/Project-MONAI repos, the MONAI Model Zoo/bundles, NVIDIA blogs, papers and HuggingFace. Report:
   - the exact name, repo and URL;
   - the licence (commercial use?) and the release date;
   - what it does: physiological signals? physics-informed modelling? medical image→physiology?
   - which models/weights exist;
   - its dependencies (MONAI version, CUDA, GPU memory);
   - the input and output formats.
   If there are several candidates, list them all, together with how sure you are about each.
2. **Map it onto BodyTwin:** which of our nodes/modules could it feed or replace? Read `notes/RESULTS_INDEX.md` (the last 150 lines), `bodytwin_core/`, `results/HYPERREALISM_ROADMAP_20260924/FUNCTIONALITY_SEARCHSPACE.md`, `results/HYPERREALISM_ROADMAP_20260924/SEED_EXPANSION_20260925/DEVELOPMENT_EXPANSION.md`, and the graph (`./graph working rank --query '<word>'`). Name concrete touch points, for example:
   - CT/MRI segmentation → geometry/attachments;
   - physiological time series → the models' parameters;
   - image-based registration (our missing implant↔marker registration, A1785);
   - the surrogate/speed track.
3. **Try it for real** (the licence is not a gate): install it in a separate venv under /media/anton/sdc1-tmp/bodytwin/monai_physio/ (NOT on /, the disk is tight). Run the smallest example or bundle on public test data. Measure the time and the GPU memory on an RTX 5070: a single process, NEVER run `nvidia-smi -q`; use only field queries such as `nvidia-smi --query-gpu=memory.used --format=csv`. Do not install anything system-wide.
4. **Recommendation per touch point:** use it directly / use it as a component / do not use it, with the reason, the integration cost, and the first concrete step.

Deliver `results/CX-MONAIPHYSIO/RESULTS.md` starting with `# CX-MONAIPHYSIO`, results.json, and any run log. Internal data (the collaborator, GC, restricted model data) must NEVER go to any external service. lane runner has full permissions in the workspace; `~/projects/bodytwin` is read-only.

## Addition (Anton): less focus on licensing
Licence: just one line per candidate (what it is), not an obstacle to trying or recommending it. Put the effort into what it does, how it connects and running it.
