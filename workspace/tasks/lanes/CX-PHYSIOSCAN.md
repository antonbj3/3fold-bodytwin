# CX-PHYSIOSCAN — find projects SIMILAR to NVIDIA's monai-physio: open physiology/digital-twin/biomechanics frameworks we can use

Anton: "also send a Sol agent to look for similar projects, because it's the first time I've heard of it." CX-MONAIPHYSIO is investigating monai-physio itself. This lane searches BROADLY for what is similar or better.

## Tasks
1. **Search the web and GitHub broadly (2023–2026, weighted to the newest).** Categories:
   - AI/physics frameworks for physiology and digital twins (NVIDIA: MONAI, Clara, Holoscan, PhysicsNeMo/Modulus, Omniverse/Isaac for medicine; others: Google, Microsoft, Meta, DeepMind);
   - open human-body/organ simulators (SOFA, OpenSim/Moco, MuJoCo MyoSuite, Simbody, FEBio, Artisynth, Chaste, OpenCOR/CellML/Physiome, OpenCMISS, SimVascular, BioDigitalTwin/VPH projects, the Living Heart);
   - learned surrogates/foundation models for medical imaging → physiology/biomechanics (segmentation to meshes: TotalSegmentator, nnU-Net, MedSAM; muscle/bone models from images);
   - datasets + toolkits that combine anatomy and physiology;
   - newer 2025–2026 projects and papers that we might have missed.
2. **For each candidate** (aim for 30–60, deduplicated): the name, URL, what it does in one sentence, activity (latest commit/release), GPU/CPU, input/output, and ONE line on the licence (not a focus).
   Then state the concrete touch point in BodyTwin. Read `notes/RESULTS_INDEX.md` (the last 150 lines), `bodytwin_core/`, `results/HYPERREALISM_ROADMAP_20260924/FUNCTIONALITY_SEARCHSPACE.md` and `SEED_EXPANSION_20260925/DEVELOPMENT_EXPANSION.md`, and grep `results/HYPERREALISM_ROADMAP_20260924/PUBLIC_MODULE_INVENTORY.tsv` and `EXTERNAL_CANDIDATES.json` for what is ALREADY inventoried. Mark each candidate as new or known.
3. **Top 10 by value to BodyTwin / integration cost:**
   - one line on why each;
   - for the top 3, a runnable first step: install in a venv under external_media (not on /), and run the smallest example on public data if it takes under 20 min. Single process on the GPU, NEVER `nvidia-smi -q`.

Deliver `results/CX-PHYSIOSCAN/RESULTS.md` starting with `# CX-PHYSIOSCAN`, and `candidates.json`. Internal data must never go to external services. lane runner has full permissions in the workspace; `~/projects/bodytwin` is read-only.
