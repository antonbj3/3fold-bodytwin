# CX-DESIGNLOOP — a real design call: shape change → attachments/mass properties → rebuilt mechanics → forces → candidate check, with no archived force curve

Background: proof_lane's night priority §4 (research/COLLABORATOR_PRESENTATION_20260924/NIGHT_PRIORITIES.md). `bodytwin_core.demo` builds geometry but solves SAVED matrices. The separate trials CX-WHATIF2 and CX-D1PARITY actually rebuild the systems:
- `results/CX-D1PARITY/`: build_direct.py, direct_geometry.py, femur_massprops.py, massprop.py, full_step.py;
- `results/CX-WHATIF2/`: its builders.
Read both, plus `bodytwin_core/` (demo.py, solver, determined), before writing any code.

## Tasks (PREREG.md + sha256 first)
1. Build ONE callable function, `design_eval(shape_params, subject, task)`, in bodytwin_core:
   - (a) shape change (e.g. femoral neck length/anteversion, or tibial slope; choose a bounded parameter that the builders already support);
   - (b) new attachment points and mass properties;
   - (c) rebuilt mechanics (A, b, F0, contact projection) from scratch;
   - (d) forces: the LP band + min-effort + the parameter-free law where relevant;
   - (e) a candidate check: the capacity boundary/feasibility, and joint contact against a baseline.
   It must NOT fetch a finished force curve from the archive.
2. Run a new parameter case (at least 5 values across one parameter) on subject z001, plus a second anatomy or load (z009 or another task) as a check.
   - Report sources, frozen assumptions, changed quantities, and the actual build and solve time per step.
   - Verify one step with a full recomputation through the existing CX-D1PARITY chain (it must be identical within tolerance).
3. A search step: find the parameter value that minimises the peak knee (or hip) contact under feasibility, and check the optimum with a full recomputation.
4. Leave out the contact layer (CART-KNEE is not activated) and say where it would plug in.

Deliver RESULTS.md starting with `# CX-DESIGNLOOP`, results.json, pytest (including a parity test against CX-D1PARITY), and the timing. Internal data stays local; heavy runs go under bigmem.lock with 2 threads. lane runner has full permissions in the workspace; `~/projects/bodytwin` is read-only. Every outcome is a node that expands.
