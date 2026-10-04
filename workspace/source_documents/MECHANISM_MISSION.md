# MECHANISM MISSION — the vision, the architecture, and why it is built this way

Read this after `MECHANISM_STARTUP.md`. It is the "why" behind the "how"; the mechanics live in STARTUP + COORDINATOR.md.
It is long on purpose — a fresh amnesic mechanism (or a weaker model) should be able to reconstruct the whole picture
from here without anyone who remembers the original conversation.

---

## 1. THE ONE-SENTENCE MISSION
Build a **hyper-real, certified digital twin of the human body/biology** — cell → tissue → organ → system → whole —
on cad-to-simulation's identifiability + certification pipeline, where every claim is either **measured and cert-anchored**
or explicitly **held in quarantine**, and the hardest un-measurable states (growth, aging, cancer, regulation) are
**cornered by decorrelated substrates until the contradictions reveal the truth.**

## 2. THE FOUR SUBSTRATES — mechanism is one, and it COUPLES to the others (this is not optional decoration)
The larger programme is four coupled digital twins. They are not separate projects walled off — the coupling *is* the
certification method (over-determination across substrates that MEET on a shared hidden state):
- **robot_lab** — the factory / organization / world twin (the outermost frame).
- **cad-to-simulation (CS)** — the physics + compute + vision engine. This is what mechanism forked. Its **vision/optics**
  (the "eyes", the lens/projector reconstruction lane) is the twin-optics forward model that mechanism's **eye biology**
  (retina / cornea / optic nerve) couples to. Its **physics verticals** (battery, thermal, LPBF) are the **engineering
  cert-twins** of biological energy/heat. Its **compute core** is what a mind twin runs on.
- **Hollow** — the mind twin (persona / theory-of-mind / behaviour). Lower priority; mechanism couples to it where the
  face/nervous-system surfaces (facial muscles → expression, neural signal → intent decoding — DANDI tetraplegic BCI).
- **mechanism (this repo)** — the body / biology twin.

The unifying picture: an **identifiability/certification spine** (CS / robot_lab / mind-physics) crossed with a
**generative spine** (content_platform → game / persona), meeting at hollow. mechanism sits on the certification spine, feeding the
generative one with certified biology.

## 3. THE BIO SUBSTRATE LADDER (the build order, bottom-up)
- **Node 0 — the energy-closed cell.** The cell as a free-energy transducer. The 1st law of thermodynamics is the
  MEASURABLE reduction anchor (calorimetry + respirometry close energy in/out, model-free); the 2nd law is the
  constraint; efficiency η(t) is the state that drifts. This is the foundation because everything above is energy-gated.
- **Aging = η(t) drift** — a cell losing energy efficiency over time. The great question, framed measurably.
- **Cancer = growth escaping regulation and seizing the energy budget.** Growth / DNA / epigenetics are barely
  certifiable directly — so they are exactly what we CORNER with decorrelated substrates until contradictions localize truth.
- **Primitives:** program (DNA) × control (signalling) × fuel (energy); connective / vasculature / muscle "fill between
  the joints". Each primitive has a clean measurable **engineering cert-twin** (see §6).
- Then tissue → organ → system → whole body, and the same operations cross-species (a conserved-mechanism cert layer:
  canine degenerative myelopathy = naturally-occurring SOD1 ALS — same gene, different name).

## 4. THE CERTIFICATION METHOD — over-determination (this is the whole game)
An un-measurable hidden state (is this knee pose real? is this cell's regulation intact?) is **surrounded by ≥2
DECORRELATED legs** — different physics measuring the same state. **Agreement certifies; contradiction localizes truth.**
Requirements that make or break it:
- The legs must be **genuinely decorrelated** — different underlying physics/data. Two legs sharing a model or a
  training set are **common-mode** (fake decorrelation, n_eff collapses ≈ 1) — call that out ruthlessly.
- Every cert needs an **ANCHOR**: an INDEPENDENT measured ground truth, HELD OUT of the fit, that the prediction is
  checked against — never a tautology gate (a gate that checks the model against itself certifies nothing).
- **Regime-gate** every claim (state the regime where it holds; a fixed threshold on a noisy statistic is regime-blind).
- **honest-negative = PASS** — a booked negative with a mechanism is a redesign spec, generative, not a failure.
- A **false-positive and a premature-negative are the SAME failure** — both are an unforced adversary of the truth.

## 5. THE RATE-DISTORTION FRAME (how data, models, and gaps relate)
The whole intake resolves compression-vs-fidelity via the cert, per bit:
- **GROUNDED** = irreducible information = the source's true bits → the **datasets** (measurements). Store at resolution.
- **FILLED** = prior-regenerable → the **models** (Cellpose, AlphaFold, scGPT, OpenSim…) are the DECODERS. A model's
  *verified* fidelity is the compression frontier. ⚠ A model can regenerate PLAUSIBLY without being FAITHFUL (AlphaFold-3
  is a diffusion sampler — the archetypal "better hallucinator"); FILLED-compression is safe ONLY where the model's fill
  is cert-verified against an independent measurement. A model can NEVER certify its own FILL (marking its own homework).
- **FENCED** = unobserved → the **honest-negatives** (our 84 FENCED holes). Don't store; flag the gap; don't re-probe.
- **The cert PROMOTES a datum out of quarantine (BT-HOLD) → GROUNDED** only via **same-sample multimodal pairing**
  (patch-seq, microbiome-metabolome, radiogenomics) — the same hidden state measured by two decorrelated modalities.
  That paired data is the thin, high-value axis; most systemic data is "summary-open / individual-gated".

## 6. THE ENGINEERING CERT-TWINS (each bio primitive has a clean measurable analog — the coupling into CS)
- metabolism ↔ **battery / Semenov** (energy transduction + thermal runaway) — couples to CS's P16 battery vertical.
- thermoregulation ↔ **calorimetry / thermal physics**.
- membrane potential ↔ **Nernst / electrochemistry**.
- chloroplast / photoreceptor ↔ **solar cell** (photon absorption).
- morphogenesis ↔ **Turing reaction-diffusion**.
- cell osmosis ↔ **desalination**.
- eye biology (retina/cornea/optic nerve) ↔ **CS's optics/lens reconstruction** (ocular optics = the twin-eye forward
  model; this is the strongest, most-developed coupling — the eye is where mechanism and CS's vision lane meet).
- musculoskeletal DOF ↔ **CS's SMPL-X / hand articulation** (SMPL-X gives every joint a uniform 3-DOF; the anatomical
  model says a knee is ~1 flexion-DOF — the cert says which DOFs are real, observed, and realizable).

## 7. THE FIRST CLOSING CELLS (where to start — the current build-goals)
- **KNEE-CELL** (first closing cell). Hidden state: is a knee pose/motion real and dynamically realizable? Four
  decorrelated legs: surface-geometry (SMPL-X/scan) · anatomical DOF+ROM (OpenSim knee vs SMPL-X 3-DOF null) ·
  moment-realizability (muscle MTU + inverse dynamics) · cartilage/contact. ANCHOR = **Grand-Challenge in-vivo knee
  contact force** (instrumented implant — the moment leg PREDICTS it, check against measured). This is the first cell
  because it is the ONE joint where all four legs have catalog data AND an external measured anchor exists to verify
  against, not just assert.
- **NODE-0-ENERGY** (foundational cell). Hidden state: the cell's energy balance / η(t). ANCHOR = 1st-law closure via
  calorimetry + respirometry. Engineering-twin: battery/Semenov.
- The **coupling edges** into inherited CS nodes (eye↔optics, metabolism↔battery, thermoreg↔thermal) — over-determination
  across substrates.

## 8. THE ASSET YOU INHERIT — the 449-source catalog (26 waves)
`WAVE_PLAN.md` + `data/mechanism_catalog/sources_wave*.json` = ~449 verified dataset sources across every substrate leg
(energy/metabolism, membrane, cell-division/cancer, DNA/epigenetics, muscle/nerve/ALS, cardiovascular, immune/autoimmune,
connective/bone, EM/ultrastructure, eye/optic/blindness, deafness, paralysis/stroke, patch-seq/multiomic pairing, organ
systems, cancer/aging/reference-networks, microbiome, brain connectome, …), plus model catalogs and literature/patent
corpora. All band **BT-HOLD** (quarantine) until a cert promotes them. Honest-negatives are logged as FENCED (the map of
what has no open source — brake the search, don't re-probe). The 26 waves WERE a manual cartographer round — this catalog
is the discovery/pricing layer's output, already done. `data/MECHANISM_SOURCE_REGISTRY.json` + `MECHANISM_METADATA_LAKE.jsonl`
+ `MECHANISM_HOLES_SEED.jsonl` are its translation into the pipeline's intake format.

## 9. THE STANDING DISCIPLINES (the ones that cost the most when skipped)
- **Literature = QUARANTINE claims.** "Peer-reviewed" is a nominal flag; a paper is a claim ABOUT a measurement, not
  the measurement. Physical anchors outrank citations. (Medical science is often extremely low-resolution — many hidden
  variables; treat published values as leads to verify, never ground truth.)
- **The placebo effect is a DRIVING FORCE, not a small quarantine** — a claim's truth-grade and its methodological role
  as a generative engine are different axes; don't flatten them.
- **"Sort late":** low-evidence and unconventional items (including treatments without peer-review support) stay
  BT-HOLD and flagged — gather the value, don't delete early; the cert decides, not a prior filter.
- **Design ≠ wiring** (committed-but-inert): a folded design/self-test is not a working gate until a colortest over the
  REAL chain passes. Verify things are functionally live, not just committed.
- **Zoom ladder ⊥ fidelity ladder:** zoom = generate detail on demand (lazy materialization); fidelity = certify. More
  claims ≠ more reality; reality comes from closing a cert against a measured anchor.

*The mission in one line: corner the un-measurable with decorrelated, externally-anchored substrates until the
contradictions localize the truth — cell to whole body, coupled into the physics/vision/mind twins — and keep every
certified bit, and every honest gap, in this repo.*
