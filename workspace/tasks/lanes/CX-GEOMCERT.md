# CX-GEOMCERT — Instance certificate for Geometry Manager

Builds on `results/N7c/geomgr`, `results/BT-C1b` (source keys, units, covariance), `results/BT-R2b` (validation gates), `results/BT-V7`/`results/BT-LV2` (arithmetic validity), `results/CX-JWGEOM`. A scientific status must not be derived from a technical certificate alone.

## Task
1. Define a machine-readable certificate for an individual instance: stable anatomical IDs, mm/m, body frame, right/left, oriented surface, non-inverted shape, attachments within tolerance to the bone surface, wrap and contact validity, covariance/provenance and export to the reference model/FE.
2. Build certification function + adapter to N7c instance and at least three actual individuals. C1b's sealed envelope must carry certificate status and source without doubling dependent observations.
3. Adversarial mutations: unit change, mirroring, vertex order, wrong segment, copied source, missing density and uncertain SDF sampling. Review which cases must give UNKNOWN instead of FAIL.

## PREREG and gate
≥1000 mutated cases; 100% of critical errors detected and falsely rejected valid instances ≤1%; deterministic hash identical across 2 runs. Countertest: disable one check at a time and show at least one leaking error; the reference model export nodes must be reimported without coordinate deviation >0,1 mm.

## Gemensamma regler
- Skriv endast under `results/CX-GEOMCERT/`, samt stora mellanresultat i `external_media`Read it. `tasks/NIGHT_PREAMBLE.md`, and relevant A-rows in `notes/RESULTS_INDEX.md`. `~/projects/bodytwin` and other sessions' workspaces are read-only.
- Before the first computation: `PREREG.md` and `PREREG.sha256` with hypothesis, data/selection, numerical gate, strongest baseline, countertests and error definition. Document `Builds on` with graph node and source files plus `Not redone`. Preserve negative results and `UNKNOWN`.
- Deliver `results/CX-GEOMCERT/RESULTS.md` with first line `# CX-GEOMCERT`, `results.json`, executable code, provenance/hashes and meaningful checks. Report both the number of valid and excluded units.
