BT-DAT-Q077

## Built on

`AGENT_EXIT=1` with an empty `agent.log`: the earlier session left only `inputs/`, so this run reuses `inputs/Q077_model.py` (`Parameters`, lines 21-74) as sole authority on which numbers a dataset may touch, `inputs/Q077_PREREG.md` for its two declared gaps (assumed ANT kinetics, assumed effective diffusivity), and `inputs/Q077_RESULTS.md` for the failed `R_flux>=1.10` gate. `PREREG.md` + `PREREG.sha256` were frozen and checksum-verified before any byte was downloaded. No `Q077` node id exists here: `UNKNOWN`.

## Found

`DATA_SOURCES.json`, 10 candidates: 4 `verified`, 2 `metadata-only`, 4 `unreachable`/absent. Coverage of `Parameters`: **18 of 25** physical fields have a verified or metadata-only public dataset; `diffusion_coefficient` is named only by an absent one; `ant_k_adp`, `ant_k_atp`, `ant_fraction_pmf`, `membrane_potential` have none (three further misses are solver knobs, where none is expected). Geometry is the well-covered branch — exactly as `PREREG.md` predicted.

Usable anchors: **EMPIAR-11542** (cryo-ET tilt series, C. elegans whole mitochondria + released cristae, 159.01 GB, MRC, 5.4/3.576 O/px) for `crista_length`, `junction_area_total`, `junction_length`, `total_ims_volume`; **EMD-18991** for site densities; **MTBLS4307** (mitochondria-vs-cytoplasm fractionation, 38 samples) for the adenylate pool split; **BIOMD0000000090** and **MODEL1603150001** (SBML) for turnover numbers — both labelled as models, never as data.

## Sample

`samples/emd_18991.map.gz`, 2 688 407 B, sha256 `948fa9b9…c55e1a`, decompressing to 10 977 024 B. Form: 3-D scalar density volume, 2 744 000 voxels. Units: angstrom, arbitrary density scale. Coordinate frame: origin `(0,0,0)`, col/row/sec ordering = X/Y/Z fast-to-slow, 3.58 O/voxel, box 501.2 Oh, space group 1, C2 applied. **Found a defect**: the deposited MRC `mode` word says int32 while the EMDB XML and the bytes say float32 — the float32 reading reproduces the depositor's own min/max to 1e-8, the int32 reading spans the whole dtype range including the `-2147483648` sentinel.

## What fell

Four probes, reported not hidden: EMPIAR raw bytes (FTP 404; PDBj mirror TLS mismatch), SABIO-RK REST (serves a 404 shell), HMDB (Cloudflare 403), BRENDA export (click-through, fails C1). EMPIAR-11542 therefore stays `metadata-only` and no tilt frame was parsed. Two explicit negatives: no public ADP-diffusivity-in-IMS dataset, no IMS adeninate time series.

## Next

BRENDA and SABIO-RK need one human click-through each, not new science — they would close `as_kcat` and the ANT scale. The diffusivity gap is real and must stay `UNKNOWN`; `diffusion_coefficient` drives every transport time linearly. No model was re-run; no parameter was changed.
