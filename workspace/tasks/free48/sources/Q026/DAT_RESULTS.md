BT-DAT-Q026

# Which public measured datasets constrain the Fgf8a gradient model in `inputs/`?

Job type: **data sourcing**. The model is frozen in `inputs/Q026_model.py`; it was not re-run and
not re-fitted. Preregistration: `PREREG.md` (sha256 `PREREG.sha256`). Machine-readable:
`DATA_SOURCES.json`, `results.json`, `sample_load.json`.

## Headline finding: the anchor paper deposited nothing

The model's only **VERIFIED** numeric parameter, `D = 55 um^2 s^-1`, comes from Harish et al. 2023
(`10.1242/dev.201559`, PMC10565248, CC BY). Its Data Availability statement reads, verbatim:

> "All relevant data can be found within the article and its supplementary information."

Read from `https://www.ebi.ac.uk/europepmc/webservices/rest/PMC10565248/fullTextXML`. There is no
repository deposit. `C0 = 8 nM` and the `40-70 um` / `70-140 um` windows exist only as figure values
in a PDF. So the strongest measured number in the model has **no public machine-readable dataset
behind it**, and any independent bound is doing real work here.

## The six sources, and what each one actually constrains

Full detail in `DATA_SOURCES.json`. License and size were taken from each publisher's own metadata
endpoint at build time, not from a search snippet.

| id | dataset | license | size | subjects | constrains |
|---|---|---|---|---|---|
| SRC1 | PDB **2FDB**, FGF8b·FGFR2c, 2.28 A X-ray | CC0 (PDB public domain) | **595,036 B** mmCIF | n/a, 4 protein + 4 solvent chains | `K_D_nM` structurally; ligand identity |
| SRC2 | UniProtKB **P55075** FGF8_HUMAN | CC BY 4.0 | ~10 kB JSON | n/a | source boundary condition (annotated "Secreted") |
| SRC3 | UniProtKB **P11362** FGFR1_HUMAN | CC BY 4.0 | ~10 kB JSON | n/a | `theta` readout population; 25 um `sigma` must not be confused with the ~30-40 nm ectodomain |
| SRC4 | Dryad **10.5061/dryad.qg8dt** Bajanca et al. in vivo Dystrophin FRAP, zebrafish | **CC0-1.0** | 540,209,369 B (1 zip) | zebrafish embryos, n = UNKNOWN | `D`, and the **excluded 7% slow/bound pool** |
| SRC5 | Cell Tracking Challenge **Fluo-N3DL-DRO**, Drosophila embryo | **not open** - explicit permission required | 5.8 GB zip | 1 embryo (Keller, Janelia) | `sigma_um`, `half_width_um`, time anchors (30 s steps) |
| SRC6 | EMBL-EBI **Expression Atlas** | CC BY 4.0 | 2,610,109 B index | 4562 experiments; per-experiment n = UNKNOWN | which cells are "receptor-bearing" for `theta` |

Two honest negatives worth more than the positives: **SRC5 cannot be used** (its conditions of use
require explicit permission for any non-CTC scientific use, and the file is 5.8 GB), and **SRC4 is
the only measured free-vs-bound diffusion dataset in the set** but is 515 MB, 10x over the sample
limit. SRC1 is FGF8**b**, while the model anchors on Fgf8**a** - the shared receptor-binding core
transfers, the splice-variant segment does not, and PDB stores no `KD`, so SRC1 bounds `K_D` as a
prior and not as a number.

## Sample actually downloaded and loaded

`fetch_sample.py` (enforces the 50 MB cap, aborts and deletes rather than truncating) pulled
`https://files.rcsb.org/download/2FDB.cif` into `samples/`:

- **595,036 B**, within the 52,428,800 B limit; `Content-Type: chemical/x-cif`; HTTP 200
- **sha256** `d2d91ae0114cfbebaccdd3e3c3d4974e81ff980d8348f20f874e7dec3e506b56`, recomputed from
  disk after writing and identical to the streaming digest
- **form**: mmCIF text block, first line `data_2FDB`, 8,908 lines
- **units**: angstrom - cell lengths 170.846 / 46.908 / 109.616 A, beta = 91.66 deg; 173 K,
  d_min 2.28 A
- **coordinate frame**: orthonormal Cartesian in the deposited unit cell, space group `C 1 2 1`
  (International Tables no. 5), origin at the crystallographic cell corner
- **shape**: 5,345 atoms read from the `_atom_site` loop across 8 chains
  (A-D: 1161/1506/1128/1425; E-H: 34/39/21/31 solvent atoms); x range -14.444..89.793 A,
  y -29.759..34.349 A, z -13.533..59.069 A; composition LYS 411, VAL 410, LEU 368, TYR 332

Full stdout and `sample_load.json` back this. Values above are read out of the file on disk, not
from the RCSB landing page.

## Derived from the frozen model (no new measurement, no re-fit)

From `inputs/Q026_results.json` only: `k_clear = 0.0015 s^-1` implies a residence time
**tau = 666.7 s = 11.1 min**, and `L = sqrt(2 D tau)` is **270.8 um** for `D = 55` but **73.0 um**
for `D = 4 um^2 s^-1`. That is the whole argument in one line: the model's assumed clearance sets an
11-minute clock, and whether the reported `0.4 C0` window lands at 70-140 um depends on which of the
two diffusion components dominates. The frozen run already stands at `x_0.4 = 203.81 um`, outside
that window - the pre-registered falsifier, still standing, still unresolved.

## What the hypotheses got wrong

- **H1 `k_clear` - not satisfied.** No public measured in-tissue residence time for a secreted FGF
  ligand was reachable. `k_clear = 0.0015 s^-1` is an **unmeasured assumption**, and it is the
  parameter the frozen sensitivity says moves `x_0.8` from 86.0 to 56.2 um over +/-50%
  (`inputs/Q026_results.json:sensitivity`). This is the model's weakest point.
- **H2 `D` - partial.** SRC4 is the right *kind* of measurement (measured, in vivo, free vs bound
  pool) but a sarcolemmal protein, so it bounds the scale and confirms no number.
- **H3 `K_D` - not satisfied as a number.** No freely-licensed, machine-readable FGF8-FGFR affinity
  table was found; BindingDB's REST service returned 404 on every endpoint tried. `K_D = 2 nM`
  stays an assumption, and `x_KD = 459.52 um` is therefore a distance no measured cell class can be
  shown to sample.
- **H4 `sigma` - qualitative only.** `sigma_um = 25 um` and `half_width_um = 500 um` remain
  assumptions; SRC5 would be the way to test them but is unusable.
- **H5 time anchors - satisfied for two sources** (SRC4 FRAP series, SRC5 30 s light-sheet steps),
  metadata verified, neither downloaded.

## Two identifier errors caught before they became findings

ENSG00000160789, recalled as human FGFR1 for a Human Protein Atlas route, is **LMNA**; P34443,
recalled as human FGF8, is **C. elegans RHEB-1**. Both were caught by API check, and the HPA source
was dropped rather than carried with a wrong accession. Human FGF8 is **P55075**, human FGFR1 is
**P11362** - both re-verified at build time. Recorded in `results.json:what_fell` because the honest
version of this report is one that says what it nearly got wrong.

## Next step

Get SRC4 (Dryad, CC0, ~515 MB) outside the 50 MB budget and fit a two-pool free/bound model to its
FRAP recovery curves. That is the single route that would turn `k_clear` and the excluded 7% pool
from assumptions into fitted parameters, and it is simultaneously the test of the standing
`x_0.4 > 140 um` falsifier. Second priority: a measured FGF8-FGFR affinity table to pin `K_D`, which
currently swings `x_KD` from `null` to 351.80 um across the model's own +/-50% sensitivity band.

## Reproduce

```
sha256sum -c PREREG.sha256
python3 fetch_sample.py          # downloads samples/2FDB.cif, prints form/units/frame
python3 build_data_sources.py    # re-probes every repository, writes DATA_SOURCES.json + results.json
```
