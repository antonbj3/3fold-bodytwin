# Cell portability fix — COMPLEMENT_DISCRIMINATION, CORNEA_SCATTER, MITOSTRESS

Status: **PENDING_INDEPENDENT_REVIEW**
Date: 2026-10-03
Scope: path resolution only. No physics, no parameters, no numbers changed.
Workspace: `` (nothing written to `source_repository/`).

## 1. Paths found and replaced

Line numbers are from the pre-change files (`git show HEAD:<path>`).

| Cell | File | Line | Before (absolute / depth-bound) | After (relative to the file) |
|---|---|---|---|---|
| COMPLEMENT_DISCRIMINATION | cell.py | 16 | `Path('source_repository/data/complement_cascade/complement_cascade_results.json')` | `_acquired('data/.../complement_cascade_results.json')` → `HERE/<name>` if present, else `$BODYTWIN_REPO` (default = old path) |
| COMPLEMENT_DISCRIMINATION | cell.py | 17 | `Path('source_repository/scripts/msk/complement_cascade.py')` | `_acquired('scripts/msk/complement_cascade.py')`, same resolver |
| CORNEA_SCATTER | cell.py | 16 | `WORK = Path(__file__).resolve().parents[4]` (fixed directory depth) | `_eye_dir()`: `$BODYTWIN_EYE_DIR`, then `HERE`, then any ancestor `results/LANE_EYE_OPTICAL_TWIN`, then canonical lane dir; probes for `scatter_port_r1.py` |
| CORNEA_SCATTER | cell.py | 22 | `Path('source_repository/data/corneal_transparency/corneal_transparency_results.json')` | `HERE/<name>` if present, else `$BODYTWIN_REPO/<rel>` (default = old path) |
| MITOSTRESS | oxphos_scrambled_stoichiometry_voidfloor.py | 21 | `sys.path.insert(0, "source_repository/data/body_twin/agent_scratch_preserved")` | `sys.path.insert(0, str(_dep_dir()))`: `$BODYTWIN_SCRATCH_DIR`, then `HERE`, then canonical dir; probes for the reimpl module |
| MITOSTRESS | oxphos_scrambled_stoichiometry_voidfloor.py | 128 | `open("/tmp/mt_vfoxphos_n7c4/voidfloor_oxphos_results.json", "w")` | `open(OUT_DIR / "voidfloor_oxphos_results.json", "w")`, `OUT_DIR = $CELL_OUT_DIR` or the cell's own directory |
| MITOSTRESS | oxphos_ode_emergent_downstream_propagation.py | 26 | same `sys.path.insert` absolute | `_dep_dir()` as above |
| MITOSTRESS | oxphos_ode_emergent_downstream_propagation.py | 209 | `out_path = "source_repository/data/oxphos_ode_emergent_downstream_propagation_results.json"` (a write into the read-only repo) | `out_path = OUT_DIR / "oxphos_ode_emergent_downstream_propagation_results.json"` |
| MITOSTRESS | oxphos_corrected_nhatp_resimulation.py | 41 | same `sys.path.insert` absolute | `_dep_dir()` as above |
| MITOSTRESS | oxphos_corrected_nhatp_resimulation.py | 113 | `open("/tmp/oxphos_corrected_nhatp_final.json", "w")` | `open(OUT_DIR / "oxphos_corrected_nhatp_final.json", "w")` |

10 call sites in 5 files: 9 absolute literals plus 1 fixed-depth `parents[4]`.
Base in every case is `Path(__file__).resolve().parent`; the absolute literal survives
only as the last-resort default for *acquired inputs* that do not live inside the cell,
so in-place behaviour is unchanged.

## 2. Numeric comparison, before vs after (run in place)

Compared field by field over all numeric JSON leaves (recursive flatten, key sets
checked for identity as well).

| Cell / script | numeric fields compared | max abs diff | non-numeric diffs |
|---|---|---|---|
| COMPLEMENT_DISCRIMINATION / cell.py | 21 | 0.000e+00 | 0 |
| CORNEA_SCATTER / cell.py | 15 | 0.000e+00 | 0 |
| MITOSTRESS / oxphos_scrambled_stoichiometry_voidfloor.py | 33 | 0.000e+00 | 0 |
| MITOSTRESS / oxphos_ode_emergent_downstream_propagation.py | 30 | 0.000e+00 | 0 |
| MITOSTRESS / oxphos_corrected_nhatp_resimulation.py | 30 | 0.000e+00 | 0 |

Baseline harness note: the three MITOSTRESS "before" runs could not be executed as
written. Two wrote into fixed locations (`source_repository/data/...`, i.e.
the read-only repo; and `/tmp/oxphos_corrected_nhatp_final.json`), and one targeted
`/tmp/mt_vfoxphos_n7c4/`, a directory that does not exist, so it raised
`FileNotFoundError` after computing. For the baseline the unmodified sources were copied
to a scratch directory with **only the output path string** rewritten to scratch, nothing
else; that is the harness, not a change to the delivered files.

`COMPLEMENT_DISCRIMINATION/query.py` carries no absolute path (it already used
`Path(__file__).with_name('cell.py')`) and was not modified.

## 3. Temporary-directory test

All three cells were copied to `/tmp/coordinator-1000/portcheck_20261003/` and run with the
working directory set to `/` (not the cell directory, not the workspace):

| Script run from /tmp | exit | numbers vs in-place baseline |
|---|---|---|
| COMPLEMENT_DISCRIMINATION/cell.py | 0 | max abs diff 0.000e+00 (21 fields) |
| CORNEA_SCATTER/cell.py | 0 | max abs diff 0.000e+00 (15 fields) |
| oxphos_scrambled_stoichiometry_voidfloor.py | 0 | max abs diff 0.000e+00 (33 fields) |
| oxphos_ode_emergent_downstream_propagation.py | 0 | max abs diff 0.000e+00 (30 fields) |
| oxphos_corrected_nhatp_resimulation.py | 0 | max abs diff 0.000e+00 (30 fields) |

The MITOSTRESS outputs landed next to the copied cell in
`/tmp/coordinator-1000/portcheck_20261003/MITOSTRESS/`, not at any fixed location.

Negative control, same /tmp directory, pre-change sources from `git show HEAD:`:

- `CORNEA_SCATTER/cell.py` → exit 1, `ModuleNotFoundError: No module named 'scatter_port_r1'` (the `parents[4]` depth assumption).
- `oxphos_scrambled_stoichiometry_voidfloor.py` → exit 1, `FileNotFoundError: '/tmp/mt_vfoxphos_n7c4/voidfloor_oxphos_results.json'`.
- `COMPLEMENT_DISCRIMINATION/cell.py` → exit 0 even before the change. Its two absolute
  input paths happen to exist on this machine, so this cell's defect is portability to
  another machine or a copied input tree, not to another working directory. Stated as
  measured, not inflated.

`find source_repository/ -newermt <session start>` returned nothing: the
read-only repo was not written to at any point.

## 4. Same pattern in the remaining 42 cells

Searched with: absolute literals under `/home /mnt /opt /data /var /srv /tmp`, f-string
absolute paths, `os.path.join` with an absolute first element, any `sys.path.insert|append`,
fixed-depth `parents[N]`, absolute paths inside `*.json`, and `*.sh/*.yaml/*.yml/*.toml/*.cfg`.

**6 further cells carry runtime-breaking absolute paths in Python** (13 files):
`IMMUNITY`, `SOLBENCH`, `SURG_COLLAGEN`, `SURG_HEALING`, `SURG_HEMOSTASIS`, `SURG_INCISION`.
Typical forms: `OUT_DIR = "source_repository/data/msk_smoketest/..."` followed
by `os.makedirs(OUT_DIR, exist_ok=True)`; `REPO_ROOT = "source_repository/"`;
absolute `out_path`/`RESULTS_PATH` writes; absolute cross-cell input JSON
(`WOUND_JSON`, `COAG_JSON`, `TEWL_JSON`).

**1 further cell, metadata only:** `Q090/Q090_results.json` records
`"path": "/opt/agents/jobs/BT-HX-Q090/PREREG.md"`. Provenance, not executed — counted
separately. The absolute paths inside `COMPLEMENT_DISCRIMINATION/PORT*.json` and
`CORNEA_SCATTER/PORT*.json` are likewise provenance pointers and were left alone.

Also present, same class of fragility, not fixed here: fixed-depth
`Path(__file__).resolve().parents[2]` in `MITOSTRESS/mitochondrial_oxphos.py`,
`MITOSTRESS/autophagy_mitophagy.py`, `IMMUNITY/tcell_activation_exhaustion.py`.

**Count missed by the original grep: 7 cells** — 6 executable, 1 metadata-only.
No f-string absolute path, no `os.path.join` with an absolute first element, and no
shell/YAML/TOML/cfg hit anywhere in the 45 cells.
