"""Curated files only. Never copies a source checkout."""
from pathlib import Path
import hashlib, json, shutil
R = Path(__file__).resolve().parents[1]
S = R.parent
v2 = S / 'PROOF_LANE_GENCAD_V2/gencad_bench_v2'
rows = []
files = {**{str(v2 / (n + '.py')): f'code/legacy/{n}.py' for n in ['geometry', 'checks', 'generators']}, **{str(v2 / 'vendor' / n): f'code/legacy/vendor/{n}' for n in ['obstacle.py', 'continuous.py', 'crown_fit_geometry.py']}, str(S / 'LANE_X1B_CROWN_LOOP/vendor/crown_design_geometry.py'): 'vendor/crown_design_geometry.py', str(S / 'LANE_X1B_CROWN_LOOP/vendor/crown_fit_geometry.py'): 'vendor/crown_fit_geometry.py', str(S / 'LANE_X18B_PREP_SURFACE/code/surface.py'): 'vendor/x18b_surface.py', str(S / 'LANE_X18B_PREP_SURFACE/code/implicit_prep.py'): 'vendor/x18b_implicit_prep.py'}
for (a, b) in files.items():
    p = R / b
    p.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(a, p)
    rows.append(dict(source=a, local=b, sha256=hashlib.sha256(p.read_bytes()).hexdigest()))
(R / 'code/legacy/__init__.py').touch()
(R / 'SOURCE_REUSE.json').write_text(json.dumps(dict(files=rows, review=str(S / 'LANE_XREVIEW_GENCAD_V2/XREVIEW.md'), existing_search='Targeted dental/results, cells, swarm48 dental RESULTS and bodytwin LANE RESULTS; no equivalent v3 found', old_cells_index='notes/OLD_DENTAL_CELL_INDEX.md and jsonl; no old measured preparation quantity imported'), indent=2) + '\n')
