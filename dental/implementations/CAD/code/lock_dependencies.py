from dental_release.paths import expand as _release_expand
from cadlib import *
paths = set()
for (lane, folder) in [(_release_expand('GENCAD_V2'), 'gencad_bench_v2'), (_release_expand('GENCAD_V6'), 'code'), ('PROOF_LANE_FULL_CROWN_R4', 'code'), (_release_expand('X85'), 'code'), (_release_expand('X82'), 'code'), (_release_expand('X89'), 'code'), (_release_expand('X53'), 'code'), ('LANE_X52_PREP_ERROR_LIBRARY', 'code'), ('LANE_NEXT_D_INSERTION_PROOF', 'code')]:
    paths.update((BASE / lane / folder).glob('*.py'))
paths.update((BASE / 'PROOF_LANE_GENCAD_V2/gencad_bench_v2/vendor').glob('*.py'))
paths.update((BASE / 'PROOF_LANE_GENCAD_V2/gencad_bench_v2/vendor/v1/gencad_bench/checks').glob('*.py'))
paths.update([BASE / 'LANE_X13_CEMENT_GAP/cement_port.py', BASE / 'LANE_X53_MILLING_TOOLS/NOMINAL_TIP_REFERENCE_CARD.json', BASE / 'PROOF_LANE_GENCAD_V4/code/quality.py'])
for r in read(ROOT / 'INPUT_LOCK_R3.json')['rows']:
    paths.add(Path(r['task_path']))
    if r['status'] == 'AVAILABLE':
        for k in ['public_geometry', 'die_path', 'reference_path']:
            paths.add(Path(r[k]))
        paths.add(DATA.parent / 'PROOF_LANE_GENCAD_V2/axial_R3' / (r['key'] + '_axial_flare.npz'))
for r in read(ROOT / 'raw/R1_PREDICTIONS.json'):
    if 'source' in r:
        paths.add(Path(r['source']))
paths.update([ROOT / 'INPUT_LOCK_R3.json', ROOT / 'raw/R1_PREDICTIONS.json', ROOT / 'raw/R2_PREDICTIONS.json', ROOT / 'raw/R5_PREDICTIONS.json', DATA / 'R1_references.npz'])
for row in read(ROOT / 'raw/R5_PREDICTIONS.json'):
    if row['status'] == 'GENERATED':
        paths.add(Path(row['mesh_path']))
files = [dict(path=p, bytes=p.stat().st_size, sha256=sha(p)) for p in sorted(paths)]
freeze(ROOT / 'DEPENDENCIES.json', dict(files=files, scope='Read-only execution dependencies and immutable geometry. Native archive members are hash-validated by predecessor pair loader; archives not copied or rehashed wholesale. Standard Python package versions recorded separately.'))
import platform, scipy, trimesh, networkx, matplotlib
freeze(ROOT / 'RUNTIME_VERSIONS.json', dict(python=platform.python_version(), executable=sys.executable, numpy=np.__version__, scipy=scipy.__version__, trimesh=trimesh.__version__, networkx=networkx.__version__, matplotlib=matplotlib.__version__, platform=platform.platform()))
