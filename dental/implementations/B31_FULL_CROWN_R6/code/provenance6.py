"""Read-only source integrity audit and replay dependency lock.

This is a post-experiment dependency inventory, not a retrospective preregistration.
Original prediction manifests and their freeze timestamps remain authoritative.
"""
from dental_release.paths import expand as _release_expand
from local import *
import score6, construct_f, construct_k, exact_wall6
from fractions import Fraction
sys.path.insert(0, str(BASE / _release_expand('PROOF_LANE')))
from gencad_bench.checks import exact
import importlib.metadata as metadata

def run():
    if (R6 / 'DEPENDENCIES.json').exists():
        raise RuntimeError('Dependency lock already exists; do not replace it')
    checks = []
    paths = set()

    def verify(path, expected, kind):
        path = Path(path).resolve()
        actual = sha(path)
        checks.append(dict(path=str(path), expected_sha256=expected, actual_sha256=actual, match=actual == expected, kind=kind))
        paths.add(path)
    for rec in read(BASE / 'PROOF_LANE_FULL_CROWN_R4/FROZEN_INPUTS_B.json')['records']:
        for kind in ('public', 'private'):
            verify(rec[kind + '_path'], rec[kind + '_sha256'], kind)
    for rec in read(BASE / 'PROOF_LANE_FULL_CROWN_R4/FROZEN_PREDICTIONS_C.json')['rows']:
        if rec['method'] == 'exact_margin':
            verify(rec['mesh_path'], rec['mesh_sha256'], 'unchanged_virtual_preparation')
    for rec in read(BASE / 'PROOF_LANE_FULL_CROWN_R5/FROZEN_PREDICTIONS_B.json')['rows']:
        if rec['method'] == 'distance_local_thickening':
            verify(rec['mesh_path'], rec['mesh_sha256'], 'R5_control')
    for p in sorted(R6.glob('*.sha256')):
        parts = p.read_text().split()
        expected = parts[0]
        name = parts[1] if len(parts) > 1 else p.with_suffix('.json').name
        verify(R6 / name, expected, 'original_frozen_manifest')
    for p in sorted(R6.glob('FROZEN_PREDICTIONS_*.json')):
        for rec in read(p).get('rows', []):
            if 'mesh_path' in rec:
                verify(rec['mesh_path'], rec['mesh_sha256'], 'R6_frozen_candidate')
    audit = dict(audited_utc=now(), all_pass=all((r['match'] for r in checks)), checks=checks, source_cohort_teeth=18, source_arches=6, scope='Encoded inputs and original frozen candidates; no physical measurement')
    save(R6 / 'raw/INPUT_INTEGRITY.json', audit)
    assert audit['all_pass'], 'A frozen source changed'
    for m in list(sys.modules.values()):
        p = getattr(m, '__file__', None)
        if p and str(BASE) in p and Path(p).is_file():
            paths.add(Path(p).resolve())
    paths.update((p.resolve() for p in (R6 / 'code').glob('*.py')))
    paths.update((p.resolve() for p in (R6 / 'code').glob('*.cpp')))
    paths.update((p.resolve() for p in (R6 / 'code').glob('*.sh')))
    for p in [R6 / 'run_all.sh', R6 / 'RESULTS_D.json', D6 / 'D_exterior.npz', BASE / 'PROOF_LANE_FULL_CROWN_R4/FROZEN_TEST_SET.json', BASE / 'PROOF_LANE_FULL_CROWN_R4/FROZEN_INPUTS_B.json', BASE / 'PROOF_LANE_FULL_CROWN_R4/FROZEN_PREDICTIONS_C.json', BASE / 'PROOF_LANE_FULL_CROWN_R5/FROZEN_PREDICTIONS_B.json', BASE / 'PROOF_LANE_FULL_CROWN_R5/RESULTS_B.json', BASE / 'LANE_X13_CEMENT_GAP/FACIT.csv', BASE / 'LANE_X49_DESIGN_GATE/README_DEMO.md']:
        paths.add(p.resolve())
    for name in ('intersections', 'repair', 'repair_local', 'repair_protected'):
        paths.add(D6 / name)
    paths.update(D6.glob('*.deb'))
    paths.update((D / 'deps').glob('igl*.so'))
    versions = {}
    for package in ('numpy', 'scipy', 'trimesh', 'scikit-image', 'shapely', 'rtree', 'threadpoolctl'):
        try:
            versions[package] = metadata.version(package)
        except metadata.PackageNotFoundError:
            versions[package] = 'UNKNOWN'
    save(R6 / 'DEPENDENCIES.json', dict(locked_utc=now(), scope='Post-experiment local replay lock; not a portable installer or new prediction freeze', files=[dict(path=str(p), sha256=sha(p), bytes=p.stat().st_size) for p in sorted(paths)], runtime_versions=versions, python=sys.version, source_integrity_sha256=sha(R6 / 'raw/INPUT_INTEGRITY.json')))
    print('All', len(checks), 'frozen checks passed;', len(paths), 'dependencies locked')
if __name__ == '__main__':
    run()
