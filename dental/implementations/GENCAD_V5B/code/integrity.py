from common import *

def inputs():
    paths = set()
    for r in cohort():
        paths.add(Path(r['prep_path']))
        paths.add(V4 / 'payload/whole_private' / r['key'] / 'reference.npz')
        if r['status'] == 'SCORED':
            paths.add(Path(r['mesh_path']))
    for f in [V5 / 'raw/FUNCTIONAL_ROWS.json', V5 / 'raw/R6_EXTERNAL_SUPPORT.json', V5 / 'raw/PHYSICAL_ROWS.json.gz', V5 / 'LEADERBOARD.md', V5 / 'HANDOFF.md', V5 / 'results.json', V4 / 'FROZEN_LITERATURE_PREDICTIONS.json', V4 / 'payload/literature/CROWN_CURATED.json', V4 / 'payload/literature/CEMENT.csv', D / 'LANE_X1B_CROWN_LOOP/code/calibrate.py', D / 'LANE_X1B_CROWN_LOOP/raw/CALIBRATION_R1.json', D / 'LANE_X1B_CROWN_LOOP/raw/LITERATURE_GROUPS.json', D / 'LANE_X1B_CROWN_LOOP/raw/MATCHED_R2.json', D / 'LANE_X55_METROLOGY/code/metrology.py', D / 'LANE_X13_CEMENT_GAP/cement_port.py', D / 'LANE_NEXT_E_CEMENT_SQUEEZE/squeeze.py', D / 'LANE_X23_CROWN_WEIBULL_FLOOR/HANDOFF.md', D / 'LANE_X23_CROWN_WEIBULL_FLOOR/code/controls.py']:
        if f.exists():
            paths.add(f)
    for f in (D / 'LANE_X1B_CROWN_LOOP/inputs/literature').glob('*'):
        if f.is_file():
            paths.add(f)
    return sorted(paths)

def lock():
    p = P / 'INPUT_LOCK.json'
    if not p.exists():
        freeze(p, dict(scope='Selected read-only sources, pinned after exploratory pilots and before final full demo; no claim of preregistered source blindness', files=[dict(path=str(f), sha256=sha(f), bytes=f.stat().st_size) for f in inputs()]))
    return verify()

def verify():
    frozen = P / 'FROZEN_PREDICTIONS_R7.json'
    if frozen.exists():
        if sha(frozen) != frozen.with_suffix('.sha256').read_text().strip():
            raise ValueError('Frozen prediction bytes changed')
        for (name, digest) in read(frozen)['implementation_hashes'].items():
            if sha(P / 'code' / name) != digest:
                raise ValueError('Frozen numerical code changed; create a new version before generating ' + name)
    lock = read(P / 'INPUT_LOCK.json')
    bad = [r['path'] for r in lock['files'] if not Path(r['path']).exists() or sha(r['path']) != r['sha256']]
    if bad:
        raise ValueError('SOURCE DRIFT: ' + str(bad))
    for f in P.glob('PREREG_R*.json'):
        if sha(f) != f.with_suffix('.sha256').read_text().strip():
            raise ValueError('PREREG DRIFT ' + str(f))
    supp = read(P / 'CODE_SOURCE_LOCK.json') if (P / 'CODE_SOURCE_LOCK.json').exists() else {'files': []}
    for r in supp['files']:
        if sha(r['path']) != r['sha256']:
            raise ValueError('Source code drift ' + r['path'])
    return dict(pass_gate=True, supplementary_source_files=len(supp['files']), source_files=len(lock['files']), protocols=len(list(P.glob('PREREG_R*.json'))), source_bytes=sum((r['bytes'] for r in lock['files'])))

def manifest():
    paths = sorted(DATA.glob('R[1-7]/*/*'))
    rows = [dict(path=str(f), bytes=f.stat().st_size, sha256=sha(f)) for f in paths if f.is_file()]
    total = sum((r['bytes'] for r in rows))
    if total > 3000000000:
        raise ValueError('3GB lane intermediate budget exceeded')
    dump(P / 'ARTIFACT_MANIFEST.json', dict(files=rows, total_bytes=total, all_arrays_over_50MB_location=str(DATA), budget_bytes=3000000000))
    return dict(file_count=len(rows), total_bytes=total, manifest_sha256=sha(P / 'ARTIFACT_MANIFEST.json'))
if __name__ == '__main__':
    print(lock())
