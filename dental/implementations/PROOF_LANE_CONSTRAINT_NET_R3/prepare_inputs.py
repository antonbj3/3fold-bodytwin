"""Freeze selected source records, never source trees. Re-running verifies, never refreshes."""
from dental_release.paths import expand as _release_expand
import hashlib, json, shutil
from pathlib import Path
from datetime import datetime, timezone
P = Path(__file__).resolve().parent
D = Path(_release_expand('@DENTAL_INPUT_ROOT@/workspace'))
R = D / 'results'
R2 = R / _release_expand('CONSTRAINT_R2')
R1 = R / _release_expand('CONSTRAINT_NET')

def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def main():
    if (P / 'INPUT_LOCK.json').exists():
        x = json.loads((P / 'INPUT_LOCK.json').read_text())
        for row in x['files'].values():
            assert sha(P / row['snapshot']) == row['sha256']
        print('Frozen inputs verified')
        return
    (P / 'sources').mkdir(exist_ok=True)
    files = {}
    originals = {}

    def add(logical, actual=None, expected=None):
        logical = str(logical)
        actual = Path(actual or logical)
        h = sha(actual)
        if expected:
            assert h == expected, (logical, h, expected)
        if logical in files:
            assert files[logical]['sha256'] == h, logical
            return
        assert actual.stat().st_size < 50000000
        dst = 'sources/' + h + actual.suffix
        if not (P / dst).exists():
            shutil.copyfile(actual, P / dst)
        files[logical] = {'sha256': h, 'snapshot': dst, 'bytes': actual.stat().st_size}
        if Path(logical).is_file():
            originals[logical] = sha(logical)

    def walk(x):
        if isinstance(x, dict):
            if all((k in x for k in ['source_file', 'snapshot', 'source_sha256'])):
                p = R2 / x['snapshot']
                if not p.exists():
                    p = R1 / x['snapshot']
                add(x['source_file'], p, x['source_sha256'])
            for v in x.values():
                walk(v)
        elif isinstance(x, list):
            for v in x:
                walk(v)
    net = json.loads((R2 / 'CONSTRAINT_NET_DENTAL_R2.json').read_text())
    walk(net)
    for (f, row) in json.loads((R2 / 'INPUT_LOCK.json').read_text())['files'].items():
        add(f, R2 / row['snapshot'], row['sha256'])
    props = D.parent / 'bodytwin/results/ASSEMBLY_DENTAL_GAP_PROPOSALS/PROPOSALS_V1.json'
    add(props)
    for p in json.loads(props.read_text())['proposals']:
        add(p['source_file'], expected=p['source_sha256'])
    names = ['PROOF_LANE_CONSTRAINT_NET_R2/CONSTRAINT_NET_DENTAL_R2.json', 'PROOF_LANE_CONSTRAINT_NET_R2/VALIDATION.json', 'PROOF_LANE_CONSTRAINT_NET_R2/results.json', 'LANE_XREVIEW_BATCH20/REVIEW_PROOF_LANE_CONSTRAINT_NET_R2.json', 'LANE_XREVIEW_BATCH20/patches/PROOF_LANE_CONSTRAINT_NET_R2/validate.py', 'LANE_X71_MEASUREMENT_PRIORITY/inputs/MEASUREMENTS.json', 'LANE_X71_MEASUREMENT_PRIORITY/inputs/DECISIONS.json', 'LANE_X71_MEASUREMENT_PRIORITY/MEASUREMENT_PRIORITY.md', 'LANE_NIGHT_I_DECISION_CERT/CALIBRATION_ACQUISITION_SPEC.json', 'LANE_X55_METROLOGY/results.json', 'LANE_X56_MICROMOTION/results.json', 'LANE_X64_MANUFACTURING/results.json', 'LANE_X68_INSERTION_TORQUE/results.json', 'PROOF_LANE_GENCAD_V5B/results.json', 'LANE_X58_CANAL_WALL_SPREAD/raw/PER_SITE_LOCATED.json', 'LANE_XREVIEW_BATCH15/REVIEW_LANE_X58_CANAL_WALL_SPREAD.json', 'LANE_XREVIEW_BATCH3/REVIEW_LANE_X31_CLINICAL_ANSWERS.json', 'LANE_XREVIEW_BATCH3/REVIEW_LANE_X26_DECIDABILITY.json', 'LANE_XREVIEW_BATCH8/REVIEW_LANE_X44_MISSING_CHAINS.json', 'LANE_XREVIEW_BATCH16/REVIEW_LANE_X65_ISO14801_FATIGUE.json', 'LANE_XREVIEW_BATCH17/REVIEW_LANE_X67_DRILL_HEAT.json']
    for n in names:
        add(R / n)
    for f in R2.rglob('*'):
        if f.is_file() and '__pycache__' not in str(f):
            originals[str(f)] = sha(f)
    (P / 'INPUT_LOCK.json').write_text(json.dumps({'created_utc': datetime.now(timezone.utc).isoformat(), 'files': files}, ensure_ascii=False, indent=2) + '\n')
    (P / 'READ_ONLY_MANIFEST.json').write_text(json.dumps(originals, indent=2) + '\n')
    print('Frozen', len(files), 'logical files;', sum((p.stat().st_size for p in (P / 'sources').iterdir())), 'snapshot bytes')
if __name__ == '__main__':
    main()
