"""Freeze only explicitly selected small JSON/text records, not dataset trees."""
from dental_release.paths import expand as _release_expand
import datetime
from common import *
JOBS = [_release_expand('X69'), _release_expand('X70'), _release_expand('X71'), _release_expand('X72'), _release_expand('X73'), _release_expand('X74'), _release_expand('X75'), 'LANE_X76_DATASET_COMPLETION', _release_expand('FULL_CROWN'), _release_expand('GENCAD_V6'), _release_expand('FALT_TANDLAST')]

def main():
    if (HERE / 'INPUT_LOCK.json').exists():
        for row in read(HERE / 'INPUT_LOCK.json')['files'].values():
            assert sha(HERE / row['snapshot']) == row['sha256']
        print('Using immutable INPUT_LOCK; no source refresh')
        return
    selected = [R1 / 'CONSTRAINT_NET_DENTAL.json', R1 / 'CONSTRAINT_STRESS_MAP_DENTAL.json', R1 / 'results.json', DENT / 'results/LANE_XREVIEW_BATCH18/CORRECTIONS.jsonl']
    for job in JOBS + [_release_expand('CONSTRAINT_NET')]:
        for batch in [18, 19]:
            selected.extend((DENT / f'results/LANE_XREVIEW_BATCH{batch}').glob('REVIEW_' + job + '.json'))
        p = DENT / 'results' / job / 'results.json'
        if p.exists():
            selected.append(p)
    for job in ['LANE_X24_CBCT_HU_CALIBRATION', 'LANE_X33_LASER_PULP', 'LANE_X44_MISSING_CHAINS', 'LANE_X26_DECIDABILITY', 'LANE_X31_CLINICAL_ANSWERS']:
        selected.append(DENT / 'results' / job / 'results.json')
    selected += [DENT / 'results/LANE_X2_OCCLUSION_B2B/raw/EXTERNAL_REFERENTS.json', DENT / 'results/LANE_X71_MEASUREMENT_PRIORITY/inputs/DECISIONS.json', DENT / 'results/LANE_X71_MEASUREMENT_PRIORITY/inputs/MEASUREMENTS.json', DENT / 'results/LANE_X72_FALT_FIXES/raw/X18.json', DENT / 'results/LANE_X72_FALT_FIXES/raw/X54.json', DENT / 'results/LANE_XREVIEW_BATCH18/evidence/independent_checks.json']
    net = read(R1 / 'CONSTRAINT_NET_DENTAL.json')
    for e in net['edges']:
        if e['id'] in ['D-E-SINTER-FIT', 'D-E-REPLICA-CT', 'D-E-PRELOAD-HISTORY', 'D-E-CANAL-RELEASE']:
            selected += [Path(x['source_file']) for x in e['evidence']]
    for p in list(selected):
        if p.name.startswith('REVIEW_'):
            d = read(p)
            q = Path(d['result_file'])
            if q.exists():
                selected.append(q)
    files = {}
    (HERE / 'sources').mkdir(exist_ok=True)
    for p in sorted(set(selected)):
        size = p.stat().st_size
        if size > 20000000:
            raise ValueError('Unexpectedly large record ' + str(p))
        h = sha(p)
        dst = HERE / 'sources' / (h + p.suffix)
        dst.write_bytes(p.read_bytes())
        files[str(p)] = {'sha256': h, 'bytes': size, 'snapshot': str(dst.relative_to(HERE))}
    dump('INPUT_LOCK.json', {'frozen_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(), 'files': files, 'jobs': JOBS, 'review_cutoff': 'Only completed review JSON present at freeze; SAMPLE is not a review.'})
    print('Frozen', len(files), 'records;', sum((x['bytes'] for x in files.values())), 'bytes')
if __name__ == '__main__':
    main()
