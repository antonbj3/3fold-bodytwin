from dental_release.paths import expand as _release_expand
from common import *
from experiment import metrics
import argparse, csv

def init():
    m = read(X7 / 'raw/DATA_MANIFEST.json')
    candidates = read('raw/R1_ERROR_CASES.json') + read('raw/R3_ERROR_CASES.json')
    ids = sorted({r['case_id'] for r in candidates})[:11]
    if _release_expand('@DENTAL_CASE_ID@') in m['test'] and _release_expand('@DENTAL_CASE_ID@') not in ids:
        ids.append(_release_expand('@DENTAL_CASE_ID@'))
    freeze('FROZEN_MEASUREMENT_PANEL.json', dict(frozen_utc=now(), case_ids=ids, phase='No new clinical measurement performed. Existing test-patient observations already exposed.', prediction_files={k: dict(path=str(ROOT / ('raw/' + k + '_PREDICTIONS.json')), sha256=sha(ROOT / ('raw/' + k + '_PREDICTIONS.json'))) for k in ['R1', 'R3']}, representation='Separate native presence, visible prosthetic site replacement, restoration/crown state, arch spacing and wear. Anonymous rater + repeat; unknown allowed.', prereg_sha256=sha(ROOT / 'PREREG_R1_RESTORATION_GEOMETRY.json'), material_port_sha256=sha(ROOT / 'raw/MATERIAL_APPLICABILITY_PORT_R5.json')))
    p = ROOT / 'raw/MEASUREMENT_TEMPLATE.csv'
    if not p.exists():
        with p.open('w', newline='') as fh:
            w = csv.DictWriter(fh, fieldnames=['case_id', 'fdi', 'rater_id', 'repeat', 'restoration_present', 'native_surface_present', 'prosthetic_site_present', 'crown_present', 'wear_grade', 'pair_neighbor_fdi', 'spacing_gap_mm', 'observation_method', 'acquisition_time', 'uncertainty', 'source_locator'])
            w.writeheader()
            for c in ids:
                for q in range(1, 5):
                    for pos in range(1, 8):
                        w.writerow(dict(case_id=c, fdi=q * 10 + pos))

def compare(path):
    panel = read('FROZEN_MEASUREMENT_PANEL.json')
    for r in panel['prediction_files'].values():
        assert sha(r['path']) == r['sha256']
    records = list(csv.DictReader(open(path)))
    filled = [r for r in records if r['restoration_present'] != '']
    if not filled:
        write('raw/MATCHED_MEASUREMENT_RESULT.json', dict(status='NOT_RUN_NO_FILLED_OBSERVATIONS', provided_rows=len(records)))
        print('NOT_RUN: no filled independent observations', flush=True)
        return
    for r in filled:
        assert r['case_id'] in panel['case_ids'] and r['source_locator'] and r['rater_id'] and r['observation_method']
        assert r['restoration_present'] in ['0', '1', 'UNKNOWN']
    rows = [dict(case_id=r['case_id'], fdi=int(r['fdi']), y=int(r['restoration_present'])) for r in filled if r['restoration_present'] in ['0', '1']]
    assert len({(r['case_id'], r['fdi']) for r in rows}) == len(rows), 'For prediction accuracy select one rater/repeat per tooth; analyze repeated-rater agreement separately'
    pred = {(r['case_id'], r['fdi']): r for r in read(panel['prediction_files']['R1']['path'])}
    out = metrics(rows, [pred[r['case_id'], r['fdi']]['predicted']['shape'] for r in rows])
    out.update(status='MEASURED_COMPARISON_UNREVIEWED', measurement_path=str(Path(path).resolve()), measurement_sha256=sha(path), review_state='PENDING_INDEPENDENT_REVIEW', no_refit=True)
    write('raw/MATCHED_MEASUREMENT_RESULT.json', out)
    print('Matched measurement result saved', flush=True)
if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--observations')
    args = ap.parse_args()
    init()
    compare(args.observations or ROOT / 'raw/MEASUREMENT_TEMPLATE.csv')
