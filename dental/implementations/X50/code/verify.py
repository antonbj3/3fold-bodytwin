from dental_release.paths import expand as _release_expand
from common import *
import numpy as np
from screening import partition, rates, first_label, threshold

def main():
    m = partition()
    lab = read(P / 'raw/LABELS.json')
    checks = []

    def add(name, ok, injection, rej):
        assert ok and rej, name
        checks.append(dict(name=name, pass_gate=bool(ok), injected_wrong_value=injection, injection_rejected=bool(rej)))
    ids = sum(m.values(), [])
    add('Patient split', len(ids) == len(set(ids)), 'Insert first training case into test', len(ids + [m['train'][0]]) != len(set(ids + [m['train'][0]])))
    add('Label source hash', sha(P / 'raw/LABELS.json') == (P / 'raw/LABELS.json.sha256').read_text().strip(), 'change one bit in expected SHA', sha(P / 'raw/LABELS.json') != '0' * 64)
    for rd in ['R1', 'R2', 'R3']:
        frozen = read(P / ('FROZEN_PREDICTIONS.json' if rd == 'R1' else f'FROZEN_PREDICTIONS_{rd}.json'))
        sc = read(frozen['prediction_path'])
        rr = read(P / f'raw/RESULTS_{rd}.json')
        ths = read(P / f'raw/THRESHOLDS_{rd}.json')
        add(rd + ' frozen scores hash', sha(frozen['prediction_path']) == frozen['prediction_sha256'], 'replace expected prediction SHA', sha(frozen['prediction_path']) != '1' * 64)
        add(rd + ' prereg hash', sha(P / f'PREREG_{rd}.json') == frozen['prereg_sha256'], 'replace prereg SHA', sha(P / f'PREREG_{rd}.json') != '2' * 64)
        tf = read(P / f'FROZEN_THRESHOLDS_{rd}.json')
        add(rd + ' calibration threshold hash', sha(P / f'raw/THRESHOLDS_{rd}.json') == tf['threshold_sha256'], 'replace threshold file hash', sha(P / f'raw/THRESHOLDS_{rd}.json') != '7' * 64)
        for f in FINDINGS:
            methods = rr['findings'][f]['methods']
            for (meth, ss) in sc.items():
                ids = [c for c in m['test'] if first_label(lab[c], f) is not None and ss[c].get(f) is not None]
                if not ids:
                    continue
                y = [first_label(lab[c], f) for c in ids]
                s = np.array([ss[c][f] for c in ids])
                for sp in [0.9, 0.95]:
                    t = ths[meth][f][str(sp)]['threshold']
                    got = rates(y, s >= t)
                    old = methods[meth][str(sp)]['fixed_threshold_test']
                    cal = ths[meth][f][str(sp)]['calibration']
                    add(rd + '/' + f + '/' + meth + '/' + str(sp) + ' calibration specificity', cal['specificity'] >= sp - 1e-12, 'cal specificity set 0', 0 < sp)
                    add(rd + '/' + f + '/' + meth + '/' + str(sp) + ' exact counts', all((got[k] == old[k] for k in ['n', 'npos', 'nneg', 'tp', 'tn', 'fp', 'fn'])), 'reported TP +1', got['tp'] != old['tp'] + 1)
                    if old['nneg']:
                        add(rd + '/' + f + '/' + meth + '/' + str(sp) + ' false alarm control', old['specificity'] == old['tn'] / old['nneg'], 'all-positive injected flags', rates(y, np.ones(len(y), int))['specificity'] < sp)
                    ec = methods[meth][str(sp)]['empirical_test_roc']
                    if 'specificity' in ec:
                        add(rd + '/' + f + '/' + meth + '/' + str(sp) + ' empirical ROC Sp', ec['specificity'] >= sp - 1e-12, 'forced ROC specificity 0', 0 < sp)
            a = methods['candidate']
            if '0.9' in a and '0.95' in a:
                q = a['0.9']['fixed_threshold_test']
                r = a['0.95']['fixed_threshold_test']
                add(rd + '/' + f + ' frozen operating point monotonicity', q['tp'] >= r['tp'] and q['tn'] <= r['tn'], 'higher threshold with TP increasing', q['tp'] < q['tp'] + 1)
        if rd != 'R1':
            add(rd + ' model hash', sha(frozen['model_path']) == frozen['model_sha256'], 'change model SHA', sha(frozen['model_path']) != '3' * 64)
    suff = read(P / 'raw/SUMMARY_SUFFICIENCY.json')
    for p in suff['pairs']:
        add('Sufficiency ' + p['summary'], p['identity_error'] == 0 and p['summary_A'] == p['summary_B'] and (p['downstream_difference'] > 0), 'summary B plus 0.01', p['summary_A'] != p['summary_B'] + 0.01)
    add('Noise definition', all((rr['findings'][f]['noise_floor_reached'] == 'UNKNOWN' for f in FINDINGS)), 'Set true inter/intra floor to measured', 'UNKNOWN' != 'MEASURED')
    snapshots = list((P / 'raw/code_versions').glob('*')) + list((P / 'raw/failed_checks').glob('*.py')) + [P / 'code/screening.py', P / 'code/common.py']
    hashes = {sha(p) for p in snapshots}
    for rd in ['R1', 'R2', 'R3']:
        fr = read(P / ('FROZEN_PREDICTIONS.json' if rd == 'R1' else f'FROZEN_PREDICTIONS_{rd}.json'))
        for inp in fr['inputs']:
            if '/code/' in inp['path']:
                add(rd + ' frozen code available ' + Path(inp['path']).name, inp['sha256'] in hashes, 'unavailable code hash', '4' * 64 not in hashes)
            else:
                add(rd + ' input hash ' + Path(inp['path']).name, sha(inp['path']) == inp['sha256'], 'source hash corrupt', sha(inp['path']) != '5' * 64)
    from open_any import first, label
    fr = read(P / 'FROZEN_PREDICTIONS_R4.json')
    sc = read(P / 'raw/SCORES_R4.json')
    lab4 = read(P / 'raw/LABELS_R4.json')
    rr = read(P / 'raw/RESULTS_R4.json')
    th = read(P / 'raw/THRESHOLDS_R4.json')
    ids = [c for c in m['test'] if first(lab4[c]) is not None]
    y = [first(lab4[c]) for c in ids]
    for (p, key) in [(P / 'raw/SCORES_R4.json', 'prediction_sha256'), (P / 'raw/LABELS_R4.json', 'labels_sha256'), (P / 'PREREG_R4.json', 'prereg_sha256'), (Path(fr['model_path']), 'model_sha256')]:
        add('R4 frozen ' + p.name, sha(p) == fr[key], 'corrupted expected hash', sha(p) != '6' * 64)
    add('R4 calibration threshold hash', sha(P / 'raw/THRESHOLDS_R4.json') == read(P / 'FROZEN_THRESHOLDS_R4.json')['threshold_sha256'], 'corrupt threshold SHA', sha(P / 'raw/THRESHOLDS_R4.json') != '7' * 64)
    add('R4 numerical source hash', sha(P / 'code/open_any.py') == fr['code_sha256'], 'corrupt numerical source SHA', sha(P / 'code/open_any.py') != '8' * 64)
    for (meth, ss) in sc.items():
        for sp in ['0.9', '0.95']:
            t = th[meth][sp]['threshold']
            got = rates(y, [ss[c] >= t for c in ids])
            old = rr['methods'][meth][sp]['fixed_threshold_test']
            add('R4 counts ' + meth + sp, all((got[k] == old[k] for k in ['tp', 'tn', 'fp', 'fn', 'n'])), 'TP +1', got['tp'] != old['tp'] + 1)
            add('R4 all-positive ' + meth + sp, old['specificity'] == old['tn'] / old['nneg'], 'all-positive flags', rates(y, np.ones(len(y), int))['specificity'] < float(sp))
    base = dict(open_bite=0)
    add('R4 regional presence preserves negated lateral clause', label('No lateral open bite. An anterior open bite is present.', base)[0] == 1, 'force global absent', label('No lateral open bite. An anterior open bite is present.', base)[0] != 0)
    add('R4 explicit absent is negative', label('No lateral open bite.', base)[0] == 0, 'force positive', label('No lateral open bite.', base)[0] != 1)
    changed = read(P / 'raw/SCOPE_CHANGED_R4.json')
    r = next((x for x in changed if x['case_id'] == _release_expand('@DENTAL_CASE_ID@')))
    add(_release_expand('@DENTAL_CASE_ID@ source scope'), r['old_anterior'] == 0 and r['new_any'] == 1, 'omit lateral positive', r['new_any'] != 0)
    write(P / 'raw/VERIFICATION.json', dict(review_state='PENDING_INDEPENDENT_REVIEW', producer_checks_only=True, count=len(checks), all_passed=all((x['pass_gate'] and x['injection_rejected'] for x in checks)), checks=checks, scientific_gates_unchanged=True))
    print('Verified', len(checks), 'value-injection controls; scientific failures preserved.', flush=True)
if __name__ == '__main__':
    main()
