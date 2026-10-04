from registration import *

def main():
    t = time.monotonic()
    g = extract()
    r = json.loads((ROOT / 'raw/R1_RESULTS.json').read_text())
    p = g['points']
    q = g['target']
    rng = np.random.default_rng(69)
    q = q[rng.choice(len(q), min(18000, len(q)), replace=False)]
    (T, loss, n) = icp(p, q, np.array(r['source_to_target']))
    (d, w, c) = nearest_surface(transform(p, T), g['target'], g['faces'])
    hold = g['region'] == 1
    result = {'kind': 'Conventional full-support rigid ICP same information, R1 fit seed then all crown points', 'claim_type': 'control_only_no_algorithm_claim', 'fit_loss_mm2': loss, 'iterations': n, 'source_to_target': T.tolist(), 'heldout_summary': stats(d[hold]), 'all_summary': stats(d), 'seconds': time.monotonic() - t, 'fault_control': {'injection': 'Translate recovered pose by10 mm inx', 'frozen_rejection_threshold_mm': 0.5}}
    bad = T.copy()
    bad[0, 3] += 10
    (db, _, _) = nearest_surface(transform(p[hold][::20], bad), g['target'], g['faces'])
    result['fault_control'].update(p95_mm=float(np.percentile(db, 95)), fault_rejected=bool(np.percentile(db, 95) > 0.5))
    dump(ROOT / 'raw/R1_CONTROL.json', result)
if __name__ == '__main__':
    main()
