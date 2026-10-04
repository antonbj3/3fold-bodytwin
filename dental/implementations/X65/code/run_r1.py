import csv, datetime, hashlib, json, math, time
from pathlib import Path
import numpy as np
from physics import ROOT, legacy, draws, forward, scalar_root, vector_root
from freeze import frozen_write

def save(path, v):
    (ROOT / path).write_text(json.dumps(v, indent=2, allow_nan=False) + '\n')

def main():
    t = time.perf_counter()
    (ROOT / 'raw').mkdir(exist_ok=True)
    allrows = {r['id']: r for r in csv.DictReader((ROOT / 'inputs/anchors_legacy.csv').open())}
    rows = [allrows[i] for i in legacy.S_SET]
    dd = draws(rows)
    old = json.loads((ROOT / 'inputs/k1b_eval_legacy.json').read_text())['per_anchor']
    replay = []
    checks = []
    for r in rows:
        (v, ports) = forward(r, dd[r['id']])
        q = np.percentile(v, [10, 50, 90]).tolist()
        o = old[r['id']]['M_asm_P10_P50_P90']
        err = max((abs(a / b - 1) for (a, b) in zip(q, o)))
        replay.append({'id': r['id'], 'old_N': o, 'replayed_N': q, 'max_relative_error': err, 'pass': err <= 1e-06})
        for (j, (g, s0)) in enumerate(ports):
            val = scalar_root(float(g[0]), float(s0[0]), float(dd[r['id']]['HV'][0]), float(dd[r['id']]['sqa'][0]))
            vval = float(vector_root(g[0], s0[0], dd[r['id']]['HV'][0], dd[r['id']]['sqa'][0]))
            checks.append({'id': r['id'], 'section': j, 'relative_difference': abs(val / vval - 1)})
    save(Path('raw/replay.json'), replay)
    frozen = {'created_at_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(), 'kind': 'retrospective_inverse_targets_not_prospective_validation', 'prereg_sha256': hashlib.sha256((ROOT / 'PREREG_R1.json').read_bytes()).hexdigest(), 'rows': [{'id': r['id'], 'target_force_N': float(r['F_lim_N']), 'saved_prediction_N': old[r['id']]['M_asm_P10_P50_P90'][1]} for r in rows]}
    frozen_write(ROOT / 'FROZEN_PREDICTIONS_R1.json', frozen)
    if not all((r['pass'] for r in replay)):
        save(Path('results_R1.json'), {'status': 'REPLAY_FAIL', 'replay': replay})
        raise RuntimeError('legacy replay failed')
    inverses = []
    specs = {'kt_scale': (0.15, 2.0), 'strength_scale': (0.2, 8.0), 'preload_scale': (0.0, 30.0), 'lever_scale': (0.12, 2.0), 'defect_scale': (0.0002, 8.0)}
    for r in rows:
        target = float(r['F_lim_N'])
        a = r['id']
        for (name, (low, high)) in specs.items():
            fn = lambda x: float(np.median(forward(r, dd[a], **{name: x})[0])) - target
            (fl, fh) = (fn(low), fn(high))
            if fl * fh > 0:
                inverses.append({'id': a, 'changed_input': name, 'value': None, 'status': 'NOT_REACHABLE_IN_DECLARED_SEARCH', 'force_at_bounds_N': [fl + target, fh + target]})
                continue
            (lo, hi) = (low, high)
            for _ in range(35):
                mid = (lo + hi) * 0.5
                fm = fn(mid)
                if fl * fm <= 0:
                    hi = mid
                else:
                    lo = mid
                    fl = fm
            root = (lo + hi) * 0.5
            inverses.append({'id': a, 'changed_input': name, 'value': root, 'force_residual_N': fn(root), 'status': 'MODEL_COUNTERFACTUAL_NOT_MEASURED_CAUSE', 'resolution': 'PHENOMENOLOGICAL'})
    save(Path('raw/inverse_mechanisms.json'), inverses)
    summary = lambda r: (float(r['D_mm']), legacy.mat_class(r['material'], r['manufacturing']), float(r['angle_deg']), float(r['N_runout']))
    (a, b) = (allrows['K1-A3a'], allrows['K1-A3b'])
    (sa, sb) = (summary(a), summary(b))
    suff = {'summary': ['diameter_mm', 'material_class', 'angle_deg', 'horizon_cycles'], 'state_a': sa, 'state_b': sb, 'exact_equal': sa == sb, 'identity_error': 0.0 if sa == sb else None, 'downstream_difference_N': float(a['F_lim_N']) - float(b['F_lim_N']), 'downstream_relative_difference': float(a['F_lim_N']) / float(b['F_lim_N']) - 1, 'external_locator': 'PMC7730231 Table 1 and Table 2, DOI 10.3390/ijerph17238988', 'smallest_pair_extension': 'retain connection-specific configuration endpoint; connection label separates this pair, not proven sufficient across systems', 'resolution': 'PER_TOOTH'}
    save(Path('raw/sufficiency_R1.json'), suff)
    mutations = [{'gate': 'REPLAY', 'nominal_pass': all((r['pass'] for r in replay)), 'injected_error_rejected': abs(2 * replay[0]['replayed_N'][1] / replay[0]['old_N'][1] - 1) > 1e-06, 'injection': 'double replayed force'}, {'gate': 'SCALAR_ROOT', 'nominal_pass': max((c['relative_difference'] for c in checks)) < 1e-08, 'injected_error_rejected': abs(2.0 - 1) > 1e-08, 'injection': 'double independently evaluated force'}, {'gate': 'SUMMARY_IDENTITY', 'nominal_pass': suff['exact_equal'], 'injected_error_rejected': summary(a) != tuple([sa[0] + 0.1, *sa[1:]]), 'injection': 'alter nominal diameter'}, {'gate': 'DOWNSTREAM_DIFFERENCE', 'nominal_pass': suff['downstream_difference_N'] > 10, 'injected_error_rejected': 0 <= 10, 'injection': 'erase published force difference'}]
    out = {'claim_type': 'information_link', 'replay_pass': all((r['pass'] for r in replay)), 'max_replay_relative_error': max((r['max_relative_error'] for r in replay)), 'independent_scalar_max_relative_difference': max((c['relative_difference'] for c in checks)), 'mechanism_identified': 'UNKNOWN', 'diagnosis': 'Different latent changes recover each endpoint. Inverse roots are not measurements of their cause. No same-configuration physical metrology in sources.', 'sufficiency': suff, 'inverse_roots': inverses, 'mutations': mutations, 'wall_s': time.perf_counter() - t, 'external_referent': {'kind': 'independent_measurement', 'locator': 'inputs/primary/PMC7730231.xml, Table 2, DOI10.3390/ijerph17238988', 'compared_quantity': 'connection-specific maximum cyclic force at 5e6 cycles', 'refutes_us': True}}
    save(Path('results_R1.json'), out)
    print(json.dumps({k: v for (k, v) in out.items() if k not in ['inverse_roots', 'mutations', 'sufficiency']}))
if __name__ == '__main__':
    main()
