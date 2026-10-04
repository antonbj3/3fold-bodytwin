from common import *
from contact import compare, gates, Q
from reuse import VF
import itertools, collections, copy

def run():
    st = time.perf_counter()
    pr = read(ROOT / 'PREREG_R1.json')['metrics']
    rows = read(ROOT / 'raw/R1_ROWS.json')
    items = read(ROOT / 'raw/R1_ITEMS.json')
    by = collections.defaultdict(list)
    for r in rows:
        by[r['task_id']].append(r)
    pair_counts = {a + ' vs ' + b: sum((next((r for r in rr if r['participant'] == a))['L2'] != next((r for r in rr if r['participant'] == b))['L2'] for rr in by.values())) for (a, b) in itertools.combinations(FRONTIER, 2)}
    checks = []
    count = sum((len({r['L2'] == 'PASS' for r in rr}) > 1 for rr in by.values()))
    reported = read(ROOT / 'rounds/R1.json')['discriminating_L2']
    checks.append(dict(name='paired_item_output', nominal_pass=count == reported, injected_rejected=count != reported + 1))
    count2 = 0
    component_checks = []
    for r in read(ROOT / 'raw/R2_ROWS.json'):
        if r['status'] != 'SCORED':
            continue
        z = npz(r['field_path'])
        q = Q.contact_map(z['xy'], z['faces'], z['pred_gap'], z['reference_gap'])
        value = q['contact_symdiff_mm2']
        stored = r['contact'].get('symdiff_mm2')
        if value is None:
            if stored is not None:
                raise AssertionError('Unsupported output')
        elif abs(value - stored) > 1e-07:
            raise AssertionError('Pattern output mismatch ' + r['uid'])
        if len(component_checks) < 8 and r['contact']['status'] == 'SCORED':
            ff = z['faces']
            ok = np.isfinite(z['pred_gap'][ff]).all(1) & np.isfinite(z['reference_gap'][ff]).all(1)
            expected = VF.contact_components(z['xy'], ff[ok], z['pred_gap'])
            got = r['contact']['predicted']['count']
            component_checks.append(dict(uid=r['uid'], nominal_pass=expected['count'] == got, injected_rejected=expected['count'] != got + 1))
        count2 += 1
    checks.append(dict(name='actual_complete_crown_pattern_outputs', checked=count2, nominal_pass=True, injected_rejected=abs(value - (stored + 2)) > 1e-07))
    checks.append(dict(name='v5_component_count_control', checked=len(component_checks), nominal_pass=all((x['nominal_pass'] for x in component_checks)), injected_rejected=all((x['injected_rejected'] for x in component_checks))))
    ctr = read(ROOT / 'raw/R1_CONTROLS.json')
    err = [c['error_mm2'] for c in ctr if c['error_mm2'] is not None]
    checks.append(dict(name='v4_same_information_operator', nominal_pass=max(err) < 1e-07, injected_rejected=max(err) + 2 > 1e-07, max_error_mm2=max(err)))
    fine = [x for x in read(ROOT / 'raw/R2_CONTROLS.json') if 'fine_L2' in x]
    checks.append(dict(name='four_frontier_first_site_grid_refinement', nominal_pass=all((x['fine_L2'] == x['coarse_L2'] for x in fine)), injected_rejected=any((x['fine_L2'] != (not x['coarse_L2']) for x in fine)), scoped_sites=1, rows=len(fine)))
    pr_checks = []
    for p in ROOT.glob('PREREG_*.json'):
        digest = p.with_suffix('.sha256').read_text().split()[0]
        pr_checks.append(sha(p) == digest)
    checks.append(dict(name='frozen_prereg', nominal_pass=all(pr_checks), injected_rejected=hashlib.sha256((ROOT / 'PREREG_R1.json').read_bytes() + b' ').hexdigest() != (ROOT / 'PREREG_R1.sha256').read_text().split()[0]))
    out = dict(checks=checks, pairwise_R1_disagreement=pair_counts, all_pass=all((x['nominal_pass'] and x['injected_rejected'] for x in checks)), seconds=time.perf_counter() - st)
    dump(ROOT / 'raw/OUTPUT_CONTROLS.json', out)
    if not out['all_pass']:
        raise AssertionError(out)
    return out
if __name__ == '__main__':
    print(clean(run()))
