"""Expand the incidence set using fine-raster named regions, then exact gaps."""
import copy, json, time, datetime
import numpy as np
from contact_model import P, D, sha, write, state
from run_r3 import continuous, difference
from run_r2 import force_control

def main():
    start = time.perf_counter()
    rows = []
    manifest = []
    for case in [1, 2, 3]:
        state('R4_RUNNING', 'R3 gap minima checked; additional fine-raster regions remain', 'Expand contact incidence and re-solve case' + str(case))
        file = P / 'rounds/R4' / f'case{case:03d}.json'
        if file.exists():
            rows.append(json.loads(file.read_text()))
            continue
        fine = json.loads((P / 'raw' / f'geometry_{case:03d}_h01.json').read_text())
        coarse = json.loads((P / 'raw' / f'geometry_{case:03d}_h02.json').read_text())
        orig = {(p['upper_fdi'], p['lower_fdi'], p['region']) for p in coarse['patches']}
        known = {(p['upper_fdi'], p['lower_fdi'], p['region']) for p in fine['patches']}
        union = copy.deepcopy(fine)
        for patch in coarse['patches']:
            if (patch['upper_fdi'], patch['lower_fdi'], patch['region']) not in known:
                union['patches'].append(copy.deepcopy(patch))
        (g, checks) = continuous(union)
        write(P / 'raw' / f'geometry_{case:03d}_continuous_union.json', g)
        (s, a) = force_control(g, 100.0)
        if a is None:
            raise ValueError('R4_100N_ROOT_REJECTED:' + str(case))
        ap = D / f'R4_case{case:03d}.npz'
        np.savez_compressed(ap, **a)
        old = json.loads((P / 'rounds/R3' / f'case{case:03d}.json').read_text())['force_control_prediction']
        new = [i for (i, p) in enumerate(g['patches']) if (p['upper_fdi'], p['lower_fdi'], p['region']) not in orig]
        edit = copy.deepcopy(g)
        high = int(np.argmax([p['lambda_N'][0] for p in s['rows']]))
        edit['patches'][high]['gap_mm'] += 0.03 * edit['patches'][high]['basis'][0][2]
        edit['patches'][high]['vertical_gap_mm'] += 0.03
        (relieved, ra) = force_control(edit, 100.0)
        reduction = s['rows'][high]['lambda_N'][0] - relieved['rows'][high]['lambda_N'][0] if ra is not None else None
        delta = difference(old, s)
        row = dict(case=case, prediction=s, checks=checks, arrays=dict(path=str(ap), sha256=sha(ap)), original_patch_count=len(coarse['patches']), expanded_patch_count=len(g['patches']), new_patch_indices=new, new_contact_force_N=[s['rows'][i]['lambda_N'][0] for i in new], new_active_count=sum((s['rows'][i]['lambda_N'][0] > 1e-07 for i in new)), max_tooth_share_difference_pp=delta, whole_set_gate=delta <= 2.0, argmax_before=old['argmax'], argmax_after=s['argmax'], relief_reduction_N=reduction, relief_gate=reduction is not None and reduction >= 1.0, scope='continuous gaps over union of two sampled region inventories; topology completeness not enclosed')
        write(file, row)
        rows.append(row)
        manifest.append(dict(path=str(file), sha256=sha(file), arrays=row['arrays']))
        print('R4', case, 'newactive', row['new_active_count'], 'differencepp', delta, 'argmax', old['argmax'], s['argmax'], flush=True)
    fp = P / 'FROZEN_PREDICTIONS_R4.json'
    if not fp.exists():
        write(fp, dict(frozen_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), prereg_sha256=sha(P / 'PREREG_R4.json'), manifest=manifest, physical_measurement_status='NOT_RUN'))
        (P / 'FROZEN_PREDICTIONS_R4.json.sha256').write_text(sha(fp) + '\n')
    write(P / 'rounds/R4/results.json', dict(round='R4', claim_type='capability', answer='New named contact regions are included in the same force-controlled compliance query', rows=rows, all_whole_set_gates_pass=all((r['whole_set_gate'] for r in rows)), all_LP_controls_pass=all((c['pass_gate'] for r in rows for c in r['checks']['controls'])), all_NCP_pass=all((r['prediction']['numerics']['all_pass'] for r in rows)), external_referent=dict(kind='closed_form', locator='barycentric affine triangle gap minimum LP', compared_quantity='regional projected gapmm; force differences conditionalsimulation', refutes_us=True), physical_force_validation='UNKNOWN; R1/R2 external physical improvement gatesFAIL', cost=dict(wall_s=time.perf_counter() - start), rigorous_enclosure='MISSING: sampled union is not an exhaustive continuous contact candidate set'))
    state('R4_COMPLETE', 'Contact incidence expansion gates saved; physical support/preload/sensor still UNKNOWN', 'Independent physical calibration and geometric/elastic model validation')
if __name__ == '__main__':
    main()
