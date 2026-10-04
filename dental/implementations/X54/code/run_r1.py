import argparse, datetime, json, time, shutil
import numpy as np
from contact_model import P, D, BASE, ENGINE, LABELS, sha, write, state, geometry, solve

def profile_gate(prediction, reference, area):
    mae = float(np.mean(np.abs(np.asarray(prediction) - np.asarray(reference))))
    area_mae = float(np.mean(np.abs(np.asarray(area) - np.asarray(reference))))
    return dict(MAE_pp=mae, rigid_MAE_pp=area_mae, improvement_pp=area_mae - mae, gate=bool(mae <= 5 and area_mae - mae >= 2))

def area_shares(geom, jaw):
    x = geom['rigid_area_mm2'][jaw]
    total = sum(x.values())
    return {t: 100 * v / max(total, 1e-30) for (t, v) in x.items()}

def group(shares):
    return np.array([sum((v for (k, v) in shares.items() if int(k) % 10 == t)) for t in range(1, 8)])

def external(predictions, geoms):
    refpath = BASE / 'LANE_NEXT_P_OCCLUSION_VALIDATION/EXTERNAL_REFERENTS.json'
    refs = json.loads(refpath.read_text())['studies']
    (f, h) = refs[:2]
    fm = group({str(t): v for (t, v) in zip(f['fdi'], f['mean_percent'])})
    hm = 2 * np.asarray(h['single_side_mean_percent'])
    out = {}
    for (jaw, target, study) in [('upper', fm, f), ('lower', hm, h)]:
        pred = np.mean([group(x['force_shares_pp'][jaw]) for x in predictions], axis=0)
        area = np.mean([group(area_shares(g, jaw)) for g in geoms], axis=0)
        tested = profile_gate(pred, target, area)
        mae = tested['MAE_pp']
        am = tested['rigid_MAE_pp']
        out[jaw] = dict(source_id=study['id'], source_locator=study['pmcid'] if 'pmcid' in study else study['doi'], compared_quantity='bilateral I1/I2/C/P1/P2/M1/M2 force percent; jaw-matched', reference_pp=target, candidate_pp=pred, rigid_area_pp=area, MAE_pp=mae, rigid_MAE_pp=am, improvement_pp=am - mae, gate=tested['gate'], resolution='POPULATION_AGGREGATE_OF_PER_TOOTH', independence='reference profile not fitted; support and closures not selected against these numbers', physical_validation='UNKNOWN: unmatched cohort, unmeasured preload and sensor transfer', reference_sum_percent=float(target.sum()), predicted_sum_percent=float(pred.sum()))
        corrupt = pred + 50
        out[jaw]['injected_50pp_error_rejected'] = not profile_gate(corrupt, target, area)['gate']
    out['reference_file'] = str(refpath)
    out['reference_sha256'] = sha(refpath)
    return out

def main(args):
    prereg = json.loads((P / 'PREREG_R1.json').read_text())
    start = time.perf_counter()
    D.mkdir(parents=True, exist_ok=True)
    if sum((p.stat().st_size for p in D.rglob('*') if p.is_file())) > 3000000000:
        raise RuntimeError('DISK_LIMIT')
    geoms = []
    central = []
    manifest = []
    errors = []
    source_dir = P / 'rounds/R1/source_v1'
    source_dir.mkdir(parents=True, exist_ok=True)
    for fn in ['contact_model.py', 'run_r1.py', 'preregister.py']:
        dest = source_dir / fn
        if not dest.exists():
            shutil.copyfile(P / 'code' / fn, dest)
    for case in prereg['cases'][:args.cases]:
        state('R1_RUNNING', 'PREREG_R1 frozen; mechanical/physical acceptance not tested', f'case{case}: build and solve conditional0/.05/.1mm contact', completed=len(central))
        try:
            gf = P / 'raw' / f'geometry_{case:03d}_h02.json'
            if gf.exists():
                g = json.loads(gf.read_text())
            else:
                g = geometry(case)
                write(gf, g)
            geoms.append(g)
            for closure in prereg['loading']['closure_mm']:
                of = P / 'rounds/R1' / f'case{case:03d}_d{closure:.2f}.json'
                if of.exists():
                    s = json.loads(of.read_text())
                else:
                    (s, a) = solve(g, closure, controls=case == 1 and closure == 0.05)
                    arraypath = D / f'R1_case{case:03d}_d{closure:.2f}.npz'
                    np.savez_compressed(arraypath, **a)
                    s['arrays'] = dict(path=str(arraypath), sha256=sha(arraypath))
                    s['geometry_sha256'] = sha(gf)
                    write(of, s)
                manifest.append(dict(path=str(of), sha256=sha(of), arrays=s['arrays']))
                print('case', case, 'closure', closure, 'contacts', len(s['rows']), 'N', round(s['total_normal_N'], 3), 'residual', s['numerics']['natural_residual_N'], 's', round(s['cost']['solve_s'], 2), flush=True)
                if closure == 0.05:
                    central.append(s)
        except Exception as e:
            err = dict(case=case, type=type(e).__name__, error=str(e))
            errors.append(err)
            write(P / 'raw' / f'R1_failure_{case:03d}.json', err)
            print(err, flush=True)
    frozen = P / 'FROZEN_PREDICTIONS.json'
    if not frozen.exists():
        write(frozen, dict(frozen_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), prereg_sha256=sha(P / 'PREREG_R1.json'), code_sha256={x: sha(P / 'code' / x) for x in ['contact_model.py', 'run_r1.py']}, manifest=manifest, scope='All R1 predicted regional force rows hashed before external comparison; no future physical measurement read', physical_measurement_status='NOT_RUN'))
        (P / 'FROZEN_PREDICTIONS.json.sha256').write_text(sha(frozen) + '\n')
    cmp = external(central, geoms[:len(central)]) if central else {}
    decisions = []
    for (s, g) in zip(central, geoms):
        entry = dict(case=s['case'], rigid_nearest={}, rigid_area={}, compliant=s['argmax'], changed={}, max_share_difference_pp={})
        for jaw in ('upper', 'lower'):
            areas = area_shares(g, jaw)
            entry['rigid_area'][jaw] = max(areas, key=areas.get)
            entry['rigid_nearest'][jaw] = min(g['nearest_gap_mm'][jaw], key=g['nearest_gap_mm'][jaw].get)
            entry['changed'][jaw] = entry['rigid_area'][jaw] != s['argmax'][jaw]
            entry['max_share_difference_pp'][jaw] = max((abs(s['force_shares_pp'][jaw].get(t, 0) - v) for (t, v) in areas.items()))
        decisions.append(entry)
    change = {jaw: dict(changed=sum((x['changed'][jaw] for x in decisions)), n=len(decisions), fraction=sum((x['changed'][jaw] for x in decisions)) / max(len(decisions), 1), large_differences=sum((x['max_share_difference_pp'][jaw] >= 10 for x in decisions))) for jaw in ('upper', 'lower')}
    result = dict(round='R1', claim_type=['information_link', 'capability'], answer='Conditional regional reactions and friction branches now executable; individual force remains unidentified', external_referent=dict(kind='independent_measurement', locator='PMCID:PMC5726858 Table3; doi:10.7144/sgf.2.111 Fig1', compared_quantity='jaw-matched bilateral per-tooth-type relative force/sensor percentage', refutes_us=True), cases=len(central), preregistered_cases=len(prereg['cases']), numerical_all_pass=all((s['numerics']['all_pass'] for s in central)), controls=[s['controls'] for s in central if s['controls']], decisions=decisions, decision_change=change, external_comparison=cmp, errors=errors, unloaded_pose_compatible=sum((g['census']['unloaded_pose_compatible'] for g in geoms)), dropout=dict(case_fraction=len(errors) / max(args.cases, 1), cases_rejected=errors, patch_normal_filter_mean_fraction=float(np.mean([g['census']['candidate_dropout_fraction'] for g in geoms])) if geoms else None), resolution='PER_SURFACE_REGION; external per-tooth/POPULATION aggregation', rigorous_enclosure='MISSING, no affine sensitivity claim', full_cost=dict(this_run_wall_s=time.perf_counter() - start, solve_wall_s=sum((x['cost']['total_s'] for x in central)), inherited_fit='X11 training/inference reused, original cumulative cost UNKNOWN', human_provider_cost='UNKNOWN'), manifest=manifest)
    write(P / 'rounds/R1/results.json', result)
    state('R1_FORWARD_COMPLETE', 'External numerical profile gates recorded; physical validation UNKNOWN', 'Write R1 handoff; freeze support/pose sufficiency and local crown edit construction', cases=len(central), result_file='rounds/R1/results.json')
if __name__ == '__main__':
    a = argparse.ArgumentParser()
    a.add_argument('--cases', type=int, default=12)
    main(a.parse_args())
