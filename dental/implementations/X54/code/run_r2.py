"""Force-controlled queries, real-geometry sufficiency test and local crown edit."""
import copy, datetime, json, time
import numpy as np
from scipy.optimize import brentq
from contact_model import P, D, sha, write, state, geometry, solve
from run_r1 import external, area_shares

def force_control(g, target=100.0, **kwargs):
    evaluations = []
    answers = {}

    def fun(delta):
        (s, a) = solve(g, closure=delta, **kwargs)
        answers[delta] = (s, a)
        evaluations.append(dict(closure_mm=delta, vertical_force_N=s['total_vertical_N'], normal_force_N=s['total_normal_N'], solver_gate=s['numerics']['all_pass'], cost=s['cost']))
        return s['total_vertical_N'] - target
    lo = fun(-0.3)
    hi = fun(0.3)
    if lo > 0 or hi < 0:
        return (dict(status='REJECTED_FORCE_BRACKET', force_target_N=target, range_N=[lo + target, hi + target], evaluations=evaluations, reason='No supported root in frozen[-.3,.3]mm; fixed patch set may miss later contacts'), None)
    try:
        delta = brentq(fun, -0.3, 0.3, xtol=1e-11, rtol=1e-13, maxiter=60)
        if delta not in answers:
            fun(delta)
        (s, a) = answers[delta]
        s['force_control'] = dict(status='ROOT_FOUND', force_target_N=target, force_error_N=abs(s['total_vertical_N'] - target), closure_mm=delta, evaluations=evaluations, gate=abs(s['total_vertical_N'] - target) <= 1e-05)
        if s['total_normal_N'] <= 1e-07:
            s['argmax'] = {j: None for j in ('upper', 'lower')}
        return (s, a)
    except Exception as e:
        return (dict(status='REJECTED_ROOT', force_target_N=target, reason=str(e), evaluations=evaluations), None)

def persist(case, name, s, a):
    path = P / 'rounds/R2' / f'case{case:03d}_{name}.json'
    if a is not None:
        fp = D / f'R2_case{case:03d}_{name}.npz'
        np.savez_compressed(fp, **a)
        s['arrays'] = dict(path=str(fp), sha256=sha(fp))
    write(path, s)
    return dict(path=str(path), sha256=sha(path), arrays=s.get('arrays'))

def load_or_run(g, name, target=100.0, **kw):
    of = P / 'rounds/R2' / f"case{g['case']:03d}_{name}.json"
    if of.exists():
        s = json.loads(of.read_text())
        return (s, dict(path=str(of), sha256=sha(of), arrays=s.get('arrays')))
    (s, a) = force_control(g, target, **kw)
    return (s, persist(g['case'], name, s, a))

def main():
    spec = json.loads((P / 'PREREG_R2.json').read_text())
    tick = time.perf_counter()
    manifest = []
    allrows = []
    central = []
    geoms = []
    suff = []
    refine = []
    for case in spec['cases']:
        state('R2_RUNNING', 'R1 external improvementFAIL; force-control contract frozen', 'Force control, support permutation and local regional relief case' + str(case))
        g = json.loads((P / 'raw' / f'geometry_{case:03d}_h02.json').read_text())
        row = dict(case=case, force_targets={}, scenarios={})
        for target in spec['force_targets_N']:
            (s, m) = load_or_run(g, 'force' + str(int(target)), target)
            manifest.append(m)
            row['force_targets'][str(int(target))] = s
        base = row['force_targets']['100']
        if 'force_control' not in base:
            allrows.append(row)
            print('R2 case', case, 'no100N root', flush=True)
            continue
        central.append(base)
        geoms.append(g)
        for (name, kw) in [('mu0', dict(mu=0.0)), ('mu04', dict(mu=0.4)), ('no_crown', dict(crown=False)), ('lateral02_mu02', dict(lateral=0.02)), ('lateral02_mu0', dict(lateral=0.02, mu=0.0)), ('offset_minus05', dict(offset=-0.05)), ('offset_plus05', dict(offset=0.05))]:
            (s, m) = load_or_run(g, name, **kw)
            manifest.append(m)
            row['scenarios'][name] = s
        if all((t in g['centers']['upper'] for t in ('16', '26'))):
            scaleA = {('upper', 16): 0.5, ('upper', 26): 2.0}
            scaleB = {('upper', 16): 2.0, ('upper', 26): 0.5}
            (aa, ma) = load_or_run(g, 'support_A', support_scale=scaleA)
            (bb, mb) = load_or_run(g, 'support_B', support_scale=scaleB)
            manifest.extend([ma, mb])
            row['scenarios']['support_A'] = aa
            row['scenarios']['support_B'] = bb
            if 'force_control' in aa and 'force_control' in bb:
                delta = {jaw: max((abs(aa['force_shares_pp'][jaw].get(t, 0) - bb['force_shares_pp'][jaw].get(t, 0)) for t in aa['force_shares_pp'][jaw])) for jaw in ('upper', 'lower')}
                bgeom = json.loads(json.dumps(g))
                identity = float(np.max(np.abs(np.array([p['gap_mm'] for p in g['patches']]) - np.array([p['gap_mm'] for p in bgeom['patches']]))))
                multisetA = sorted([1 / 0.5, 1 / 2.0])
                multisetB = sorted([1 / 2.0, 1 / 0.5])
                suff.append(dict(case=case, summary='entire rigid gap array and per-tooth near-area array; total support reciprocal multiset', identity_error_mm=identity, area_identity_error_mm2=0.0, reciprocal_sum_identity_error=float(sum(multisetA) - sum(multisetB)), rigid_summary_sha256_A=sha(P / 'raw' / f'geometry_{case:03d}_h02.json'), rigid_summary_sha256_B=sha(P / 'raw' / f'geometry_{case:03d}_h02.json'), downstream_max_share_difference_pp=delta, argmaxA=aa['argmax'], argmaxB=bb['argmax'], summary_sufficient=bool(max(delta.values()) < 1e-12), smallest_extension_in_this_model='named tooth-specific6DOF compliance plus load/preload; root/PDL measurement replaces PHENOMENOLOGICAL tensor', force_target_identity='same exact100N input; solved output has separate finite root residual, not claimed machine-identical'))
        high = int(np.argmax([p['lambda_N'][0] for p in base['rows']]))
        edit = copy.deepcopy(g)
        p = edit['patches'][high]
        p['gap_mm'] += 0.03 * p['basis'][0][2]
        p['vertical_gap_mm'] += 0.03
        (s, m) = load_or_run(edit, 'relief30um')
        manifest.append(m)
        row['scenarios']['relief30um'] = s
        if 'force_control' in s:
            before = base['rows'][high]['lambda_N'][0]
            after = s['rows'][high]['lambda_N'][0]
            row['crown_edit'] = dict(patch_index=high, upper_fdi=p['upper_fdi'], region=p['region'], source_face=p['upper_source_face'], relief_mm=0.03, before_N=before, after_N=after, reduction_N=before - after, gate=before - after >= 1.0, argmax_before=base['argmax'], argmax_after=s['argmax'], scope='regional gap modification, fixed normals/area; not exported fullmesh or manufactured/measured surface')

        def diff(s):
            if 'force_control' not in s:
                return None
            return max((abs(s['force_shares_pp'][jaw].get(t, 0) - base['force_shares_pp'][jaw].get(t, 0)) for jaw in ('upper', 'lower') for t in base['force_shares_pp'][jaw]))
        row['scenario_max_share_difference_pp'] = {k: diff(v) for (k, v) in row['scenarios'].items()}
        if case in [1, 2, 3]:
            gf = P / 'raw' / f'geometry_{case:03d}_h01.json'
            if gf.exists():
                fine = json.loads(gf.read_text())
            else:
                fine = geometry(case, 0.1)
                write(gf, fine)
            ff = P / 'rounds/R2' / f'case{case:03d}_refine_R1.json'
            if ff.exists():
                f = json.loads(ff.read_text())
            else:
                (f, a) = solve(fine, closure=0.05)
                persist(case, 'refine_R1', f, a)
            coarse = json.loads((P / 'rounds/R1' / f'case{case:03d}_d0.05.json').read_text())
            dd = {jaw: max((abs(f['force_shares_pp'][jaw].get(t, 0) - coarse['force_shares_pp'][jaw].get(t, 0)) for t in set(f['force_shares_pp'][jaw]) | set(coarse['force_shares_pp'][jaw]))) for jaw in ('upper', 'lower')}
            refine.append(dict(case=case, grid_coarse_mm=0.2, grid_fine_mm=0.1, max_share_difference_pp=dd, gate=max(dd.values()) <= 2.0, coarse_regions=len(g['patches']), fine_regions=len(fine['patches']), interpretation='sampled regional convergence diagnostic; not rigorous geometry enclosure'))
        allrows.append(row)
        write(P / 'rounds/R2' / f'case{case:03d}_summary.json', row)
        print('R2 case', case, '100Nclosure', round(base['closure_mm'], 4), 'relief', row.get('crown_edit', {}).get('reduction_N'), 'refinement', refine[-1] if case in [1, 2, 3] else '', flush=True)
    fp = P / 'FROZEN_PREDICTIONS_R2.json'
    if not fp.exists():
        write(fp, dict(frozen_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), prereg_sha256=sha(P / 'PREREG_R2.json'), code_sha256=sha(P / 'code/run_r2.py'), manifest=manifest, physical_measurement_status='NOT_RUN', scope='Frozen force-controlled scenarios before external comparison; no parameter calibrated to comparison'))
        (P / 'FROZEN_PREDICTIONS_R2.json.sha256').write_text(sha(fp) + '\n')
    cmp = external(central, geoms) if central else {}
    regional = dict(summary='tooth total normal forceN', state_A_N=[32.0, 96.0], state_B_N=[96.0, 32.0], area_mm2=[0.25, 1.0], total_A_N=128.0, total_B_N=128.0, identity_error_N=0.0, mean_pressure_A_MPa=[128.0, 96.0], mean_pressure_B_MPa=[384.0, 32.0], max_pressure_difference_MPa=256.0, smallest_extension='regional force and contact area; peak pressure needs a resolved elastic pressure field', external_referent=dict(kind='our_own_fixture', locator='code/run_r2.py: exact binary two-region algebra', compared_quantity='force-sum sufficiency for regional traction', refutes_us=True))
    result = dict(round='R2', claim_type=['information_link', 'capability'], answer='Force-controlled conditional regional contact and crown relief available; rigid geometry summary is tested for sufficiency', external_referent=dict(kind='independent_measurement', locator='PMCIDPMC5726858 Table3; doi10.7144/sgf.2.111 Fig1', compared_quantity='jaw-matched bilateral tooth-type relative force/sensor signal', refutes_us=True), central_accepted=len(central), attempted=len(allrows), force_root_counts={str(int(t)): sum(('force_control' in row['force_targets'][str(int(t))] for row in allrows)) for t in spec['force_targets_N']}, dropout=dict(central_fraction=1 - len(central) / len(allrows), reasons={str(row['case']): row['force_targets']['100'].get('reason') for row in allrows if 'force_control' not in row['force_targets']['100']}), external_comparison=cmp, sufficiency=suff, regional_sufficiency=regional, refinement=refine, rows=allrows, cost=dict(this_run_wall_s=time.perf_counter() - tick, total_evaluations=sum((len(s.get('force_control', s).get('evaluations', [])) for row in allrows for s in list(row['force_targets'].values()) + list(row['scenarios'].values()))), inherited_fit='X11 original training total UNKNOWN, reused read-only'), rigorous_enclosure='MISSING; finite support/offset/mu scenarios are sensitivity diagnostics, not box enclosure; no linear sensitivity reported', manifest=manifest)
    write(P / 'rounds/R2/results.json', result)
    state('R2_COMPLETE', 'Force-controlled, support-sufficiency, relief and refinement gates saved', 'Freeze observation-driven regional force calibration/validation port', accepted=len(central))
if __name__ == '__main__':
    main()
