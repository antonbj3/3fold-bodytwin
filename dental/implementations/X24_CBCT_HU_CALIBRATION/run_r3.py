from dental_release.paths import expand as _release_expand
import os
from pathlib import Path
import json, hashlib, gzip, sys, time, math
import numpy as np
from calibration import bone_branch_b, apply
P = Path(__file__).resolve().parent
ROOT = P.parents[1]
pre = json.loads((P / 'PREREG_R3.json').read_text())
sys.path.insert(0, str(ROOT / 'cells/design'))
import design_primary_stability as DP
import implant_insertion as II
from decision_certificate import certify_decision
INSERT = Path(_release_expand('@DENTAL_EXTERNAL_ROOT@/media/sdc1-tmp/datasets/dental_3fold_extra/scratch_INSERT'))
SITES = Path(_release_expand('@DENTAL_EXTERNAL_ROOT@/media/sdc1-tmp/datasets/dental_3fold_extra/scratch_DESIGN/sites'))

def sha(p):
    h = hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda : f.read(1048576), b''):
            h.update(b)
    return h.hexdigest()

def clean(d):
    if isinstance(d, np.generic):
        return clean(d.item())
    if isinstance(d, dict):
        return {k: clean(v) for (k, v) in d.items()}
    if isinstance(d, list):
        return [clean(v) for v in d]
    if isinstance(d, float) and (not math.isfinite(d)):
        return None
    return d
OUT = Path(os.environ.get('X24_OUTPUT_DIR', str(P)))
OUT.mkdir(parents=True, exist_ok=True)

def run():
    t0 = time.perf_counter()
    assert sha(P / 'PREREG_R3.json') == (P / 'PREREG_R3.sha256').read_text().strip()
    bas = json.loads((P / 'PHANTOM_BASELINE.json').read_text())
    R1 = json.loads((P / 'results_R1.json').read_text())
    V = {r['material']: r['gray_median'] for r in next((d for d in bas if d['device'] == 'Varian'))['inserts']}
    obstruction = []
    for d in bas:
        G = {r['material']: r['gray_median'] for r in d['inserts']}
        if G['PMP'] >= G['LDPE'] and V['PMP'] < V['LDPE']:
            obstruction.append({'device': d['device'], 'material_pair': ['PMP', 'LDPE'], 'gray_pair': [G['PMP'], G['LDPE']], 'reference_HU_pair': [V['PMP'], V['LDPE']], 'monotone_minimax_lower_bound_HU': (V['LDPE'] - V['PMP']) / 2, 'refutes_40HU_goal': (V['LDPE'] - V['PMP']) / 2 > 40, 'resolution': 'PER_SURFACE_REGION'})
    rhos = []
    rawrows = []
    source = []
    for f in sorted((INSERT.parent / 'scratch_DESIGN2/main').glob('v2_*.jsonl')):
        source.append({'path': str(f), 'sha256': sha(f)})
        rawrows.extend((json.loads(s) for s in f.read_text().splitlines() if s))
    assert len({(r['case'], r['site']) for r in rawrows}) == len(rawrows)
    materials = {}
    for f in sorted(INSERT.glob('B_picks_*.jsonl')):
        source.append({'path': str(f), 'sha256': sha(f)})
        for s in f.read_text().splitlines():
            r = json.loads(s)
            materials[r['case'], r['site']] = r.get('material')
    drops = []
    answers = []
    certs = []
    fault = []
    candidate_total = 0
    candidate_drop = 0
    maxinverse = 0
    for r in rawrows:
        case = r['case']
        tooth = r['site']
        key = (case, tooth)
        m = materials.get(key)
        if not m:
            drops.append({'case': case, 'site': tooth, 'reason': 'no material'})
            continue
        rho = m['rho_app']
        H = (rho / 2 - 0.19) / 0.001
        Eref = float(bone_branch_b(H)['E_MPa'])
        if not bool(bone_branch_b(H)['Keller_ash_range']) or not bool(bone_branch_b(H)['published_NewTom_range']):
            drops.append({'case': case, 'site': tooth, 'reason': 'outside_joint_Keller_NewTom_domain', 'rho_app': rho, 'H': H})
            continue
        meta = INSERT / ('A_' + case + '.json')
        npz = INSERT / ('A_' + case + '.npz')
        cf = SITES / (case + '.json.gz')
        if not all((p.exists() for p in [meta, npz, cf])):
            drops.append({'case': case, 'site': tooth, 'reason': 'missing profiles/candidates'})
            continue
        source.extend(({'path': str(p), 'sha256': sha(p)} for p in [meta, npz, cf]))
        A = json.loads(meta.read_text())['sites'][str(tooth)]
        groups = {(tuple(g[0]), g[1]): i for (i, g) in enumerate(A['groups'])}
        with gzip.open(cf, 'rt') as f:
            allcase = json.load(f)
        site = next((s for s in allcase['sites'] if s['site'] == tooth))
        shortlist = []
        seen = set()
        for nm in pre['shortlist']:
            x = r[nm] if nm == 'B1' else r.get(nm, {}).get('x')
            if x is None:
                candidate_drop += 1
                candidate_total += 1
                continue
            k = (tuple(x[0]), x[1], x[2])
            if k in seen:
                continue
            seen.add(k)
            candidate_total += 1
            c = next((c for c in site['candidates'] if (tuple(c['pose']), c['D'], c['L']) == k), None)
            if c is None or (tuple(x[0]), x[1]) not in groups:
                candidate_drop += 1
                continue
            shortlist.append((nm, x, c))
        if not shortlist:
            drops.append({'case': case, 'site': tooth, 'reason': 'empty matched shortlist'})
            continue

        def eval_E(E):
            cc = []
            with np.load(npz) as z:
                for (nm, x, c) in shortlist:
                    gi = groups[tuple(x[0]), x[1]]
                    can = II.bone_cancellous(m['sigma_c'], E)
                    a = DP.anchorage(x[1], x[2], x[0][4], z['nb_' + str(tooth)][gi], z['nc1.0_' + str(tooth)][gi], can, x[1] - 0.5, tc=1.0, crest_band=True)
                    geo = c['wall'] >= 1 and c['adj'] >= 1.5 and (c['J1']['base'] <= 0.1)
                    cc.append({'name': nm, 'x': x, 'u_um': a['u'], 'k_N_per_um': a['k'], 'necro': a['necro'], 'geometry_eligible': geo, 'eligible': geo and (not a['necro']), 'converged': math.isfinite(a['u']), 'E_MPa': E})
            feas = [c for c in cc if c['eligible'] and c['converged']]
            choice = max(feas, key=lambda c: c['k_N_per_um'])['x'] if feas else None
            return (cc, choice)
        (cref, xref) = eval_E(Eref)
        hbudget = (H + 190) * min(1.05 ** (1 / 2.57) - 1, 1 - 0.95 ** (1 / 2.57))
        for sc in R1['scanner_calibration']:
            c = np.array(sc['a_inverse_Href'])
            g = (H - c[0]) / c[1]
            bad = bone_branch_b(g)
            Ebad = float(bad['E_MPa'])
            (cbad, xbad) = eval_E(Ebad)
            Hcorrect = float(apply(c, g))
            Ecorrect = float(bone_branch_b(Hcorrect)['E_MPa'])
            maxinverse = max(maxinverse, abs(Ecorrect / Eref - 1))
            chosenref = next((c for c in cref if c['x'] == xref), None)
            chosenbadref = next((c for c in cbad if c['x'] == xref), None)
            u = chosenref['u_um'] if chosenref else float('nan')
            margins = [{'id': 'u50', 'unit': 'um', 'definition': '50 minus modeled crest micromotion at30N/7mm', 'reports': [{'margin': 50 - u if math.isfinite(u) else None, 'sigma': None, 'node': 'scanner', 'unit': 'um'}]}]
            cert = certify_decision(margins, {'required_margin_ids': ['u50'], 'alpha': 0.05})
            certs.append({'case': case, 'site': tooth, 'device': sc['device'], 'status': cert.status, 'missing': cert.missing, 'resolution': 'PER_TOOTH', 'scope': 'partial u50 predicate; full K42 also needs anatomy/material/load/contact law'})
            if sc['device'] != 'Varian':
                aa = -c[0] / c[1]
                bb = 1 / c[1]
                gs = 165.0
                gd = 165.0 + 1048.8599348534203 * 1.535
                gt = gs + 1048.8599348534203 * 0.598 * rho
                normalized = (aa + bb * gt - (aa + bb * gs)) / (aa + bb * gd - (aa + bb * gs)) * 1.535 / 0.598
                corrupted = (aa + bb * gt - (aa + bb * gs + 100)) / (aa + bb * gd - (aa + bb * gs + 100)) * 1.535 / 0.598
                fault.append({'device': sc['device'], 'case': case, 'site': tooth, 'normalization_relative_error': abs(normalized / rho - 1), 'anchor_plus100gray_relative_error': abs(corrupted / rho - 1), 'injected_rejected_5percent': abs(corrupted / rho - 1) > 0.05})
            answers.append({'case': case, 'site': tooth, 'device': sc['device'], 'rho_app_reference_closure': rho, 'H_reference_closure': H, 'gray_affine_scenario': float(g), 'E_reference_closure_MPa': Eref, 'E_naive_MPa': Ebad, 'E_naive_ratio': Ebad / Eref, 'naive_density_valid': bool(bad['density_valid']), 'reference_shortlist': cref, 'naive_shortlist': cbad, 'reference_choice': xref, 'naive_clipped_choice': xbad, 'choice_flip_including_abstention': xref != xbad, 'qualified_naive_choice': xbad if bool(bad['density_valid']) else None, 'u_ratio_at_reference_choice': chosenbadref['u_um'] / chosenref['u_um'] if chosenbadref and chosenref and math.isfinite(chosenbadref['u_um']) and (chosenref['u_um'] > 0) else None, 'u50_reference': bool(chosenref and chosenref['u_um'] <= 50), 'u50_naive_at_reference_choice': bool(chosenbadref and chosenbadref['u_um'] <= 50), 'HU_symmetric_budget_5percent_E': hbudget, 'resolution': 'PER_TOOTH', 'model_quantity_resolution': 'PHENOMENOLOGICAL', 'replacement_measurement': 'site-matched HA calibration and directional mechanical test; interface/micromotion experiment', 'status': 'CONDITIONAL_SCENARIO_NOT_PATIENT_TRUTH'})
        print(len(answers) // 4, case, tooth, flush=True)
    if not answers:
        raise RuntimeError('no valid domain sites')
    ele = [a for a in answers if a['device'] != 'Varian']
    summary = {}
    for d in sorted({a['device'] for a in ele}):
        a = [a for a in ele if a['device'] == d]
        rat = [r['E_naive_ratio'] for r in a]
        ur = [r['u_ratio_at_reference_choice'] for r in a if r['u_ratio_at_reference_choice'] is not None]
        summary[d] = {'n_sites': len(a), 'E_naive_ratio_min_median_max': [min(rat), float(np.median(rat)), max(rat)], 'u_ratio_min_median_max': [min(ur), float(np.median(ur)), max(ur)] if ur else None, 'restricted_shortlist_choice_flips': sum((r['choice_flip_including_abstention'] for r in a)), 'naive_density_invalid': sum((not r['naive_density_valid'] for r in a)), 'u50_predicate_flips': sum((r['u50_reference'] != r['u50_naive_at_reference_choice'] for r in a)), 'resolution': 'PER_TOOTH', 'all_are_conditional': True}
    rr = json.loads((P / 'results_R2.json').read_text())['rows']
    fault_ref = all((abs(r['error_HU'] + 1000) > 40 for r in rr))
    gate = {'affine_inverse_replay': maxinverse <= 1e-09, 'normalization_invariance': max((r['normalization_relative_error'] for r in fault)) <= 1e-09, 'wrong_reference_rejected': fault_ref, 'wrong_anchor_rejected': all((r['injected_rejected_5percent'] for r in fault)), 'K42_unknown_abstains': all((c['status'] == 'ABSTAIN_UNKNOWN' for c in certs))}
    tolerances = [a['HU_symmetric_budget_5percent_E'] for a in answers]
    out = {'round': 'R3', 'claim_type': 'capability', 'outcome': 'SCOPED_ANSWER_WITH_CONDITIONAL_CONSUMER_AND_ABSOLUTE_UNKNOWN', 'external_referent': pre['external_referent'], 'prereg_sha256': sha(P / 'PREREG_R3.json'), 'gates': gate, 'monotone_order_obstruction': obstruction, 'max_affine_inverse_relative_error': maxinverse, 'HU_symmetric_budget_5percent_E_min_median_max': [min(tolerances), float(np.median(tolerances)), max(tolerances)], 'summary_by_device': summary, 'answers': clean(answers), 'K42_certificates': certs, 'fault_injections': fault, 'dropout': {'source_sites': len(rawrows), 'retained_sites': len(answers) // 4, 'rejected_sites': len(drops), 'rejected_fraction': len(drops) / len(rawrows), 'reasons': drops, 'shortlist_candidates_considered': candidate_total, 'shortlist_candidates_rejected': candidate_drop, 'note': 'identical shortlist choices deduplicated; all rejected sites retained as rows; no full optimization rerun'}, 'cost': {'wall_s': time.perf_counter() - t0, 'patient_volumes': 0, 'new_physical_measurements': 0, 'calls': sum((len(a['reference_shortlist']) + len(a['naive_shortlist']) for a in answers))}, 'limitations': ['no measured E or implant outcome at these sites', 'radiotherapy CBCT response cannot be presumed dental-machine response', 'original K37 uses internally normalized gray/ranks: pure affine scanner effect cancels; raw HU result is explicitly counterfactual practice', 'K37 shortlist3 not full multigoal catalog', 'clipped1MPa for impossible density reported numerically but rejected physically', 'single phantom composition and field spatial variation do not calibrate bone law']}
    (OUT / 'results_R3.json').write_text(json.dumps(clean(out), indent=2, allow_nan=False) + '\n')
    unique = {r['path']: r for r in source}
    (OUT / 'SOURCE_MANIFEST_K37.json').write_text(json.dumps(list(unique.values()), indent=2) + '\n')
    (OUT / 'HANDOFF_R3.md').write_text('# R3 : material fields and consumer bounded\n\n' + json.dumps(summary, indent=2) + '\n\n' + json.dumps(clean(gate)) + '\n\nAbsolute patient-E and complete implant decision UNKNOWN. Next physical construction: water + two certified HA levels in the bone area and a separate HA intermediate step in the same image, repeated position/protocol validation and paired mandibulatory mechanics test.\n')
    s = json.loads((P / 'CURRENT_WORK_STATE.json').read_text())
    s.update(phase='R3_COMPLETE', latest_gate=clean(gate), next_operation='package one-command demo, independent referents, minimum measurement contract and figure')
    (OUT / 'CURRENT_WORK_STATE.json').write_text(json.dumps(s, indent=2) + '\n')
    print(summary, gate)
if __name__ == '__main__':
    run()
