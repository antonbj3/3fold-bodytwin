from common import *
import time, csv
from scipy.stats import spearmanr

def valid_port(r, p, j, f, t):
    return r.get('patient') == p and r.get('jaw') == j and (r.get('fdi') == f) and (r.get('week') == t) and (r.get('unit') == 'N') and (r.get('moment_unit') == 'Nmm') and (r.get('moment_origin') == 'baseline crown centroid') and (r.get('timescale') == 'HANDOVER_WEEK') and (r.get('resolution') == 'PER_TOOTH')

def correlation(a, b):
    if len(a) < 3 or np.ptp(a) == 0 or np.ptp(b) == 0:
        return None
    return float(spearmanr(a, b).statistic)

def run():
    start = time.monotonic()
    freeze = read(ROOT / 'FROZEN_PREDICTIONS_R3.json')
    assert sha(ROOT / freeze['prediction_file']) == freeze['sha256']
    pred = read(ROOT / freeze['prediction_file'])
    motion = read(ROOT / 'raw/R2_RELATIVE_POSE.json')
    out = []
    rejections = []
    faults = {}
    keys = sorted({(r['patient'], r['jaw'], r['fdi']) for r in pred['rows']})
    for (p, j, f) in keys:
        q = [r for r in motion if r['patient'] == p and r['jaw'] == j and (r['fdi'] == f) and (r['week'] == 9) and (r['interval'] == 'cumulative')]
        fs = [r for r in pred['rows'] if r['patient'] == p and r['jaw'] == j and (r['fdi'] == f)]
        reason = None
        if len(q) != 2:
            reason = 'no two non-self incisor references (central incisors have only one)'
        elif any((r['planned_angular_deg'] <= 4 for r in q)):
            reason = 'planned angular magnitude <= 4deg conditional radius in at least one reference'
        elif len(fs) != 9 or any((not valid_port(r, p, j, f, r['week']) for r in fs)):
            reason = 'missing/mismatched frozen force rows'
        if reason:
            rejections.append(dict(patient=p, jaw=j, fdi=f, reason=reason))
            continue
        alpha = float(np.mean([r['signed_angular_projection_fraction'] for r in q]))
        lag = 1 - alpha
        norm_ratios = [r['angular_magnitude_ratio'] for r in q]
        row = dict(patient=p, jaw=j, fdi=f, resolution='PER_TOOTH', weeks=9, lag_signed=lag, signed_projection_fraction=alpha, angular_norm_fraction_mean=float(np.mean(norm_ratios)), lag_is_diagnostic_not_measured_biological_efficiency=True, conditional_angular_lag_both_refs=all((r['conditional_angular_lag'] for r in q)), conditional_plan_tracking_rejected_both_refs=all((r['conditional_angular_tracking_rejected'] for r in q)), reference_alpha_spread=float(np.ptp([r['signed_angular_projection_fraction'] for r in q])), sum_force_proxy_N=float(sum((r['force_N'] for r in fs))), sum_force_proxy_note='sum of weekly instantaneous scenario magnitudes; NOT Nweek impulse nor measured exposure', geometry_stiffness_N_per_mm=fs[0]['geometry_stiffness_N_per_mm'], planned_angular_mean_deg=float(np.mean([r['planned_angular_deg'] for r in q])), material_status='PHENOMENOLOGICAL_PETG_TRANSFER_TO_NATURALIGNER', physical_force='UNKNOWN', absolute_motion='UNKNOWN')
        out.append(row)
    probe = next((r for r in pred['rows'] if valid_port(r, r['patient'], r['jaw'], r['fdi'], r['week'])))
    args = [probe['patient'], probe['jaw'], probe['fdi'], probe['week']]
    faults['valid_port_accepts'] = valid_port(probe, *args)
    for (k, v) in [('patient', 'WRONG'), ('fdi', 99), ('week', 10), ('unit', 'mN'), ('moment_origin', 'root CoR'), ('timescale', 'SIMULTANEOUS')]:
        bad = {**probe, k: v}
        faults[f'wrong_{k}_rejected'] = not valid_port(bad, *args)
    by = {}
    for p in ['3485', '6457']:
        rr = [r for r in out if r['patient'] == p]
        y = [r['lag_signed'] for r in rr]
        force = correlation([r['sum_force_proxy_N'] for r in rr], y)
        stiff = correlation([r['geometry_stiffness_N_per_mm'] for r in rr], y)
        plan = correlation([r['planned_angular_mean_deg'] for r in rr], y)
        passed = force is not None and plan is not None and (force >= 0.3) and (force >= plan + 0.1)
        by[p] = dict(teeth=len(rr), rho_force_lag=force, rho_geometric_stiffness_lag=stiff, rho_plan_angular_lag=plan, descriptive_gain_gate='PASS' if passed else 'FAIL', inference_p_value=None, reason='two clustered patient series; no population inference or causal prediction')
    allpass = all((r['descriptive_gain_gate'] == 'PASS' for r in by.values()))
    res = dict(round='R3', claim_type='information_link', external_referent=read(ROOT / 'PREREG_R3.json')['external_referent'], patients=2, eligible_teeth=len(out), candidate_teeth=len(keys), rejected_teeth=len(rejections), rejected_fraction=len(rejections) / len(keys), rejections=rejections, by_patient=by, descriptive_gain_gate='PASS' if allpass else 'FAIL', physical_force_to_movement='UNKNOWN_UNMEASURED_FORCE_MATERIAL_CONTACT_EXPOSURE_AND_BIOLOGICAL_MOBILITY', same_info_control=pred['control'], port_fault_tests=faults, wall_s=time.monotonic() - start, negative_result=not allpass, minimum_new_information='Same patient/tooth/stage fixed-reference repeat IOS plus independently calibrated six-axis wrench history on same finished aligner; independent patient holdout before population prediction')
    write(ROOT / 'raw/R3_RESULTS.json', res)
    write(ROOT / 'raw/R3_TOOTH_LINKS.json', out)
    with (ROOT / 'raw/x19_force_lag.csv').open('w') as o:
        w = csv.DictWriter(o, fieldnames=list(out[0]))
        w.writeheader()
        w.writerows(out)
    txt = f"R3: descriptive X19 force-rank gain {res['descriptive_gain_gate']} on {len(out)}/{len(keys)} teeth, separately per patient. {by}. Physical force-to-movement UNKNOWN; PETG proxy is not Naturaligner material. No population p-value, no force-power-law fit.\nNext: independent repeat reference and matched seated force/time-history measurement; freeze heldout patient outcome before acquisition.\n"
    (ROOT / 'HANDOFF_R3.md').write_text(txt)
    (ROOT / 'HANDOFF.md').write_text(txt)
    state('R3_ADJUDICATED', res['descriptive_gain_gate'], 'Audit source STL/registration, package one-command demo and source comparison with PHENOMENOLOGICAL debts')
    print(res)
    return res
if __name__ == '__main__':
    run()
