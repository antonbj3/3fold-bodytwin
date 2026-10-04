import sys, time, json, hashlib, resource
from pathlib import Path
from datetime import datetime, timezone
import numpy as np
from contact_model import make_system, solve, insertion
R = Path(__file__).resolve().parents[1]
DENT = R.parents[1]
sys.path.insert(0, str(DENT / 'results/LANE_X19_ALIGNER_FORCE/code'))
from model import response

def put(n, x):
    (R / n).write_text(json.dumps(x, indent=2, allow_nan=False) + '\n')

def main():
    start = time.perf_counter()
    pr = json.loads((R / 'PREREG_R1.json').read_text())
    assert hashlib.sha256((R / 'PREREG_R1.json').read_bytes()).hexdigest() == json.loads((R / 'PREREG_R1.sha256.json').read_text())['sha256']
    arches = json.loads((R / 'inputs/ARCHES.json').read_text())
    regions = json.loads((R / 'inputs/PARK_REGIONS.json').read_text())
    rows = []
    history = []
    fails = []
    pathmax = 0.0
    controls = []
    relax = []
    for arch in arches:
        for fdi in [11, 13, 15]:
            try:
                system = make_system(arch, fdi, regions)
            except ValueError as err:
                fails.append(dict(case=arch['case'], active_fdi=fdi, reason=str(err)))
                continue
            out = dict(case=arch['case'], active_fdi=fdi, resolution_level='PER_TOOTH', claim_type='information_link')
            for kind in ['regional', 'nominal_contact', 'tooth_type_only']:
                s = system if kind == 'regional' else make_system(arch, fdi, regions, kind)
                out[kind] = solve(s)
            out['practice'] = response(arch, fdi, 0.2, 2746.0, 0.75)
            out['regional_iqr_scenarios'] = [solve(make_system(arch, fdi, regions, quantile=q)) for q in ['q1_mm', 'q3_mm']]
            out['iqr_note'] = 'Changing h and g together at marginal quartile endpoints is a sensitivity scenario, not a joint interval guarantee'
            for load in [0.0, 30.0]:
                finals = []
                for path in ['front_first', 'back_first', 'simultaneous']:
                    tr = insertion(system, path, load)
                    history.extend((dict(case=arch['case'], **x) for x in tr))
                    finals.append(tr[-1]['wrenches'][str(fdi)])
                pathmax = max(pathmax, float(np.max(np.ptp(np.array(finals), axis=0))))
                out['held_30N' if load else 'released_0N'] = solve(system, seating_N=load)
            control = solve(system, method='primal')
            parity = max((float(np.max(np.abs(np.array(control['wrenches'][k]) - out['regional']['wrenches'][k]))) for k in control['wrenches']))
            controls.append(dict(case=arch['case'], fdi=fdi, max_wrench_difference=parity, control_residuals=control['residuals']))
            for (i, factor, label) in [(122, 0.7962, 'tension_37C_3h'), (101, 1 / 24, 'wet_bending_37C_24h')]:
                for load in [0.0, 30.0]:
                    a = solve(system, seating_N=load, relaxation_factor=factor)
                    base = out['held_30N' if load else 'released_0N']
                    naive = np.array(base['wrenches'][str(fdi)]) * factor
                    relax.append(dict(case=arch['case'], fdi=fdi, coupon_record=i, load_case=label, seating_N=load, time_scale='HANDOVER', resolution_level='PER_TOOTH', coupon_resolution_level='PHENOMENOLOGICAL', transfer_status='UNVALIDATED_CONDITIONAL_COUPON_MODULUS_SCALING', retained_factor=factor, reequilibrated_wrench=a['wrenches'][str(fdi)], naive_scaled_wrench=naive.tolist(), wrench_difference=float(np.linalg.norm(np.array(a['wrenches'][str(fdi)]) - naive)), contacts_before=[(r['fdi'], r['region']) for r in base['regional'] if r['contact']], contacts_after=[(r['fdi'], r['region']) for r in a['regional'] if r['contact']]))
            rows.append(out)
    pred = dict(rows=rows, failures=fails, controls=controls, insertion_endpoint_path_difference=pathmax, relaxation=relax)
    put('raw/PREDICTIONS_R1.json', pred)
    put('raw/INSERTION_TRACE.json', history)
    frozen = dict(frozen_utc=datetime.now(timezone.utc).isoformat(), prereg_sha256=hashlib.sha256((R / 'PREREG_R1.json').read_bytes()).hexdigest(), predictions_sha256=hashlib.sha256((R / 'raw/PREDICTIONS_R1.json').read_bytes()).hexdigest(), source_visibility='Published endpoints seen during source read; no fit and no post-score changes; not blinded', predictions=pred)
    if (R / 'FROZEN_PREDICTIONS_R1.json').exists():
        old = json.loads((R / 'FROZEN_PREDICTIONS_R1.json').read_text())
        if old['predictions_sha256'] != frozen['predictions_sha256']:
            raise RuntimeError('Frozen R1 prediction drift; preserve output and investigate')
    else:
        put('FROZEN_PREDICTIONS_R1.json', frozen)
    put('raw/COST_R1.json', dict(compute_s=time.perf_counter() - start, peak_RSS_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, source_arches=len(arches), tooth_cases=len(rows), trace_steps=len(history), fit_s=0, threads_max=4, heavier_than_4GB=False, discovery_tokens=None, questions=0, inherited_preparation_s=6.73))
    print('R1 predictions frozen:', len(rows), 'tooth cases;', len(history), 'insertion steps; peak RSS', resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, 'MiB')
if __name__ == '__main__':
    main()
