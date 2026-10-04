"""Reviewer fault probes using the delivered scorer and its original software fixture.
Run with /usr/bin/python3 -s review_control.py from a reviewed X1c copy.
No physical observations, new model or refrozen scientific predictions.
"""
import sys, copy, math, json, datetime
from pathlib import Path
R = Path(__file__).resolve().parent
sys.path.insert(0, str(R / 'code'))
import compare_lab as cl

def main():
    rr = []
    now = datetime.datetime.now(datetime.timezone.utc)
    for (prot, p) in cl.PROTOCOLS.items():
        for d in ['D1', 'M2', 'M1']:
            base = {'D1': 1000.0, 'M2': 3000.0, 'M1': 1000 ** (1 - cl.W) * 3000 ** cl.W}[d]
            for i in range(12):
                r = {k: '' for k in cl.REQUIRED}
                r.update(specimen_id=f'SYNTHETIC_{prot}_{d}_{i}', protocol_id=prot, design_id=d, force_unit='N', fracture_force_N=str(base * (0.9 + 0.2 * i / 11)), die_E_MPa=str(p['E']), indenter_diameter_mm=str(p['diameter']), crosshead_mm_min=str(p['speed']), load_angle_deg=str(p['angle']), material_batch='SYNTHETIC', crown_product=p['crown'], cement_product=p['cement'], cement_batch='SYNTHETIC', die_product=p['die'], die_batch='SYNTHETIC', surface_protocol=p['surface'], storage_days='7', interlayer_spec='SYNTHETIC_fixed_spec', manufacturing_spec_sha256='1' * 64, failure_mode='crown_tensile_fracture', fracture_origin='intaglio_tensile_zone', force_trace_sha256='2' * 64, measured_utc=(now + datetime.timedelta(seconds=60 if d == 'M1' else -60)).isoformat(), status='measured')
                rr.append(r)
    cal = cl.calibration([r for r in rr if r['design_id'] != 'M1'])
    good = cl.compare(rr, cal)
    control = {p['protocol']: math.exp(math.log(p['group_summaries']['D1']['mean_N']) + cl.W * (math.log(p['group_summaries']['M2']['mean_N']) - math.log(p['group_summaries']['D1']['mean_N']))) for p in good['protocols']}
    bad = copy.deepcopy(rr)
    for r in bad:
        if r['design_id'] == 'M1':
            r['fracture_force_N'] = str(float(r['fracture_force_N']) / 1.6)
    wrong = cl.compare(bad, cal)
    mode = copy.deepcopy(rr)
    for r in mode:
        if r['protocol_id'] == 'C0_PROTT' and r['design_id'] == 'M1':
            r['failure_mode'] = 'die_fracture'
            r['fracture_origin'] = 'die'
    mech = cl.compare(mode, cal)
    identity = max((abs(a['residual_log'] - b['residual_log']) for (a, b) in zip(good['protocols'], mech['protocols'])))
    original = cl.compare
    try:
        cl.compare = lambda *a, **k: good
        blind_mutation_detected = not all((p['force_verdict'] == 'REJECTED' for p in cl.compare(bad, cal)['protocols']))
    finally:
        cl.compare = original
    out = dict(scope='existing software fixtures only', conventional_prediction_N=control, baseline_verdicts=[p['force_verdict'] for p in good['protocols']], corrupted_control_multiplier=1.6, wrong_control_verdicts=[p['force_verdict'] for p in wrong['protocols']], actual_scorer_called=True, always_support_scorer_mutation_detected=blind_mutation_detected, sufficiency=dict(summary='all force values and log residuals', identity_error=identity, exact_identity=identity == 0, baseline=good['protocols'][0]['mechanism_gates'], other=mech['protocols'][0]['mechanism_gates'], minimum_extension='observed failure mode and origin alongside force; these do not identify a constitutive mechanism by themselves'), conclusion_changed=False, physical_validation='UNKNOWN')
    assert all((p['force_verdict'] == 'REJECTED' for p in wrong['protocols'])) and identity == 0 and blind_mutation_detected
    (R / 'raw/REVIEW_CONTROL.json').write_text(json.dumps(out, indent=2) + '\n')
    print(json.dumps(out))
if __name__ == '__main__':
    main()
