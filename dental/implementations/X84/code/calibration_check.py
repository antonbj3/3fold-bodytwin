from patient_geometry import P, write, sha
import json, time, resource
import numpy as np

def main():
    start = time.perf_counter()
    d = json.load(open(P / 'inputs/hattori1994_calibration.json'))
    a = d['published_gain']
    b = d['published_offset_N']
    rows = []
    for (fdi, v) in d['channels'].items():
        for (x, y) in zip(v['applied_force_N'], v['readout_N']):
            pred = (y - b) / a
            err = abs(pred - x) / x
            rows.append(dict(fdi=int(fdi), applied_N=x, raw_readout_N=y, corrected_N=pred, relative_error=err, pass_gate=err <= 0.1, identity_relative_error=abs(y - x) / x, injected100N_rejected=abs((y + 100 - b) / a - x) / x > 0.1, resolution='PER_TOOTH'))
    out = dict(round='R2', claim_type='information_link', outcome='PASS' if all((r['pass_gate'] for r in rows)) else 'FAIL_UNIFORM_GAIN_ACCURACY', rows=rows, retained_rows=len(rows), dropout_fraction=0.0, n_pass=sum((r['pass_gate'] for r in rows)), n_fail=sum((not r['pass_gate'] for r in rows)), raw_identity_pass=sum((r['identity_relative_error'] <= 0.1 for r in rows)), max_relative_error=max((r['relative_error'] for r in rows)), all_injected_errors_rejected=all((r['injected100N_rejected'] for r in rows)), patient_gain_transfer='REJECTED_DIFFERENT_PATIENT_AND_HISTORICAL_SENSOR', external_referent=dict(kind='independent_measurement', locator='https://doi.org/10.2186/jjps.38.835; Table1 p837 and equation p838', compared_quantity='48 externally applied reference loads N on eight tooth channels', refutes_us=True), same_information_control='Direct affine inverse gives identical values; no algorithm claim', source_sha256=sha(P / 'inputs/hattori1994_calibration.json'), cost=dict(wall_s=time.perf_counter() - start, maxrss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, source_transcription='UNKNOWN', physical_patient_acquisition='NOT_RUN'))
    write(P / 'rounds/R2/results.json', out)
    write(P / 'CURRENT_WORK_STATE.json', dict(status='R2_COMPLETE', latest_gate=out['outcome'], next_operation='R3: direct force interval + matched acquisition state; never transfer historical gain', updated_utc=__import__('datetime').datetime.now(__import__('datetime').timezone.utc).isoformat()))
    print(json.dumps({k: out[k] for k in ['outcome', 'n_pass', 'n_fail', 'max_relative_error', 'raw_identity_pass', 'all_injected_errors_rejected']}, indent=2))
if __name__ == '__main__':
    main()
