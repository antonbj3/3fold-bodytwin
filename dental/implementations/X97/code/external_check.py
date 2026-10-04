"""Compare inherited external facts to the freshly checked primary numeric table."""
from pathlib import Path
import json, hashlib
import numpy as np
R = Path(__file__).resolve().parents[1]

def run():
    p = R / 'raw/EXTERNAL_SOURCE_FACTS.json'
    src = json.loads(p.read_text())
    old = json.loads((R / 'inputs/X94_EXTERNAL_VALIDATION.json').read_text())['source_quantile_table']
    checks = {}
    for (key, label) in [('R', 'right_pp'), ('posterior', 'posterior_pp')]:
        observed = np.array(old['values_share'][key]) * 100
        truth = np.array(src[label])
        err = float(abs(observed - truth).max())
        bad = observed.copy()
        bad[4] += 1.0
        baderror = float(abs(bad - truth).max())
        checks[key] = dict(max_error_pp=err, match=err <= 0.05, actual_injected_1pp_rejected=baderror > 0.05, resolution='POPULATION')
    out = dict(external_referent=src, checks=checks, tolerance_pp=0.05, compared_quantity='X94 retained primary table transcription, independently rechecked from publisher TABLE1; no patient force validation', same_specimen_force_rows_accepted=0, source_table_rows_rejected_for_individual_force=18, rejection_fraction=1.0, rejection_reason='Wrong quantity/resolution; no absolute force or specimen binding', source_sha256=hashlib.sha256(p.read_bytes()).hexdigest())
    (R / 'raw/EXTERNAL_FACIT_CHECK.json').write_text(json.dumps(out, indent=2) + '\n')
    assert all((c['match'] and c['actual_injected_1pp_rejected'] for c in checks.values()))
    print(json.dumps({'external_table_match': True, 'individual_force_calibration_rows_accepted': 0}))
if __name__ == '__main__':
    run()
