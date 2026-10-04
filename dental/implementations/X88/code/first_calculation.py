"""Frozen, source-only calculation. Executed before external study search."""
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from scipy.stats import beta
ROOT = Path(__file__).resolve().parents[1]

def cp(k, n, alpha=0.05):
    return [0.0 if k == 0 else float(beta.ppf(alpha / 2, k, n - k + 1)), 1.0 if k == n else float(beta.ppf(1 - alpha / 2, k + 1, n - k))]

def main():
    prereg = ROOT / 'PREREG_R1.json'
    assert hashlib.sha256(prereg.read_bytes()).hexdigest() == (ROOT / 'PREREG_R1.sha256').read_text().split()[0]
    s = json.loads((ROOT / 'raw/PRIMARY_INPUTS.json').read_text())
    rows = []
    for g in ['normal', 'widening', 'lesion']:
        k = s['CBCT_completed_counts']['fallback_RCT'][g]
        n = k + s['CBCT_completed_counts']['P'][g]
        rows.append(dict(category=g, k=k, n=n, p=k / n, CI95=cp(k, n), simultaneous_CI95=cp(k, n, 0.05 / 3), resolution='POPULATION'))
    result = dict(timestamp_utc=datetime.now(timezone.utc).isoformat(), claim_type='information_link', own_calculation_before_web=True, rows=rows, pooled=dict(k=25, n=86, p=25 / 86, CI95=cp(25, 86)), missing=dict(n=1, endpoint='conversion', fraction=1 / 86), external_transport_prediction=dict(endpoint='unsuccessful haemostasis causing conversion/exclusion from assigned/attempted full pulpotomy', pooled_probability=25 / 86, test='Exact two-sided binomial, alpha .05; protocol-matched only. No CBCT-conditional external calibration without categories.', conditional_predictions={r['category']: r['p'] for r in rows}))
    out = ROOT / 'raw/FIRST_CALCULATION.json'
    assert not out.exists(), 'First calculation is immutable; use demo for replay.'
    out.write_text(json.dumps(result, indent=2) + '\n')
    frozen = dict(timestamp_utc=result['timestamp_utc'], prediction=result['external_transport_prediction'], already_known='Source Table 4 cells and prior R6 negative predictive gate; these are NOT blind predictions of the primary cohort.', not_yet_read='New independent primary study selected by endpoint after this freeze.', sha256_of_first_calculation=hashlib.sha256(out.read_bytes()).hexdigest(), sha256_of_prereg=hashlib.sha256(prereg.read_bytes()).hexdigest())
    p = ROOT / 'FROZEN_PREDICTIONS.json'
    p.write_text(json.dumps(frozen, indent=2) + '\n')
    (ROOT / 'FROZEN_PREDICTIONS.sha256').write_text(hashlib.sha256(p.read_bytes()).hexdigest() + '  FROZEN_PREDICTIONS.json\n')
    print(json.dumps(result, indent=2))
if __name__ == '__main__':
    main()
