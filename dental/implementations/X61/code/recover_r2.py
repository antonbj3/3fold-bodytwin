"""Recover already completed numerical R2 after numpy-bool JSON serialization failed.
No simulation, threshold, seeds or metric changed. Full one-command replay does
not use this recovery; it reruns all trials and measures its own cost.
"""
import json, numpy as np
import evaluate as ev
from lab_alarm import R
raw = json.loads((R / 'raw/TRIALS_R2.json').read_text())
cov = json.loads((R / 'raw/POSTERIOR_COVERAGE_R2.json').read_text())
forms = json.loads((R / 'raw/MODEL_LOCALIZATION_R2.json').read_text())
from scipy.stats import beta
null = [x for x in raw if x['kind'] == 'nominal']
n = len(null)
k = sum((x['any_alarm'] for x in null))
tr = dict(null=dict(replicates=n, alarms=k, rate=k / n, binomial_95CI=[float(beta.ppf(0.025, k, n - k + 1)) if k else 0.0, float(beta.ppf(0.975, k + 1, n - k)) if k < n else 1.0], horizon_specimens=72, global_alpha=0.05), faults={}, wall_seconds=None)
for (fault, channel) in ev.FAULTS.items():
    rr = [x for x in raw if x['kind'] == fault]
    ns = [x['first_alarm_at_specimen'].get(channel) for x in rr]
    seen = [v for v in ns if v is not None]
    frac = len(seen) / len(ns)
    tr['faults'][fault] = dict(replicates=len(ns), detected=len(seen), detection_fraction=frac, channel=channel, median_measurements_among_detected=float(np.median(seen)) if seen else None, p10_p90_measurements_among_detected=np.quantile(seen, [0.1, 0.9]).tolist() if seen else None, not_detected_by_72=len(ns) - len(seen), right_observation_link_fraction=frac, causal_localization='Observation link; physical cause vs metrology remains candidate set')
    if fault == 'one_bad_record':
        tr['faults'][fault]['measurement_events_after_fault'] = 1
post = dict(replicates=len(cov), coverage={key: sum((x['coverage'][key] for x in cov)) / len(cov) for key in cov[0]['coverage']}, model_localization_fraction={'nominal': sum((x['nominal'] >= 0.95 for x in forms)) / len(forms), 'separable': sum((x['model_form'] <= 0.05 for x in forms)) / len(forms)}, wall_seconds=None)
ev.trials = lambda *args: (tr, {})
ev.posterior_checks = lambda *args: post
if __name__ == '__main__':
    ev.main()
