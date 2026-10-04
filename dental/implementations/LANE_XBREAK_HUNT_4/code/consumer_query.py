"""Source-cohort pathway query. Never a treatment recommendation."""
import argparse, json, hashlib
from pathlib import Path

def main():
    p = argparse.ArgumentParser()
    p.add_argument('state')
    p.add_argument('--cbct', required=True, choices=['normal', 'widening', 'lesion', 'missing'])
    p.add_argument('--out')
    a = p.parse_args()
    state = json.loads(Path(a.state).read_text())
    category = a.cbct
    risk = next((r for r in state['observed_categories'] if r['CBCT'] == category), None)
    timing = next((r for r in state['timing_intervals'] if r['query'] == category), None)
    response = {'observation': {'CBCT': category}, 'consumer': 'Research protocol and scheduling specification', 'source_cohort_fallback_count': None if risk is None else {'k': risk['fallback_RCT'], 'n': risk['n']}, 'source_cohort_fallback_fraction': None if risk is None else risk['sample_risk'], 'conditional_IID_simultaneous95_interval': None if risk is None else risk['IID_simultaneous95_interval'], 'nominal_moment_compatible_category_mean_first_visit_min': None if timing is None else timing['interval_min'], 'external_calibrated_forecast': 'UNKNOWN', 'treatment_choice': 'UNKNOWN', 'full_treatment_time': 'UNKNOWN', 'review_state': state['review_state'], 'assumptions': state['assumptions'], 'empirical_validity': 'Published trial sample only. Time intervals use nominal moments and column-N timing closure; not patient-time or population confidence intervals.', 'next_measurement': state['next_measurement']}
    text = json.dumps(response, indent=2) + '\n'
    if a.out:
        Path(a.out).write_text(text)
    print(text)
if __name__ == '__main__':
    main()
