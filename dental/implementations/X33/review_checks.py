"""Validate an actual R3 result against separately executed support LPs."""
import copy
import math

def validate_support(row, independent_peak, threshold=5.5, tolerance=1e-07):
    value = row['worst_peak_pulp_C']
    return isinstance(value, (int, float)) and math.isfinite(value) and (abs(value - independent_peak) < tolerance) and (row['coarse_model_below_threshold'] is (independent_peak < threshold)) and (row['LP_pass'] is True)

def validate_finer(row, threshold=5.5):
    return math.isfinite(row['maximum_C']) and row['below_threshold'] is (row['maximum_C'] < threshold)

def validate_interval_label(candidate, lower, upper, threshold):
    if not all((math.isfinite(x) for x in [lower, upper, threshold])) or lower > upper:
        return False
    independent = 'BELOW' if upper < threshold else 'ABOVE' if lower > threshold else 'ABSTAIN'
    return candidate == independent

def interval_probes(actual, lower, upper, threshold):
    positive = validate_interval_label(actual, lower, upper, threshold)
    poisoned = 'BELOW' if actual != 'BELOW' else 'ABOVE'
    return {'positive_interval_label_accepted': positive, 'claiming_interval_below_rejected': not validate_interval_label(poisoned, lower, upper, threshold)}

def result_probes(rows, independent_peaks, fine):
    assert len(rows) == len(independent_peaks) and rows and fine
    positive = all((validate_support(r, p) for (r, p) in zip(rows, independent_peaks)))
    positive = positive and all((validate_finer(r) for r in fine))
    poisoned = copy.deepcopy(rows[-1])
    poisoned['worst_peak_pulp_C'] += 1.0
    numerical = not validate_support(poisoned, independent_peaks[-1])
    mislabeled = copy.deepcopy(rows[-1])
    mislabeled['coarse_model_below_threshold'] = not mislabeled['coarse_model_below_threshold']
    native_label = not validate_support(mislabeled, independent_peaks[-1])
    exceeding = next((r for r in fine if r['maximum_C'] >= 5.5))
    false_fine = copy.deepcopy(exceeding)
    false_fine['below_threshold'] = True
    finer_label = not validate_finer(false_fine)
    return {'positive_rows_accepted': positive, 'support_plus_1C_rejected': numerical, 'native_wrong_label_rejected': native_label, 'false_below_label_on_finer_history_rejected': finer_label}
