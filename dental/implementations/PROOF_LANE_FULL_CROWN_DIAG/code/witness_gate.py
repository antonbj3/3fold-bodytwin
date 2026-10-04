"""Check existing sufficiency witnesses against their saved triangle states."""
import copy
import numpy as np

def observed(data_file):
    with np.load(data_file, allow_pickle=False) as z:
        states = {k: z[k] for k in z.files}

    def band(t):
        area = np.linalg.norm(np.cross(t[:, 1] - t[:, 0], t[:, 2] - t[:, 0]), axis=1) / 2
        gap = 0.25 - t.mean(axis=1)[:, 2]
        mask = (gap >= 0) & (gap <= 0.1)
        total = float(area[mask].sum())
        centroid = (t.mean(axis=1)[mask, :2] * area[mask, None]).sum(axis=0) / total if total else None
        return (total, centroid, mask, area)
    (a, b) = (band(states['a']), band(states['b']))
    (a2, b2) = (band(states['a2']), band(states['b2']))
    assert np.array_equal(states['a2'][:, :, :2], states['b2'][:, :, :2])
    return {'sign_areas': [a[0], b[0]], 'sign_difference': a[0] - b[0], 'sign_gaps': [0.25 - float(states[k][:, :, 2].mean()) for k in ('a', 'b')], 'location_areas': [a2[0], b2[0]], 'centroid_distance': float(np.linalg.norm(a2[1] - b2[1])), 'symmetric_difference': float(a2[3][a2[2] != b2[2]].sum())}

def report_matches(report, data_file):
    obs = observed(data_file)
    (s, l) = (report['sign_witness'], report['location_witness'])
    keys = ('p95_mm', 'rms_equal_surface_mm', 'sampled_max_mm')
    identity = max((abs(s['state_a'][k] - s['state_b'][k]) for k in keys))
    location_identity = max((abs(l['state_a'][k] - l['state_b'][k]) for k in keys))
    return identity == s['identity_error_mm'] == 0 and location_identity == l['identity_error_mm'] == 0 and (s['contact_areas_mm2'] == obs['sign_areas']) and (s['gaps_mm'] == obs['sign_gaps']) and (s['downstream_difference_mm2'] == obs['sign_difference'] > 0) and (obs['location_areas'] == [l['contact_area_both_mm2']] * 2) and (l['contact_centroid_distance_mm'] == obs['centroid_distance'] > 0.5) and (l['symmetric_difference_mm2'] == obs['symmetric_difference'])

def corrected_controls(report, data_file):
    assert report_matches(report, data_file)
    bad_identity = copy.deepcopy(report)
    bad_identity['sign_witness']['identity_error_mm'] = 0.001
    bad_location = copy.deepcopy(report)
    bad_location['location_witness']['contact_centroid_distance_mm'] = 30.0
    bad_gap = copy.deepcopy(report)
    bad_gap['sign_witness']['gaps_mm'][0] += 1.0
    return {'nonzero_identity_mutation_rejected': not report_matches(bad_identity, data_file), 'swapped_contact_location_rejected': not report_matches(bad_location, data_file), 'wrong_gap_value_rejected': not report_matches(bad_gap, data_file)}
