from common import *
from local_pairs import metrics, distance, sample
from witness_gate import corrected_controls

def tile(x, z):
    v = np.array([[x, 0, z], [x + 1, 0, z], [x + 1, 1, z], [x, 1, z]], float)
    return v[np.array([[0, 1, 2], [0, 2, 3]])]

def run():
    ref = tile(0, 0)
    a = tile(0, 0.25)
    b = tile(0, -0.25)
    ma = metrics(a, ref, 128)
    mb = metrics(b, ref, 128)
    ident = max((abs(ma[k] - mb[k]) for k in ['p95_mm', 'rms_equal_surface_mm', 'sampled_max_mm']))
    gap_a = 0.25 - a[:, :, 2].mean()
    gap_b = 0.25 - b[:, :, 2].mean()
    area_a = float(0 <= gap_a <= 0.1)
    area_b = float(0 <= gap_b <= 0.1)
    ref2 = np.r_[tile(0, 0), tile(3, 0)]
    a2 = np.r_[tile(0, 0.25), tile(3, -0.25)]
    b2 = np.r_[tile(0, -0.25), tile(3, 0.25)]
    aa = metrics(a2, ref2, 128)
    bb = metrics(b2, ref2, 128)
    controls = dict(identity_pass=float(distance(ref, sample(ref, 128)).max()) <= 1e-10, translated_identity_rejected=float(distance(ref + np.array([0, 0, 1.0]), sample(ref, 128)).max()) > 1e-10)
    assert ident == 0 and area_a - area_b >= 0.5 and all(controls.values())
    result = dict(claim_type='capability', resolution='PER_SURFACE_REGION', source='constructed counterexamples, not clinical reference', external_referent=dict(kind='closed_form', locator='Euclidean distance between parallel congruent unit squares; code/sufficiency.py', compared_quantity='distance and intersection area under an explicit rigid planar geometry', refutes_us=True), sign_witness=dict(state_a=ma, state_b=mb, identity_error_mm=ident, gaps_mm=[gap_a, gap_b], contact_areas_mm2=[area_a, area_b], downstream_difference_mm2=area_a - area_b, minimal_extension_for_this_pair='one signed gap to antagonist'), location_witness=dict(state_a=aa, state_b=bb, identity_error_mm=max((abs(aa[k] - bb[k]) for k in ['p95_mm', 'rms_equal_surface_mm', 'sampled_max_mm'])), contact_area_both_mm2=1, component_count_both=1, contact_centroid_distance_mm=3.0, symmetric_difference_mm2=2.0, minimal_extension_for_this_pair='one contact location coordinate; scalar area alone still fails'), controls=controls, interpretation='p95, RMS and Hausdorff are not sufficient for contact. A region vector is not claimed sufficient for all clinical outcomes.', rigorous_scope='binary exact coordinates; analytic squares; floating triangle calculation agrees exactly. No general mesh interval certificate.')
    np.savez_compressed(DATA / 'sufficiency.npz', ref=ref, a=a, b=b, ref2=ref2, a2=a2, b2=b2)
    controls.update(corrected_controls(result, DATA / 'sufficiency.npz'))
    assert all(controls.values())
    dump(ROOT / 'raw/SUFFICIENCY.json', result)
    return result
if __name__ == '__main__':
    run()
