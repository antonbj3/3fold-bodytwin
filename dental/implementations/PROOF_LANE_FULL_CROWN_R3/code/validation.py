from common_r3 import *
sys.path.insert(0, str(V6 / 'code'))
from contact import compare, gates

def run():
    pr = read(ROOT / 'PREREG_A.json')['metrics']['v6']
    xy = np.array([[0, 0], [1, 0], [0, 1], [1, 1], [3, 0], [4, 0], [3, 1], [4, 1]], float)
    faces = np.array([[0, 1, 2], [1, 3, 2], [4, 5, 6], [5, 7, 6]])
    a = np.array([0.0625] * 4 + [0.25] * 4)
    b = np.roll(a, 4)
    aa = compare(xy, faces, a, a)
    bb = compare(xy, faces, b, a)
    summary_a = [aa['predicted']['area_mm2'], aa['predicted']['count'], aa['negative_gap_area_mm2'], float(a.mean()), float(a.min())]
    summary_b = [bb['predicted']['area_mm2'], bb['predicted']['count'], bb['negative_gap_area_mm2'], float(b.mean()), float(b.min())]
    assert summary_a == summary_b
    witness = dict(claim_type='capability', external_referent=dict(kind='our_own_fixture', locator='code/validation.py', compared_quantity='Exact-summary counterexample only; not empirical dental facit', refutes_us=True), resolution='PER_SURFACE_REGION', summary_a=summary_a, summary_b=summary_b, summary_fields=['area_mm2', 'count', 'negative_area_mm2', 'mean_gap_mm', 'minimum_gap_mm'], identity_error=0.0, downstream_centroid_difference_mm=bb['centroid_distance_mm'], downstream_symdiff_mm2=bb['symdiff_mm2'], same_summary_gate_a=gates(aa, pr), same_summary_gate_b=gates(bb, pr), smallest_extension_for_this_pair='one spatial centroid coordinate distinguishing the two locations; this does not prove a full field minimally necessary', time_scale='SIMULTANEOUS')
    controls = []
    for (name, c) in [('native', aa), ('height_plus_1mm', compare(xy, faces, a - 1, a)), ('missing_support', compare(xy, faces, np.full_like(a, np.nan), a)), ('wrong_location', bb)]:
        result = gates(c, pr)
        controls.append(dict(name=name, gates=result, pass_=all(result.values())))
    assert controls[0]['pass_'] and (not any((r['pass_'] for r in controls[1:])))
    corrupt = []
    import copy
    for (name, field, value) in [('support', 'coverage', 0.9), ('pattern', 'symdiff_mm2', 2.0), ('count', 'count_error', 2), ('location', 'centroid_distance_mm', 1.0), ('interference', 'negative_gap_area_mm2', 0.2)]:
        c = copy.deepcopy(aa)
        c[field] = value
        gs = gates(c, pr)
        assert not gs[name]
        corrupt.append(dict(injected_field=field, injected_value=value, expected_failed_gate=name, observed_gates=gs))
    dump(ROOT / 'raw/SUFFICIENCY.json', witness)
    dump(ROOT / 'raw/NEGATIVE_CONTROLS.json', dict(geometry_controls=controls, gate_field_controls=corrupt, scope='Geometry perturbations execute original spatial evaluator; field injections additionally exercise each threshold. These do not validate source labels or clinical function.'))
    return witness
if __name__ == '__main__':
    run()
