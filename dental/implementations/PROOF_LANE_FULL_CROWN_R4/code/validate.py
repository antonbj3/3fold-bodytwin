from score import *

def tile(x, z):
    v = np.array([[x, 0, z], [x + 1, 0, z], [x + 1, 1, z], [x, 1, z]], float)
    return v[np.array([[0, 1, 2], [0, 2, 3]])]

def witnesses():
    ref = np.r_[tile(0, 0), tile(3, 0)]
    a = np.r_[tile(0, 0.25), tile(3, -0.25)]
    b = np.r_[tile(0, -0.25), tile(3, 0.25)]
    ma = metrics(a, ref, 128)
    mb = metrics(b, ref, 128)
    (xy, f) = grid(ref, 0.25)
    ceiling = np.full(len(xy), 0.25)
    ga = ceiling - height(a, xy)
    gb = ceiling - height(b, xy)
    (sa, _) = contact.spatial(xy, f, ga)
    (sb, _) = contact.spatial(xy, f, gb)
    cc = contact.compare(xy, f, ga, gb)
    equal_fields = ['p95_mm', 'rms_mm', 'sampled_max_mm']
    identity = max((abs(ma[k] - mb[k]) for k in equal_fields))
    result = dict(identity_error_mm=identity, summaries=[ma, mb], contact_areas_mm2=[sa['area_mm2'], sb['area_mm2']], component_counts=[sa['count'], sb['count']], centroid_difference_mm=cc['centroid_distance_mm'], symmetric_difference_mm2=cc['symdiff_mm2'], minimal_extension='For this pair, one spatial contact-location coordinate; scalar p95/RMS/max/area/count is insufficient.', source_status='OUR_OWN_FIXTURE; validation witness, not external dental truth', resolution='PER_SURFACE_REGION')
    assert identity == 0 and result['contact_areas_mm2'] == [1.0, 1.0] and (result['component_counts'] == [1, 1]) and (abs(result['centroid_difference_mm'] - 3.0) <= 1e-12) and (result['symmetric_difference_mm2'] == 2.0)
    out = DATA / 'sufficiency.npz'
    np.savez_compressed(out, reference=ref, a=a, b=b, xy=xy, faces=f, ga=ga, gb=gb)
    result['raw_path'] = out
    result['raw_sha256'] = sha(out)
    wrong = dict(result, centroid_difference_mm=30.0)
    result['report_mutation_rejected'] = abs(wrong['centroid_difference_mm'] - cc['centroid_distance_mm']) > 1e-12
    return result

def run():
    st = time.perf_counter()
    controls = {}
    pr = read(ROOT / 'PREREG_B.json')['metrics']['v6']
    ref = tile(0, 0)
    (x, f) = grid(ref, 0.25)
    g = np.full(len(x), 0.05)
    base = contact.compare(x, f, g, g)
    mut = contact.compare(x, f, g + 1, g)
    controls['contact_actual_field_mutation'] = dict(baseline=all(contact.gates(base, pr).values()), injected_rejected=not all(contact.gates(mut, pr).values()), baseline_measurement=base, mutant_measurement=mut)
    inner = tile(0, -0.75)
    bw = wall(ref, inner)
    mw = wall(ref, tile(0, -0.25))
    controls['wall_actual_surface_mutation'] = dict(baseline=bw['continuous_lower_mm'] >= 0.5, injected_rejected=mw['sampled_min_mm'] < 0.5, baseline_measurement=bw, mutant_measurement=mw)
    m = trimesh.creation.box(extents=[1, 1, 1])
    mutm = trimesh.Trimesh(m.vertices, m.faces[:-1], process=False)
    controls['topology_missing_face'] = dict(baseline=bool(m.is_watertight), injected_rejected=not mutm.is_watertight)
    curve = np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0.25], [0, 1, 0]], float)
    correct = curve_error(curve, curve)
    wrong = curve_error(curve + [0, 0, 1.0], curve)
    controls['margin_coordinate_mutation'] = dict(baseline=correct <= 0.025, injected_rejected=wrong > 0.025, baseline_mm=correct, mutant_mm=wrong)
    controls['shape_actual_surface_mutation'] = dict(baseline=metrics(ref, ref, 128)['p95_mm'] <= 0.35, injected_rejected=metrics(ref + [0, 0, 1], ref, 128)['p95_mm'] > 0.35)
    dp = float(distance(tile(0, 0), np.array([[0.5, 0.5, 0.25]]))[0])
    controls['infeasible_wall_inequality'] = dict(baseline=dp + 0.025 < 0.5, injected_rejected=not dp + 1 + 0.025 < 0.5, measured_distance_mm=dp)
    assert all((c['baseline'] and c['injected_rejected'] for c in controls.values()))
    w = witnesses()
    assert w['report_mutation_rejected']
    out = dict(controls=controls, sufficiency=w, all_pass=True, seconds=time.perf_counter() - st, peak_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, scope='Checks geometry and numeric report integrity; no physical calibration')
    dump(ROOT / 'raw/VALIDATION.json', out)
    return out
if __name__ == '__main__':
    run()
