"""A real measurement-state consumer and future lab artifact, not a crown FE law."""
import datetime as dt, hashlib, json, math, struct
from pathlib import Path
import ports
from sources import load
from make_demo import freeze_predictions, write, sufficiency
ROOT = Path(__file__).resolve().parent

def coef(r1=6.05, r2=0.8, r3=7.5, b=1.0, nu=0.3):
    X = (1 + nu) * math.log((r2 / r3) ** 2) + (1 - nu) / 2 * (r2 / r3) ** 2
    Y = (1 + nu) * (1 + math.log((r1 / r3) ** 2)) + (1 - nu) * (r1 / r3) ** 2
    return -0.2387 * (X - Y) / b ** 2

def control_coef(r1=6.05, r2=0.8, r3=7.5, b=1.0, nu=0.3):
    diff = (1 + nu) * (2 * math.log(r2 / r1) - 1) + (1 - nu) * (0.5 * r2 * r2 - r1 * r1) / (r3 * r3)
    return -0.2387 * diff / b ** 2

def coefficient_enclosure():
    import mpmath as mp
    mp.iv.dps = 40
    nu = mp.iv.mpf('0.3')
    r1 = mp.iv.mpf('6.05')
    r2 = mp.iv.mpf('0.8')
    r3 = mp.iv.mpf('7.5')
    diff = (1 + nu) * (2 * mp.iv.ln(r2 / r1) - 1) + (1 - nu) * (mp.iv.mpf('0.5') * r2 * r2 - r1 * r1) / (r3 * r3)
    k = -mp.iv.mpf('0.2387') * diff
    return [math.nextafter(float(k.a), -math.inf), math.nextafter(float(k.b), math.inf)]

def disc_sdf(x, y, z, r=7.5, b=1.0):
    dr = math.hypot(x, y) - r
    dz = abs(z) - b / 2
    return math.hypot(max(dr, 0), max(dz, 0)) + min(max(dr, dz), 0)

def export_disc(path, n=128, r=7.5, b=1.0):
    verts = [(r * math.cos(2 * math.pi * i / n), r * math.sin(2 * math.pi * i / n)) for i in range(n)]
    facets = []
    for i in range(n):
        j = (i + 1) % n
        (ai, aj) = (verts[i], verts[j])
        lowi = (*ai, -b / 2)
        lowj = (*aj, -b / 2)
        highi = (*ai, b / 2)
        highj = (*aj, b / 2)
        facets.extend([((0, 0, -b / 2), lowj, lowi), ((0, 0, b / 2), highi, highj), (lowi, lowj, highj), (lowi, highj, highi)])
    payload = b'X45 source-protocol coupon, mm; conditional lab predictions'.ljust(80, b' ') + struct.pack('<I', len(facets))
    for tri in facets:
        payload += struct.pack('<12fH', 0.0, 0.0, 0.0, *(v for p in tri for v in p), 0)
    path.write_bytes(payload)
    vol = n * r * r * math.sin(2 * math.pi / n) * b / 2
    true = math.pi * r * r * b
    sag = r * (1 - math.cos(math.pi / n))
    return {'path': str(path), 'sha256': hashlib.sha256(payload).hexdigest(), 'triangles': len(facets), 'analytic_volume_mm3': true, 'polygon_volume_mm3': vol, 'relative_polygon_volume_deficit': 1 - vol / true, 'radial_chord_error_mm': sag, 'STL_max_radial_error_enclosure_mm': sag + 5e-07, 'resolution': 'PER_SURFACE_REGION', 'units': 'mm', 'geometry_kind': 'our_own_fixture exported to reproduce source coupon dimensions'}

def main():
    prereg = ROOT / 'PREREG_COUPON_CONSUMER.json'
    receipt = ROOT / 'PREREG_COUPON_CONSUMER_FREEZE.json'
    if not receipt.exists():
        write(receipt.name, {'frozen_utc': dt.datetime.now(dt.timezone.utc).isoformat(), 'prereg_sha256': hashlib.sha256(prereg.read_bytes()).hexdigest(), 'before': 'any coupon numerical calculation or prediction'})
    if json.loads(receipt.read_text())['prereg_sha256'] != hashlib.sha256(prereg.read_bytes()).hexdigest():
        raise ValueError('coupon prereg drift')
    version = 'K34_LTD_COUPON_CONSUMER_R2'
    prereg_wrapper = ROOT / ('PREREG_' + version + '.json')
    if not prereg_wrapper.exists():
        write(prereg_wrapper.name, {**json.loads(prereg.read_text()), 'construction': version, 'parent': 'K34_LTD_COUPON_CONSUMER', 'frozen_utc': dt.datetime.now(dt.timezone.utc).isoformat(), 'delta': 'Enclose stress-coefficient arithmetic too; preserve R1 load predictions and R1 parameter-only intervals. No empirical gate or physical parameter changed.', 'data_visibility': 'Published fit summaries known; no physical future observations read here.', 'freeze_receipt': json.loads(receipt.read_text())})
        if (ROOT / 'raw/COUPON_CONSUMER.json').exists():
            (ROOT / 'raw/COUPON_CONSUMER_R1.json').write_bytes((ROOT / 'raw/COUPON_CONSUMER.json').read_bytes())
    data = load()
    sufficiency(data)
    k = coef()
    kc = control_coef()
    ki = coefficient_enclosure()
    records = {}
    queries = {}
    for (grade, states) in data['LTD'].items():
        records[grade] = {}
        queries[grade] = {}
        for (state, v) in states.items():
            q = ports.material_query(v, 0.05)
            loadN = q['quantile_MPa'] / k
            bounds = q['printed_rounding_enclosure_MPa']
            queries[grade][state] = {'Q05_load_N': loadN, 'printed_rounding_load_enclosure_N': [math.nextafter(bounds[0] / ki[1], -math.inf), math.nextafter(bounds[1] / ki[0], math.inf)], 'tail_probability_at_load': -math.expm1(-(loadN * k / v['sigma0_MPa']) ** v['m']), 'sampling_CI': 'UNKNOWN', 'resolution': 'POPULATION', 'actual_new_bench_measurement': 'NOT_RUN'}
            records[grade][state] = {**v, 'law': 'coupon_Weibull_strength_after_declared_exposure', 'unit_sigma0': 'MPa', 'resolution': 'POPULATION', 'valid_geometry': 'source15mm_disc_1mm_thick', 'exposure_state': state, 'hours_at122C': 8 if state == 'aged' else 0, 'cycles_at50N': 240000, 'frequency_Hz': 1.1, 'new_lot_transfer': 'UNKNOWN', 'clinical_years': 'UNKNOWN', 'law_leaf_status': 'CONSTITUTIVE_CLOSURE', 'parameter_leaf_status': 'EXTERNALLY_MEASURED'}
    freeze_predictions(version, queries, records)
    write('MATERIAL_STATE_REGISTRY.json', records)
    exp = ROOT / 'exports'
    exp.mkdir(exist_ok=True)
    mesh = export_disc(exp / 'source_coupon_mm.stl')
    probes = [((0, 0, 0), -0.5), ((7.5, 0, 0), 0), ((0, 0, 0.5), 0), ((8, 0, 0), 0.5), ((8, 0, 1), math.sqrt(0.5))]
    sdf_error = max((abs(disc_sdf(*p) - v) for (p, v) in probes))
    controls = [{'name': 'independent biaxial stress/load formula', 'relative_error': abs(k - kc) / kc, 'pass': abs(k - kc) / kc <= 1e-12}, {'name': 'coupon inverse CDF0.05', 'pass': all((abs(v['tail_probability_at_load'] - 0.05) <= 1e-12 for st in queries.values() for v in st.values()))}, {'name': 'SDF exact probes', 'max_error_mm': sdf_error, 'pass': sdf_error <= 1e-12}, {'name': 'injected stress coefficient sign', 'detected': -k <= 0}, {'name': 'injected thickness-squared omission at b2', 'detected': abs(k - coef(b=2)) / coef(b=2) > 1e-12}]
    write('raw/COUPON_CONSUMER.json', {'claim_type': 'capability', 'stress_per_load_MPa_N': k, 'stress_coefficient_arithmetic_enclosure_MPa_N': ki, 'stress_coefficient_resolution': 'PER_SURFACE_REGION', 'poisson_ratio_leaf': 'CONSTITUTIVE_CLOSURE', 'queries': queries, 'export': mesh, 'controls': controls, 'external_referent': json.loads(prereg.read_text())['external_referent'], 'outcome': 'FROZEN_FUTURE_COUPON_PREDICTIONS', 'physical_validation': 'NOT_RUN', 'uncertainty': 'printed parameter rounding only; statistical parameter and fixture/dimension uncertainty UNKNOWN'})
    print('Coupon consumer:8 conditional load predictions, geometry export, controls', all((c.get('pass', c.get('detected')) for c in controls)))
if __name__ == '__main__':
    main()
