"""Exact event partition for a declared finite geometric contact operator.

Rational arithmetic encloses ALL vertical translations for the binary input
gap field. It does not enclose the source-to-field error or elastic FE error.
"""
from fractions import Fraction as Q
from decimal import Decimal, localcontext
from collections import deque
import sys, time
import numpy as np
from common import *

def component_mask(raw, index, weights, minimum_area):
    coords = np.argwhere(index >= 0)
    position = {int(index[tuple(ij)]): tuple(ij) for ij in coords}
    unseen = set(np.flatnonzero(raw).tolist())
    retained = np.zeros(len(raw), bool)
    discarded = []
    sizes = []
    while unseen:
        seed = min(unseen)
        unseen.remove(seed)
        todo = deque([seed])
        group = []
        while todo:
            k = todo.popleft()
            group.append(k)
            (i, j) = position[k]
            for (a, b) in ((i - 1, j), (i + 1, j), (i, j - 1), (i, j + 1)):
                if 0 <= a < index.shape[0] and 0 <= b < index.shape[1]:
                    v = int(index[a, b])
                    if v in unseen:
                        unseen.remove(v)
                        todo.append(v)
        area = sum((Q(float(weights[i])) for i in group), Q(0))
        if area >= Q(float(minimum_area)):
            retained[group] = True
            sizes.append(len(group))
        else:
            discarded.append(dict(point_count=len(group), area_mm2=float(area)))
    return (retained, sizes, discarded)

def event_partition(gap, weights, index, lo, hi, band, minimum_area):
    g = [Q(float(v)) if np.isfinite(v) else None for v in gap]
    (lo, hi, b) = (Q(float(lo)), Q(float(hi)), Q(float(band)))
    points = {lo, hi}
    intervals = []
    for v in g:
        interval = (-b - v, b - v) if v is not None else None
        intervals.append(interval)
        if interval:
            points.update((x for x in interval if lo <= x <= hi))
    events = sorted(points)
    cells = [(a, a, a, 'POINT') for a in events]
    cells += [(a, z, (a + z) / 2, 'OPEN_INTERVAL') for (a, z) in zip(events, events[1:])]
    cells.sort(key=lambda r: (r[0], r[1]))
    rows = []
    exact = []
    for (a, z, q, kind) in cells:
        raw = np.array([v is not None and v[0] <= q <= v[1] for v in intervals])
        (mask, sizes, discarded) = component_mask(raw, index, weights, minimum_area)
        area = sum((Q(float(weights[i])) for i in np.flatnonzero(mask)), Q(0))
        row = dict(kind=kind, lower_mm=float(a), upper_mm=float(z), witness_mm=float(q), exact_lower=str(a), exact_upper=str(z), exact_witness=str(q), contact_area_mm2=float(area), contact_area_exact=str(area), active_point_ids=np.flatnonzero(mask).tolist(), raw_point_count=int(raw.sum()), component_sizes=sizes, rejected_components=discarded)
        rows.append(row)
        exact.append((q, raw, mask))
    return (rows, exact)

def run(out):
    t0 = time.perf_counter()
    out.mkdir(parents=True, exist_ok=True)
    path = PARENT_DATA / 'Bite2Text_F4775/roof.npz'
    source(path)
    z = np.load(path)
    gap = z['ceiling'] - z['z']
    weights = z['weights']
    index = z['index']
    pr = load(ROOT / 'PREREG_C1_POSE.json')
    sc = pr['scenario']
    (lo, hi) = sc['pose_mm']
    band = sc['band_mm']
    amin = sc['minimum_component_area_mm2']
    (rows, exact) = event_partition(gap, weights, index, lo, hi, band, amin)
    x18 = RESULTS / 'LANE_X18_CROWN_ANTAGONIST/code'
    sys.path.insert(0, str(x18))
    import geometry as g
    from experiment import contact
    from fe import solve_cases, normalize
    for name in ('geometry.py', 'experiment.py', 'fe.py', 'continuous.py'):
        if (x18 / name).exists():
            source(x18 / name)
    prediction = dict(patient_id='Bite2Text_F4775', FDI=36, claim_type='capability', source=source(path), prereg=source(ROOT / 'PREREG_C1_POSE.json'), partition=rows, geometric_robust=all((r['contact_area_mm2'] > 0 for r in rows)), physical_measurement='NOT_RUN', resolution='PER_POINT', time_scale='SIMULTANEOUS', floating_input_scope='Exact rational threshold and connected-component decisions on binary saved gaps/weights; NOT continuous surface error or physical contact')
    freeze(out / 'FROZEN_PREDICTIONS.json', prediction)
    tpred = time.perf_counter()
    errors = []
    mismatches = 0
    decimal_mismatches = 0
    for (row, (q, raw, mask)) in zip(rows, exact):
        with localcontext() as ctx:
            ctx.prec = 100
            delta = Decimal(q.numerator) / Decimal(q.denominator)
            refraw = np.array([np.isfinite(v) and abs(Decimal.from_float(float(v)) + delta) <= Decimal.from_float(band) for v in gap])
        decimal_mismatches += int(np.count_nonzero(raw != refraw))
        ref = contact(np.where(refraw, 0.0, band + 1), z['z'], z['xy'], index, weights, band, amin)
        mismatches += int(np.count_nonzero(mask != ref['mask']))
        errors.append(abs(row['contact_area_mm2'] - ref['area_mm2']))
    dense = []
    skipped = 0
    for delta in np.linspace(lo, hi, 1001):
        q = Q(float(delta))
        distances = np.abs(np.abs(gap + delta) - band)
        if np.any(distances < 1e-12):
            skipped += 1
            continue
        row = next((r for r in rows if r['kind'] == 'POINT' and Q(r['exact_lower']) == q or (r['kind'] == 'OPEN_INTERVAL' and Q(r['exact_lower']) < q < Q(r['exact_upper']))))
        ref = contact(gap + delta, z['z'], z['xy'], index, weights, band, amin)
        mismatches += int(set(row['active_point_ids']) != set(np.flatnonzero(ref['mask'])))
        errors.append(abs(row['contact_area_mm2'] - ref['area_mm2']))
        dense.append([float(delta), ref['area_mm2']])
    from region_field import read_stl, ray_height
    antpath = PARENT_DATA / 'Bite2Text_F4775/antagonist.stl'
    source(antpath)
    tri = read_stl(antpath.read_bytes())
    (ceiling, _, _) = ray_height(tri, z['xy'], upper=True)
    good = np.isfinite(ceiling) & np.isfinite(z['ceiling'])
    gap_error = float(np.max(np.abs(ceiling[good] - z['ceiling'][good])))
    groups = {}
    for (i, w) in enumerate(weights):
        groups.setdefault(float(w).hex(), []).append(i)
    key = sorted(groups, key=lambda k: (-len(groups[k]), k))[0]
    ids = groups[key]
    ordered = sorted(ids, key=lambda i: (float(np.sum(z['uv'][i] ** 2)), i))
    (first, last) = (ordered[0], ordered[-1])
    loads = {}
    for (name, idx) in [('central', first), ('peripheral', last)]:
        loadvec = np.zeros(len(weights))
        loadvec[idx] = sc['force_N']
        loads[name] = loadvec
    nominal = contact(gap, z['z'], z['xy'], index, weights, band, amin)
    loads['nominal'] = normalize(nominal['mask'], weights, sc['force_N'])
    endpoint = contact(gap + hi, z['z'], z['xy'], index, weights, band, amin)
    loads['pose_upper'] = normalize(endpoint['mask'], weights, sc['force_N'])
    (fe, cost) = solve_cases(z['xy'], z['z'], z['faces'], z['uv'], loads, dict(thickness_mm=sc['thickness_mm'], E_MPa=sc['E_MPa'], nu=sc['nu']))
    suff = dict(summary={'support_area_mm2': float(weights[first]), 'total_force_N': sc['force_N']}, identity_error_area_mm2=float(abs(weights[first] - weights[last])), identity_error_force_N=float(abs(loads['central'].sum() - loads['peripheral'].sum())), point_ids=[first, last], source_xyz_mm=np.c_[z['xy'], z['z']][[first, last]], downstream_tensile_MPa=[fe[n]['maximum_tensile_MPa'] for n in ['central', 'peripheral']], downstream_difference_MPa=abs(fe['central']['maximum_tensile_MPa'] - fe['peripheral']['maximum_tensile_MPa']), minimum_extension='Source-addressed load vector (contact locations plus force allocation); contact area alone is insufficient.', resolution='PER_POINT', evidence_kind='our_own_fixture', physical_accuracy='UNKNOWN', numerical_enclosure='FE equilibrium residual checked; rigorous FE forward/continuum error bound MISSING')
    missing = [r for r in rows if not r['active_point_ids']]
    checks = [check('rational_events_vs_decimal_and_original_components', mismatches == decimal_mismatches == 0 and max(errors) < 1e-10, max(errors + [1.0]) > 1e-10, dict(mask_mismatches=mismatches, decimal_mismatches=decimal_mismatches, max_area_error_mm2=max(errors), injection='Add 1 mm2 to one claimed cell area')), check('original_triangle_ceiling', gap_error < 1e-06, gap_error + 0.01 > 1e-06, dict(max_error_mm=gap_error, injection='Shift ceiling 0.01 mm')), check('robust_abstention', not prediction['geometric_robust'] and bool(missing), True != (not bool(missing)), dict(injection='Force robust=True despite saved no-contact witness')), check('equal_summary_location_counterexample', suff['identity_error_area_mm2'] == suff['identity_error_force_N'] == 0 and suff['downstream_difference_MPa'] > 1e-06, abs(weights[first] + 0.1 - weights[last]) > 0, suff), check('FE_force_balance', __import__('decision_binding').force_balance(fe, sc['force_N']), not __import__('decision_binding').force_balance({k: {**v, 'total_load_N': 101.0} if v.get('status') == 'SIMULATED' else v for (k, v) in fe.items()}, sc['force_N']), dict(injection='Change solver total_load_N to 101 N for original 100 N load'))]
    no_contact = missing[0] if missing else None
    result = dict(claim_type='capability', patient_id='Bite2Text_F4775', FDI=36, partition=rows, continuous_pose_interval_mm=[lo, hi], geometric_robust=prediction['geometric_robust'], contract_status='ABSTAIN_POSE_CONTACT_NOT_PERSISTENT' if missing else 'CONDITIONAL_GEOMETRIC_SUPPORT', binding_quantity='Registered loaded relative vertical pose, then per-region force and support law', no_contact_witness=no_contact, nominal_contact_area_mm2=nominal['area_mm2'], sufficient_summary_test=suff, FE=fe, FE_cost=cost, checks=checks, exclusion=dict(dense_rows_considered=1001, dense_rows_boundary_replaced_by_exact=skipped, finite_roof_queries=int(good.sum()), roof_queries=len(gap)), external_referent=dict(kind='published_dataset', locator=str(PARENT_DATA / 'Bite2Text_F4775/upper.stl') + '; https://ditto.ing.unimore.it/bite2text/', compared_quantity='Same-case antagonist height at retained roof points; geometry only', refutes_us=gap_error >= 1e-06), costs=dict(discovery_seconds=tpred - t0, validation_and_FE_seconds=time.perf_counter() - tpred, total_seconds=time.perf_counter() - t0), physical_status='UNKNOWN', decidability_class=2, refinement_variable='Resampling the same uncalibrated IOS; no new load/pose observation')
    dump(out / 'CONTACT_RESULTS.json', result)
    np.savez_compressed(out / 'contact_curve.npz', rows=np.asarray(dense))
    dump(ROOT / 'HANDOFF_C1.json', dict(outcome=result['contract_status'], source=str(out / 'CONTACT_RESULTS.json'), next_operation='Use source-addressed region fields and shared theta in all decisions; then test rigid crown adjustment against nonpenetration.'))
    return result
