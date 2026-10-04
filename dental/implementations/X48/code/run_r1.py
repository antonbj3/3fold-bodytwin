from common import *
import time, csv, resource
from sufficiency import run as sufficiency

def run():
    start = time.monotonic()
    (rows, ix) = pose_index()
    reg = read(ROOT / 'PREREG_R1.json')
    out = []
    keys = sorted({(r['patient'], r['jaw'], r['fdi']) for r in rows})
    for (p, j, f) in keys:
        ax = axes(ix, p, j, f)
        c = np.array(ix[p, j, 'Sirona', 0, f]['centroid_base_mm'])
        for t in range(1, 10):
            for (interval, s) in [('weekly', t - 1), ('cumulative', 0)]:
                rr = [ix[p, j, src, k, f] for (src, k) in [('Sirona', t), ('Sirona', s), ('Bottmedical', t), ('Bottmedical', s), ('Bottmedical', 0)]]
                (st, sp, pt, pp, p0) = rr
                vals = {}
                for name in REFS:
                    (A, B) = (np.array(st['H'][name]), np.array(sp['H'][name]))
                    P = np.array(pt['H'][name]) @ np.linalg.inv(p0['H'][name])
                    Q = np.array(pp['H'][name]) @ np.linalg.inv(p0['H'][name])
                    (av, at, aw) = delta(A, B, c, ax)
                    (pv, ptv, pw) = delta(P, Q, c, ax)
                    vals[name] = {k: dict(achieved=av[k], planned=pv[k], signed_ratio=av[k] / pv[k] if abs(pv[k]) > 1e-15 else None) for k in av}
                for kind in vals[REFS[0]]:
                    unit = 'deg' if kind in ['buccolingual_tip', 'mesiodistal_tip', 'axial_rotation'] else 'mm'
                    if unit == 'mm':
                        ua = 0.2 + st['sampling_spread_mm'] + sp['sampling_spread_mm']
                        up = 0.2 + pt['sampling_spread_mm'] + pp['sampling_spread_mm'] + 2 * p0['sampling_spread_mm']
                    else:
                        ua = 2.0
                        up = 2.0
                    perref = {k: {**v[kind], 'ratio_interval': ratio_box(v[kind]['achieved'], v[kind]['planned'], ua, up)} for (k, v) in vals.items()}
                    denom = all((abs(v['planned']) > 5 * up for v in perref.values()))
                    bounds = [v['ratio_interval'] for v in perref.values()]
                    envelope = [min((b[0] for b in bounds)), max((b[1] for b in bounds))] if denom and all((b is not None for b in bounds)) else None
                    out.append(dict(patient=p, jaw=j, fdi=f, week=t, interval=interval, motion=kind, unit=unit, resolution='PER_TOOTH', timescale='HANDOVER_WEEK', reference_variants=perref, achieved_radius=ua, planned_radius=up, radii_resolution='PHENOMENOLOGICAL', radius_is_empirical_error_bound=False, conditional_ratio_interval=envelope, ratio_status='CONDITIONAL_REFERENCE_ENVELOPE' if envelope is not None else 'UNIDENTIFIED_DENOMINATOR_OR_REFERENCE', absolute_motion_status='UNKNOWN_NO_FIXED_ANATOMICAL_REFERENCE', surface_gate_pass=all((r['surface_fit_pass'] for r in rr)), plan_proxy_excluded_in_conditional_envelope=bool(envelope and (envelope[1] < 1 or envelope[0] > 1)), bodily_status='UNKNOWN_ROOT_NOT_OBSERVED'))
    suff = sufficiency()
    paper = [r for r in out if r['jaw'] == 'OK' and r['fdi'] in [17, 27] and (r['week'] == 8) and (r['interval'] == 'cumulative') and (r['motion'] == 'buccolingual_tip')]
    write(ROOT / 'raw/R1_SIGNED_MOTION.json', out)
    write(ROOT / 'raw/R1_PAPER_MOLAR_COMPARISON.json', paper)
    with (ROOT / 'raw/signed_motion.csv').open('w') as o:
        fields = ['patient', 'jaw', 'fdi', 'week', 'interval', 'motion', 'unit', 'resolution', 'absolute_motion_status', 'ratio_status', 'achieved', 'planned', 'signed_ratio', 'reference', 'achieved_radius', 'planned_radius']
        w = csv.DictWriter(o, fieldnames=fields)
        w.writeheader()
        for r in out:
            for (ref, v) in r['reference_variants'].items():
                w.writerow({**{k: r[k] for k in fields if k in r}, **{k: v[k] for k in ['achieved', 'planned', 'signed_ratio']}, 'reference': ref})
    r = dict(round='R1', claim_type='information_link', external_referent=reg['external_referent'], patients=2, teeth=len(keys), weekly_tooth_intervals=len(keys) * 9, signed_component_rows=len(out), absolute_identified_rows=0, conditional_ratio_rows=sum((r['conditional_ratio_interval'] is not None for r in out)), conditional_plan_proxy_exclusions=sum((r['plan_proxy_excluded_in_conditional_envelope'] for r in out)), denominator_or_reference_unidentified=sum((r['conditional_ratio_interval'] is None for r in out)), absolute_gate='FAIL_FIXED_ANATOMICAL_REFERENCE_ABSENT', paper_molar_comparison='FOUR_TARGETS_EXPORTED; different crown gauges, no independent quantitative six-DOF target', bodily='UNKNOWN_ROOT_NOT_MEASURED', source_scanner_floor=dict(value_um=50.0, sd_um=25.0, quantity='Euclidean distance to crown-like PEEK features', resolution='PER_SURFACE_REGION', locator='10.3390/oral4040039 Table 1, p494', patient_repeat_scan_floor='UNKNOWN_NO_REPEATED_SCANS; cannot turn source surface SD into a clinical pose bound'), numerical_sensitivity_debt='0.2 mm /2 deg per two-pose component plus sampling spread; PHENOMENOLOGICAL; no rigorous anatomical enclosure', sufficiency=suff, wall_s=time.monotonic() - start, peak_rss_mib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024)
    write(ROOT / 'raw/R1_RESULTS.json', r)
    text = f"R1: {len(keys)} teeth, {len(keys) * 9} weekly intervals. Absolute gate FAIL: no fixed skull reference. {r['conditional_ratio_rows']}/{len(out)} component rows have a conditional reference envelope; these are not empirical bounds. Body motion UNKNOWN. Directional sufficiency fails exactly (identity error 0, downstream difference 0.5 mm).\nNext construction: freeze a gauge-invariant relative SE3 crown field; compare signed rotation to paper molar-vs-incisor finding and plan proxy.\n"
    (ROOT / 'HANDOFF_R1.md').write_text(text)
    (ROOT / 'HANDOFF.md').write_text(text)
    state('R1_ADJUDICATED', r['absolute_gate'], 'Freeze R2 relative crown pose and exact plan-direction completion; preserve R1 denominator failures')
    print({k: v for (k, v) in r.items() if k not in ['sufficiency']})
    return r
if __name__ == '__main__':
    run()
