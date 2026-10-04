"""Synthetic specimen series, explicitly distinct from independent external facit."""
import csv, copy, hashlib, json
from pathlib import Path
import numpy as np
from lab_alarm import R, model_groups, verified, csv_rows
FIELDS = list(csv_rows(R / 'inputs/LAB_MEASUREMENTS_TEMPLATE.csv')[0])
SIDE = ['specimen_id', 'pre_sinter_ref_length_mm', 'post_sinter_ref_length_mm', 'sinter_ratio_sd', 'assembled_film_mean_um', 'assembled_film_sd_um', 'trace_peak_force_N', 'source_trace_sha256']

def series(seed, p, frozen, fault='nominal'):
    rng = np.random.default_rng(seed)
    n = p['null']
    groups = model_groups(frozen)
    ids = csv_rows(R / 'inputs/LAB_MEASUREMENTS_TEMPLATE.csv')
    cells = [(d, a) for d in ('D1', 'M1', 'M2') for a in (0, 30)]
    ordered = []
    lookup = {(r['design_id'], int(r['load_angle_deg']), int(r['specimen_id'].split('_')[-1])): r for r in ids}
    for rep in range(1, 13):
        for k in rng.permutation(6):
            ordered.append(lookup[(*cells[k], rep)].copy())
    shape = 3.0 if fault == 'shifted_weibull' else n['m']
    scale = n['lambda_reference_N'] * (0.65 if fault == 'shifted_weibull' else 1.0)
    model = 'separable' if fault == 'model_form' else n['model']
    smean = 0.806 if fault == 'wrong_sinter_factor' else n['sinter_factor']
    cmean = 85.0 if fault == 'cement_film' else n['cement_film_mean_um']
    rows = []
    for (i, row) in enumerate(ordered):
        g = groups[row['design_id'], int(row['load_angle_deg'])]
        factor = g['base_factor'] * (g['Q'] if model == 'FE' and g['interaction'] else 1.0)
        f = float(scale * factor * rng.weibull(shape) * np.exp(rng.normal(0, n['force_log_noise_sd'])))
        gap = abs(float(g['gap_um'] + n['dry_gap_bias_um'] + rng.normal(0, n['dry_gap_sd_um'])))
        scan = hashlib.sha256(('SYNTHETIC_SCAN_' + row['specimen_id']).encode()).hexdigest()
        trace = hashlib.sha256(('SYNTHETIC_TRACE_' + row['specimen_id']).encode()).hexdigest()
        row.update(material_batch='SYNTHETIC_3Y', cement_batch='SYNTHETIC_CEMENT', die_material_batch='SYNTHETIC_DIE', die_E_MPa='18000', indenter_diameter_mm='5', crosshead_mm_min='.5', MG_distance_mean_um=str(gap), MG_vertical_mean_um='', fracture_force_N=str(f), initial_stiffness_N_per_mm='', failure_mode='crown_tensile_fracture', fracture_origin='intaglio_tensile_zone', source_scan_sha256=scan, force_trace_sha256=trace)
        sval = float(smean + rng.normal(0, n['sinter_sd']))
        cval = float(cmean + rng.normal(0, n['cement_film_sd_um']))
        s = dict(specimen_id=row['specimen_id'], pre_sinter_ref_length_mm='10', post_sinter_ref_length_mm=str(10 * sval), sinter_ratio_sd=str(n['sinter_sd']), assembled_film_mean_um=str(cval), assembled_film_sd_um=str(n['cement_film_sd_um']), trace_peak_force_N=str(f), source_trace_sha256=trace)
        if fault == 'one_bad_record' and i == 8:
            row['fracture_force_N'] = str(3 * f)
        row['_sidecar'] = s
        rows.append(row)
    return rows

def save_series(path, rows):
    path = Path(path)
    with path.open('w') as h:
        w = csv.DictWriter(h, fieldnames=FIELDS)
        w.writeheader()
        w.writerows(({k: r.get(k, '') for k in FIELDS} for r in rows))
    side = path.with_name(path.stem + '_SIDECAR.csv')
    with side.open('w') as h:
        w = csv.DictWriter(h, fieldnames=SIDE)
        w.writeheader()
        w.writerows((r['_sidecar'] for r in rows))
    return (path, side)
if __name__ == '__main__':
    p = verified(R / 'inputs/DEMO_PROFILE.json')
    f = json.loads((R / 'inputs/X1B_FROZEN_PREDICTIONS.json').read_text())
    for fault in ['nominal', 'shifted_weibull', 'wrong_sinter_factor', 'one_bad_record', 'model_form', 'cement_film']:
        save_series(R / 'raw' / ('SYNTHETIC_' + fault + '.csv'), series(6100, p, f, fault))
