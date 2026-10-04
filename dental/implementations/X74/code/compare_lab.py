"""Source-replication scorer. A software fixture never establishes measurement."""
import argparse, math, json
import numpy as np
from common import ROOT, check_frozen, dump, sha
from optics import de00

def score(group, predictions=None):
    p = check_frozen('PREREG_R4.json')
    pred = predictions or check_frozen('FROZEN_LAB_PREDICTIONS.json')['predictions']
    errors = []

    def require(condition, reason):
        if not condition:
            errors.append(reason)
    fixture = group.get('evidence_type') == 'our_own_fixture'
    require(group.get('evidence_type') in ('physical_measurement', 'our_own_fixture'), 'missing evidence type')
    rows = group.get('specimens', [])
    require(len(rows) >= p['metrics']['future_source_replication_n_min'], 'fewer than frozen n=10 specimens')
    require(len({r.get('specimen_id') for r in rows}) == len(rows), 'missing/duplicate specimen IDs')
    if group.get('kind') == 'bond_coupon':
        key = group.get('system')
        covered = group.get('with_cover')
        ref = next((r for r in pred['bond'] if r['system'] == key and r['with_cover'] == covered), None)
        require(ref is not None, 'unknown paired adhesive/cement/cover state')
        expected = pred['bond_protocol']
        for (k, v) in expected.items():
            require(group.get(k) == v, 'protocol mismatch: ' + k)
        adhesives = {'RXU': 'Scotchbond Universal', 'SC2': 'Prime & Bond Universal', 'MEC': 'Optibond Universal'}
        require(group.get('adhesive_product') == adhesives.get(key), 'missing/wrong paired adhesive product')
        if covered:
            cover = group.get('measured_cover_thickness_bounds_mm')
            require(isinstance(cover, list) and len(cover) == 2 and (1.0 - 0.01 <= cover[0] <= cover[1] <= 1.0 + 0.01), 'missing/unsupported local cover thickness bounds')
            space = group.get('measured_cover_internal_space_mm')
            require(isinstance(space, (float, int)) and abs(space - 0.3) <= 0.01, 'missing/unsupported cover internal space')
            require(bool(group.get('cover_shape_scan_locator')), 'missing actual capsule/cover dome scan')
        values = []
        for row in rows:
            require(row.get('failure_mode') in ('ADHESIVE', 'MIXED'), 'wrong/unknown failure mode')
            force = row.get('peak_force_N')
            area = row.get('measured_bonded_area_mm2')
            require(isinstance(force, (float, int)) and math.isfinite(force) and (force > 0), 'missing/invalid peak force [N]')
            require(isinstance(area, (float, int)) and math.isfinite(area) and (area > 0), 'missing/invalid measured bonded area [mm2]')
            require(bool(row.get('failure_origin_record')), 'missing fracture origin record')
            require(bool(row.get('raw_observation_locator')), 'missing raw observation locator')
            if isinstance(force, (float, int)) and isinstance(area, (float, int)) and (force > 0) and (area > 0):
                values.append(force / area)
        mu = float(np.mean(values)) if values else None
        sd = float(np.std(values, ddof=1)) if len(values) > 1 else None
        err = abs(mu - ref['mean_MPa']) / ref['mean_MPa'] if mu is not None and ref else None
        require(err is not None and err <= p['metrics']['future_source_replication_relative_SBS_mean_error_max'], 'SBS mean outside frozen source-replication tolerance')
        numeric = {'mean_MPa': mu, 'sample_SD_MPa': sd, 'relative_source_mean_error': err}
    elif group.get('kind') == 'optical_stack':
        ref = next((r for r in pred['optical'] if r['substrate'] == group.get('substrate') and r['cement'] == group.get('cement')), None)
        require(ref is not None, 'unknown substrate/cement stack')
        require(group.get('cement') != 'CG', 'uncemented control has no qualified replication exposure contract in this scorer')
        for (k, v) in pred['optical_protocol'].items():
            require(group.get(k) == v, 'protocol mismatch: ' + k)
        values = []
        for row in rows:
            lab = row.get('Lab')
            h = row.get('measured_tile_thickness_mm')
            film = row.get('measured_cement_thickness_um')
            require(isinstance(lab, list) and len(lab) == 3 and all((isinstance(x, (float, int)) and math.isfinite(x) for x in lab)), 'missing/invalid Lab')
            require(isinstance(h, (float, int)) and math.isfinite(h) and (abs(h - 1.0) <= p['metrics']['protocol_thickness_error_mm_max']), 'unsupported tile thickness')
            require(isinstance(film, (float, int)) and math.isfinite(film) and (film >= 0), 'missing/invalid cement film measurement')
            require(bool(row.get('raw_observation_locator')), 'missing raw observation locator')
            if isinstance(lab, list) and len(lab) == 3 and all((isinstance(x, (float, int)) and math.isfinite(x) for x in lab)):
                values.append(lab)
        mu = np.mean(values, axis=0).tolist() if values else None
        err = de00(mu, ref['mean_Lab']) if mu is not None and mu[2] > 0 and ref else None
        require(err is not None and err <= p['metrics']['future_source_replication_Lab_de00_error_max'], 'Lab outside frozen source-replication tolerance')
        numeric = {'mean_Lab': mu, 'de00_error_to_source_mean': err}
    else:
        errors.append('unknown specimen kind')
        numeric = {}
    require(bool(group.get('material_batch_id')), 'missing material batch ID')
    return {'gate': 'FAIL' if errors else 'PASS_SOFTWARE_FIXTURE' if fixture else 'PASS_SOURCE_REPLICATION', 'errors': errors, **numeric, 'claim_type': 'capability', 'resolution_level': 'POPULATION', 'evidence_type': group.get('evidence_type'), 'patient_transfer': 'UNKNOWN: source replication does not validate curved crown transfer', 'scientific_admission': False}

def main():
    a = argparse.ArgumentParser()
    a.add_argument('observation_json')
    a.add_argument('--output')
    args = a.parse_args()
    group = json.loads(open(args.observation_json).read())
    out = score(group)
    out['observation_sha256'] = sha(args.observation_json)
    if args.output:
        dump(args.output, out)
    print(json.dumps(out, ensure_ascii=False))
    raise SystemExit(1 if out['gate'] == 'FAIL' else 0)
if __name__ == '__main__':
    main()
