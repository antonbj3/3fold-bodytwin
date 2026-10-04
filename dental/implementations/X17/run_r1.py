from dental_release.paths import expand as _release_expand
import sys, os, json, time, hashlib, resource, platform
from datetime import datetime, timezone
from pathlib import Path
import importlib.util
import numpy as np
from airway import profile, landmarks, crop_pharynx, scenarios, pose, gate
ROOT = Path(__file__).resolve().parent
DATA = Path(_release_expand('@DENTAL_WORK_ROOT@/X17-sleep-apnea'))
DATA.mkdir(parents=True, exist_ok=True)
S = _release_expand('@DENTAL_INPUT_ROOT@/workspace/cells/geometry/tf2_io.py')
spec = importlib.util.spec_from_file_location('tf2_io', S)
tf = importlib.util.module_from_spec(spec)
spec.loader.exec_module(tf)
prereg = json.loads((ROOT / 'PREREG_R1.json').read_text())
assert hashlib.sha256((ROOT / 'PREREG_R1.json').read_bytes()).hexdigest() == (ROOT / 'PREREG_R1.json.sha256').read_text().split()[0]

def write(name, obj):
    (ROOT / name).write_text(json.dumps(obj, indent=2, ensure_ascii=False, allow_nan=False) + '\n')

def state(op, gate_value, next_op):
    write('CURRENT_WORK_STATE.json', dict(lane='X17-sleep-apnea', updated_utc=datetime.now(timezone.utc).isoformat(), current_operation=op, latest_gate=gate_value, next_operation=next_op))

def freeze(name, payload):
    digest = hashlib.sha256(json.dumps(payload, sort_keys=True, allow_nan=False).encode()).hexdigest()
    p = ROOT / name
    if p.exists():
        old = json.loads(p.read_text())
        assert old['prediction_sha256'] == digest, 'Frozen prediction drift'
    else:
        write(name, {'frozen_at_utc': datetime.now(timezone.utc).isoformat(), 'prediction_sha256': digest, 'predictions': payload, 'source_aware_retrospective': True, 'note': 'No new physical measurement. Published target values were visible during source selection; this is a reproducibility freeze, not prospective or blinded validation.'})
    return digest

def source_inputs():
    from xml.etree import ElementTree as E
    path = Path(_release_expand('@DENTAL_CORPUS_ROOT@/europepmc/fulltext/PMC10492877.xml'))
    r = E.parse(path).getroot()
    tables = []
    for t in r.findall('.//table-wrap'):
        table = {'label': ''.join(t.find('label').itertext()), 'rows': []}
        for tr in t.findall('.//tr'):
            table['rows'].append([' '.join(''.join(c.itertext()).split()) for c in tr])
        tables.append(table)
    write('raw/source_tables.json', {'path': str(path), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest(), 'doi': '10.1007/s00784-023-05186-w', 'tables': [t for t in tables if t['label'] in ['Table 2', 'Table 4', 'Table 5']]})
    t2 = next((t for t in tables if t['label'] == 'Table 2'))
    t4 = next((t for t in tables if t['label'] == 'Table 4'))

    def row(t, label):
        return next((x for x in t['rows'] if x[0].strip().replace(' ', '') == label.replace(' ', '')))

    def val(s):
        return float(s.split('±')[0].strip())
    dose = row(t2, 'Advancement (mm)')
    ap = row(t4, 'A-P (mm)')
    lat = row(t4, 'La (mm)')
    v = row(t4, 'V (cm3)')
    return ([dict(group=g, n=n, dose_mm=val(dose[j + 1]), ap0_mm=val(ap[1 + 2 * j]), ap1_mm=val(ap[2 + 2 * j]), la0_mm=val(lat[1 + 2 * j]), la1_mm=val(lat[2 + 2 * j]), v0_cm3=val(v[1 + 2 * j]), v1_cm3=val(v[2 + 2 * j])) for (j, (g, n)) in enumerate([('responders', 15), ('nonresponders', 16)])], t4)

def source_area(t4):
    x = next((x for x in t4['rows'] if x[0].startswith('CSAmin')))
    return [[float(x[1 + 2 * j].split('±')[0]), float(x[2 + 2 * j].split('±')[0])] for j in range(2)]

def main():
    (ROOT / 'raw').mkdir(exist_ok=True)
    start = time.monotonic()
    allcase = []
    raw = {}
    state('R1 measured profile and motion scenarios', 'PREREG frozen', 'gate R1 then change operation')
    for case in prereg['cases']:
        t = time.monotonic()
        print('case', case, flush=True)
        record = {'case': case}
        try:
            (arr, spacing, hdr) = tf.load(case)
            record['spacing_zyx_mm'] = list(spacing)
            record['shape_zyx'] = list(arr.shape)
            record['label_payload_sha256'] = hashlib.sha256(arr).hexdigest()
            record['header'] = hdr
            (crop, offset, frac, ncomp) = crop_pharynx(arr)
            np.savez_compressed(DATA / (case + '_airway.npz'), mask=crop, spacing=np.array(spacing), offset=np.array(offset))
            f = DATA / (case + '_airway.npz')
            record['crop_file'] = {'path': str(f), 'bytes': f.stat().st_size, 'sha256': hashlib.sha256(f.read_bytes()).hexdigest()}
            pr = profile(crop, spacing, offset[0])
            ids = np.where(pr['eligible'])[0]
            ix = int(ids[np.argmin(pr['area_mm2'][ids])])
            record.update(volume_mm3=pr['volume_mm3'], raw_annotation_minimum_mm2=float(pr['area_mm2'][pr['area_mm2'] > 0].min()), interior_minimum_mm2=float(pr['area_mm2'][ix]), minimum_array_z_mm=float(pr['z_mm'][ix]), margin_mm=5.0, largest_component_retained_fraction=frac, n_components=ncomp, cap_sensitivity={str(m): float(profile(crop, spacing, offset[0], m)['area_mm2'][profile(crop, spacing, offset[0], m)['eligible']].min()) for m in [3, 5, 8]})
            expected = crop.sum() * np.prod(spacing)
            record['volume_relative_error'] = float(abs(expected - pr['volume_mm3']) / expected)
            raw[case] = {k: v.tolist() if isinstance(v, np.ndarray) else v for (k, v) in pr.items()}
            try:
                lm = landmarks(arr, spacing)
                record['landmarks'] = lm
                rows = scenarios(pr, lm, spacing, offset, crop, [0, 2, 4, 6, 8, 10], [0, 2, 4, 6], [0, 0.25, 0.5, 0.75, 1])
                record['scenarios'] = rows
                record['zero_motion_error_mm2'] = max((abs(s['minimum_mm2'] - record['interior_minimum_mm2']) for s in rows if s['advancement_mm'] == 0 and s['rotation_deg'] == 0))
                err = []
                for m in [0, 2, 4, 6, 8, 10]:
                    for a in [0, 2, 4, 6]:
                        inc = np.array(lm['incisor_zyx_mm'])
                        (out, _) = pose(inc, lm['pivot_zyx_mm'], inc, m, a, lm['ap_sign'], lm['si_sign'])
                        err.append(abs((out[1] - inc[1]) * lm['ap_sign'] - m))
                record['pose_error_mm'] = max(err)
                mm = [r['minimum_mm2'] for r in rows if r['advancement_mm'] == 6]
                record['response_interval_at_6mm_mm2'] = [min(mm), max(mm)]
                record['relative_width_at_6mm'] = float((max(mm) - min(mm)) / record['interior_minimum_mm2'])
                record['geometric_status'] = 'CONDITIONAL_SCENARIOS'
            except ValueError as e:
                record['geometric_status'] = 'UNKNOWN'
                record['pose_exclusion'] = str(e)
            del arr, crop
        except Exception as e:
            record.update(geometric_status='UNKNOWN', exclusion=type(e).__name__ + ': ' + str(e))
        record['seconds'] = time.monotonic() - t
        allcase.append(record)
        write('raw/cases_R1.json', allcase)
        write('raw/profiles.json', raw)
    (source, t4) = source_inputs()
    areas = source_area(t4)
    preds = []
    for (inp, (a0, a1)) in zip(source, areas):
        preds.append(dict(inp, area0_mm2=a0, ap_only_predicted_area1_mm2=a0 * inp['ap1_mm'] / inp['ap0_mm'], proportional_volume_predicted_area1_mm2=a0 * inp['v1_cm3'] / inp['v0_cm3']))
    digest = freeze('FROZEN_PREDICTIONS_R1.json', {'external_predictions': preds, 'geometry_scenarios': [{k: c[k] for k in ['case', 'scenarios'] if k in c} for c in allcase]})
    refs = []
    for (p, (_, target)) in zip(preds, areas):
        refs.append(dict(p, area1_observed_mm2=target, observed_change_per_mm_mm2=(target - p['area0_mm2']) / p['dose_mm'], ap_only_relative_error=abs(p['ap_only_predicted_area1_mm2'] - target) / target, volume_control_relative_error=abs(p['proportional_volume_predicted_area1_mm2'] - target) / target))
    eligible = [c for c in allcase if 'relative_width_at_6mm' in c]
    controls = {'volume_integration': all((gate(c['volume_relative_error'], 1e-10) for c in allcase if 'volume_relative_error' in c)), 'zero_motion': all((gate(c['zero_motion_error_mm2'], 1e-09) for c in eligible)), 'pose_command': all((gate(c['pose_error_mm'], 1e-09) for c in eligible)), 'usable_width': len(eligible) > 0 and sum((gate(c['relative_width_at_6mm'], 0.5) for c in eligible)) / len(eligible) >= 0.8, 'published_ap_only': all((gate(r['ap_only_relative_error'], 0.2) for r in refs)), 'published_volume_control': all((gate(r['volume_control_relative_error'], 0.2) for r in refs))}
    mutations = {'volume_integration': not gate(0.01, 1e-10), 'zero_motion': not gate(1.0, 1e-09), 'pose_command': not gate(1.0, 1e-09), 'usable_width': not gate(1.0, 0.5), 'published_ap_only': not gate(1.0, 0.2), 'published_volume_control': not gate(1.0, 0.2)}
    out = {'round': 'R1', 'outcome': 'FAIL_IDENTIFIABILITY' if not controls['usable_width'] else 'CONDITIONAL_ONLY', 'cases_total': len(allcase), 'pose_eligible': len(eligible), 'controls': controls, 'injected_error_rejected': mutations, 'external_comparisons': refs, 'external_referent': prereg['external_referent'], 'prediction_sha256': digest, 'cost': {'execution_seconds': time.monotonic() - start, 'peak_rss_MiB': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, 'CPU_seconds': resource.getrusage(resource.RUSAGE_SELF).ru_utime + resource.getrusage(resource.RUSAGE_SELF).ru_stime, 'disk_bytes': sum((p.stat().st_size for p in DATA.glob('*') if p.is_file())), 'preparation_seconds': 'UNMEASURED', 'external_measurement_acquisition_cost': 'UNMEASURED; existing expert annotations and literature'}, 'software': {'python': platform.python_version(), 'numpy': np.__version__, 'tf2_reader_path': S, 'tf2_reader_sha256': hashlib.sha256(Path(S).read_bytes()).hexdigest()}}
    write('raw/round_R1.json', out)
    state('R1 adjudicated', out['outcome'], 'R2 add independently measured lateral shape response; test held-out CSA and preserve failed gates')
    print(json.dumps(out, indent=2), flush=True)
if __name__ == '__main__':
    main()
