import json, hashlib, time
from pathlib import Path
from datetime import datetime, timezone
import numpy as np
R = Path(__file__).resolve().parents[1]

def put(n, x):
    (R / n).write_text(json.dumps(x, indent=2, allow_nan=False) + '\n')

def main():
    start = time.perf_counter()
    pr = json.loads((R / 'PREREG_R6.json').read_text())
    ge = pr['gates']
    assert hashlib.sha256((R / 'PREREG_R6.json').read_bytes()).hexdigest() == json.loads((R / 'PREREG_R6.sha256.json').read_text())['sha256']
    source = json.loads((R / 'inputs/THESIS_EXTERNAL_ADDENDUM.json').read_text())
    r5 = json.loads((R / 'round5/results.json').read_text())
    r1 = json.loads((R / 'round1/results.json').read_text())
    raw = json.loads((R / 'raw/PREDICTIONS_R1.json').read_text())['rows']
    arches = json.loads((R / 'inputs/ARCHES.json').read_text())
    gamma = r5['gain']
    central = next((r for r in r1['summary'] if r['fdi'] == 11))
    neighbours = []
    materials = []
    for ref in source['neighbours']:
        values = []
        nominal = []
        for row in raw:
            if row['active_fdi'] != ref['active_fdi']:
                continue
            arch = next((a for a in arches if a['case'] == row['case']))
            e = np.array(arch['teeth'][str(ref['response_fdi'])]['buccal_unit'])
            values.append(float(gamma * np.array(row['regional']['wrenches'][str(ref['response_fdi'])][:3]) @ e))
            w = row['practice']['wrenches'][str(ref['response_fdi'])][:3]
            nominal.append(float(-np.array(w) @ e * central['external_N'] / central['practice_N']))
        neighbours.append(dict(active_fdi=ref['active_fdi'], response_fdi=ref['response_fdi'], predicted_N=float(np.median(values)), prediction_range_N=[min(values), max(values)], same_calibrated_nominal_N=float(np.median(nominal))))
    for mat in source['materials']:
        gain = mat['incisor_N'] / central['regional_N']
        nominal_gain = mat['incisor_N'] / central['practice_N']
        for (fdi, key) in [(13, 'canine_N'), (15, 'premolar_N')]:
            r = next((x for x in r1['summary'] if x['fdi'] == fdi))
            materials.append(dict(material=mat['material'], fdi=fdi, predicted_N=gain * r['regional_N'], same_calibrated_nominal_N=nominal_gain * r['practice_N'], calibration_incisor_N=mat['incisor_N']))
    payload = dict(neighbours=neighbours, materials=materials)
    h = hashlib.sha256(json.dumps(payload, sort_keys=True).encode()).hexdigest()
    frozen = dict(frozen_utc=datetime.now(timezone.utc).isoformat(), prereg_sha256=hashlib.sha256((R / 'PREREG_R6.json').read_bytes()).hexdigest(), predictions_sha256=h, predictions=payload, source_visibility=pr['source_visibility'])
    if (R / 'FROZEN_PREDICTIONS_R6.json').exists():
        if json.loads((R / 'FROZEN_PREDICTIONS_R6.json').read_text())['predictions_sha256'] != h:
            raise RuntimeError('R6 drift')
    else:
        put('FROZEN_PREDICTIONS_R6.json', frozen)
    for (row, ref) in zip(neighbours, source['neighbours']):
        truth = ref['force_N']
        row.update(external_N=truth, external_SD_N=ref['sd_N'], relative_error=abs(row['predicted_N'] - truth) / abs(truth), nominal_relative_error=abs(row['same_calibrated_nominal_N'] - truth) / abs(truth), unit='N', resolution_level='PER_TOOTH', locator=source['neighbour_force_table'])
        row['gate'] = row['relative_error'] <= ge['neighbour_signed_relative_error_max']
        row['injected_plus10N_rejected'] = abs(row['predicted_N'] + 10 - truth) / abs(truth) > ge['neighbour_signed_relative_error_max']
    for row in materials:
        ref = next((r for r in source['materials'] if r['material'] == row['material']))
        key = 'canine_N' if row['fdi'] == 13 else 'premolar_N'
        truth = ref[key]
        row.update(external_N=truth, relative_error=abs(row['predicted_N'] - truth) / truth, nominal_relative_error=abs(row['same_calibrated_nominal_N'] - truth) / truth, unit='N', resolution_level='PER_TOOTH', locator=source['active_force_table'])
        row['gate'] = row['relative_error'] <= ge['material_unfit_relative_error_max']
        row['injected_plus10N_rejected'] = abs(row['predicted_N'] + 10 - truth) / truth > ge['material_unfit_relative_error_max']
    gates = dict(neighbour_transfer=all((r['gate'] for r in neighbours)), new_material_transfer=all((r['gate'] for r in materials)), corruptions_refuted=all((r['injected_plus10N_rejected'] for r in neighbours + materials)))
    result = dict(round='R6', claim_type='information_link', outcome='CALIBRATED_ACTIVE_TYPE_FORCE_DOES_NOT_VALIDATE_NEIGHBOUR_OR_GENERAL_MATERIAL_TRANSFER' if not gates['neighbour_transfer'] or not gates['new_material_transfer'] else 'NEIGHBOUR_AND_NEW_MATERIAL_FORCE_TRANSFER_PASS', external_referent=pr['external_referent'], neighbour_rows=neighbours, material_rows=materials, gates=gates, neighbour_groups_pass=sum((r['gate'] for r in neighbours)), new_material_groups_pass=sum((r['gate'] for r in materials)), corroboration=dict(brand_now_verified='Taglus PETG; NOT Park Duran', Kaur_active_force_table='Table5.1 PETG row repeats source mean forces', Kaur_moment_origin='Approximate CoR via Faro Jacobian; coordinates remain unavailable'), wall_s=time.perf_counter() - start, physical_moments_validated=False, source_visibility=pr['source_visibility'], next_operation='Replace artificial housing ground supports with measured full-arch shell/contact support topology, while preserving fitted active-tooth datum and refuting with all6 adjacent force means; exact specimen gap/shape is the missing prerequisite')
    put('round6/results.json', result)
    (R / 'round6/HANDOFF.md').write_text('# R6 primary-table reaction test\n\n' + result['outcome'] + '\n\n' + json.dumps(gates) + '\n\nPrimary thesis Table4.4 PDFpage97 supplies unfit neighbouring forces; Table5.1 PDFpage111 supplies new PET/PU means, with only each central mean calibrated. No material/neighbor endpoint fit or moment-origin fabrication.\n\n' + result['next_operation'] + '\n')
    put('CURRENT_WORK_STATE.json', dict(lane='X27-aligner-regional', state='R6_COMPLETE_PACKAGING', latest_gate=result['outcome'], next_operation=result['next_operation']))
    print(json.dumps(dict(outcome=result['outcome'], gates=gates, neighbours=neighbours, materials=materials), indent=2))
if __name__ == '__main__':
    main()
