"""Consume a prospective registered force-fraction JSON, without fitting tensors."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from force_budget import upper_stress
ROOT = Path(__file__).resolve().parent

def query(path):
    prereg = ROOT / 'PREREG_R3.json'
    if hashlib.sha256(prereg.read_bytes()).hexdigest() != (ROOT / 'PREREG_R3.sha256').read_text().strip():
        raise ValueError('frozen prereg hash drift')
    frozen = json.loads(prereg.read_text())
    manifest_path = ROOT / 'SOURCE_MANIFEST_R3.json'
    if hashlib.sha256(manifest_path.read_bytes()).hexdigest() != frozen['source_manifest_sha256']:
        raise ValueError('frozen basis-manifest hash drift')
    value = json.loads(Path(path).read_text())
    required = ['case', 'fdi', 'basis_sha256', 'force_fractions', 'epsilon_fraction', 'total_force_N', 'provenance_locator']
    if any((value.get(k) is None or value.get(k) == '' for k in required)):
        raise ValueError('NOT_RUN_MISSING_MEASUREMENT: every measured field and provenance locator required')
    manifest = json.loads(manifest_path.read_text())
    entry = next((x for x in manifest if x.get('fdi') == value['fdi'] and x.get('case') == value['case']), None)
    if not entry or value['basis_sha256'] != entry['sha256']:
        raise ValueError('measurement/basis case, tooth or hash mismatch')
    p = ROOT / entry['local']
    if hashlib.sha256(p.read_bytes()).hexdigest() != entry['sha256']:
        raise ValueError('basis hash drift')
    if value['total_force_N'] != 100:
        raise ValueError('this frozen screen is for 100 N; total-force uncertainty and rescaling need a separate contract')
    basis = np.load(p)['stress_tensors_MPa']
    result = upper_stress(basis, value['force_fractions'], value['epsilon_fraction'])
    result.update(case=value['case'], fdi=value['fdi'], measurement_input_sha256=hashlib.sha256(Path(path).read_bytes()).hexdigest(), provenance_locator=value['provenance_locator'], basis_sha256=entry['sha256'], engineering_threshold_MPa=100, conditional_below_100MPa=result['upper_MPa'] < 100, physical_accuracy='UNKNOWN until independent strain/load-displacement and mesh-convergence validation')
    return result
if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('measurement_json')
    a = p.parse_args()
    print(json.dumps(query(a.measurement_json), indent=2))
