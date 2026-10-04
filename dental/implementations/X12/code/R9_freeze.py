import datetime
import hashlib
import json
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]

def freeze(name, obj):
    p = ROOT / name
    if p.exists():
        raise FileExistsError('Frozen file exists: ' + str(p))
    obj['frozen_utc'] = datetime.datetime.now(datetime.timezone.utc).isoformat()
    payload = json.dumps(obj, indent=2, sort_keys=True) + '\n'
    p.write_text(payload)
    p.with_suffix('.sha256').write_text(hashlib.sha256(payload.encode()).hexdigest() + '  ' + p.name + '\n')
if __name__ == '__main__':
    freeze('PREREG_R1_IDENTITY.json', {'capability': 'Population pulp-to-surface measurements on externally annotated paired scans', 'obstacle': 'Pulpy3D lacks whole-tooth masks and archive central directory', 'changed_operation': 'CRC-verified local record reading, then independent same-image identity matching to local ToothFairy2', 'consumer': 'K23/K41 preparation and heat-distance ports, O material-feasibility necessary conditions, K2 witness inputs', 'metrics': {'source_member_gate': 'exact decompressed length and CRC32', 'identity_gate': 'same shape and physical affine within 1e-5 mm, and all voxel values identical; any crop/resampling requires separately frozen registration trial', 'population_gate': 'at least 20 verified same-image cases for direct atlas; otherwise FAIL_INPUT_COVERAGE', 'mutation_gate': 'flip one CT voxel or one CRC bit must reject'}, 'strongest_equally_informed_control': 'Standard ZIP reader for intact archives; exhaustive voxel equality for pairing; both get same bytes', 'falsifiers': ['CRC mismatch', 'non-identical source images', 'no compatible dimensions/physical frames', 'surface annotations absent'], 'cost': {'preparation': 'metadata and header census wall/CPU/RSS measured', 'fit': 0, 'discovery': 'all considered cases reported', 'validation': 'CRC and CT equality separately charged', 'queries': 'unknown before execution', 'fallback': 'no pairing -> freeze a changed image-boundary measurement with explicit closures'}, 'external_referent': {'kind': 'published_dataset', 'locator': 'https://ditto.ing.unimore.it/pulpy3d/', 'compared_quantity': 'supplied pulp label semantics and scan count', 'refutes_us': False}, 'disk_limit_bytes': 3000000000, 'threads_max': 4})
