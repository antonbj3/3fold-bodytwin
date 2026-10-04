"""Apply the unchanged packaged X55 CLI to a matched measured specimen."""
from dental_release.paths import expand as _release_expand
import argparse, datetime, hashlib, json, subprocess
from pathlib import Path
R = Path(__file__).resolve().parent
DATA = Path(_release_expand('@DENTAL_WORK_ROOT@/X70_end_to_end'))
P = R.parent / 'DEMO48_PACKAGE/batch14/demos/X55/code/metrology.py'

def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--scans', nargs='+', required=True)
    ap.add_argument('--metadata', required=True)
    ap.add_argument('--output')
    a = ap.parse_args()
    frozen = R / 'FROZEN_PREDICTIONS.json'
    if sha(frozen) != (R / 'FROZEN_PREDICTIONS.json.sha256').read_text().strip():
        raise ValueError('Frozen patient prediction drift')
    f = json.loads(frozen.read_text())['payload']
    latest = json.loads((R / 'LATEST_RUN.json').read_text())
    source = Path(latest['runs'][0]['directory'])
    for name in ['crown.stl', 'REGIONS.json']:
        if sha(source / name) != f['files'][name]['sha256']:
            raise ValueError('Frozen geometry/region changed')
    meta = json.loads(Path(a.metadata).read_text())
    if meta.get('design_sha256') != f['files']['crown.stl']['sha256']:
        raise ValueError('Metadata must identify the frozen X70 design_sha256')
    for name in a.scans:
        row = meta.get('scans', {}).get(Path(name).name, {})
        if not row.get('specimen_id', '').startswith('X70_PATIENT27_'):
            raise ValueError('Same X70 patient/design specimen ID required; no D1 mapping')
        if row.get('physical_measurement') is not True:
            raise ValueError('Declare actual physical measurement; synthetic controls use run_end_to_end.sh')
    dest = Path(a.output).resolve() if a.output else DATA / 'lab_measurements' / datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%S_%fZ')
    if dest.exists():
        raise ValueError('Use a new output directory; prior measurement gates preserved')
    dest.parent.mkdir(parents=True, exist_ok=True)
    cmd = ['/usr/bin/bwrap', '--die-with-parent', '--unshare-net', '--ro-bind', '/', '/', '--bind', str(dest.parent), str(dest.parent), '--tmpfs', '/tmp', '--proc', '/proc', '--dev', '/dev', _release_expand('@DENTAL_PYTHON@'), '-B', str(P), '--design', str(source / 'crown.stl'), '--regions', str(source / 'REGIONS.json'), '--units', 'mm', '--scans', *[str(Path(p).resolve()) for p in a.scans], '--metadata', str(Path(a.metadata).resolve()), '--registration', 'datum', '--output', str(dest)]
    code = subprocess.run(cmd).returncode
    if code:
        return code
    report = json.loads((dest / 'report.json').read_text())
    report['patient_frozen_predictions_sha256'] = sha(frozen)
    report['physical_accuracy_interpretation'] = 'CAD discrepancy only unless independent reference/datum/noise supplied; not seated cement film or calibrated force'
    (dest / 'X70_MEASUREMENT_BINDING.json').write_text(json.dumps(report, indent=2) + '\n')
    print('X55 physical scan comparison saved; X70 force/film likelihood remains UNKNOWN:', dest)
    return 0
if __name__ == '__main__':
    try:
        raise SystemExit(main())
    except (ValueError, KeyError) as e:
        print('X70_METROLOGY_REJECT:', e)
        raise SystemExit(2)
