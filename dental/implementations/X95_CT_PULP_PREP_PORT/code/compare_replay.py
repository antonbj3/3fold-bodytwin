import argparse, json
from pathlib import Path
from identity import ROOT, write

def normalized(obj):
    if isinstance(obj, dict):
        return {k: normalized(v) for (k, v) in obj.items() if k not in ['wall_s', 'peak_RSS_MiB']}
    if isinstance(obj, list):
        return [normalized(v) for v in obj]
    if isinstance(obj, str) and obj.startswith('/') and ('/X95-ct-pulp-prep-port/' in obj or '/LANE_X95_CT_PULP_PREP_PORT/' in obj):
        return Path(obj).name
    return obj

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--replay', type=Path, required=True)
    a = parser.parse_args()
    rows = []
    for name in ['R1_ROWS.json', 'R1_WITNESSES.json', 'R1_INTERVAL_ROWS.json', 'R1_FRAMES.json', 'R1_CONTROLS.json', 'R2_ROWS.json', 'R2_CONTROLS.json', 'R3_ROWS.json', 'R3_CONTROLS.json', 'INDEPENDENT_FACE_CHECKS.json', 'R1_OUTCOME.json', 'R2_OUTCOME.json', 'R3_OUTCOME.json', 'EXPORT_OUTCOME.json', 'EXPORT_CONTROLS.json']:
        left = json.loads((ROOT / 'raw' / name).read_text())
        right = json.loads((a.replay / 'raw' / name).read_text())
        rows.append({'file': name, 'all_numeric_and_scientific_fields_exact': normalized(left) == normalized(right)})
    artifacts = []
    for name in ['R1_MANIFEST.json', 'R2_MANIFEST.json', 'R3_MANIFEST.json', 'EXPORT_RECEIPTS.json']:
        left = json.loads((ROOT / 'raw' / name).read_text())
        right = json.loads((a.replay / 'raw' / name).read_text())
        l = {Path(r['path']).name: r['sha256'] for r in left}
        rr = {Path(r['path']).name: r['sha256'] for r in right}
        artifacts.append({'manifest': name, 'artifacts': len(l), 'all_bytes_SHA256_identical': l == rr})
    assert all((r['all_numeric_and_scientific_fields_exact'] for r in rows))
    assert all((r['all_bytes_SHA256_identical'] for r in artifacts))
    out = {'replay_root': str(a.replay), 'kind': 'fresh complete producer replay, not independent scientist review', 'raw_table_comparisons': rows, 'artifact_comparisons': artifacts, 'passed': True, 'normalization': 'Only current output paths and volatile wall/RSS measurements excluded; all numerical fields and artifact bytes retained'}
    write(ROOT / 'raw/REPLAY_COMPARISON.json', out)
    print(json.dumps(out, indent=2))
if __name__ == '__main__':
    main()
