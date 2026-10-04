"""One-time small-file snapshot. Never copy trees or refresh a frozen source."""
from dental_release.paths import expand as _release_expand
from pathlib import Path
import datetime
import hashlib
import json
import shutil
ROOT = Path(__file__).resolve().parent
SOURCES = ROOT.parent

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    if (ROOT / 'SOURCE_MANIFEST.json').exists():
        raise SystemExit('Snapshot already frozen; use a new version for source changes')
    selections = {_release_expand('X5'): ['results.json', 'R3_CONNECTED_COUNTEREXAMPLES.json', 'SUMMARY_R1.json', 'SUMMARY_R2.json', 'SUMMARY_R3.json', 'PREREG_R3.json', 'FROZEN_PREDICTIONS_R3.json', 'DATA_MANIFEST.json', 'DATA_MANIFEST_R3.json', 'MODEL_EVALUATION_EXAMPLE.json', 'assess_model.py', 'voxel_decisions.py', 'decision_error.py', 'SOURCE_DESIGN_MANIFEST.json'], 'LANE_X26_DECIDABILITY': ['results.json', 'RAW_R2_CANAL.json', 'PREREG_R2.json', 'FROZEN_PREDICTIONS.json', 'inputs/X5_margins.csv', 'inputs/GEOM_field.json', 'decidability.py'], 'LANE_X31_CLINICAL_ANSWERS': ['results.json', 'tables/GUIDE_MARGIN.csv', 'raw/GUIDE_PER_SITE.csv', 'raw/R1_IDENTIFIABILITY.json', 'inputs/guide_parameters.json'], _release_expand('X35_STATISTICS'): ['results.json', 'raw/ENDPOINT_RECORDS.json'], 'PROOF_LANE_XREVIEW_CLINICAL': ['REVIEW_LANE_X5_DECISION_SEG_ERROR.json', 'CORRECTIONS.jsonl']}
    manifest = []
    omitted_optional = []
    for (lane, names) in selections.items():
        for name in names:
            src = SOURCES / lane / name
            if not src.exists():
                if name == 'inputs/guide_parameters.json':
                    omitted_optional.append(str(src))
                    continue
                raise FileNotFoundError(src)
            dst = ROOT / 'inputs' / lane / name
            if src.stat().st_size > 20000000:
                raise ValueError(f'Selected file too large: {src}')
            dst.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(src, dst)
            manifest.append({'source': str(src), 'local': str(dst.relative_to(ROOT)), 'bytes': dst.stat().st_size, 'sha256': sha(dst)})
    probes = json.loads((ROOT / 'inputs/LANE_X5_DECISION_SEG_ERROR/results.json').read_text())
    for entry in probes['data_artifacts']:
        if entry['path'].endswith('_local.npz'):
            src = Path(entry['path'])
            dst = ROOT / 'inputs/geometry' / src.name
            dst.parent.mkdir(parents=True, exist_ok=True)
            if sha(src) != entry['sha256']:
                raise ValueError(f'Source NPZ hash drift: {src}')
            shutil.copyfile(src, dst)
            manifest.append({'source': str(src), 'local': str(dst.relative_to(ROOT)), 'bytes': dst.stat().st_size, 'sha256': sha(dst)})
    now = datetime.datetime.now(datetime.timezone.utc).isoformat()
    (ROOT / 'SOURCE_MANIFEST.json').write_text(json.dumps({'frozen_utc': now, 'files': manifest, 'optional_not_located': omitted_optional, 'no_whole_tree_copy': True}, indent=2) + '\n')
    for name in ['voxel_decisions.py', 'assess_model.py', 'decision_error.py']:
        shutil.copyfile(ROOT / 'inputs/LANE_X5_DECISION_SEG_ERROR' / name, ROOT / name)
    shutil.copyfile(ROOT / 'inputs/LANE_X26_DECIDABILITY/decidability.py', ROOT / 'decidability.py')
    pred = {'frozen_utc': now, 'design': 'Retrospective known-output replay, no new lab predictions', 'prereg_sha256': sha(ROOT / 'PREREG_R1_MANUSCRIPT.json'), 'source_manifest_sha256': sha(ROOT / 'SOURCE_MANIFEST.json'), 'predictions': {'connected_probe_count': 6, 'minimum_added_voxels': [11, 4, 8, 3, 9, 2], 'conditional_class_counts': {'1': 3648, '2': 484}, 'physical_certificates': 0}, 'future_physical_prediction': 'NONE; no independent physical measurement available'}
    (ROOT / 'FROZEN_PREDICTIONS.json').write_text(json.dumps(pred, indent=2) + '\n')
    (ROOT / 'FROZEN_PREDICTIONS.sha256').write_text(sha(ROOT / 'FROZEN_PREDICTIONS.json') + '\n')
    print(json.dumps({'files': len(manifest), 'bytes': sum((x['bytes'] for x in manifest)), 'frozen_utc': now, 'optional_not_located': omitted_optional}))
if __name__ == '__main__':
    main()
