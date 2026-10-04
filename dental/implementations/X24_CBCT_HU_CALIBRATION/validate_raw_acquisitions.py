"""Read the 22 original local acquisitions without changing frozen evidence."""
import datetime
import hashlib
import json
import os
import time
import zipfile
from pathlib import Path
import numpy as np
from measure_phantom import P, SOURCE, old, prepare_rois, sample, sha
from run_r2 import cached_sample, parse_nii

def main():
    started = time.perf_counter()
    baseline = json.loads((P / 'PHANTOM_BASELINE.json').read_text())
    new = json.loads((P / 'PHANTOM_NEW_R2.json').read_text())
    rows = []

    def check(record, image, affine, actual_sha, expected_sha, mapping):
        assert actual_sha == expected_sha, 'Acquisition bytes changed'
        assert list(image.shape) == record['image_shape_zyx']
        assert np.array_equal(affine, np.array(record['image_affine']))
        if mapping == 'original_masks':
            (rois, _) = prepare_rois(SOURCE / 'extracted' / record['device'])
            (observations, raw) = sample(image, affine, rois)
        else:
            source = next((r for r in baseline if r['device'] == record['device']))
            with np.load(source['raw_samples']['path']) as cache:
                (observations, raw) = cached_sample(image, cache)
        expected = {r['material']: r for r in record['inserts']}
        max_error = max((abs(r['gray_median'] - expected[r['material']]['gray_median']) for r in observations))
        with np.load(record['raw_samples']['path']) as reference:
            arrays_match = set(raw) == set(reference.files) and all((np.array_equal(array, reference[key]) for (key, array) in raw.items()))
        assert max_error == 0 and arrays_match, 'Measurement replay changed'
        rows.append({'device': record['device'], 'source': record.get('zip_member', record.get('source')), 'sha256': actual_sha, 'mapping': mapping, 'ROIs': len(observations), 'max_median_error_gray': max_error, 'all_raw_arrays_match': arrays_match})
        print(len(rows), record['device'], 'original acquisition PASS', flush=True)
    for record in baseline:
        path = Path(record['source'])
        (image, affine) = old.read_nii(path)
        check(record, image, affine, sha(path), record['sha256'], 'original_masks')
        del image
    with zipfile.ZipFile(SOURCE / 'cbct_study.zip') as archive:
        for record in new:
            blob = archive.read(record['zip_member'])
            (image, affine) = parse_nii(blob)
            mapping = 'baseline_point_mapping' if record['affine_matches_baseline'] else 'original_masks'
            check(record, image, affine, hashlib.sha256(blob).hexdigest(), record['compressed_member_sha256'], mapping)
            del blob, image
    report = {'status': 'PASS', 'original_local_acquisitions_read': len(rows), 'ROIs_remeasured': sum((r['ROIs'] for r in rows)), 'observations': rows, 'new_physical_acquisitions': 0, 'wall_s': time.perf_counter() - started, 'completed_at': datetime.datetime.now(datetime.timezone.utc).isoformat()}
    out = Path(os.environ.get('X24_OUTPUT_DIR', str(P / 'demo_outputs')))
    out.mkdir(parents=True, exist_ok=True)
    (out / 'RAW_ACQUISITION_REPLAY.json').write_text(json.dumps(report, indent=2) + '\n')
    print('PASS', len(rows), 'original acquisitions /', report['ROIs_remeasured'], 'ROIs')
if __name__ == '__main__':
    main()
