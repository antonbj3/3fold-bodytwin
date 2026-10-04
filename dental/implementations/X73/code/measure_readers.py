from common import *
from blinded_port import compare_readers
import argparse

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('submission')
    ap.add_argument('--out', default=str(ROOT / 'raw/READERS_RESULT.json'))
    a = ap.parse_args()
    s = json.loads(Path(a.submission).read_text())
    fp = ROOT / 'FROZEN_PREDICTIONS_R5_LOCAL_CONTOURS.json'
    assert sha(fp) == fp.with_suffix('.sha256').read_text().split()[0] == s['frozen_predictions_sha256']
    rois = {r['roi_id']: r for r in json.loads((ROOT / 'raw/BLINDED_ROI_INPUT.json').read_text())['rois']}
    out = compare_readers(s['records'], rois)
    dump(a.out, out)
    print(json.dumps(out, indent=2))
if __name__ == '__main__':
    main()
