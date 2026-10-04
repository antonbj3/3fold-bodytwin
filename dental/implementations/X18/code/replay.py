"""Preserve freeze files before intentional design reproduction, no source cache."""
import shutil, datetime, sys
from geometry import H, D, dump
folder = H / 'raw/replays' / datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%SZ')
folder.mkdir(parents=True)
for name in ['FROZEN_PREDICTIONS_R2.json', 'FROZEN_PREDICTIONS_R3.json', 'FROZEN_PREDICTIONS_R4.json', 'raw/PREDICTIONS_R2.json', 'raw/PREDICTIONS_R4.json', 'raw/PREDICT_R2_COST.json', 'rounds/R1.json', 'rounds/R2.json', 'rounds/R3.json', 'rounds/R4.json', 'results.json', 'GRAPH_FEEDBACK.json', 'raw/PREDICT_COST.json', 'raw/EVALUATE_R1.log']:
    p = H / name
    if p.exists():
        shutil.copy2(p, folder / p.name)
dump(folder / 'README.json', dict(scope='Reproduction from previously observed sources, not a new blind trial. Original design/prediction files hashes preserved; generated arrays will be deterministically overwritten in data directory.', data_root=str(D)))
if len(sys.argv) > 1 and sys.argv[1] == '--design':
    for name in ['raw/PREDICTIONS_R2_checkpoint.json', 'raw/PREDICTIONS_R4_checkpoint.json']:
        p = H / name
        if p.exists():
            shutil.copy2(p, folder / p.name)
            p.unlink()
