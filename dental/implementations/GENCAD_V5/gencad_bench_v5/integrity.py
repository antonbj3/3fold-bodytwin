from common import *

def verify_blob(blob, expected):
    return hashlib.sha256(blob).hexdigest() == expected

def input_files():
    out = [V4 / 'results.json', V4 / 'raw/QUALITY_ROWS.csv.gz', V4 / 'FROZEN_LITERATURE_PREDICTIONS.json', V4 / 'payload/literature/CEMENT.csv', V4 / 'payload/literature/CROWN_CURATED.json', V4 / 'payload/private/COHORT.json', V4 / 'payload/whole_inputs/RECORDS.json', ROOT.parent / 'LANE_XREVIEW_BATCH10/CORRECTIONS.jsonl', ROOT.parent / 'LANE_XREVIEW_BATCH10/patches/PROOF_LANE_GENCAD_V4/review_controls.py', ROOT.parent / 'LANE_X60_CROWN_OPTIMIZER/raw/SCORED_ROWS.json', ROOT.parent / 'LANE_X1B_CROWN_LOOP/raw/MATCHED_R2.json', ROOT.parent / 'LANE_X1B_CROWN_LOOP/raw/CALIBRATION_R1.json']
    for folder in ['whole_inputs', 'whole_private', 'whole_predictions', 'prep_predictions', 'public/tasks', 'public/scenes']:
        out.extend((p for p in (V4 / 'payload' / folder).rglob('*') if p.is_file() and p.suffix in ['.json', '.npz']))
    out.extend((p for p in (V4 / 'payload/participant_code').rglob('*.py')))
    out.append(V4 / 'code/quality.py')
    out.extend((DATA / 'external_inputs' / p.name for p in (DATA / 'external_inputs').glob('*.npz')))
    out.extend((DATA / 'external_predictions' / p.name for p in (DATA / 'external_predictions').glob('*.npz')))
    return sorted(set(out))

def seal():
    freeze(ROOT / 'INPUT_LOCK.json', dict(scope='Post-development source pin; original predictions separately frozen before scoring. No prospective input-discovery claim', files={str(p): dict(sha256=sha(p), bytes=p.stat().st_size) for p in input_files()}))

def verify():
    files = read(ROOT / 'INPUT_LOCK.json')['files']
    bad = [p for (p, m) in files.items() if not Path(p).is_file() or sha(p) != m['sha256']]
    if bad:
        raise ValueError('Input hash mismatch ' + repr(bad[:5]))
    frozen_files = sorted({p.name for pattern in ['PREREG_*.json', 'FROZEN_*.json', 'EXTERNAL_INPUTS.json', 'EXTERNAL_GENERATOR_R3.json', 'INPUT_LOCK.json'] for p in ROOT.glob(pattern)})
    for f in frozen_files:
        expected = (ROOT / Path(f).with_suffix('.sha256')).read_text().split()[0]
        if sha(ROOT / f) != expected:
            raise ValueError('Frozen contract drift ' + f)
    model = read(ROOT / 'raw/EXTERNAL_MODEL_LOCK.json')
    model_files = model['files']
    for m in model_files:
        if sha(m['path']) != m['sha256']:
            raise ValueError('External model/source drift ' + m['path'])
    old = read(ROOT / 'history/external_100_steps/FROZEN_PREDICTIONS.json')
    for (name, m) in old['files'].items():
        if sha(DATA / 'external_predictions_100' / name) != m['sha256']:
            raise ValueError('Retained100-step prediction drift ' + name)
    blob = Path(next(iter(files))).read_bytes()
    valid = verify_blob(blob, hashlib.sha256(blob).hexdigest())
    mutated = verify_blob(blob + b'\nINJECTED', hashlib.sha256(blob).hexdigest())
    assert valid and (not mutated)
    return dict(status='PASS', files=len(files), bytes=sum((m['bytes'] for m in files.values())), frozen_contracts=len(frozen_files), model_source_files=len(model_files), model_source_bytes=sum((m['bytes'] for m in model_files)), retained100_predictions=len(old['files']), injected_bad_bytes_rejected=not mutated)
if __name__ == '__main__':
    if '--seal' in sys.argv:
        seal()
    else:
        print(verify())
