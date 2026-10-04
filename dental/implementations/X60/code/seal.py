from common import *

def run():
    files = {}
    for pattern in ['code/*.py', 'run_all.sh', 'PREREG*.json', 'INPUT_LOCK.json', 'PANEL_SELECTION.json', 'FROZEN*.json', 'exports/whole*/*']:
        for p in ROOT.glob(pattern):
            if p.is_file():
                files[str(p.relative_to(ROOT))] = dict(sha256=sha(p), bytes=p.stat().st_size)
    dump(ROOT / 'RELEASE_LOCK.json', dict(schema='X60-executable-and-prediction-seal-v1', sealed_utc=now(), scope='Final code/input/prediction reproducibility seal, after numerical results; not a new blind experiment', files=files))
    (ROOT / 'RELEASE_LOCK.sha256').write_text(sha(ROOT / 'RELEASE_LOCK.json') + '  RELEASE_LOCK.json\n')
    print('sealed', len(files), 'files')
if __name__ == '__main__':
    run()
