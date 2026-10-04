"""Author-side release creation. Never run automatically by the demo."""
from common import *

def run():
    dirs = ['code', 'payload/public', 'payload/private', 'payload/participant_code', 'payload/models', 'payload/runtime', 'payload/runner', 'payload/predictions', 'payload/whole_inputs', 'payload/whole_private', 'payload/whole_runner', 'payload/whole_predictions', 'payload/prep_runner', 'payload/prep_predictions', 'payload/literature']
    files = {}
    for d in dirs:
        for p in sorted((ROOT / d).rglob('*')):
            if p.is_file():
                files[str(p.relative_to(ROOT))] = {'sha256': sha(p), 'bytes': p.stat().st_size}
    for p in sorted(ROOT.glob('PREREG_*.json')) + sorted(ROOT.glob('FROZEN_*.json')) + [ROOT / 'SOURCE_REUSE.json', ROOT / 'DECOMPOSITION.json', ROOT / 'run_all.sh', ROOT / 'revisions/EXECUTION_ATTEMPT_1/rounds/R2.json', ROOT / 'raw/GENERATION_PROCESS.json', ROOT / 'raw/WHOLE_GENERATION_PROCESS.json', ROOT / 'raw/PREP_GENERATION_PROCESS.json']:
        files[str(p.relative_to(ROOT))] = {'sha256': sha(p), 'bytes': p.stat().st_size}
    freeze(ROOT / 'RELEASE_LOCK.json', dict(files=files, closed_directories=dirs, threat_model='Inputs, prediction bytes, scorer/config/runtime drift and private-reference isolation. Replacing the trusted release anchor is a different release, not an accepted submission. No defense against a malicious host/root kernel.'))
    (ROOT / 'TRUSTED_RELEASE.sha256').write_text(sha(ROOT / 'RELEASE_LOCK.json') + '  RELEASE_LOCK.json\n')
    (ROOT / 'EXECUTABLES.sha256').write_text(sha(DATA / 'runtime/bin/python3') + '  payload/runtime/bin/python3\n' + sha(ROOT / 'code/integrity.py') + '  code/integrity.py\n')
    print('sealed', len(files), 'files', sum((x['bytes'] for x in files.values())))
if __name__ == '__main__':
    run()
