from common import *
import platform, sys
out = []
for p in sorted(ROOT.rglob('*')):
    if not p.is_file() or '__pycache__' in p.parts:
        continue
    if p.name in ['ARTIFACT_MANIFEST.json', 'CURRENT_WORK_STATE.json'] or p.suffix == '.log' or 'RESOURCES' in p.name:
        continue
    out.append(dict(file=str(p.relative_to(ROOT)), bytes=p.stat().st_size, sha256=sha(p)))
write(ROOT / 'ARTIFACT_MANIFEST.json', dict(sealed_utc=now(), files=out, excluded='live state; manifest itself; runtime logs/resource records that finish after seal', environment=dict(python=sys.version, platform=platform.platform(), numpy=np.__version__), review_state='PENDING_INDEPENDENT_REVIEW'))
s = read(ROOT / 'CURRENT_WORK_STATE.json')
s.update(artifact_manifest_sha256=sha(ROOT / 'ARTIFACT_MANIFEST.json'), one_command='PASS', verification_file='raw/VERIFICATION.json')
write(ROOT / 'CURRENT_WORK_STATE.json', s)
print('Sealed', len(out), 'artifacts')
