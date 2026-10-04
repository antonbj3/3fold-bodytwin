"""Check the consumed local artifacts when a package is moved; leave external data read-only."""
from dental_release.paths import expand as _release_expand
import argparse, hashlib, json
from pathlib import Path
SOURCE = Path(_release_expand('@DENTAL_IMPLEMENTATIONS@/X12'))

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def resolve(p, package):
    p = Path(p)
    if p.is_relative_to(SOURCE):
        return package / p.relative_to(SOURCE)
    return p if p.is_absolute() else package / p

def verify(package):
    r = json.loads((package / 'results.json').read_text())
    rows = r['numerical_artifact_manifest'] + json.loads((package / 'raw/R8_manifest.json').read_text())
    failures = []
    for item in rows:
        p = resolve(item['path'], package)
        if not p.is_file() or sha(p) != item['sha256']:
            failures.append(str(p))
    return failures
if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--package', type=Path, required=True)
    args = p.parse_args()
    bad = verify(args.package.resolve())
    print(json.dumps({'pass': not bad, 'failed': bad}, indent=2))
    raise SystemExit(bool(bad))
