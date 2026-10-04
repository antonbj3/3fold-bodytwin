"""Replace only our own byte-identical outputs with hard links; keep every locator/hash."""
import json, os
from pathlib import Path
from prepare_replay import HERE, DATA, sha, write
seen = {}
saved = 0
n = 0
for manifest in sorted((HERE / 'raw/output_manifests').glob('*.json')):
    for r in json.loads(manifest.read_text()):
        if r.get('kind') != 'new_run_output':
            continue
        p = Path(r['locator'])
        assert p.is_relative_to(DATA)
        if not p.exists():
            raise RuntimeError(str(p))
        key = r['sha256']
        if key in seen:
            canonical = seen[key]
            assert sha(p) == key and sha(canonical) == key
            if p.stat().st_ino != canonical.stat().st_ino:
                size = p.stat().st_size
                temp = p.with_name(p.name + '.dedup_pending')
                os.link(canonical, temp)
                os.replace(temp, p)
                saved += size
                n += 1
        else:
            seen[key] = p
write('raw/OWN_OUTPUT_DEDUPLICATION.json', {'relinked_files': n, 'logical_duplicate_bytes_removed': saved, 'every_original_locator_preserved': True, 'sha256_unchanged': True, 'source_files_modified': 0})
print(n, saved)
