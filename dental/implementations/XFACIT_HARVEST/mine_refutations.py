from dental_release.paths import expand as _release_expand
from pathlib import Path
import json, re, hashlib, datetime, collections, time
STORE = Path(_release_expand('@DENTAL_EXTERNAL_ROOT@/storage/research/bunny48_20260926/dental'))
OUT = Path(__file__).resolve().parent
NUM = re.compile('(?<![A-Za-z])[-+]?\\d+(?:\\.\\d+)?(?:[eE][-+]?\\d+)?')
start = time.monotonic()
rows = []
manifest = []
bad = []
counts = collections.Counter()
paths = sorted(STORE.glob('*/results.json'))
paths += sorted(Path(_release_expand('@DENTAL_EXTERNAL_ROOT@/media/sdc1-tmp/dental_bunny48_archive')).glob('*/results.json'))
seen = set()
for f in paths:
    resolved = str(f.resolve())
    if resolved in seen:
        counts['duplicate_paths'] += 1
        continue
    seen.add(resolved)
    try:
        b = f.read_bytes()
        d = json.loads(b)
    except Exception as e:
        bad.append({'path': str(f), 'error': str(e)})
        continue
    counts['readable_reports'] += 1
    ext = d.get('external_referent')
    if not ext:
        counts['without_external_referent'] += 1
        continue
    if not isinstance(ext, (list, dict)):
        counts['malformed_external_referent'] += 1
        continue
    exts = ext if isinstance(ext, list) else [ext]
    for (j, x) in enumerate(exts):
        if not isinstance(x, dict):
            counts['malformed_external_referent'] += 1
            continue
        counts['referents'] += 1
        counts['refutes_us_true'] += x.get('refutes_us') is True
        s = json.dumps(x, ensure_ascii=False)
        loc = str(x.get('locator', ''))
        score = 20 * bool(re.search('10\\.\\d{4,9}/|PMID|PMC\\d+|pubmed', loc, re.I)) + 10 * bool(re.search('table|table|fig|abstract|section', loc, re.I)) + 8 * bool(re.search('MPa|GPa|µm|μm|mm|Ncm|N cm|°C|rpm|Hz|\\bN\\b|%', s)) + 5 * bool(x.get('refutes_us') is True) + 3 * bool(x.get('our_value'))
        row = {'candidate_id': f.parent.name + (':' + str(j) if len(exts) > 1 else ''), 'job': f.parent.name, 'path': str(f), 'resolved_path': resolved, 'source_sha256': hashlib.sha256(b).hexdigest(), 'target_id': d.get('target_id'), 'claim_type': d.get('claim_type'), 'decision': d.get('decision'), 'external_referent': x, 'rank_score': score, 'fields_missing': [k for k in ['value', 'unit', 'population', 'protocol'] if k not in x]}
        rows.append(row)
        manifest.append({'path': str(f), 'sha256': row['source_sha256'], 'bytes': len(b)})
rows.sort(key=lambda r: (-r['rank_score'], r['candidate_id']))
for (rank, r) in enumerate(rows, 1):
    r['rank'] = rank
(OUT / 'ALL_EXTERNAL_REFERENTS.jsonl').write_text(''.join((json.dumps(r, ensure_ascii=False) + '\n' for r in rows)))
(OUT / 'CORPUS_MANIFEST.json').write_text(json.dumps({'at': datetime.datetime.now(datetime.timezone.utc).isoformat(), 'counts': dict(counts), 'bad': bad, 'files': manifest, 'scan_wall_s': time.monotonic() - start}, indent=2))
print(json.dumps({'counts': dict(counts), 'bad': len(bad), 'scan_s': time.monotonic() - start}))
for r in rows[:70]:
    print(r['rank'], r['job'][-14:], r['rank_score'], str(r['target_id'])[:40], json.dumps(r['external_referent'], ensure_ascii=False)[:1250])
