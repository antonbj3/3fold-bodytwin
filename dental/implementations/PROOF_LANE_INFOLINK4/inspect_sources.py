from dental_release.paths import expand as _release_expand
from pathlib import Path
from lxml import etree as E
import json, re, subprocess, time, hashlib
H = Path(__file__).resolve().parent
R = H.parent.parent
P = Path(_release_expand('@DENTAL_CORPUS_ROOT@/europepmc/fulltext'))
cs = json.loads((H / 'CANDIDATES_FROZEN.json').read_text())
ids = [x['source'] for x in cs if x['source'].startswith('PMC')]
(H / 'source_review').mkdir(exist_ok=True)
for sid in ids:
    p = P / (sid + '.xml')
    root = E.parse(str(p))
    lines = []
    for xp in ['//article-title', '//article-id', '//permissions', '//abstract', '//sec[title[contains(.,"Method") or contains(.,"Material")]]', '//table-wrap']:
        for e in root.xpath(xp):
            if e.tag == 'table-wrap':
                lines.append('TABLE ' + str(e.get('id')) + ' ' + re.sub('\\s+', ' ', ' '.join(e.xpath('./caption//text()'))))
                for (i, tr) in enumerate(e.xpath('.//tr'), 1):
                    lines.append(f'  tr[{i}] ' + ' | '.join((re.sub('\\s+', ' ', ' '.join(x.itertext())).strip() for x in tr)))
                lines.append('Footnotes: ' + ' '.join(e.xpath('./table-wrap-foot//text()')))
            else:
                lines.append(E.ElementTree(root.getroot()).getpath(e) + ' ' + re.sub('\\s+', ' ', ' '.join(e.itertext())))
    (H / 'source_review' / f'{sid}.txt').write_text('\n\n'.join(lines))
paths = []
for base in [R / 'results', Path(_release_expand('@DENTAL_EXTERNAL_ROOT@/coupled-model/results')), Path(_release_expand('@DENTAL_EXTERNAL_ROOT@/storage/research/bunny48_20260926/dental'))]:
    pattern = 'LANE_*/RESULTS.md' if base.parent.name == 'bodytwin' else '*/RESULTS.md'
    paths.extend((str(p) for p in base.glob(pattern) if p.is_file()))
paths.extend((str(p) for p in (R / 'cells').glob('**/*.py')))
pat = '|'.join(ids + ['g1335', 'ohbinding', 'retention.py', 'tail_family_cert'])
t = time.perf_counter()
rs = []
for i in range(0, len(paths), 400):
    a = subprocess.run(['rg', '-n', '-i', pat, *paths[i:i + 400]], capture_output=True, text=True)
    rs.append(a.stdout)
    if a.returncode not in (0, 1):
        rs.append(a.stderr)
(H / 'PRIOR_WORK_SEARCH.txt').write_text(''.join(rs))
(H / 'SEARCH_MANIFEST.json').write_text(json.dumps({'paths_scanned': len(paths), 'search_pattern': pat, 'seconds': time.perf_counter() - t, 'roots': 'targeted one-level result globs, cells/**/*.py, LANE_*/RESULTS.md; no archive recursion or copies', 'output_bytes': sum((len(x) for x in rs))}, indent=2))
print('Primary dumps:', len(ids), 'Search files:', len(paths), 'search result bytes:', sum((len(x) for x in rs)))
