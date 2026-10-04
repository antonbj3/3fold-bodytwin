"""Reproduce the local title screen; this is not a fulltext evidence admission."""
from dental_release.paths import expand as _release_expand
from pathlib import Path
from lxml import etree
from common import ROOT, dump
base = Path(_release_expand('@DENTAL_CORPUS_ROOT@/europepmc/fulltext'))
selected = []
count = 0
for path in base.glob('*.xml'):
    count += 1
    try:
        tree = etree.parse(str(path))
        title = ' '.join(tree.xpath('//front//article-title')[0].itertext())
        low = title.lower()
        if any((k in low for k in ['color', 'colour', 'translucen', 'bond strength', 'surface treatment', 'abrasion', 'primer'])) and any((k in low for k in ['zircon', 'ceram', 'cement', 'crown', 'silicate'])):
            selected.append({'path': str(path), 'title': title, 'doi': tree.xpath('//article-id[@pub-id-type="doi"]/text()')})
    except Exception:
        pass
dump(ROOT / 'raw/TITLE_SCREEN_REPLAY.json', {'files_screened': count, 'title_selected': selected})
print('Title screen only:', count, 'files,', len(selected), 'hits; no physical evidence accepted.')
