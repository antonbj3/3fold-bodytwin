from dental_release.paths import expand as _release_expand
from pathlib import Path
from bs4 import BeautifulSoup
import xml.etree.ElementTree as E, json, re
S = Path(_release_expand('@DENTAL_WORK_ROOT@/XFACIT-harvest/sources'))

def clean(el):
    if el is None:
        return ''
    el = E.fromstring(E.tostring(el))
    for parent in el.iter():
        for child in list(parent):
            if child.tag == 'tex-math':
                parent.remove(child)
    return re.sub('\\s+', ' ', ''.join(el.itertext())).strip()
for n in ['PMC11833108', 'PMC8144379', 'PMC10291793', 'PMC3134825']:
    root = E.parse(S / (n + '.xml')).getroot()
    print('\nSOURCE', n)
    for t in root.findall('.//table-wrap'):
        print('TABLE', t.attrib, clean(t.find('caption'))[:200])
        if n == 'PMC11833108' and t.get('id') != 'Tab1':
            continue
        for row in t.findall('.//tr'):
            cells = [clean(c) for c in row]
            line = ' | '.join(cells)
            if n == 'PMC8144379' and (not any((s in line.lower() for s in ['good', 'soballe', 'mean', 'micromotion', 'human']))):
                continue
            if n == 'PMC3134825' and (not any((s in line.lower() for s in ['cattaneo', 'meyer', 'yoshida']))):
                continue
            print(line[:1300])
    for p in root.findall('.//body//p'):
        line = clean(p)
        if n == 'PMC11833108' and any((s in line.lower() for s in ['observed stress', 'test set', 'tested', 'specimen'])):
            print('PAR', p.attrib, line[:2300])
        if n == 'PMC8144379' and any((s in line.lower() for s in ['15 to 750', '0.41', 'mean micromotion'])):
            print('PAR', p.attrib, line[:2100])
