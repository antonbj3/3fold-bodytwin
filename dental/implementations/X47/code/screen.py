"""Offline, one-pass targeted screening; no source tree copy."""
from dental_release.paths import expand as _release_expand
from pathlib import Path
import json, re, subprocess, hashlib, xml.etree.ElementTree as ET
ROOT = Path(__file__).resolve().parents[1]
CORPUS = Path(_release_expand('@DENTAL_CORPUS_ROOT@/europepmc/fulltext'))

def txt(e):
    return ' '.join(''.join(e.itertext()).split()) if e is not None else ''

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    keyword = 'retention|dislodg[e]?ment|pull.off|tensile bond|thermocycling|aging|convergence angle|preparation height'
    p = subprocess.run(['rg', '-l', '-i', keyword, str(CORPUS)], check=True, capture_output=True, text=True)
    rows = []
    for f in sorted(p.stdout.splitlines()):
        path = Path(f)
        try:
            tree = ET.parse(path)
        except ET.ParseError:
            rows.append(dict(pmc=path.stem, path=f, decision='REJECT', reason='xml_parse_error'))
            continue
        title = txt(tree.find('.//article-title'))
        abstract = ' '.join((txt(x) for x in tree.findall('.//abstract')))
        ids = {x.get('pub-id-type'): txt(x) for x in tree.findall('./front/article-meta/article-id')}
        crown = bool(re.search('crown|coping|abutment|fixed dental prosthes', title, re.I))
        endpoint = bool(re.search('retent|pull.off|dislodg|debond|tensile bond', title + ' ' + abstract, re.I))
        if not crown:
            reason = 'title_not_whole_crown_or_abutment'
        elif not endpoint:
            reason = 'no_retention_endpoint_in_title_or_abstract'
        elif re.search('review|meta.analysis|clinical|patients|\\bFEA\\b|finite.element|correction|retrospective|case.series|shielding|application of', title, re.I):
            reason = 'review_clinical_FE_or_correction'
        elif re.search('bracket|veneer|denture|endocrown|splinted', title, re.I):
            reason = 'different_restoration_endpoint'
        else:
            reason = 'READ_PRIMARY_TABLES'
        rows.append(dict(pmc=path.stem, path=f, title=title, doi=ids.get('doi'), decision='READ' if reason == 'READ_PRIMARY_TABLES' else 'REJECT', reason=reason))
    (ROOT / 'raw/SCREENING.json').write_text(json.dumps({'all_xml_files': len(list(CORPUS.glob('*.xml'))), 'keyword_files': len(rows), 'records': rows}, indent=2) + '\n')
    keep = [x for x in rows if x['decision'] == 'READ']
    for r in keep:
        tree = ET.parse(r['path'])
        sections = []
        for tab in tree.findall('.//table-wrap'):
            lines = [tab.get('id', ''), txt(tab.find('label')), txt(tab.find('caption'))]
            for row in tab.findall('.//tr'):
                lines.append(' | '.join((txt(x) for x in row if x.tag in ['td', 'th'])))
            lines += [txt(x) for x in tab.findall('.//table-wrap-foot')]
            sections.append('\n'.join(lines))
        methods = []
        for para in tree.findall('./body//p'):
            s = txt(para)
            if re.search('\\bmm\\b|height|taper|convergence|thermocycl|aging|ageing|storage|stored|water|cement|sample|specimen|crosshead', s, re.I):
                methods.append(s)
        license_text = ' '.join((txt(x) for x in tree.findall('.//permissions')))
        r.update(sha256=sha(Path(r['path'])), license=license_text)
        (ROOT / 'raw' / f"{r['pmc']}_READ.txt").write_text(r['title'] + '\nDOI ' + str(r['doi']) + '\nLICENSE ' + license_text + '\n\nTABLES\n' + '\n\n'.join(sections) + '\n\nBODY PARAGRAPHS\n' + '\n\n'.join(methods))
    (ROOT / 'raw/SOURCE_MANIFEST.json').write_text(json.dumps(keep, indent=2) + '\n')
    print(f"XML={len(list(CORPUS.glob('*.xml')))} keyword={len(rows)} primary_table_reads={len(keep)}")
    for r in keep:
        print(r['pmc'], r['title'])
if __name__ == '__main__':
    main()
