"""Re-read original table cells: a modified extracted numerical value must fail."""
import json, hashlib, copy, re, os, xml.etree.ElementTree as ET
from pathlib import Path
LANE = Path(__file__).resolve().parents[1]

def check(d):
    source = Path(d['source_path'])
    if source.exists() and os.environ.get('X56_PORTABLE') != '1':
        if hashlib.sha256(source.read_bytes()).hexdigest() != d['source_sha256']:
            return False
    else:
        snapshots = json.loads((LANE / 'raw/table_snapshot_manifest.json').read_text())
        snap = next((s for s in snapshots if s['pmcid'] == d['study']))
        source = LANE / snap['snapshot_path']
        if snap['source_sha256'] != d['source_sha256'] or hashlib.sha256(source.read_bytes()).hexdigest() != snap['snapshot_sha256']:
            return False
    root = ET.parse(source).getroot()
    tb = next((t for t in root.findall('.//table-wrap') if t.get('id') == d['table_id']))
    if not any((a.text == d['study'] for a in root.findall('.//article-id') if a.get('pub-id-type') == 'pmcid')):
        return False
    cell = tb.findall('.//tr')[d['source_cell']['tr'] - 1][d['source_cell']['cell'] - 1]
    s = ''.join(cell.itertext())
    values = [float(v.replace('−', '-')) for v in re.findall('[−-]?\\d+(?:\\.\\d+)?', s)]
    return values[0] == d['motion_um'] and d['unit'] == 'um' and (not d['interface_relative_direct'])

def main():
    rows = json.loads((LANE / 'raw/measurements.json').read_text())
    results = [check(d) for d in rows]
    bad = copy.deepcopy(rows[0])
    bad['motion_um'] *= 1000
    badunit = copy.deepcopy(rows[0])
    badunit['unit'] = 'mm'
    badquantity = copy.deepcopy(rows[0])
    badquantity['interface_relative_direct'] = True
    badlink = copy.deepcopy(rows[0])
    badlink['study'] = 'PMC7150554'
    out = {'rows_checked': len(rows), 'source_matches': sum(results), 'source_mismatches': len(results) - sum(results), 'injected_x1000_rejected': not check(bad), 'injected_wrong_unit_rejected': not check(badunit), 'injected_interface_label_rejected': not check(badquantity), 'permuted_study_link_rejected': not check(badlink), 'scope': 'independently reread exact original table-cell number/hash, or portable verbatim numeric-table snapshot when original absent; label injection rejected by source-defined proxy contract'}
    (LANE / 'raw/source_checks.json').write_text(json.dumps(out, indent=2) + '\n')
    assert all(results) and all((out[k] for k in ['injected_x1000_rejected', 'injected_wrong_unit_rejected', 'injected_interface_label_rejected', 'permuted_study_link_rejected']))
    print(json.dumps(out, indent=2))
if __name__ == '__main__':
    main()
