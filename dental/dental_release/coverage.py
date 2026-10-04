"""Validate per-delivery scope identity against the independent frozen census."""
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]

def validate(catalogue=None):
    catalogue=catalogue or json.loads((ROOT/'demos.json').read_text())['demos']
    expected=json.loads((ROOT/'provenance/EXPECTED_COVERAGE.json').read_text())
    entries={r['id']:r for r in catalogue}
    if len(entries)!=len(catalogue):raise ValueError('Duplicate release identity')
    kept=0
    for row in expected['deliveries']:
        if row.get('superseded'):
            if row['replacement_entry'] not in entries:raise ValueError('Missing replacement')
            if entries[row['replacement_entry']]['reviewed_result_sha256']!=row['replacement_result_sha256']:raise ValueError('Changed replacement hash')
            continue
        valid=False
        for ident in row['release_entries']:
            if ident not in entries:continue
            e=entries[ident]
            if e.get('reviewed_result_sha256') not in row['accepted_result_sha256']:continue
            if 'document' in e:
                valid=(ROOT/e['document']).is_file()
            else:
                valid=(ROOT/'implementations'/ident/'README.md').is_file() and bool(e['code_files'])
            if valid:break
        if not valid:raise ValueError('Missing bound delivery: '+row['job'])
        kept+=1
    return dict(deliveries=len(expected['deliveries']),direct=kept,superseded=len(expected['deliveries'])-kept,missing=0,scope='Reviewed source/document coverage, not physical validity')
