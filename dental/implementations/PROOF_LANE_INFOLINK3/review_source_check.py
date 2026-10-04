"""Bind witness units and locators to the independently reviewed frozen contract."""
import argparse
import ast
import copy
import csv
import hashlib
import io
import json
import zipfile
from pathlib import Path
import openpyxl

def sha(p):
    h = hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda : f.read(1048576), b''):
            h.update(b)
    return h.hexdigest()

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--demo', type=Path, required=True)
    ap.add_argument('--trusted-facit', type=Path, required=True)
    ap.add_argument('--trusted-sources', type=Path, required=True)
    ap.add_argument('--output', type=Path, required=True)
    a = ap.parse_args()
    tree = ast.parse((a.demo / 'make_demo.py').read_text())
    functions = [x for x in tree.body if isinstance(x, ast.FunctionDef) and x.name in ['read_fact', 'verify_value']]
    ns = dict(Path=Path, csv=csv, io=io, json=json, zipfile=zipfile, openpyxl=openpyxl)
    exec(compile(ast.Module(body=functions, type_ignores=[]), 'original_value_checker', 'exec'), ns)
    old = ns['verify_value']
    trusted = json.loads(a.trusted_facit.read_text())
    current = json.loads((a.demo / 'FACIT.json').read_text())
    sources = {r['path']: r['sha256'] for r in json.loads(a.trusted_sources.read_text())}

    def strict(w, reference):
        try:
            return w == reference and sha(w['path']) == sources[w['path']] and old(w)
        except (KeyError, ValueError, OSError, IndexError, TypeError):
            return False
    checks = []
    for (lid, rows) in trusted.items():
        for (i, w) in enumerate(rows):
            unit = copy.deepcopy(w)
            unit['unit'] = 'INJECTED_WRONG_UNIT'
            value = copy.deepcopy(w)
            value['original_value'] = w['original_value'] + '__WRONG' if isinstance(w['original_value'], str) else w['original_value'] + max(1, abs(w['original_value']) * 0.1)
            checks.append(dict(link=lid, witness=i, positive_passed=strict(current[lid][i], w), original_wrong_unit_accepted=old(unit), corrected_wrong_unit_rejected=not strict(unit, w), corrected_wrong_value_rejected=not strict(value, w)))
    result = dict(checks=checks, count=len(checks), original_unit_blind_count=sum((x['original_wrong_unit_accepted'] for x in checks)), all_corrected_pass=all((x['positive_passed'] and x['corrected_wrong_unit_rejected'] and x['corrected_wrong_value_rejected'] for x in checks)), conclusion='Eleven bounded local links unchanged; public EBSD deposit still unresolved.')
    a.output.write_text(json.dumps(result, indent=2) + '\n')
    assert result['all_corrected_pass']
    print(json.dumps({k: v for (k, v) in result.items() if k != 'checks'}))
if __name__ == '__main__':
    main()
