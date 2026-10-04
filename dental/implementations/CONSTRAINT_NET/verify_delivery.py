"""Check exported binding identities and immutable-source seals."""
from common import *

def main():
    n = read(HERE / 'CONSTRAINT_NET_DENTAL.json')
    b = read(HERE / 'BINDINGS_RECEIPT.json')
    errors = []
    rows = [json.loads(s) for s in Path(b['path']).read_text().splitlines() if s.strip()]
    expected = {'variable': n['variables'], 'edge': {e['id']: e for e in n['edges']}}
    for row in rows:
        if 'constraint_kind' in row:
            obj = expected[row['constraint_kind']].get(row['constraint_id'])
            if obj != row['source_contract']:
                errors.append('Binding contract mismatch: ' + row['id'])
        else:
            obj = next((c for c in n['chains'] if c['id'] == row['northstar_id']), None)
            if obj != row['source_contract']:
                errors.append('Chain contract mismatch: ' + row['id'])
        if row['status'] != 'OPEN' or row['no_scientific_admission'] is not True:
            errors.append('Improper admission')
    if len(rows) != len(n['variables']) + len(n['edges']) + 52:
        errors.append('Binding count')
    if sha(b['path']) != b['sha256']:
        errors.append('Binding file hash')
    for r in MANIFEST['files']:
        if sha(HERE / r['snapshot']) != r['sha256']:
            errors.append('Frozen source: ' + r['source_file'])
    for name in ['PREREG_R1.json', 'PREREG_R2.json', 'FROZEN_PREDICTIONS.json', 'FROZEN_PREDICTIONS_R2.json']:
        if sha(HERE / name) != (HERE / (name + '.sha256')).read_text().strip():
            errors.append('Freeze: ' + name)
    validation = read(HERE / 'VALIDATION.json')
    if not validation['pass'] or not validation['mutation_tests']['pass']:
        errors.append('Core validation')
    if not read(HERE / 'VALIDATION_LIVE_VALUES.json')['pass']:
        errors.append('Current source values')
    if not read(HERE / 'RESULTS_R2.json')['checks_pass']:
        errors.append('R2 check')
    out = {'pass': not errors, 'errors': errors, 'binding_rows': len(rows), 'frozen_source_files': len(MANIFEST['files']), 'sealed_preregs_and_predictions': 4, 'independent_scientific_review': 'PENDING', 'formal_decomposition_timing': 'RETROSPECTIVE; disclosed protocol deviation'}
    dump('DELIVERY_VERIFICATION.json', out)
    print(json.dumps(out))
    if errors:
        raise SystemExit(1)
if __name__ == '__main__':
    main()
