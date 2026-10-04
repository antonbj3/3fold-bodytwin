"""Integrity and scope checks against stored observations, no refit."""
import ast, json, hashlib, pathlib, datetime
R = pathlib.Path(__file__).resolve().parents[1]

def sha(p):
    return hashlib.sha256(pathlib.Path(p).read_bytes()).hexdigest()

def read(p):
    return json.loads(pathlib.Path(p).read_text())

def run():
    for p in (R / 'code').glob('*.py'):
        ast.parse(p.read_text())
    freezes = []
    for p in R.glob('*.sha256'):
        (h, name) = p.read_text().strip().split(None, 1)
        name = name.strip()
        valid = sha(R / name) == h
        assert valid, (p, name)
        freezes.append(dict(file=name, valid=valid))
    chain = '0' * 64
    stream = []
    for line in (R / 'FROZEN_PREDICTIONS_D_STREAM.jsonl').read_text().splitlines():
        r = json.loads(line)
        assert r['previous_record_sha256'] == chain
        chain = hashlib.sha256(line.encode()).hexdigest()
        stream.append((r['row']['key'], chain))
    D = read(R / 'RESULTS_D.json')['rows']
    assert len(D) == 18 and sum((r.get('geometric_conjunction', False) for r in D)) == 15
    for r in D:
        assert r['individual_prediction_sha256'] == dict(stream)[r['key']]
        if r['status'] == 'GENERATED':
            assert sha(r['mesh_path']) == r['mesh_sha256']
    A2 = read(R / 'RESULTS_A2.json')['rows']
    assert len(A2) == 18
    for r in A2:
        assert r['apex']['exact_minimax_residual_mm'] == '1/20'
        assert r['native_wall_witness']['exact']['distance_squared_mm2'] == '0'
        assert r['generated_wall_witness']['exact']['distance_squared_mm2'] == '0'
    for tag in ['E', 'G']:
        assert read(R / f'RESULTS_{tag}.json')['admitted'] == 0
    S = read(R / 'RESULTS_SCOPE.json')['rows']
    assert len(S) == 15
    for r in S:
        assert 0 < r['virtual_preparation_volume_mm3'] < r['occupied_cavity_volume_mm3'] < r['support_volume_mm3']
    for name in ['raw/SOURCE_HASH_CHECK.json', 'raw/VALIDATION.json', 'raw/GEOMETRY_FAULT_CONTROLS.json']:
        assert read(R / name)['all_pass']
    result = read(R / 'results.json')
    assert result['cohort']['full_clinical_qualified'] == 0 and result['cohort']['full_clinical_status'] == 'UNKNOWN'
    out = dict(all_pass=True, utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), sha_sidecars=freezes, stream_records=len(stream), generated_mesh_hashes_verified=15, round_count=len(result['rounds']), scope='Stored record integrity and consistency; independent scientific review still pending')
    (R / 'raw/PACKAGE_AUDIT.json').write_text(json.dumps(out, indent=2) + '\n')
    print('Package audit PASS:18 stream records,15 mesh hashes,all freeze sidecars and retained failures')
if __name__ == '__main__':
    run()
