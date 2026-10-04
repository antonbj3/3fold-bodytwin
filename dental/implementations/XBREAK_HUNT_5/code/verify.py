from common import *
from source_decode import decode
from columns import distance, continuous_gap
import trimesh
inputs = json.loads((ROOT / 'sources/INPUTS.json').read_text())
records = []
for source in inputs['test_files']:
    if sha(source['file']) != source['sha256']:
        raise ValueError('source hash drift')
    case = source['name'].split('/')[-1].split('_')[0]
    (v, f, lab) = decode(source['file'])
    z = np.load(DATA / (case + '_R1.npz'))
    if not (np.array_equal(v, z['v']) and np.array_equal(f, z['f']) and np.array_equal(lab, z['lab'])):
        raise ValueError('source decoding/array mismatch')
for prereg in sorted(ROOT.glob('PREREG_R[1-4].sha256')):
    for line in prereg.read_text().splitlines():
        (expected, name) = line.split(maxsplit=1)
        if sha(ROOT / name.strip()) != expected:
            raise ValueError('frozen prereg drift')
for line in (ROOT / 'PREREG_R1_ADDENDUM.sha256').read_text().splitlines():
    (expected, name) = line.split(maxsplit=1)
    if sha(ROOT / name.strip()) != expected:
        raise ValueError('pre-run addendum drift')
r1 = json.loads((ROOT / 'raw/R1.json').read_text())
r2 = json.loads((ROOT / 'raw/R2.json').read_text())
r3 = json.loads((ROOT / 'raw/R3.json').read_text())
r4 = json.loads((ROOT / 'raw/R4.json').read_text())
assert r1['available_n'] == 11 and r1['cone_pass_n'] == 0 and (r1['control_agreement_n'] == 11)
assert r2['validation_gate'] == 'PASS' and r2['sample_pass_n'] == 0
for row in r3['rows']:
    source_arrays = np.load(DATA / (row['case'] + '_R1.npz'))
    source_triangles = source_arrays['v'][source_arrays['f'][source_arrays['lab'] == 1]]
    for grid in row['grids']:
        if sha(grid['stl']) != grid['stl_sha256']:
            raise ValueError('export hash drift')
        mesh = trimesh.load(grid['stl'], force='mesh')
        if not (mesh.is_watertight and mesh.is_winding_consistent):
            raise ValueError('export solid failure')
        saved = grid['gap']
        if not (np.isfinite(saved['lower_max']) and np.isfinite(saved['upper_max']) and (0 <= saved['lower_max'] <= saved['upper_max'])):
            raise ValueError('invalid saved gap enclosure')
        measured = continuous_gap(mesh, source_triangles)
        if measured['gate'] != saved['gate']:
            raise ValueError('recomputed whole-surface gap gate differs')
        if max((abs(measured[k] - saved[k]) for k in ['lower_max', 'upper_max'])) > 1e-05:
            raise ValueError('recomputed whole-surface gap bounds differ')
        point = np.array(grid['gap']['witness'])
        d = float(distance(mesh, point[None, :])[0])
        expected = grid['gap']['lower_max']
        if abs(d - expected) > 1e-05:
            raise ValueError('actual saved witness distance differs')
        records.append({'case': row['case'], 'h': grid['h_source_unit'], 'rerun_witness_gap': d, 'saved_gap': expected, 'error': abs(d - expected)})
cover = json.loads((ROOT / 'raw/EXTERNAL_COVER.json').read_text())
assert cover['gate'] == 'PASS'
assert r3['primary_gate'] == 'CONDITIONAL_PASS' and r4['primary_gate'] == 'CONDITIONAL_PASS'
write('raw/VERIFY.json', {'gate': 'PASS', 'sources_checked': 11, 'export_witnesses': records, 'preregs': 'all original and pre-run addendum hashes match', 'independent_scientific_review': False})
print('PASS: source hashes/arrays, frozen preregs, 6 closed exported cavities and actual distance witnesses; physical gate UNKNOWN')
