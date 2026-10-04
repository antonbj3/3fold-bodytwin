"""Future measured-pressure consumer; never fills absent patch forces."""
import argparse, csv
from geometry import *
from force_simplex import peak

def compare(path):
    frozen = json.loads((H / 'FROZEN_PREDICTIONS_R3.json').read_text())
    result = json.loads((H / 'rounds/R3.json').read_text())
    assert sha(H / 'rounds/R3.json') == frozen['result_sha256']
    measured = list(csv.DictReader(Path(path).open()))
    group = {}
    for r in measured:
        group.setdefault((int(r['case']), int(r['fdi'])), []).append(r)
    outputs = []
    for (pair, rr) in group.items():
        row = next((r for r in result['rows'] if (r['case'], r['fdi']) == pair and 'basis_file' in r))
        frozen_file = next((f for f in frozen['files'] if f['path'] == row['basis_file']))
        assert sha(row['basis_file']) == frozen_file['sha256']
        assert all((r['basis_sha256'] == frozen_file['sha256'] for r in rr)), 'Measurement-to-prediction hash mismatch'
        ids = [int(r['patch_id']) for r in rr]
        assert sorted(ids) == list(range(row['patch_count'])), 'Missing/duplicate pressure patch channel'
        values = np.array([float(r['force_N']) for r in sorted(rr, key=lambda r: int(r['patch_id']))])
        assert np.isfinite(values).all() and np.all(values >= 0) and (values.sum() > 0), 'Invalid measured force'
        basis = np.load(row['basis_file'])['stress_tensors_MPa']
        tensor = np.einsum('j,jea->ea', values / 100, basis)
        outputs.append(dict(case=pair[0], fdi=pair[1], total_measured_force_N=float(values.sum()), conditional_predicted_peak_MPa=peak(tensor), force_source_sha256=sha(path), physical_stress_accuracy='UNKNOWN until independent strain/force-displacement observation'))
    assert outputs, 'Empty measurement is NOT_RUN'
    dump(H / 'raw/PRESSURE_CONSUMER_RESULT.json', dict(input_file=str(path), input_sha256=sha(path), rows=outputs, fit='none; frozen linear tensor basis', measurement_provenance='User/operator must supply independently calibrated pressure with registered patch mapping'))
if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('csv', type=Path)
    a = p.parse_args()
    compare(a.csv)
