import argparse, numpy as np, json
from common import *
from region_field import RegionField
from material_port import resolve, InvalidMaterial

def load():
    d = np.load(D / 'region_field.npz')
    r = json.loads((L / 'raw/K09.json').read_text())
    mapping = {x['pdl_component']: x['tooth_component'] for x in r['pdl_tooth_interfaces']}
    faces = [(d[f'faces_{i}'], None) for i in (1, 2, 3)]
    return RegionField(d['V'], d['T'], d['tag'], d['ec'], d['pc'], faces, mapping)

def run(P=None):
    field = load()
    if P is None:
        ids = np.array([np.flatnonzero(field.tag == r)[0] for r in (1, 2, 3)])
        P = field.V[field.T[ids]].mean(1)
    q = field.query(P)
    registry = json.loads((L / 'MATERIAL_REGISTRY.json').read_text())
    rows = []
    for (i, reg) in enumerate(q['region']):
        key = {1: 'tooth', 2: 'pdl_mandible_unknown', 3: 'bone'}.get(int(reg))
        law = registry.get(key)
        response = dict(point_mm=P[i].tolist(), region=int(reg), body=int(q['body'][i]), tet_witness=int(q['witness'][i]), material_sdf_mm=q['sdf_mm'][i].tolist(), resolution='PER_POINT', timescale='SIMULTANEOUS')
        if law:
            try:
                C = resolve(law, int(reg), {})
                response.update(material_key=key, operator_status='CONSTITUTIVE_CLOSURE', min_stiffness_eigenvalue_MPa=float(np.linalg.eigvalsh(C).min()), debt=law['debt'])
            except InvalidMaterial as e:
                response.update(material_key=key, operator_status='UNKNOWN', consumer_must_stop=True, reason=str(e), debt=law['debt'])
        else:
            response.update(operator_status='EXTERIOR_NO_MATERIAL')
        rows.append(response)
    if P is not None and len(P) == 3:
        write(L / 'raw/DEMO_QUERY.json', rows)
    return rows
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Query source-region SDF, owner and executable material validity; coordinates in mm.')
    parser.add_argument('--point', type=float, nargs=3)
    args = parser.parse_args()
    print(json.dumps(run(np.array([args.point]) if args.point else None), indent=2))
