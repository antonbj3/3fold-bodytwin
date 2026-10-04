from common import *
from model_fit import loadmodel
from field_adapter import generate_mesh, export_research_mesh
import trimesh

def run():
    key = split()['dev_round1'][0]
    config = json.load(open(ROOT / 'FINAL_CONFIG.json'))
    b = config['branch_methods']['x42_contact_branch']
    base = DATA / 'FINAL_MODELS'
    rows = []
    dest = ROOT / 'examples'
    dest.mkdir(exist_ok=True)
    for orig in tasks(key):
        if orig['status'] != 'READY' or orig['level'] != 'normal':
            continue
        t = scene(orig)
        s = loadmodel(base / b['shape'] / (t['family'] + '.npz'))
        c = loadmodel(base / b['contact'] / (t['family'] + '.npz'))
        port = generate_mesh(t, s, c, b['beta'])
        if port['status'] != 'DESIGN':
            continue
        p = dest / (t['family'] + '_research_roof.stl')
        export_research_mesh(port, p)
        m = trimesh.load(p, process=True)
        status = all_checks(t, submit(t, port['outer_vertices'][:, 2], port['inner_vertices'][:, 2]))
        rows.append(dict(task_id=t['task_id'], path=str(p.relative_to(ROOT)), sha256=sha(p), bytes=p.stat().st_size, triangles=len(port['faces']), watertight=bool(m.is_watertight), winding_consistent=bool(m.is_winding_consistent), roof_check=status['validity'], complete_restoration='UNKNOWN', field_kernel_execution='NOT_RUN'))
    dump(ROOT / 'raw/EXPORT_RECEIPT.json', dict(rows=rows, source='public dev input only', frame='each local site in mm', resolution='PER_POINT', timescale='HANDOVER', physical_specimen_validation='NOT_RUN'))
    print('research shell exports', len(rows), 'native field kernel NOT_RUN')
if __name__ == '__main__':
    run()
