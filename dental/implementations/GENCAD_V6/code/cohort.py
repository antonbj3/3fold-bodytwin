from dental_release.paths import expand as _release_expand
from common import *

def crowns():
    old = read(V5 / 'raw/FUNCTIONAL_ROWS.json') + read(V5 / 'raw/R6_EXTERNAL_SUPPORT.json')['rows']
    adapted = {r['uid']: r for r in read(V5B / 'raw/R7_CHECKPOINT.json')}
    out = []
    for r in old:
        uid = r['track'] + '__' + r['participant'] + '__' + r['key']
        a = adapted.get(uid, {})
        out.append(dict(uid='V5B__' + uid, key=r['key'], participant=r['participant'], track='V5B_' + r['track'], family=r.get('family'), status='AVAILABLE' if a.get('mesh_path') else 'NO_MESH', mesh_path=a.get('mesh_path'), parent_geometry_pass=a.get('geometry_pass', False), source_status=a.get('status', r.get('status'))))
    for tag in ['R2', 'R3']:
        for r in read(FC / 'rounds' / f'{tag}.json')['rows']:
            if r.get('kind') != 'shell':
                continue
            path = DATA.parent / _release_expand('FULL_CROWN') / (tag + '_predictions') / r['participant'] / r['key'] / 'mesh.npz'
            out.append(dict(uid='FULL__' + tag + '__' + r['participant'] + '__' + r['key'], key=r['key'], participant=r['participant'], track='FULL_' + tag, family=r['family'], status='AVAILABLE' if r['status'] == 'SCORED' and path.exists() else 'NO_MESH', mesh_path=str(path) if path.exists() else None, parent_geometry_pass=r.get('digital_closed_shell', False), source_status=r['status']))
    return out
