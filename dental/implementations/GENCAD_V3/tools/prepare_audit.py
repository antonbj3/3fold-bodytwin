"""Capture small source-linked triangle witnesses for independent replay."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'code'))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'code/legacy/vendor'))
from util import *
from source_ops import pair, make_site
sel = read(PAYLOAD / 'private/PREP_SELECTION.json')['payload']
lookup = {r['case_key']: r for r in read(PAYLOAD / 'private/COHORT.json')['payload']['cases']}
out = []
for s in sel:
    row = lookup[s['case_key']]
    dest = PAYLOAD / 'private/audit'
    dest.mkdir(exist_ok=True)
    try:
        (p, sources) = pair(row)
        site = make_site(p, 'molar_crown', compute_obstacles=False)
        a = p['lower']
        idx = np.flatnonzero(a['owner'] == 36)
        original = a['v'][a['f'][idx]]
        np.savez_compressed(dest / (row['case_key'] + '.npz'), original_world_triangles=original, source_face_ids=idx, base=site['base'], R=site['R'], xy=site['xy'], reference=site['reference'])
        out.append(dict(case_key=row['case_key'], status='CAPTURED', sources=sources))
    except (KeyError, ValueError) as e:
        out.append(dict(case_key=row['case_key'], status='UNAVAILABLE', reason=str(e)))
dump(PAYLOAD / 'private/audit/PROVENANCE.json', out)
