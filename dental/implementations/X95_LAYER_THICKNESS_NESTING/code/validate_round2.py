from pose_join import *
import numpy as np

def run():
    rows = []
    for e in json.loads((ROOT / 'raw/R2_RESULTS.json').read_text())['rows']:
        a = dict(np.load(ROOT / e['array_path'], allow_pickle=False))
        k = a['kept_yml'].copy()
        g = a['guard_min_mm']
        h = float(a['cube_h_mm'])
        alt = k.copy()
        keep = int(np.flatnonzero(k)[0])
        bad = int(np.flatnonzero(~k)[0])
        alt[keep] = False
        alt[bad] = True
        identity = (int(k.sum()) - int(alt.sum())) * h ** 3
        ga = float(g[k].min())
        gb = float(g[alt].min())
        assert identity == 0.0 and ga >= 0 and (gb < 0)
        raw = dict(np.load(next((x['mesh_path'] for x in json.loads((ROOT / 'inputs/R4_FROZEN_EXPORTS.json').read_text())['exports'] if x['family'] == e['family'])), allow_pickle=False))
        ov = raw['vertices'][np.unique(raw['faces'][raw['roles'] == 0])]
        identity_outer = float(abs(ov - a['locked_outer_vertices']).max())
        assert identity_outer == 0
        rows.append({'family': e['family'], 'same_volume_identity_error_mm3': identity, 'same_cube_count': int(k.sum()), 'state_A_guard_min_mm': ga, 'state_B_guard_min_mm': gb, 'downstream_guard_difference_mm': gb - ga, 'state_A_rule_pass': True, 'state_B_rule_pass': False, 'swapped_kept_cube': keep, 'swapped_forbidden_cube': bad, 'outer_vertex_identity_error_mm': identity_outer, 'minimal_extension': 'spatial cube membership plus rule/pose incidence; volume alone is insufficient', 'witness_origin': 'our_own_fixture based on R2 core, not empirical reference', 'resolution': 'PER_POINT'})
    dump(ROOT / 'raw/R2_SUFFICIENCY.json', rows)
    return rows
if __name__ == '__main__':
    run()
