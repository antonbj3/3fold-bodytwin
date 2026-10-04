from dental_release.paths import expand as _release_expand
from common import *
from stats import summarize
check()
p = ROOT / 'FROZEN_PREDICTIONS.json'
assert sha(p) == p.with_suffix('.json.sha256').read_text().strip()
assert load(ROOT / 'rounds/R1.json')['gates']['coverage']
import math
r1 = load(ROOT / 'rounds/R1.json')
assert r1['gates']['control'], 'R1 numerical control gate failed'
r1_controls = load(ROOT / 'raw/R1_CONTROLS.json')
assert len(r1_controls) == 8, 'R1 expected eight fixed controls'
for control in r1_controls:
    (low, high) = control['candidate']
    distance = control['control_mm']
    assert all((math.isfinite(v) for v in [low, high, distance])), 'R1 non-finite control'
    assert low <= high, 'R1 inverted distance bracket'
    error = max(0.0, low - distance, distance - high)
    assert error <= 1e-06, 'R1 independent distance outside bracket'
    assert abs(error - control['error_mm']) <= 1e-12, 'R1 recorded error mismatch'
    assert control['injected_plus1mm_rejected'], 'R1 negative control failed'
assert all(load(ROOT / 'rounds/R2.json')['gates'].values())
r3 = load(ROOT / 'rounds/R3.json')
assert r3['gates']['coverage'] and r3['gates']['controls'] and r3['gates']['containment'] and r3['gates']['decisive_shape_ambiguity']
assert load(ROOT / 'raw/SUFFICIENCY_R1.json')['identity_error_mm'] == 0
assert load(ROOT / 'raw/SUFFICIENCY_R3.json')['identity_error_mm'] == 0
assert load(ROOT / 'raw/SUFFICIENCY_R4_GUIDE.json')['identity_error_mm'] == 0
assert not r3['gates']['summary_sufficiency']
source = load(ROOT / 'PREREG_R1.json')['inputs']
faults = []
bad = dict(source)
k = next(iter(bad))
bad[k] = '0' * 64

def source_valid(expected):
    return all((sha(p) == h for (p, h) in expected.items()))
faults.append(dict(control='full_source_hash_contract', clean_pass=source_valid(source), injected_wrong_hash_rejected=not source_valid(bad)))
canonical = load(X73 / 'raw/IMAGE_GROUPS.json')['canonical']
sites = load(X73 / 'raw/LINEAGE_SITES.json')

def alias_valid(rr):
    return all((r['image_group'] == canonical[r['patient']] for r in rr))
assert alias_valid(sites)
changed = [dict(r) for r in sites]
r = next((r for r in changed if r['patient'] == _release_expand('@DENTAL_CASE_ID@')))
r['image_group'] = _release_expand('@DENTAL_CASE_ID@')
faults.append(dict(control='image_alias_contract', clean_pass=True, injected_wrong_alias_rejected=not alias_valid(changed)))
costs = [load(ROOT / 'raw' / f'{s}_COST.json') for s in ['R1', 'R2']]
for c in costs:
    assert sha(c['artifact']['path']) == c['artifact']['sha256']
assert sha(r3['artifact']['path']) == r3['artifact']['sha256']
assert all((r['rejected'] for r in load(ROOT / 'raw/R4_METROLOGY.json')['faults']))
dump(ROOT / 'raw/VERIFICATION.json', dict(status='PASS', scope='frozen sources, protocol datum, image aliases, all geometry coverage, matched controls, exact-summary witnesses, artifact hashes, prospective port; no physical validation', faults=faults, formal_floating_point_enclosure='MISSING', review_state='PENDING_INDEPENDENT_REVIEW'))
print('Verification PASS; physical guide, tool profile, anatomy and clinical transport UNKNOWN')
