from dental_release.paths import expand as _release_expand
import json, hashlib, re, datetime
from pathlib import Path
ROOT = Path(__file__).resolve().parent
DENTAL = ROOT.parent.parent

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def write(name, data):
    p = ROOT / name
    p.write_text(json.dumps(data, ensure_ascii=False, indent=2) + '\n')
    return sha(p)
now = datetime.datetime.now(datetime.timezone.utc).isoformat()
jobs = ['CROWN-MARGIN-WALL', 'TOOTH-FORCE-CLOSED-BITE', 'GUIDE-ACCURACY', 'PULPOTOMY-CONVERSION', 'MILLING-DEFLECTION', 'CONTRALATERAL-SYMMETRY']
inventory = []
for j in jobs:
    p = DENTAL / 'results' / ('BT-DW48-HUNT-SRC-' + j) / 'RESULTS.md'
    inventory.append(dict(job=j, path=str(p), exists=p.exists(), sha256=sha(p) if p.exists() else None))
write('INPUT_INVENTORY.json', inventory)
prereg = dict(round='XSRC-verify', timestamp_utc=now, claim_type='information_link', capability='Use independently checked numerical literature observations at the correct downstream measurement level.', obstacle='Harvested numbers conflate endpoint, location, denominator and protocol; two requested harvests have no delivered report.', changed_operation='Verify article identity then exact number/unit/locator and gate applicability separately; recompute only declared transferable arithmetic.', consumers=['PROOF_LANE_CROWN_FIX_PREP_R2', 'PROOF_LANE_FULL_CROWN_DIAG', _release_expand('X94'), _release_expand('X87'), _release_expand('X96'), _release_expand('X88'), _release_expand('X89')], gate={'numeric_tolerance': 'exact displayed decimal after unit conversion; no inferred digits', 'source_gate': 'Primary article body/table or indexed abstract supports numerical claim; secondary numbers labelled secondary, not independent primary evidence.', 'transfer_gate': 'same quantity, units, measurement surface/region, protocol, time and denominator; otherwise KEEP_CONTEXT or REJECT_TRANSFER', 'prediction_gate': 'Existing recipient frozen thresholds preserved; no refit or metric change.'}, current_practice='The swarm report unverified or recipient default/missing empirical input', strongest_equally_informed_control='Independent decimal extraction from saved XML and recipient arithmetic; no algorithm advantage claim.', falsifiers=['DOI/PMID points to wrong article', 'number/unit/table absent or wrong', 'different endpoint or measurement resolution falsely treated as matched', 'injected numeric/unit/protocol corruption accepted'], full_cost={'preparation': 'targeted reading and inventory', 'fit': 0, 'discovery': 'source lookup and extraction, human reasoning time unmeasured', 'validation': 'metadata/fulltext comparison and negative controls', 'queries': 'record each API response/failure', 'fallback': 'unavailable source stays UNKNOWN, no fabricated facit'}, limits={'threads': 4, 'intermediate_bytes_max': 3000000000, 'gpu': False})
ps = write('PREREG_R1.json', prereg)
(ROOT / 'PREREG_R1.sha256').write_text(ps + '\n')
write('DECOMPOSITION.json', {'idea': 'A source observation can narrow a consumer only through a compatible observation operator.', 'mechanism': 'Observation y=O(state,region,protocol,time); a population statistic is not the state of a new individual.', 'equations': ['Accept_numeric = identity AND value AND unit AND locator', 'Accept_transfer = Accept_numeric AND same_observable AND compatible_population AND compatible_protocol', 'RCT_later / treated != intraoperative_conversion / planned', 'mean distance does not determine surface p95'], 'leaves': [{'name': 'published decimal and study protocol', 'status': 'EXTERNALLY_MEASURED', 'stop': 'Original observation data unavailable: preserve reported statistic and its sampling limits.'}, {'name': 'unit conversion and arithmetic', 'status': 'DERIVED_UNDER_ASSUMPTIONS', 'stop': 'Exact decimal/rational identities suffice.'}, {'name': 'population-to-individual or mean-to-tail transfer', 'status': 'UNKNOWN', 'stop': 'No paired observations or distribution supplied.'}, {'name': 'tool cantilever or guide radial-tail law', 'status': 'CONSTITUTIVE_CLOSURE', 'stop': 'Use only original recipient model, with unmeasured calibration debt preserved.'}]})
write('FROZEN_PREDICTIONS.json', {'timestamp_utc': now, 'prereg_sha256': ps, 'predictions': ['A per-tooth transducer MVBF cannot falsify a simultaneous small height-change response.', 'Mean labial symmetry cannot calibrate whole-crown p95.', 'Late pulpotomy failures cannot replace intraoperative conversion probability.', 'Missing source-job output is delivery absence, not physical refutation.'], 'new_physical_measurement': False})
(ROOT / 'FROZEN_PREDICTIONS.sha256').write_text(sha(ROOT / 'FROZEN_PREDICTIONS.json') + '\n')
write('OWN_FIRST_CALCULATION.json', {'site_accounting': {'quoted_sites': 236 + 28750 + 2790, '993_times_32': 993 * 32, 'identity_error': 236 + 28750 + 2790 - 993 * 32}, 'conversion_partition': {'planned': [25, 86], 'normal_image': [5, 39], 'complement': [20, 47]}, 'not_external_facit': True})
write('CURRENT_WORK_STATE.json', dict(timestamp_utc=now, phase='SOURCE_VERIFICATION', last_gate='PREREG_R1 frozen and first local arithmetic saved', next_operation='Fetch primary metadata and fulltext; separately audit numeric validity and transfer.', review_state='PENDING_INDEPENDENT_REVIEW'))
print('Frozen prereg', ps, 'reports present', sum((x['exists'] for x in inventory)))
