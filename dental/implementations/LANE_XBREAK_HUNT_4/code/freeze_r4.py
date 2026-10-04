import json, datetime, hashlib
from pathlib import Path
p = json.loads(Path('PREREG_R3.json').read_text())
p.update(round='R4', claim_type='capability', frozen_at_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), parent_artifact='raw/R3_RESULTS.json', capability='Decide whether empirical joint-reference precision survives patient sampling uncertainty.', obstacle='A unique reconstructed count is not a known probability for the next patient.', changed_operation='Union exact binomial intervals over all25 compatible tables, using Bonferroni alpha=0.05/4. Preserve selected scheduled-RCT cohort and mark transport/cluster independence UNKNOWN.', input='Every R3 admissible integer table; no new data or pooled test assumptions.')
p['metrics'] = {'coverage': 'Conditional simultaneous>=95% under independent exchangeable patients within a fixed observation stratum; union over source rounding/missing allocations. Dentist clustering is not accounted and hence empirical coverage UNKNOWN.', 'tolerance': 1e-08, 'decision_gate': 'All4 widths<=0.10 AND both C1 patterns exclude threshold0.8. Failure if either criterion fails.', 'thresholds': 'Same0.1,...,0.9 mathematical loss scenarios; no clinical cutoff.'}
p['strongest_equally_informed_control'] = 'Independent binomial-CDF inversion by brentq, compared with beta-quantile Clopper-Pearson endpoints. Same counts and Bonferroni budget.'
p['falsifier'] = 'Endpoint mismatch>1e-8; width>0.10;0.8 falls in either C1 interval; extrapolated clinic probability presented without selection model.'
p['full_cost']['validation'] = '8 endpoint comparisons over every compatible count pair; no bootstrap claim, no external clinical validation.'
Path('PREREG_R4.json').write_text(json.dumps(p, indent=2) + '\n')
Path('PREREG_R4.sha256').write_text(hashlib.sha256(Path('PREREG_R4.json').read_bytes()).hexdigest() + '\n')
d = json.loads(Path('DECOMPOSITION_R3.json').read_text())
d['round'] = 'R4'
d['changed_operation'] = p['changed_operation']
d['leaves'] += [{'leaf': 'iid_stratum_sampling', 'status': 'CONSTITUTIVE_CLOSURE', 'basis': 'If patients are IID within fixed stratum then k|n~Binomial(n,p).', 'stop_argument': '62 dentist operators share practice/protocol influences; no cluster IDs provided so actual coverage remains UNKNOWN.'}, {'leaf': 'exact_sampling_bounds', 'status': 'DERIVED_UNDER_ASSUMPTIONS', 'basis': 'Invert equal-tailed binomial tests; Bonferroni4 intervals; union over compatible count tables.', 'stop_argument': 'Finite sample coverage follows model, not target population transport.'}]
Path('DECOMPOSITION_R4.json').write_text(json.dumps(d, indent=2) + '\n')
