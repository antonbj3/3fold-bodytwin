import json, hashlib
from datetime import datetime, timezone
from run_r1 import ROOT, write
r = json.loads((ROOT / 'results.json').read_text())
verify = json.loads((ROOT / 'raw/verification.json').read_text())
r['verification'] = {'status': verify['status'], 'path': 'raw/verification.json', 'checks': len(verify['checks']), 'sha256': hashlib.sha256((ROOT / 'raw/verification.json').read_bytes()).hexdigest()}
write('results.json', r)
feedback = {'lane': 'X17-sleep-apnea', 'review_state': 'PENDING_INDEPENDENT_REVIEW', 'status': 'PENDING_INDEPENDENT_REVIEW', 'target_id': None, 'coverage_status': 'MISSING_OR_INACCESSIBLE_WORKING_COVERAGE', 'source_locator': 'notes/chains_brainstorm/03_ortho_occlusion.md O16', 'result_path': str(ROOT / 'results.json'), 'result_sha256': hashlib.sha256((ROOT / 'results.json').read_bytes()).hexdigest(), 'measured_quantity': ['expert annotated interior pharyngeal minimum area', 'minimum array-plane location', 'conditional uncertainty width', 'published group mean area reconstruction error'], 'units': ['mm2', 'mm', 'mm2', 'fraction'], 'uncertainty': 'geometric scenario interval only; anatomical ROI, orientation and tissue motion UNKNOWN', 'population_regime': '12 TF2 open-bite baseline scans; 2 groups in one published supine awake CBCT cohort. No sleep/AHI prediction.', 'frozen_gate': ['PREREG_R1.json', 'PREREG_R2.json', 'PREREG_R3.json', 'PREREG_R4.json'], 'baseline': ['proportional volume with same scenario volume', 'affine determinant with measured AP/lateral', 'standard full LP', 'standard 2-output linear dose regression'], 'outcome': 'Measured geometry capability; baseline-to-MAD response UNKNOWN; method TIE; failed shared response and scalar-volume gates', 'negative_result': True, 'dispatch_receipt': None, 'feedback_receipt': None, 'graph_failure': 'working rank --query unsupported; working rank --limit5 JSONDecodeError. No graph file edited.', 'coverage_proposal': {'goal': 'known jaw pose to registered pharyngeal area profile and argmin uncertainty', 'dependencies': ['paired movement fields', 'clinical ROI', 'posture / breathing repeat measurement'], 'proposal_only': True}}
write('GRAPH_FEEDBACK.json', feedback)
manifest = []
for p in sorted(ROOT.rglob('*')):
    if p.is_file() and '__pycache__' not in str(p) and (p.name not in ['ARTIFACT_MANIFEST.json', 'CURRENT_WORK_STATE.json', 'run_all_validation.log']):
        manifest.append({'path': str(p.relative_to(ROOT)), 'bytes': p.stat().st_size, 'sha256': hashlib.sha256(p.read_bytes()).hexdigest()})
write('ARTIFACT_MANIFEST.json', manifest)
print(r['overall_outcome'])
print('Artifact and corruption checks:', verify['status'])
print('Reader entry: README_DEMO.md; figure: figure.png; exact outcomes: results.json')
