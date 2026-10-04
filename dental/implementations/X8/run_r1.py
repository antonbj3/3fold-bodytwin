import hashlib, importlib.util, json, sys, time
from pathlib import Path
import xml.etree.ElementTree as ET
from scipy.stats import binomtest, beta
ROOT = Path(__file__).resolve().parent

def save(name, obj):
    (ROOT / name).write_text(json.dumps(obj, indent=2, ensure_ascii=False) + '\n')
sys.path.insert(0, str(Path(__import__('os').environ.get('DENTAL_PROJECT_ROOT', str(ROOT.parent.parent))) / 'cells/procedure'))
import guide_error
start = time.perf_counter()
pr = ROOT / 'PREREG_R1.json'
assert hashlib.sha256(pr.read_bytes()).hexdigest() == (ROOT / 'PREREG_R1.sha256').read_text().strip()
src = ET.parse(ROOT / 'sources/PMC12221155.xml').getroot()
table = src.find(".//table-wrap[@id='T2']")
raw = ''.join(table.itertext())
assert '316 patients' in raw and '4 hyperesthesia' in raw
assert '14 patients' in raw and '1 hypoesthesia (7.14%)' in raw
rows = [{'study': 'Karabit2018', 'bin_mm': [1, 2], 'events': 4, 'n': 316, 'unit': 'patients', 'source': 'doi:10.4317/medoral.27125 Table2'}, {'study': 'Kutuk2014', 'bin_mm': [1, 2], 'events': 0, 'n': 34, 'unit': 'patients as represented in review', 'source': 'doi:10.4317/medoral.27125 Table2'}, {'study': 'Kutuk2014', 'bin_mm': [0, 1], 'events': 1, 'n': 14, 'unit': 'patients as represented in review; implant subgroup primary lineage', 'source': 'doi:10.4317/medoral.27125 Table2; primary PMID25365650'}]
controls = []
for r in rows:
    (k, n) = (r['events'], r['n'])
    ci = [0 if k == 0 else float(beta.ppf(0.025, k, n - k + 1)), 1 if k == n else float(beta.ppf(0.975, k + 1, n - k))]
    controls.append({**r, 'rate': k / n, 'exact_95pct_ci': ci})
p = float(binomtest(1, 14, 0.68).pvalue)
clin = guide_error.clinical_model()
ang = clin['pred_angle_mean_deg']
gates = {'literal_zero_above_1mm': {'pass': all((r['events'] == 0 for r in rows if r['bin_mm'] == [1, 2])), 'actual_events': 4}, 'transfer_of_68percent': {'pass': p >= 0.01, 'two_sided_binomial_p': p, 'observed': 1 / 14, 'source_exposed': True, 'not_proof_of_primary_study_error': True}, 'K3_angle_external': {'pass': abs(ang - 3.5) <= 1.0, 'prediction_deg': ang, 'published_mean_deg': 3.5, 'abs_error_deg': abs(ang - 3.5), 'external_source': 'doi:10.1111/clr.13346'}}
out = {'construction': 'R1', 'outcome': 'FAIL' if not all((v['pass'] for v in gates.values())) else 'PASS', 'gates': gates, 'clinical_model': clin, 'equally_informed_control': controls, 'conclusion': 'Universal distance-incidence transport and inherited pose closure are rejected for this contract. Clinical per-site risk remains UNKNOWN.', 'external_referent': {'kind': 'independent_measurement', 'locator': 'doi:10.1097/ID.0000000000000160; PMID25365650; subgroup from doi:10.4317/medoral.27125 Table2', 'compared_quantity': 'Neurosensory incidence under 1mm and mean angular deviation', 'refutes_us': True}, 'cost': {'wall_seconds': time.perf_counter() - start, 'manual_discovery_seconds': None}, 'injected_faults': {'false_zero_with_one_event_rejected': not all((x == 0 for x in [0, 1])), 'bad_angle_10deg_rejected': abs(10 - 3.5) > 1, 'impossible_universal_rate_rejected': float(binomtest(1, 14, 0.99).pvalue) < 0.01}}
save('RAW_R1.json', out)
save('CURRENT_WORK_STATE.json', {'lane': 'X8-guide-nerve-risk', 'status': 'R1_COMPLETE_R2_DESIGN', 'latest_gate': out['outcome'], 'next_operation': 'Replace universal clinical curve by geometric ambiguity sets and acquire signed clearance-loss measurement contract', 'R1': 'RAW_R1.json'})
(ROOT / 'HANDOFF_R1.md').write_text('R1 FAIL. Exact gates in RAW_R1.json. Published Table2 has 4/316 events at 1–2mm; universal 68% conflicts with 1/14 subgroup. K3 angle is independently challenged. The subgroup comparison is source exposed. No sitewise clinical risk calibrated. Next: R2 conditional geometry + published magnitude moments, then independent signed clearance-loss calibration.\n')
print(json.dumps(out, indent=2))
