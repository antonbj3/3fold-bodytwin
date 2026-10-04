from dental_release.paths import expand as _release_expand
import sys, zipfile
from shared import *
out = Path(sys.argv[1])
zpath = sys.argv[2]
x7 = RESULTS / _release_expand('X7')
sys.path.insert(0, str(x7 / 'code'))
import parser
pred = load(out / 'X7_PREDICTION.json')
reports = []
x21 = RESULTS / _release_expand('X21')
sys.path.insert(0, str(x21 / 'code'))
from repair_reports import repair_text
scoped = load(x21 / 'raw/SCOPED_REPORTS_R6.json')[_release_expand('@DENTAL_SURFACE_ID_C@')]
with zipfile.ZipFile(zpath) as z:
    for n in sorted(z.namelist()):
        if n.startswith(_release_expand('@DENTAL_SURFACE_ID_C@/')) and '/reports_ios_en/' in n and n.endswith('.txt'):
            b = z.read(n)
            labels = parser.parse(b.decode('utf8'))
            original_labels = dict(labels)
            (corrected, clauses, _) = repair_text(b.decode('utf8'), labels.get('crossbite'))
            labels['crossbite'] = corrected
            comparison = {}
            for (k, v) in pred['fields'].items():
                truth = labels.get(k)
                if truth is not None:
                    comparison[k] = dict(reference=truth, prediction=v['category_set'], covered=truth in v['category_set'], forced=v['single_forced_baseline'], forced_matches=v['single_forced_baseline'] == truth)
            reports.append(dict(member=n, sha256=__import__('hashlib').sha256(b).hexdigest(), categories=labels, original_X7_labels=original_labels, X21_negation_clauses=clauses, regional_crossbite=next((r for r in scoped if r['member'] == n), None), comparison=comparison))
for r in reports:
    if 'crossbite' in r['comparison']:
        r['comparison']['crossbite']['scope'] = 'X7 categories compared against X21 negation-corrected reference; regional propositions preserved separately'
dump(out / 'REPORT_COMPARISON.json', dict(reports=reports, external_referent=dict(kind='published_dataset', locator=zpath + _release_expand('::@DENTAL_SURFACE_ID_C@/reports_ios_en/'), compared_quantity='Clinician report categories, not geometric mm/force', refutes_us=any((not a['covered'] for r in reports for a in r['comparison'].values()))), population_accuracy='Not estimated from one selected, previously evaluated case'))
