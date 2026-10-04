"""Preserve literal UNKNOWNs; resolve nested metadata and reviewed source aliases."""
from collections import Counter
import re
from planning import *
occ = read('raw/UNKNOWN_OCCURRENCES.json')
ds = read('inputs/DECISIONS.json')
jobids = {j: [d['id'] for d in ds if j in d['source_jobs']] for j in {j for d in ds for j in d['source_jobs']}}
man = read('SOURCE_MANIFEST.json')
x46 = next((x for x in man if x['job'] == 'X46' and x['path'].endswith('results.json')))
rows = json.loads(Path(x46['path']).read_text())['rows']
final = []
extras = []
for (i, original) in enumerate(occ):
    x = dict(original)
    x['occurrence_id'] = i
    job = x['job']
    ptr = x['pointer'].lower()
    segments = ptr.split('/')
    if job in ['REVIEW', 'CONTRACT'] or any(('cost' in s or 'license' in s or 'licence' in s or (s in ['model', 'model_provenance', 'tokens', 'environment', 'permissions', 'producer_claim_type_original', 'route_metadata']) for s in segments)):
        x.update(disposition='EXCLUDED_NONPHYSICAL', reason='Nested cost/model-route/licence metadata or duplicated review wording; no physical assay credit', decision_obligations=[])
    else:
        if job == 'X46' and re.match('^/rows/\\d+/', ptr):
            ix = int(ptr.split('/')[2])
            sourcejob = rows[ix]['job_id']
            x['original_source_job'] = sourcejob
            linked = jobids.get(sourcejob, [])
            if not linked:
                id_ = 'UNRANKED:' + sourcejob
                linked = [id_]
                if not any((y['id'] == id_ for y in extras)):
                    extras.append(dict(id=id_, source_job=sourcejob, decision_blocked=rows[ix].get('finding', 'See source field'), quantity='Unidentified physical closure/consumer from selected legacy review', unit='UNKNOWN', status='UNKNOWN', reason_no_credit='No matched quantity/operator/regime and prospective physical acquisition identified; law/scope review required', locator=x46['path'] + '#/rows/' + str(ix), claim_type='information_link', resolution_level='PHENOMENOLOGICAL'))
        elif job == 'X38':
            linked = ['D01', 'D08', 'D27', 'D42']
        else:
            linked = jobids.get(job, [])
        x.update(disposition='LINKED_SOURCE_OBLIGATION' if linked else 'UNRANKED_SOURCE_OBLIGATION', reason='Source-level obligation set; repeated numerical and historical UNKNOWN outputs are not independent leaves or decisions', decision_obligations=linked)
    final.append(x)
write('raw/UNKNOWN_FINAL_INVENTORY.json', final)
write('raw/UNRANKED_PHYSICAL_OBLIGATIONS.json', extras)
counts = Counter((x['disposition'] for x in final))
write('raw/INVENTORY_FINAL_SCOPE.json', dict(literal_occurrences=len(final), dispositions=dict(counts), excluded_fraction=counts['EXCLUDED_NONPHYSICAL'] / len(final), unranked_legacy_obligations=len(extras), literal_unassigned_occurrences=sum((not x['decision_obligations'] and x['disposition'] != 'EXCLUDED_NONPHYSICAL' for x in final)), ranked_consumer_questions=46, legacy_review_obligations=18, scope='44 review documents,46 index entries,5 additional X40 package consumers,18 reviewed legacy sources. Every literal UNKNOWN retained; source-level set links are not a one-to-one verified physical leaf ontology.', attrition='Metadata/cost/licence and unranked legacy leaves excluded from value; repeated sites and aliases never increase value; no corpus prevalence claimed'))
print(dict(counts), 'unranked legacy', len(extras))
