"""Retain prefixed UNKNOWN statuses as well as bare words; no new value credit."""
import hashlib, json, re, shutil
from pathlib import Path
from planning import *
manifest = read('SOURCE_MANIFEST.json')
occ = []
drift = []
old = read('raw/UNKNOWN_OCCURRENCES.json')
for x in manifest:
    p = Path(x['path'])
    if p.suffix != '.json':
        continue
    if sha(p) != x['sha256']:
        drift.append(dict(path=str(p), frozen_sha256=x['sha256'], current_sha256=sha(p), disposition='REJECT_CURRENT_VERSION_RETAIN_ORIGINAL_WORD_CENSUS'))
        occ.extend((y for y in old if y['path'] == str(p)))
        continue

    def walk(o, ptr=''):
        if isinstance(o, dict):
            for (k, v) in o.items():
                walk(v, ptr + '/' + str(k).replace('~', '~0').replace('/', '~1'))
        elif isinstance(o, list):
            for (i, v) in enumerate(o):
                walk(v, ptr + '/' + str(i))
        elif isinstance(o, str) and re.search('UNKNOWN|NOT_MEASURED|NOT_PERFORMED|NOT_RUN', o, re.I):
            occ.append(dict(job=x['job'], path=x['path'], pointer=ptr, text=o, source_sha256=x['sha256']))
    walk(json.loads(p.read_text()))
initial = HERE / 'raw/UNKNOWN_OCCURRENCES_INITIAL.json'
if not initial.exists():
    shutil.copy2(HERE / 'raw/UNKNOWN_OCCURRENCES.json', initial)
write('raw/UNKNOWN_OCCURRENCES.json', occ)
write('raw/INVENTORY_V2_CHANGE.json', dict(reason='Capture prefixed UNKNOWN_SITE/UNKNOWN_NO_CONTACT/etc.; initial bare-word census preserved', new_occurrences=len(occ), initial_occurrences=len(read('raw/UNKNOWN_OCCURRENCES_INITIAL.json')), value_or_decision_criteria_changed=False, ranked_questions_unchanged=46, drift_rejections=drift, scope='Status strings only; unknown-valued field names with known true/false controls do not become new physical unknowns'))
print('UNKNOWN status scan', len(occ))
