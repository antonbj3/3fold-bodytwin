"""Scope-aware rule parser. No LLM, no default diagnosis for unreadable clauses."""
import re, sys
from pathlib import Path
P = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(P / 'sources/upstream/src'))
from bite2text.report.parse import parse_report as published_parse
FIELDS = ['overbite', 'overjet', 'crossbite', 'midlines', 'spee', 'molar_right', 'molar_left', 'canine_right', 'canine_left']
CL = re.compile('\\bclass\\s*(iii|ii|i|3|2|1)\\b', re.I)
SIDE = re.compile('\\b(?:on the )?(right|left)(?: side)?\\b', re.I)
EDGE = re.compile('(?:edge|end|head)[ -]to[ -](?:edge|end|head)', re.I)

def tooth_class(t, tooth):
    t = re.split(',?\\s*tending (?:toward|towards|to)', t, flags=re.I)[0]
    kws = list(re.finditer('\\b' + tooth + '\\b', t, re.I))
    classes = list(CL.finditer(t))
    if not kws:
        return None
    out = []
    for k in kws:
        tail = t[k.end():]
        tail = re.split('\\b(?:molar|canine|class)\\b', tail, flags=re.I)[0]
        if re.search('not assessable|cannot be (?:assessed|evaluated)', tail, re.I):
            out.append('not assessable')
            continue
        if not classes:
            continue
        before = [c for c in classes if c.end() <= k.start()]
        after = [c for c in classes if c.start() >= k.end()]
        candidates = []
        if before:
            c = before[-1]
            between = t[c.end():k.start()]
            d = k.start() - c.end()
            candidates.append((d, c))
        if after:
            candidates.append((after[0].start() - k.end(), after[0]))
        c = min(candidates, key=lambda v: v[0])[1]
        val = {'1': 'I', '2': 'II', '3': 'III'}.get(c.group(1), c.group(1).upper())
        if val == 'II':
            ci = classes.index(c)
            lo = classes[ci - 1].end() if ci else 0
            hi = classes[ci + 1].start() if ci + 1 < len(classes) else len(t)
            prefix = t[lo:c.start()]
            prefix = re.split('[,;]|\\band\\b', prefix, flags=re.I)[-1]
            suffix = t[c.end():hi]
            suffix = re.split('[,;]|\\band\\b', suffix, flags=re.I)[0]
            val = 'II edge-to-edge' if EDGE.search(prefix + ' ' + suffix) else 'II full'
        out.append(val)
    return out[0] if out and len(set(out)) == 1 else None

def sagittal_sides(t):
    marks = list(SIDE.finditer(t))
    if not marks:
        return {'right': t, 'left': t}
    chunks = [t[:marks[0].start()]] + [t[marks[i].end():marks[i + 1].start()] for i in range(len(marks) - 1)] + [t[marks[-1].end():]]

    def info(v):
        return bool(CL.search(v) or re.search('not assessable|cannot be assessed', v, re.I))
    prefix = info(chunks[0])
    out = {}
    for (i, m) in enumerate(marks):
        v = chunks[i] if prefix else chunks[i + 1]
        out[m.group(1).lower()] = v
    return out

def parse(text):
    base = published_parse(text).as_dict()
    t = text.lower()
    sents = [s.strip() for s in re.split('(?<=[.!?])\\s+', t) if s.strip()]
    sag = ' '.join((s for s in sents if CL.search(s)))
    if sag:
        sides = sagittal_sides(sag)
        for side in ('right', 'left'):
            for tooth in ('molar', 'canine'):
                base[f'{tooth}_{side}'] = tooth_class(sides.get(side, ''), tooth)
    v = re.sub('(?:a |an |and |with |presents )?(?:right |left |bilateral )?lateral open bite[^.;]*', '', t)
    if re.search('(?:anterior )?deep bite|overbite (?:is )?increased|increased overbite', v):
        base['overbite'] = 'increased'
    elif re.search('(?:anterior )?open bite', v):
        base['overbite'] = 'open'
    elif re.search('overbite[^.;]{0,35}(?:reduced|decreased)|reduced overbite', v):
        base['overbite'] = 'reduced'
    elif re.search('correct (?:transverse and )?vertical|vertical relationships? (?:are|is) (?:correct|normal)|overbite[^.;]{0,35}(?:normal|within normal)', v):
        base['overbite'] = 'normal'
    else:
        base['overbite'] = None
    if re.search('midlines?[^.]{0,60}(?:not coincident|not centered|not centred|deviat|off)', t):
        base['midlines'] = 'deviated'
    elif re.search('midlines?[^.]{0,60}(?:coincident|centered|centred)', t):
        base['midlines'] = 'centered'
    else:
        base['midlines'] = None
    crosss = [s for s in sents if re.search('cross-?bite', s)]
    if crosss:
        positives = [s for s in crosss if not re.search('absence of (?:posterior |anterior )?cross|without (?:posterior )?cross|no (?:posterior )?cross', s)]
        base['crossbite'] = 'present' if positives else 'absent'
    elif re.search('correct transverse|transverse(?: skeletal)? relationships? (?:are )?correct|normal transverse', t):
        base['crossbite'] = 'absent'
    else:
        base['crossbite'] = None
    if re.search('(?:increased|accentuated)\\s+curve of spee', t):
        base['spee'] = 'increased'
    return base
