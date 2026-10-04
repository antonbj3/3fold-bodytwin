"""Conservative treatment extraction. Silence and all ambiguities are UNKNOWN."""
import re
TARGET = re.compile('restor(?:ation\\w*|ative\\s+(?:dental\\s+|conservative\\s+)?(?:treatment\\w*|filling\\w*))|\\bfilling\\w*|(?:prosthetic\\s+|\\ba\\s+)crown\\w*|prosthetic\\s+element\\w*|prosthetically\\s+replaced', re.I)
FDI = re.compile('\\b([1-8][1-8])\\b')
HEDGE = re.compile('\\bmay\\b|\\bmight\\b|\\bseem\\w*\\b|\\bappear\\w*\\b|\\bpossibl\\w*\\b|\\bprobably\\b|cannot.*certainty|not.*(?:visible|quality)|not of sufficient quality', re.I)
NEG = re.compile('\\bno\\b|\\bwithout\\b|\\bnot present\\b|absence of|\\bdo not\\b|\\bdoes not\\b|\\bnor\\b', re.I)

def permanent(f):
    return f // 10 in (1, 2, 3, 4) and f % 10 in range(1, 9)

def parse(text):
    emitted = []
    excluded = []
    for sentence in re.split('(?<=[.!?])\\s+', text.strip()):
        if not TARGET.search(sentence):
            continue
        if re.search('not.*(?:visible|quality)|photo.*missing', sentence, re.I):
            excluded.append(dict(clause=sentence, reason='VISIBILITY_QUALITY_UNKNOWN'))
            continue
        clauses = re.split(';|\\bwhile\\b|\\bwhereas\\b|\\bbut\\b|\\bhowever\\b|\\band\\s+(?=a prosthetic crown)', sentence, flags=re.I)
        for clause in clauses:
            if not TARGET.search(clause):
                continue
            clause = re.split(',?\\s+and\\s+(?:active |ongoing |evident |there |caries|tooth|teeth)', clause, flags=re.I)[0]
            if HEDGE.search(clause):
                excluded.append(dict(clause=clause, reason='HEDGED'))
                continue
            if re.search('sealants?\\s*/\\s*(?:conservative\\s+)?restor|sealants?\\s+or\\s+restor|restor.*?/.*?sealant', clause, re.I) and (not NEG.search(clause)):
                excluded.append(dict(clause=clause, reason='SEALANT_RESTORATION_AMBIGUOUS'))
                continue
            teeth = sorted(set((int(f) for f in FDI.findall(clause))))
            neg = bool(NEG.search(clause))
            if neg and re.search('\\bother\\b', clause, re.I):
                excluded.append(dict(clause=clause, reason='OTHER_NOT_GLOBAL_NEGATIVE'))
                continue
            kind = 'crown' if re.search('prosthetic|prosthetically', clause, re.I) else 'restoration'
            if neg and (not teeth):
                if re.search('first molars|upper|lower|maxillary|mandibular|\\bother\\b', clause, re.I):
                    excluded.append(dict(clause=clause, reason='UNLOCALIZED_RESTRICTED_NEGATIVE'))
                    continue
                emitted.append(dict(teeth=[], polarity=0, scope='ALL_REPORTED_TEETH', kind='restoration', clause=clause))
                continue
            if not teeth:
                excluded.append(dict(clause=clause, reason='UNLOCALIZED_POSITIVE'))
                continue
            good = [f for f in teeth if permanent(f)]
            bad = [f for f in teeth if not permanent(f)]
            if bad:
                excluded.append(dict(clause=clause, reason='PRIMARY_FDI_OUT_OF_SCOPE', teeth=bad))
            if good:
                emitted.append(dict(teeth=good, polarity=int(not neg), scope='NAMED_TEETH', kind=kind, clause=clause))
    pos = {f for r in emitted if r['polarity'] for f in r['teeth']}
    neg = {f for r in emitted if not r['polarity'] for f in r['teeth']}
    globalneg = any((r['scope'] == 'ALL_REPORTED_TEETH' for r in emitted))
    conflict = bool(pos & neg or (globalneg and pos))
    if conflict:
        excluded.append(dict(reason='CONFLICTING_SAME_REPORT', clause=text))
        emitted = []
    return dict(assertions=emitted, excluded=excluded, conflict=conflict)
