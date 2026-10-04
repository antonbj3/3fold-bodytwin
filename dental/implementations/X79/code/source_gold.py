from common import *
from report_parser import parse
import zipfile
EXPECTED = [[16, 17, 25, 26, 27, 36, 46, 47], [16, 26, 27, 36, 37, 47], [14, 15, 16, 24, 25, 36], [14, 15, 16, 34, 35, 36, 37, 44], [16, 26, 46], [46], [15, 16, 24, 25, 26, 35, 36, 45, 46], [16, 17, 47], [16, 35, 37, 46], [27], [36], [27, 46, 47], [16, 26, 46], [17, 26, 27], [15, 16, 17, 25, 26, 27, 37, 46], [14, 15, 16, 26, 34, 35, 37, 44, 45, 46, 47], [14, 37], [14, 16, 26, 35, 36, 46, 47], [14, 15, 16, 47], [16, 26, 36, 46], 'GLOBAL', 'GLOBAL', 'GLOBAL', 'GLOBAL', 'GLOBAL', 'GLOBAL', 'GLOBAL', 'GLOBAL', 'GLOBAL', 'GLOBAL', [], [], [], [], [], [], [], [], [], []]

def main():
    panel = read('raw/PRODUCER_GOLD_PANEL_SOURCE.json')
    rows = []
    checks = []
    source_n = 0
    emitted_n = 0
    correct = 0
    with zipfile.ZipFile(ZIP) as z:
        for (r, expected) in zip(panel, EXPECTED):
            b = z.read(r['member'])
            p = parse(b.decode())
            a = p['assertions']
            emitted = {(f, x['polarity']) for x in a for f in x['teeth']}
            globalneg = any((x['scope'] == 'ALL_REPORTED_TEETH' for x in a))
            wanted = {(f, 1) for f in expected} if isinstance(expected, list) else set()
            nn = len(emitted) + int(globalneg)
            cc = len(emitted & wanted) + int(globalneg and expected == 'GLOBAL')
            emitted_n += nn
            correct += cc
            source_n += len(wanted) + int(expected == 'GLOBAL')
            rows.append(dict(member=r['member'], sha256=digest(b), gold_positive_FDI=expected if isinstance(expected, list) else [], gold_global_negative=expected == 'GLOBAL', sampling_stratum=r['sampling_stratum']))
            checks.append(dict(member=r['member'], emitted=sorted(emitted), emitted_global_negative=globalneg, correct=cc, n_emitted=nn, false_assertions=sorted(emitted - wanted), gold=expected))
    freeze('PRODUCER_SOURCE_GOLD.json', dict(frozen_utc=now(), parser_sha256=sha(ROOT / 'code/report_parser.py'), independent_review=False, method='Producer read 40 fixed source panels sampled by hash within parser strata. Development/exposure disclosed; not blinded or random clinical audit.', rows=rows))
    precision = correct / emitted_n if emitted_n else 0
    write('raw/PARSER_GOLD_VALIDATION.json', dict(n_reports=len(rows), emitted=emitted_n, correct=correct, precision=precision, source_assertions=source_n, recall=correct / source_n if source_n else None, pass_gate=len(rows) >= 40 and precision >= 0.95, independent_review=False, checks=checks))
    print('Producer parser QA', precision, 'recall', correct / source_n, flush=True)
if __name__ == '__main__':
    main()
