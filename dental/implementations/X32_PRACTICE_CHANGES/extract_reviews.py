"""Table cell transcription checked against rendered primary PDF pages.
Facts are not pooled, denominators are retained, no source clinical endpoint is merged.
"""
from pathlib import Path
import csv, json, re, hashlib
R = Path(__file__).resolve().parent
S = [('secondary_caries', 'Caries on abutments', 2908, 1.0, 0.8, 1.4, 1685, 0.5, 0.2, 1.6), ('caries_loss', 'SCs lost due to caries', 4303, 0.5, 0.2, 1.7, 1797, 0.06, 0.007, 0.4), ('abutment_tooth_fracture', 'SCs lost due to abutment tooth fracture', 5276, 1.2, 0.7, 2.0, 1628, 0.2, 0.1, 0.8), ('endodontic', 'Loss of abutment tooth vitality', 1684, 1.7, 1.6, 1.8, 395, 0.7, 0.3, 1.4), ('discoloration', 'Marginal discoloration', 345, 1.8, 0.7, 4.4, 975, 2.3, 0.6, 8.4), ('framework_fracture', 'Framework fracture', 3075, 0.03, 0.002, 0.3, 1952, 2.3, 1.0, 5.5), ('ceramic_fracture_loss', 'SCs lost due to ceramic fractures', 4457, 0.3, 0.1, 0.6, 1155, 1.1, 0.7, 1.8), ('chipping', 'Ceramic chipping', 1146, 2.6, 1.3, 5.2, 1408, 1.5, 0.9, 2.4), ('retention_loss', 'Loss of retention', 2971, 0.6, 0.4, 1.0, 1583, 1.0, 0.4, 2.4), ('aesthetic', 'Esthetic failures', 2806, 0.5, 0.4, 0.8, 1006, 0, None, None)]
P = [('any_complication', 'Total number of SCs with complications', 1300, 13.3, 9, 19.3, 76, 16.2, 6.2, 38.4), ('soft_tissue', 'Soft tissue complications', 2118, 5.1, 2.3, 11, 234, 5.3, 2, 13.7), ('bone_loss', 'Significant marginal bone loss', 3254, 3.3, 1.6, 6.8, 670, 4.3, 3.4, 5.5), ('aesthetic', 'Aesthetic failures', 627, 1.7, 0.5, 5.6, 224, 0, 0, 2.2), ('abutment_fracture', 'Abutment fracture', 3998, 0.2, 0.06, 0.5, 790, 0.4, 0.2, 0.6), ('screw_fracture', 'Abutment or occlusal screw fracture', 3788, 0.05, 0.01, 0.2, 814, 0.1, 0.07, 0.2), ('screw_loosening', 'Abutment or occlusal screw loosening', 3954, 3.6, 1.7, 7.5, 694, 1.0, 0.5, 2.1), ('ceramic_fracture_or_chipping', 'Ceramic fracture or chipping', 4090, 2.9, 1.7, 4.7, 694, 2.8, 0.6, 12.2), ('restoration_fracture_loss', 'Failure due to fracture of the restoration', 2592, 0.2, 0.07, 0.67, 371, 2.1, 0.9, 5.1), ('retention_loss', 'Loss of retention of cemented SCs', 2211, 2.0, 1.3, 3.1, 115, 0, 0, 4.8)]

def extract():
    out = []
    checks = []
    for (tag, table, doi, groups, data) in [('sailer2015', 7, '10.1016/j.dental.2015.02.011', ['tooth_supported_metal_ceramic', 'tooth_supported_leucite_or_lithium_disilicate'], S), ('pjetursson2018', 5, '10.1111/clr.13306', ['implant_supported_metal_ceramic', 'implant_supported_veneered_zirconia'], P)]:
        text = (R / 'sources' / f'{tag}.txt').read_text()
        block = text[text.index('Table 7' if tag == 'sailer2015' else 'TA B L E 5'):]
        for (mode, q, *values) in data:
            for (i, g) in enumerate(groups):
                (n, r, lo, hi) = values[i * 4:(i + 1) * 4]
                lines = block.replace(',', '').splitlines()
                candidates = []
                for (ix, line) in enumerate(lines):
                    m = re.search('(?<!\\d)' + str(n) + '(?!\\d)', line)
                    if m and q.split()[0] in line:
                        candidates.append((ix, line, m))
                if len(candidates) > 1:
                    candidates = [c for c in candidates if q.split()[0] in c[1] and ('loosening' in c[1] if mode == 'screw_loosening' else 'chipping' in c[1] if mode == 'ceramic_fracture_or_chipping' else True)]
                assert len(candidates) == 1, (tag, q, n, 'ambiguous source row', len(candidates))
                (ix, line, m) = candidates[0]
                col = m.start()
                if tag == 'sailer2015':
                    starts = [v.start() for v in re.finditer('\\b\\d{2,4}\\s+\\d+(?:\\.\\d+)?\\s+\\d+(?:\\.\\d+)?%', line)]
                    end = min([v for v in starts if v > col], default=len(line))
                    frag = line[col:end] + '\n' + lines[ix + 1][col:end]
                else:
                    other = [v.start() for v in re.finditer('\\b\\d{2,4}\\s+\\d+(?:\\.\\d+)?\\*?', line)]
                    end = min([v for v in other if v > col], default=len(line))
                    frag = line[col:end]
                numbers = [float(x) for x in re.findall('(?<![A-Za-z])\\d+(?:\\.\\d+)?', frag)]
                assert any((abs(x - r) < 1e-08 for x in numbers)), (tag, q, 'rate', frag)
                if lo is not None:
                    assert any((abs(x - lo) < 1e-08 for x in numbers)), (tag, q, 'lo', frag)
                if hi is not None:
                    assert any((abs(x - hi) < 1e-08 for x in numbers)), (tag, q, 'hi', frag)
                out.append(dict(source=tag, doi=doi, table=table, source_quantity=q, failure_mode=mode, population=g, n_crowns_or_abutments=n, rate_5year_per100=r, ci95_low=lo, ci95_high=hi, resolution_level='POPULATION', verification='PRIMARY_TABLE_VERIFIED', source_path=f'sources/{tag}.pdf', coverage='NO_VALIDATED_CLINICAL_ENDPOINT_MODEL', window='5 years', note='Each row has its own study subset/denominator; no unique-failure share inferred'))
                checks.append(dict(source=tag, quantity=q, population=g, n=n, fragment=frag))
    with (R / 'raw/review_frequencies.csv').open('w') as f:
        w = csv.DictWriter(f, fieldnames=list(out[0]))
        w.writeheader()
        w.writerows(out)
    (R / 'raw/SOURCE_CELL_CHECKS.json').write_text(json.dumps(checks, indent=2) + '\n')
    return out
if __name__ == '__main__':
    print('Verified frequency cells:', len(extract()))
