import re
from lxml import etree
from common import ROOT

def mean_sd(s):
    s = s.replace('−', '-').replace('\xa0', ' ')
    numbers = re.findall('[+-]?\\d+(?:\\.\\d+)?', s)
    assert len(numbers) == 2, s
    return tuple(map(float, numbers))

def text(n):
    return ''.join(n.itertext())

def read_sources():
    bond = []
    tree = etree.parse(str(ROOT / 'raw/bond.xml'))
    trs = tree.xpath('//table-wrap[@id="tbl5"]//tr')
    for (i, tr) in enumerate(trs[1:], 2):
        system = text(tr[0])
        for (j, cover) in [(1, True), (2, False)]:
            (mu, sd) = mean_sd(text(tr[j]))
            bond.append({'system': system, 'with_cover': cover, 'mean_MPa': mu, 'SD_MPa': sd, 'n': 10, 'locator': f'PMC10829558: //table-wrap[@id="tbl5"]//tr[{i}]/*[{j + 1}]', 'resolution_level': 'POPULATION', 'time_scale': 'HANDOVER'})
    tree = etree.parse(str(ROOT / 'raw/optical.xml'))
    optical = []
    substrate = None
    trs = tree.xpath('//table-wrap[@id="tbl2"]//tr')
    target = list(map(float, [text(c) for c in trs[1][1:]]))
    for (i, tr) in enumerate(trs[2:], 3):
        cells = list(tr)
        if len(cells) == 5:
            substrate = text(cells.pop(0))
        cement = text(cells[0])
        vals = [mean_sd(text(c)) for c in cells[1:]]
        optical.append({'substrate': substrate, 'cement': cement, 'mean_Lab': [r[0] for r in vals], 'SD_Lab': [r[1] for r in vals], 'n': 10, 'locator': f'PMC10703855: //table-wrap[@id="tbl2"]//tr[{i}]', 'resolution_level': 'POPULATION', 'time_scale': 'SIMULTANEOUS'})
    return (bond, optical, target)
BOND_FACIT = [(8.35, 1.17), (8.57, 0.71), (5.48, 0.5), (5.57, 0.58), (3.37, 0.85), (4.04, 0.71)]
OPTICAL_FACIT = [[79.59, 0.14, 9.73], [77.03, -0.45, 8.62], [78.06, -1.46, 8.49], [77.09, -0.53, 8.63], [77.06, -1.32, 8.43], [76.54, -0.01, 12.0], [73.35, -0.73, 11.87], [74.52, -0.93, 11.84], [73.57, -0.48, 11.69], [73.48, -0.64, 11.61], [72.92, 0.08, 8.24], [63.63, -1.2, 3.76], [65.46, -1.35, 4.15], [65.27, -1.21, 4.44], [63.63, -1.07, 3.93]]
SD_FACIT = [[0.17, 0.09, 0.31], [0.19, 0.19, 0.29], [0.22, 0.11, 0.37], [0.11, 0.25, 0.29], [0.57, 0.14, 0.34], [0.3, 0.11, 0.26], [0.29, 0.15, 0.24], [0.38, 0.2, 0.17], [0.27, 0.15, 0.16], [0.49, 0.17, 0.12], [0.24, 0.18, 0.25], [0.31, 0.08, 0.32], [0.21, 0.09, 0.23], [0.43, 0.09, 0.25], [0.48, 0.12, 0.29]]

def source_error(bond, optical, target):
    errors = [abs(a - b) for (row, ref) in zip(bond, BOND_FACIT) for (a, b) in zip([row['mean_MPa'], row['SD_MPa']], ref)]
    errors += [abs(a - b) for (row, ref) in zip(optical, OPTICAL_FACIT) for (a, b) in zip(row['mean_Lab'], ref)]
    errors += [abs(a - b) for (row, ref) in zip(optical, SD_FACIT) for (a, b) in zip(row['SD_Lab'], ref)]
    errors += [abs(a - b) for (a, b) in zip(target, [76.7, 1.1, 14.7])]
    assert len(bond) == 6 and len(optical) == 15
    return max(errors)
