"""Read exact local primary XML tables. Do not fit synthetic scans to aggregates."""
from dental_release.paths import expand as _release_expand
import json
import os
import re
import xml.etree.ElementTree as ET
from pathlib import Path
import numpy as np
from scipy.optimize import minimize_scalar
from metrology import sha, write_json
R = Path(__file__).resolve().parents[1]
CORPUS = Path(_release_expand('@DENTAL_CORPUS_ROOT@/europepmc/fulltext'))

def rows(table):
    return [[' '.join(''.join(cell.itertext()).split()) for cell in row if cell.tag in ('td', 'th')] for row in table.findall('.//tr')]

def extract():
    if os.environ.get('X55_PORTABLE') == '1' or not (CORPUS / 'PMC10314363.xml').exists():
        lock = json.loads((R / 'DEMO_INPUT_LOCK.json').read_text())['files']
        if sha(R / 'raw/LITERATURE_CELLS.json') != lock['raw/LITERATURE_CELLS.json']:
            raise ValueError('Bundled primary-table transcription hash mismatch')
        saved = json.loads((R / 'raw/LITERATURE_CELLS.json').read_text())
        return (saved['cells'], saved['rejected'])
    cells = []
    rejected = []
    source = []
    for (pmc, table_id, kind) in [('PMC10314363', 'TAB1', 'crown'), ('PMC10333096', 'T15', 'crown'), ('PMC10756847', 'TB2272275-4', 'inlay'), ('PMC10756847', 'TB2272275-5', 'inlay_precision')]:
        p = CORPUS / (pmc + '.xml')
        root = ET.parse(p).getroot()
        doi = next((a.text for a in root.findall('.//article-id') if a.get('pub-id-type') == 'doi'))
        table = next((t for t in root.findall('.//table-wrap') if t.get('id') == table_id))
        raw = rows(table)
        write_json(R / 'raw' / f'{pmc}_{table_id}.json', dict(source=str(p), sha256=sha(p), table_id=table_id, rows=raw))
        source.append(dict(pmc=pmc, path=str(p), sha256=sha(p), table_id=table_id, doi=doi, licence='Individual article licence applies; no blanket corpus redistribution permission'))
        if pmc == 'PMC10314363':
            for row in raw[1:]:
                if len(row) < 4:
                    continue
                if 'casted' in row[0].lower():
                    rejected.append(dict(locator=f'{pmc}#{table_id}', row=row, reason='casting outside milled/printed target'))
                    continue
                cells.append(dict(study=pmc, doi=doi, table_id=table_id, row=row[0], column='Mean', mean_um=float(row[2]), sd_um=float(row[3]), n=int(row[1]), process='printed' if 'printed' in row[0] else 'milled', material='CoCr', object='crown_coping', region='intaglio', estimand='manufacturing_per_crown_RMS', resolution='POPULATION', precision_is='between_crown_RMS_SD; not repeated-scan precision'))
        elif pmc == 'PMC10333096':
            for row in raw:
                if len(row) == 5 and row[0] == 'Design Software':
                    row = row[1:]
                if not row or row[0] not in ['EZIS', '3Shape', 'Exocad']:
                    continue
                for (col, group) in enumerate(['Aegis HM', 'Trione Z', 'Motion 2'], 1):
                    numbers = re.findall('\\d+(?:\\.\\d+)?', row[col])
                    cells.append(dict(study=pmc, doi=doi, table_id=table_id, row=row[0], column=group, mean_um=float(numbers[0]), sd_um=float(numbers[1]), n=40, process='milled', material='zirconia', object='crown', region='whole_compared_surface', estimand='manufacturing_per_crown_RMS', resolution='POPULATION', precision_is='between_crown_RMS_SD; not repeated-scan precision'))
        else:
            for row in raw:
                if not row or row[0] not in ['TS', 'LU', 'ZR', '3D']:
                    continue
                numbers = re.findall('\\d+(?:\\.\\d+)?', row[1])
                estimand = 'manufacturing_between_inlay_pairwise_RMS' if kind == 'inlay_precision' else 'manufacturing_per_inlay_RMS'
                cells.append(dict(study=pmc, doi=doi, table_id=table_id, row=row[0], column='RMS ± SD', mean_um=float(numbers[0]), sd_um=float(numbers[1]), n=78 if kind == 'inlay_precision' else 13, independent_specimens=13, process='printed' if row[0] == '3D' else 'other', material=row[0], object='inlay', region='intaglio', estimand=estimand, resolution='POPULATION', precision_is='78 dependent pairs from 13 separately manufactured inlays; not scanner repeatability'))
    write_json(R / 'raw/LITERATURE_CELLS.json', dict(cells=cells, rejected=rejected, sources=source, crown_repeat_scan_precision='UNKNOWN_NO_ELIGIBLE_TABLE_IN_SELECTED_LOCAL_STUDIES', screened_title_candidates=len((R / 'raw/source_candidates.tsv').read_text().splitlines())))
    return (cells, rejected)

def study_effect(cells):
    studies = []
    for s in sorted({c['study'] for c in cells if c['object'] in ['crown', 'crown_coping']}):
        x = [c for c in cells if c['study'] == s and c['object'] in ['crown', 'crown_coping']]
        y = np.log([c['mean_um'] for c in x])
        v = max(((c['sd_um'] / c['mean_um']) ** 2 / c['n'] for c in x))
        studies.append(dict(study=s, log_RMS_mean=float(y.mean()), variance_proxy=float(v), groups=len(x)))
    y = np.array([s['log_RMS_mean'] for s in studies])
    v = np.array([s['variance_proxy'] for s in studies])

    def reml(tau2):
        w = 1 / (v + tau2)
        mu = np.sum(w * y) / w.sum()
        return np.log(v + tau2).sum() + np.log(w.sum()) + np.sum(w * (y - mu) ** 2)
    fit = minimize_scalar(reml, bounds=(0, 4), method='bounded')
    tau2 = 0 if reml(0) <= fit.fun else float(fit.x)
    w = 1 / (v + tau2)
    mu = float(np.sum(w * y) / w.sum())
    return dict(model='log manufacturing RMS = descriptive intercept + study random intercept', studies=studies, number_independent_studies=len(studies), tau_log_RMS=float(np.sqrt(tau2)), descriptive_geometric_mean_um=float(np.exp(mu)), uncertainty='UNKNOWN: two studies, material/region confounded with study; unknown shared-patient covariance', prospective_crown_prediction='UNKNOWN', scanner_precision_prediction='UNKNOWN', resolution='POPULATION', not_a_surface_algorithm_validation=True)

def run():
    (cells, rejected) = extract()
    model = study_effect(cells)
    write_json(R / 'raw/STUDY_RANDOM_EFFECT.json', model)
    return dict(cells=len(cells), rejected=len(rejected), random_effect=model)
if __name__ == '__main__':
    run()
