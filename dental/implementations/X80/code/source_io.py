"""Read original local JATS tables; preserve their observation operators and hashes."""
from dental_release.paths import expand as _release_expand
from pathlib import Path
import xml.etree.ElementTree as ET
import hashlib, re, json
ROOT = Path(__file__).resolve().parents[1]
CORPUS = Path(_release_expand('@DENTAL_CORPUS_ROOT@/europepmc/fulltext'))
IDS = ['PMC10223376', 'PMC10415425', 'PMC10087269', 'PMC10381089', 'PMC10415439', 'PMC10004271', 'PMC10892052', 'PMC10780983', 'PMC10182118', 'PMC10829619', 'PMC11010625', 'PMC11012777', 'PMC10097162', 'PMC10896305', 'PMC10916199', 'PMC10670842', 'PMC10767728']

def text(e):
    return ' '.join(' '.join(e.itertext()).split()) if e is not None else ''

def nums(s):
    return [float(x.replace('−', '-')) for x in re.findall('[−-]?\\d+(?:\\.\\d+)?', s)]

def grid(table):
    """Expand JATS HTML rowspan/colspan without copying source trees."""
    out = []
    pending = {}
    for tr in table.findall('.//tr'):
        row = {k: v[0] for (k, v) in pending.items()}
        new = {k: (v[0], v[1] - 1) for (k, v) in pending.items() if v[1] > 1}
        col = 0
        for td in tr:
            while col in row:
                col += 1
            val = text(td)
            cs = int(td.get('colspan', '1'))
            rs = int(td.get('rowspan', '1'))
            for c in range(col, col + cs):
                row[c] = val
                if rs > 1:
                    new[c] = (val, rs - 1)
            col += cs
        out.append([row.get(i, '') for i in range(max(row) + 1)] if row else [])
        pending = new
    return out

class Sources:

    def __init__(self):
        self.roots = {}
        self.manifest = []
        for id in IDS:
            p = CORPUS / (id + '.xml')
            b = p.read_bytes()
            r = ET.fromstring(b)
            self.roots[id] = r
            doi = next((x.text for x in r.findall('.//article-id') if x.get('pub-id-type') == 'doi'), None)
            self.manifest.append(dict(id=id, path=str(p), sha256=hashlib.sha256(b).hexdigest(), bytes=len(b), doi=doi, title=text(r.find('.//article-title')), article_type=r.get('article-type'), license=text(r.find('.//permissions')), tables=[dict(id=t.get('id'), caption=text(t.find('caption')), rows=grid(t), footnote=text(t.find('table-wrap-foot'))) for t in r.findall('.//table-wrap')]))

    def table(self, id, n):
        return grid(self.roots[id].findall('.//table-wrap')[n - 1])

    def body(self, id):
        return [text(p) for p in self.roots[id].findall('.//body//p')]

    def record(self, id, quantity, locator, values, unit, region, operator, **kwargs):
        m = next((m for m in self.manifest if m['id'] == id))
        return dict(source=id, doi=m['doi'], locator=locator, source_path=m['path'], source_sha256=m['sha256'], quantity=quantity, values=values, unit=unit, resolution='POPULATION', measurement_region=region, observation_operator=operator, **kwargs)

    def extract(self):
        records = []
        t = self.table('PMC10415425', 2)
        assert t[1][0] == 'Mean value'
        masses = {c: {'mean': float(t[1][j]), 'sd': float(t[2][j])} for (j, c) in enumerate(t[0][1:], 1)}
        records.append(self.record('PMC10415425', 'retained Ti mass', 'Table2 Mean value/SD rows', masses, 'µg', 'collected silicone matrix model', 'AAS after digestion; bulk element; not TiO2 dose'))
        t = self.table('PMC10415425', 3)
        recover = [float(r[-1].replace(',', '')) for r in t[1:]]
        records.append(self.record('PMC10415425', 'measurement-recovery', 'Table3 precision balance 19mg control', recover, 'µg', 'bulk recovery control', 'AAS after digestion'))
        t = self.table('PMC10223376', 5)
        curr = {r[0]: {'mean': nums(r[1])[0], 'sd': nums(r[1])[1]} for r in t[1:]}
        records.append(self.record('PMC10223376', 'corrosion current', 'Table5 jCORR', curr, 'µA/cm²', 'whole implant exposed area', 'potentiodynamic current; abstract uses µA/mm²: factor100 conflict', unit_conflict=True))
        organs = {}
        for r in self.table('PMC10087269', 3):
            if r[0] in ['Brain', 'Spleen', 'Lungs', 'Liver']:
                organ = r[0]
            elif r[0].startswith('Titanium'):
                organs[organ] = {'experimental_median': nums(r[1])[0], 'experimental_iqr': nums(r[1])[1], 'control_median': nums(r[2])[0], 'control_iqr': nums(r[2])[1], 'p': nums(r[3])[0]}
        records.append(self.record('PMC10087269', 'organ Ti at30d', 'Table3 Titanium rows', organs, 'ng/g', 'rat organs', 'microwave acid digestion plus ICP; total elemental Ti; source calls ions'))
        rows = [r for r in self.table('PMC10087269', 2) if r[0].startswith('Titanium')]
        blood = {'before': dict(zip(['experimental_median', 'experimental_iqr', 'control_median', 'control_iqr'], nums(rows[0][1]) + nums(rows[0][2]))), 'after': dict(zip(['experimental_median', 'experimental_iqr', 'control_median', 'control_iqr'], nums(rows[1][1]) + nums(rows[1][2])))}
        records.append(self.record('PMC10087269', 'blood Ti at30d', 'Table2 Titanium pre/post rows', blood, 'ng/g', 'rat blood', 'acid digestion plus ICP; total elemental Ti'))
        abstract = text(self.roots['PMC10223376'].find('.//abstract'))
        assert 'not exceeding 6 ppb in 30 days' in abstract
        body = self.body('PMC10223376')
        assert any(('release more ions are the rough' in p for p in body))
        records.append(self.record('PMC10223376', 'measured Ti release', 'Results Figure5 description and Abstract', {'maximum_reported_ppb': 6, 'order': 'R>H and R>L', 'bath_volume_ml': 100, 'temperature_C': 37, 'times_days': [1, 7, 14, 30]}, 'µg/L', '100mL Hanks bath', 'ICP source stated ion release; inequality only, no Figure5 point digitization'))
        t = self.table('PMC10182118', 2)
        species = {r[0]: [float(x) for x in r[1:]] for r in t if r[0] in ['TEGDMA', 'Bis-GMA', 'Bis-EMA', 'Total']}
        assert len(species) == 4
        records.append(self.record('PMC10182118', 'normalized monomer elution', 'Table2 individual species and Total rows', species, 'wt.ppm', '14x14x1.9mm specimen, 10d extraction', 'HPLC; individual species denominator=initial species mass; Total is composition weighted', column_order=['self_ethanol', 'self_ethanol_water', 'print_composite_ethanol', 'print_composite_ethanol_water', 'print_resin_ethanol', 'print_resin_ethanol_water']))
        depth = []
        for r in self.table('PMC10829619', 1):
            if len(r) == 9 and r[0] in ['0', '15', '30', '60', '90', '120'] and (r[1] in ['RT', '40', '60', '80']):
                vals = [nums(c)[:2] for c in r[3:9]]
                assert len(vals) == 6 and all((len(v) == 2 for v in vals))
                depth.append({'time_min': int(r[0]), 'temperature_C': None if r[1] == 'RT' else int(r[1]), 'temperature_label': r[1], 'mean': [v[0] for v in vals], 'sd': [v[1] for v in vals], 'depth_mm': [0, 1, 2, 3, 4, 5]})
        assert len(depth) == 24, len(depth)
        records.append(self.record('PMC10829619', 'internal hardness', 'Table1 all24 recipes, depth0..5mm', depth, 'Vickers hardness', 'cross-section z=0..5mm', 'source mean/SD; discrete depth indentations', local_resolution='PER_POINT', unit_provenance='Methods Vickers microhardness; dimensionless HV'))
        t = self.table('PMC11010625', 1)
        bilateral = {}
        for r in t:
            if r[0] in ['Top', 'Bottom']:
                thickness = '1mm' if r[0] not in bilateral else '2mm'
                if thickness == '2mm':
                    key = r[0] + '_2mm'
                else:
                    key = r[0]
                bilateral[key] = dict(zip(['NC', 'PC', 'T2', 'T4', 'T8', 'T1B1', 'T2B2', 'T4B4'], [nums(c)[:2] for c in r[1:9]]))
        records.append(self.record('PMC11010625', 'directional light response', 'Table1 A/B Top/Bottom', bilateral, 'Vickers hardness', '1/2mm discs top and bottom', 'indentation; equal total LCU cycles, different direction', local_resolution='PER_POINT'))
        paragraphs = self.body('PMC10892052')
        assert any(('(63.1%)' in p for p in paragraphs))
        assert any(('(66.2%)' in p and '(68.2%)' in p for p in paragraphs))
        records.append(self.record('PMC10892052', 'post-cure DC', 'Results: Triad lowest at15 and tied-lowest at45; Discussion: baseline63.1', {'time_min': [0, 15, 45], 'dc_percent': [63.1, 66.2, 68.2], 'oven': 'Triad2000', 'temperature_C': None, 'sd_percent': None, 'irradiance_mW_cm2': None}, '%', 'model resin ATR FTIR', 'published prose means; exact raw SD not locally tabulated'))
        clinical = {}
        for (i, quantity) in enumerate(['SBI', 'PD', 'PES', 'CBL'], 1):
            t = self.table('PMC10896305', i)
            rows = [r for r in t if r[0] in ['Group 1', 'Group 2']]
            assert len(rows) == 6
            clinical[quantity] = {'n_pairs': 15, 'unit': {'SBI': 'index', 'PD': 'mm', 'PES': 'score', 'CBL': 'mm'}[quantity], 'Zr_mean': float(rows[-2][2]), 'Zr_sd': float(rows[-2][3]), 'Ti_mean': float(rows[-1][2]), 'Ti_sd': float(rows[-1][3]), 'published_p': rows[-2][4]}
        records.append(self.record('PMC10896305', 'six-month paired abutment endpoints', 'Tables1–4 last two group rows', clinical, 'mixed units: each endpoint carries unit', '15 patients, 30 implants; customized abutments', 'paired means/SD; correlation not published; titanium fixtures in both arms'))
        self.records = records
        (ROOT / 'raw/source_manifest.json').write_text(json.dumps(self.manifest, indent=2, ensure_ascii=False) + '\n')
        (ROOT / 'raw/measurements.json').write_text(json.dumps(records, indent=2, ensure_ascii=False) + '\n')
        return {r['quantity']: r for r in records}
