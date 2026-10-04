import re, shutil
import numpy as np
from scipy.stats import t
import xml.etree.ElementTree as ET
from common import R, CORPUS, dump, read, sha, freeze, now
from compare_r2 import ratio

def text(e):
    return ' '.join(' '.join(e.itertext()).split())

def extract():
    source = CORPUS / 'PMC10788321.xml'
    local = R / 'raw/sources/PMC10788321.xml'
    if local.exists() and sha(local) != sha(source):
        raise ValueError('Source snapshot drift')
    shutil.copyfile(source, local)
    root = ET.parse(local).getroot()
    doi = root.find('.//article-id[@pub-id-type="doi"]').text
    rows = root.find('.//table-wrap[@id="Tab2"]').findall('.//tr')
    geom = root.find('.//table-wrap[@id="Tab1"]').findall('.//tr')[1]
    areas = [float(text(x)) for x in list(geom)[1:]]
    out = []
    for (i, angle) in [(2, 0), (3, 30)]:
        cells = list(rows[i])
        raw = text(cells[6])
        match = re.search('(\\d+)\\s*[A-Za-z]*\\s*\\((\\d+)\\)', raw)
        if not match:
            raise ValueError('Unparsed external failure cell ' + raw)
        out.append(dict(id=f'PMC10788321_3Y_angle{angle}', study='PMC10788321', doi=doi, mean_N=float(match[1]), sd_N=float(match[2]), n=8, material='3Y', product='3M Lava Plus Multi L', bridge_class='3unit_posterior_cantilever', cantilever_area_mm2=areas[1], abutment_area_mm2=areas[0], angle_deg=angle, setup='resilient CoCr premolar dies; RelyX Unicem2;6mm steelball;0.5mm/min', aging='10000thermalcycles6.5-60C+1.2millionchewingcycles108N; matched aging/test direction', endpoint='catastrophic_system_fracture_after_aging', origin='not_a_per_specimen_connector_dataset', table_locator=f'/article/body//table-wrap[@id="Tab2"]//tr[{i + 1}]/td[7]', area_locator='/article/body//table-wrap[@id="Tab1"]//tr[2]/td[2:3]', raw_cell=raw, source_sha256=sha(source), group_resolution='POPULATION', specimen_resolution='PER_TOOTH', geometry_resolution='PER_SURFACE_REGION', load_timescale='SIMULTANEOUS', aging_timescale='HANDOVER'))
    dump('raw/R3_MEASUREMENTS.json', out)
    return out

def run():
    d = extract()
    p = read('PREREG_R3.json')
    correct = 'https://doi.org/' + d[0]['doi']
    correction = read('PREREG_R3_SOURCE_CORRECTION.json')
    if correction['correct_locator'] != correct:
        raise ValueError('R3 locator correction differs from primary article')
    groups = []
    for a in d:
        err = t.ppf(0.975, 7) * a['sd_N'] / np.sqrt(8)
        lo = a['mean_N'] - err
        hi = a['mean_N'] + err
        groups.append(dict(id=a['id'], angle_deg=a['angle_deg'], mean_N=a['mean_N'], mean_CI95_N=[lo, hi], cantilever_area_mm2=a['cantilever_area_mm2'], abutment_area_mm2=a['abutment_area_mm2'], area_min_cantilever_mm2=12.0, area_min_abutment_mm2=9.0, area_compliant=a['cantilever_area_mm2'] >= 12 and a['abutment_area_mm2'] >= 9, full_indication_compliance='UNKNOWN: pontic replaces molar, shaped as premolar; no wall-thickness/full-IFU audit', bench_mean_shortfall=hi < p['metrics']['benchmark_required_mean_N'], mean_resolution='POPULATION', mean_margin_vs1000N=a['mean_N'] / 1000.0, mean_margin_CI95_vs1000N=[lo / 1000.0, hi / 1000.0]))
    contrast = ratio(d[1], d[0])
    contrast['upper_CI_below_0.8'] = contrast['CI95'][1] < 0.8
    out = dict(claim_type='information_link', groups=groups, directional_contrast=contrast, gates=dict(area_compliant_mean_shortfall=all((a['area_compliant'] and a['bench_mean_shortfall'] for a in groups)), direction_information=contrast['upper_CI_below_0.8']), practice_control=dict(force_ratio_for_same_area=1.0, status='Area-only rule supplies no orientation-conditioned force;1 is a force-invariance null, not a claimed manufacturer force prediction'), product_rule=dict(locator='https://multimedia.3m.com/mws/media/2122575O/3m-lava-plus-high-translucency-instructions-for-use.pdf', compared_quantity='posterior die-cantilever12mm2, die-die9mm2; indication restricted to premolar/incisor and no bruxism'), equally_informed_control='Conventional summary intervals/contrast agree exactly; no method novelty', external_referent=dict(kind='independent_measurement', locator=correct + '#Tab2', compared_quantity='aged3Y cantilever system break loads, axial vs30deg same geometry', refutes_us=True), interpretation='Area compliance cannot guarantee the frozen1500N mean-capacity bench requirement. Loading-regime information distinguishes same-area capacities; no general manufacturer safety refutation.', limitations=['only one product/batch/study', 'pointwise95%intervals; no familywise multiplicity correction', 'different aging direction as well as fracture direction: cannot isolate instantaneous angle effect', 'not end-supported bridge: cannot pool with R1', 'not measured connector-specific capacity', 'summary CI assumes independent approximately normal group means'])
    dump('raw/R3_COMPARISON.json', out)
    print('R3 matched loading-regime ratio', contrast['ratio'], 'CI', contrast['CI95'], 'gates', out['gates'])
    return out
if __name__ == '__main__':
    run()
