import numpy as np
from scipy.stats import t
from common import read, dump

def ratio(a, b):
    va = (a['sd_N'] / a['mean_N']) ** 2 / a['n']
    vb = (b['sd_N'] / b['mean_N']) ** 2 / b['n']
    df = (va + vb) ** 2 / (va ** 2 / (a['n'] - 1) + vb ** 2 / (b['n'] - 1))
    se = np.sqrt(va + vb)
    c = t.ppf(0.975, df) * se
    r = a['mean_N'] / b['mean_N']
    return dict(numerator=a['id'], denominator=b['id'], ratio=r, CI95=[float(r * np.exp(-c)), float(r * np.exp(c))], log_SE=float(se), welch_df=float(df), resolution='POPULATION', method='delta log-mean ratio with Welch t; independent groups')

def run():
    d = read('raw/MEASUREMENTS.json')
    p = read('PREREG_R2.json')
    groups = []
    rules = {'3Y': 9.0, '4Y': 16.0, '5Y': 16.0}
    for a in d:
        se = a['sd_N'] / np.sqrt(a['n'])
        c = t.ppf(0.975, a['n'] - 1) * se
        lo = a['mean_N'] - c
        hi = a['mean_N'] + c
        indication = 'within_area_indication_geometry' if a['material'] == '3Y' else 'area_below_rule'
        if a['material'] == '5Y':
            indication += '; posterior_molar_span_outside_UTML_premolar_only_indication'
        groups.append(dict(id=a['id'], material=a['material'], product=a['product'], area_mm2=a['area_mm2'], rule_min_area_mm2=rules[a['material']], below_rule_area_fraction=1 - a['area_mm2'] / rules[a['material']], manufacturer_compliant=a['area_mm2'] >= rules[a['material']], indication=indication, mean_N=a['mean_N'], mean_CI95_N=[lo, hi], mean_margin_vs1000N=a['mean_N'] / 1000.0, mean_margin_CI95_vs1000N=[lo / 1000.0, hi / 1000.0], mean_resolution='POPULATION', bench_excess_evidence=a['area_mm2'] < rules[a['material']] and lo > p['gates']['mean_lower_bound_excess_threshold_N'], bench_deficit_evidence=a['area_mm2'] >= rules[a['material']] and hi < p['gates']['bench_required_mean_N'], individual_survival='UNKNOWN_NO_SPECIMEN_DATA', connector_capacity='UNKNOWN_NONCONNECTOR_EVENT', conditional_bound='For a specimen with confirmed pontic origin under this exact setup, F_connector > F_observed. Per-specimen origin unavailable: no group connector bound claimed.'))
    shapes = [ratio(*[r for r in d if r['material'] == m]) for m in ['3Y', '4Y', '5Y']]
    material = [ratio(d[2 * j], d[2 * j + 2 * k]) for j in [0] for k in [1, 2]]
    out = dict(claim_type='information_link', groups=groups, shape_contrasts=shapes, material_contrasts=material, comparator=dict(locator='https://www.kuraraynoritake.com/world/product/cad_materials/pdf/katana_zircinia_technical_guide.pdf', page='printed6/PDF7; minimum cross-section table', quantity='connector cross-sectional area mm2; ML9/STML16/UTML16 for posterior2-3units; UTML premolar only', status='official specification verified on web; not measured force', numeric_force_prediction=None), control='Direct conventional summary-statistic t intervals are identical; information-link gain is versus fixed area rule, not a numerical algorithm.', external_referent=dict(kind='independent_measurement', locator='https://doi.org/10.4047/jap.2023.15.4.171#T2', compared_quantity='system break load mean N/SD by material at9mm2 plus documented pontic origins', refutes_us=True), conclusion='Bench mean margins can be queried at 9mm2; calibrated connector capacity and clinical excess/deficit remain UNKNOWN.', limitations=['small independent groups n6', 'pointwise95%intervals; no familywise multiplicity correction', 'mean CI is not an individual load floor', 'unaged metal die axial test only', 'mode-censored connector response cannot calibrate absolute connector strength', 'no lithium-disilicate target observations'])
    dump('raw/R2_COMPARISON.json', out)
    for a in groups:
        print(a['id'], 'CI', *[round(v, 1) for v in a['mean_CI95_N']], 'below_rule_excess', a['bench_excess_evidence'])
    return out
if __name__ == '__main__':
    run()
