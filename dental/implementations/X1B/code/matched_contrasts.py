"""Pre-registered paired-thickness predictions; no prospective force is observed here."""
import re, copy, zipfile, xml.etree.ElementTree as ET
import numpy as np
from scipy.special import gamma, digamma
from scipy.stats import t
import trimesh
from common import R, read, dump, sha, now, state

def derived_groups():
    source = R / 'inputs/literature/PMC8558575.txt'
    s = source.read_text().split('Table 1 Characteristic strength')[1].split('Fig. 4')[0]

    def values(start, end):
        return [float(x) for x in re.findall('\\d+(?:\\.\\d+)?', s.rsplit(start, 1)[1].split(end)[0])]
    F0 = values('Characteristic strength (N)', 'lower limit')
    Flo = [float(x) for x in re.findall('\\d+(?:\\.\\d+)?', s.split('Characteristic strength (N)')[1].split('lower limit')[1].split('upper limit')[0])]
    Fhi = [float(x) for x in re.findall('\\d+(?:\\.\\d+)?', s.split('upper limit')[1].split('Characteristic strength (N)')[0])]
    ms = s.split('lower limit', 1)[1]
    mhi = [float(x) for x in re.findall('\\d+(?:\\.\\d+)?', ms.split('upper limit')[1].split('Weibull modulus (m)')[0])]
    mm = values('Weibull modulus (m)', 'lower limit')
    mlo = [float(x) for x in re.findall('\\d+(?:\\.\\d+)?', s.rsplit('Weibull modulus (m)', 1)[1].split('lower limit')[1])]
    if not all((len(x) == 10 for x in [F0, Flo, Fhi, mm, mlo, mhi])):
        raise ValueError('Unexpected Weibull primary table structure')
    groups = []
    for (i, thick) in [(2, 1.5), (4, 1.0), (6, 0.8), (8, 0.5)]:
        mu = F0[i] * gamma(1 + 1 / mm[i])
        vF = (np.log(Fhi[i] / Flo[i]) / (2 * 1.6448536269514722)) ** 2
        vm = (np.log(mhi[i] / mlo[i]) / (2 * 1.6448536269514722)) ** 2
        derivative = -digamma(1 + 1 / mm[i]) / mm[i]
        variance = vF + derivative ** 2 * vm
        groups.append(dict(study='PMC8558575', row_id=f'PMC8558575:derived:G{thick}', thickness_mm=thick, mean_N=float(mu), log_mean_variance=float(variance), n=14, characteristic_load_N=F0[i], Weibull_m=mm[i], F0_90CI_N=[Flo[i], Fhi[i]], m_90CI=[mlo[i], mhi[i]], angle_deg=0, material='3Y', quantity='Weibull-derived arithmetic mean; NOT independently measured mean', locator=f'doi:10.4047/jap.2021.13.5.269 Table 1 non-fatigued G{thick}, characteristic force/m/90CI', source_sha256=sha(source), closure='reported two-parameter Weibull law; parameter covariance not published and taken zero as explicit approximation', fracture_origin='occlusal loading point per primary results; does not validate intaglio mechanism'))
    for r in read('raw/LITERATURE_GROUPS.json'):
        if r['study'] == 'PMC10817558' and r['material'] == '3Y':
            groups.append(r | dict(log_mean_variance=float(np.log1p((r['sd_N'] / r['mean_N']) ** 2) / r['n']), quantity='measured group arithmetic mean'))
    return groups

def slope(rows):
    x = np.log([r['thickness_mm'] for r in rows])
    y = np.log([r['mean_N'] for r in rows])
    w = 1 / np.array([r['log_mean_variance'] for r in rows])
    xx = x - np.sum(w * x) / w.sum()
    yy = y - np.sum(w * y) / w.sum()
    b = np.sum(w * xx * yy) / np.sum(w * xx * xx)
    v = 1 / np.sum(w * xx * xx)
    residual = y - (np.sum(w * y) / w.sum() + b * xx)
    sig = float(np.mean(residual ** 2)) if len(rows) > 2 else 0.0
    return dict(slope=float(b), variance=float(v), residual_log_variance=sig, n_groups=len(rows))

def meta(slopes):
    b = np.array([s['slope'] for s in slopes])
    v = np.array([s['variance'] for s in slopes])
    w = 1 / v
    fixed = np.sum(w * b) / w.sum()
    Q = np.sum(w * (b - fixed) ** 2)
    c = w.sum() - (w * w).sum() / w.sum()
    tau = max(0.0, float((Q - (len(b) - 1)) / c)) if len(b) > 1 else 0.0
    rw = 1 / (v + tau)
    return dict(slope=float(np.sum(rw * b) / rw.sum()), variance=float(1 / rw.sum()), between_study_slope_variance=tau, residual_log_variance=float(np.mean([s['residual_log_variance'] for s in slopes])), n_studies=len(b))

def band(m, ratio_t, nref=12, ntest=12):
    if ratio_t == 1:
        return dict(point_ratio=1.0, interval_ratio=[1.0, 1.0], interval_level=0.9, upper_lower_ratio=1.0, reference_n_min=nref, test_n_min=ntest, validity='IDENTITY', quantity='reference group divided by itself', interval_assumptions='Exact identity; the same reference observation occurs in numerator and denominator')
    x = np.log(ratio_t)
    v = x * x * (m['variance'] + m['between_study_slope_variance']) + 2 * m['residual_log_variance']
    q = t.ppf(0.95, max(1, m['n_studies'] - 1))
    h = q * np.sqrt(v)
    mu = m['slope'] * x
    (lo, hi) = np.exp([mu - h, mu + h])
    return dict(point_ratio=float(np.exp(mu)), interval_ratio=[float(lo), float(hi)], interval_level=0.9, upper_lower_ratio=float(hi / lo), reference_n_min=nref, test_n_min=ntest, validity='PROVISIONAL_TRANSFER', quantity='same-batch group arithmetic mean force ratio', interval_assumptions='Study-slope exchangeability, same batch/support/contact, log-normal summary and zero covariance in reported Weibull parameters; physical mechanism UNKNOWN')

def three_mf(mesh, path):
    ns = 'http://schemas.microsoft.com/3dmanufacturing/core/2015/02'
    ET.register_namespace('', ns)
    root = ET.Element('{' + ns + '}model', unit='millimeter')
    res = ET.SubElement(root, '{' + ns + '}resources')
    o = ET.SubElement(res, '{' + ns + '}object', id='1', type='model')
    m = ET.SubElement(o, '{' + ns + '}mesh')
    v = ET.SubElement(m, '{' + ns + '}vertices')
    f = ET.SubElement(m, '{' + ns + '}triangles')
    for p in mesh.vertices:
        ET.SubElement(v, '{' + ns + '}vertex', dict(zip(['x', 'y', 'z'], [format(float(z), '.10g') for z in p])))
    for p in mesh.faces:
        ET.SubElement(f, '{' + ns + '}triangle', dict(zip(['v1', 'v2', 'v3'], map(str, p))))
    b = ET.SubElement(root, '{' + ns + '}build')
    ET.SubElement(b, '{' + ns + '}item', objectid='1')
    with zipfile.ZipFile(path, 'w', zipfile.ZIP_DEFLATED) as z:
        for (name, data) in [('[Content_Types].xml', b'<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types"><Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/><Default Extension="model" ContentType="application/vnd.ms-package.3dmanufacturing-3dmodel+xml"/></Types>'), ('_rels/.rels', b'<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships"><Relationship Target="/3D/3dmodel.model" Id="rel0" Type="http://schemas.microsoft.com/3dmanufacturing/2013/01/3dmodel"/></Relationships>'), ('3D/3dmodel.model', ET.tostring(root))]:
            info = zipfile.ZipInfo(name, date_time=(2026, 10, 2, 0, 0, 0))
            info.compress_type = zipfile.ZIP_DEFLATED
            z.writestr(info, data)

def run():
    groups = derived_groups()
    studies = sorted({r['study'] for r in groups})
    slopes = {s: slope([r for r in groups if r['study'] == s]) for s in studies}
    m = meta(list(slopes.values()))
    folds = []
    for study in studies:
        others = [v for (s, v) in slopes.items() if s != study]
        fit = meta(others)
        rs = [r for r in groups if r['study'] == study]
        ref = min(rs, key=lambda r: abs(r['thickness_mm'] - 0.8))
        observed = []
        for r in rs:
            if r == ref:
                continue
            pred = band(fit, r['thickness_mm'] / ref['thickness_mm'])
            ratio = r['mean_N'] / ref['mean_N']
            observed.append(dict(source=r['locator'], observed_ratio=ratio, predicted=pred, log_error=float(np.log(pred['point_ratio'] / ratio)), covered=pred['interval_ratio'][0] <= ratio <= pred['interval_ratio'][1]))
        folds.append(dict(held_out_study=study, reference=ref['row_id'], training_study=next((s for s in studies if s != study)), rows=observed))
    X = np.zeros((len(groups), len(studies) + 1))
    y = np.log([r['mean_N'] for r in groups])
    V = np.diag([r['log_mean_variance'] for r in groups])
    for (j, study) in enumerate(studies):
        ix = [i for (i, r) in enumerate(groups) if r['study'] == study]
        x = np.log([groups[i]['thickness_mm'] for i in ix])
        X[ix, j] = 1
        X[ix, -1] = x
        V[np.ix_(ix, ix)] += m['between_study_slope_variance'] * np.outer(x, x)
    beta = np.linalg.solve(X.T @ np.linalg.solve(V, X), X.T @ np.linalg.solve(V, y))
    delta = abs(beta[-1] - m['slope'])
    source = read('inputs/X1_ROWS.json')
    family = []
    old = read('history/X1/FROZEN_PREDICTIONS.json')
    d1 = old['designs'][0]
    for (design, variant, thickness) in [('D1', 'thin', 0.8), ('M1', 'medium', 1.0), ('M2', 'thick', 1.5)]:
        row = next((r for r in source if r['id'] == f'L005_lo_M1_k6_{variant}_3Y'))
        b = band(m, thickness / 0.8)
        fe = {str(a): row['failure_proxy_axial_N' if a == 0 else 'failure_proxy_offaxis30_N'] / d1['fracture_load_proxy_axial_N' if a == 0 else 'fracture_load_proxy_offaxis30_N'] for a in [0, 30]}
        family.append(dict(design=design, source_id=row['id'], material='3Y', nominal_thickness_mm=thickness, reference='D1', literature_conditional_prediction=b, FE_conditional_ratios=fe, FE_absolute_diagnostic_N={'0': row['failure_proxy_axial_N'], '30': row['failure_proxy_offaxis30_N']}, gap=row['fit'], pulp_guard=row['digital_pulp_guard_pass']))
        if design != 'D1':
            mesh = trimesh.load_mesh(R / f'exports/{design}/crown.stl')
            three_mf(mesh, R / f'exports/{design}/crown.3mf')
    rmse = float(np.sqrt(np.mean([np.mean([r['log_error'] ** 2 for r in f['rows']]) for f in folds])))
    coverage = float(np.mean([np.mean([r['covered'] for r in f['rows']]) for f in folds]))
    gates = dict(log_ratio_RMSE=rmse <= 0.3, coverage=coverage >= 0.8, min_two_studies=len(studies) >= 2, upper_lower_width=all((r['literature_conditional_prediction']['upper_lower_ratio'] <= 2.5 for r in family)), independent_GLS_control_agreement=bool(delta < 1e-08))
    result = dict(round='R2', groups=groups, slopes=slopes, meta=m, folds=folds, study_balanced_log_ratio_RMSE=rmse, study_balanced_coverage=coverage, gates=gates, same_information_control=dict(kind='joint GLS, free study intercepts + random slope covariance', slope=float(beta[-1]), difference=delta, outcome='TIE' if delta < 1e-08 else 'DIFFERENT'), family=family, physical_mechanism='UNKNOWN; published forces include occlusal-origin or die failure', scale_measurement='No lab force measured; prospective bands are conditional on a future D1 reference mean')
    dump('raw/MATCHED_R2.json', result)
    state('R2_MATCHED_CONTRAST_DECIDED', gates, 'R3: angle-dependent double ratio cancels common batch scale and distinguishes geometry/contact interaction')
    (R / 'HANDOFF_R2.md').write_text('R2 matched same-tooth family D1/M1/M2; conditional ratio bands frozen in next R3 prediction artifact. Ordinary joint GLS agrees with meta slope: TIE. ' + str(gates) + '\nPublished occlusal fracture origins do not validate intaglio tensile FE. Next R3 changes observable to angle double ratio and freezes concrete experiment, rather than repeating a regression.\n')
    return result
if __name__ == '__main__':
    z = run()
    print(z['gates'])
    print(z['slopes'])
    print(z['same_information_control'])
    print([(r['design'], r['literature_conditional_prediction']['interval_ratio']) for r in z['family']])
