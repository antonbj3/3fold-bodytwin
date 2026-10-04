"""Execute identifiability first; do not disguise a singular fit with a pseudoinverse."""
import numpy as np
from scipy import linalg, optimize, stats
from common import read, dump, state

def matrix(d, section=False, material_levels=None):
    materials = material_levels or sorted({r['material'] for r in d})
    names = ['intercept'] + [f'material_{m}' for m in materials[1:]] + ['log_area']
    x = np.array([[1.0] + [float(r['material'] == m) for m in materials[1:]] + [np.log(r['area_mm2'] / 9.0)] for r in d])
    if section:
        x = np.column_stack((x, [np.log(r['height_mm'] / np.sqrt(r['area_mm2'])) for r in d]))
        names += ['log_height_over_sqrt_area']
    return (x, names)

def reml(d, section=False):
    (X, names) = matrix(d, section)
    y = np.log([r['mean_N'] for r in d])
    v = np.array([(r['sd_N'] / r['mean_N']) ** 2 / r['n'] for r in d])
    studies = sorted({r['study'] for r in d})
    Z = np.array([[r['study'] == s for s in studies] for r in d], float)
    rank = int(np.linalg.matrix_rank(X))
    out = dict(model='section' if section else 'area', n_groups=len(d), n_studies=len(studies), columns=names, rank=rank, n_columns=X.shape[1], full_rank=rank == X.shape[1], study_random_effect_estimable=len(studies) >= 3)
    if rank < X.shape[1] or len(studies) < 3:
        out.update(outcome='UNKNOWN_NOT_IDENTIFIABLE', beta=None, tau2=None, loso=None)
        return out

    def fit(tau):
        V = np.diag(v) + tau * Z @ Z.T
        Vi = linalg.inv(V)
        G = X.T @ Vi @ X
        C = linalg.inv(G)
        b = C @ X.T @ Vi @ y
        e = y - X @ b
        obj = np.linalg.slogdet(V)[1] + np.linalg.slogdet(G)[1] + e @ Vi @ e
        return (obj, b, C)
    opt = optimize.minimize_scalar(lambda t: fit(t)[0], bounds=(0.0, 4.0), method='bounded')
    tau = float(opt.x) if fit(opt.x)[0] < fit(0.0)[0] else 0.0
    (_, b, C) = fit(tau)
    out.update(outcome='FIT_AVAILABLE', beta=dict(zip(names, b)), tau2=tau, covariance=C, loso=[])
    for s in studies:
        train = [r for r in d if r['study'] != s]
        test = [r for r in d if r['study'] == s]
        levels = sorted({r['material'] for r in d})
        if {r['material'] for r in test} - {r['material'] for r in train}:
            out['loso'].append(dict(study=s, outcome='UNKNOWN_UNSEEN_MATERIAL'))
            continue
        f = reml_no_loso(train, section, levels)
        if f is None:
            out['loso'].append(dict(study=s, outcome='UNKNOWN_UNSEEN_OR_RANK'))
            continue
        (Xt, _) = matrix(test, section, levels)
        if Xt.shape[1] != len(b):
            out['loso'].append(dict(study=s, outcome='UNKNOWN_UNSEEN_MATERIAL'))
            continue
        for (r, xx) in zip(test, Xt):
            mu = xx @ f[0]
            se = np.sqrt(xx @ f[1] @ xx + f[2] + (r['sd_N'] / r['mean_N']) ** 2 / r['n'])
            out['loso'].append(dict(id=r['id'], pred_N=np.exp(mu), log_error=mu - np.log(r['mean_N']), covered=abs(mu - np.log(r['mean_N'])) <= 1.96 * se))
    return out

def reml_no_loso(d, section, levels=None):
    (X, _) = matrix(d, section, levels)
    y = np.log([r['mean_N'] for r in d])
    v = np.array([(r['sd_N'] / r['mean_N']) ** 2 / r['n'] for r in d])
    if np.linalg.matrix_rank(X) < X.shape[1] or len({r['study'] for r in d}) < 3:
        return None
    Z = np.array([[r['study'] == s for s in sorted({r['study'] for r in d})] for r in d], float)

    def f(t):
        V = np.diag(v) + t * Z @ Z.T
        Vi = linalg.inv(V)
        C = linalg.inv(X.T @ Vi @ X)
        b = C @ X.T @ Vi @ y
        e = y - X @ b
        return (np.linalg.slogdet(V)[1] + np.linalg.slogdet(X.T @ Vi @ X)[1] + e @ Vi @ e, b, C)
    o = optimize.minimize_scalar(lambda t: f(t)[0], bounds=(0.0, 4.0), method='bounded')
    t = float(o.x) if f(o.x)[0] < f(0.0)[0] else 0.0
    (_, b, C) = f(t)
    return (b, C, t)

def run():
    d = read('raw/MEASUREMENTS.json')
    p = read('PREREG_R1.json')
    strict = [r for r in d if r['connector_observation'] == 'uncensored_connector_failure']
    counts = dict(system_groups=len(d), system_studies=len({r['study'] for r in d}), strict_connector_groups=len(strict), strict_connector_studies=len({r['study'] for r in strict}), distinct_areas=len({r['area_mm2'] for r in d}))
    out = dict(claim_type='information_link', counts=counts, candidate=reml(d, True), equally_informed_control=reml(d, False), strict_connector_fit='UNKNOWN_NO_UNCENSORED_EVENTS', gate_minimum_studies=counts['system_studies'] >= p['gates']['minimum_studies'], gate_area_variation=counts['distinct_areas'] >= p['gates']['minimum_distinct_areas'], outcome='ABSOLUTE_CONNECTOR_TRANSFER_UNKNOWN', external_referent=dict(kind='independent_measurement', locator='https://doi.org/10.4047/jap.2023.15.4.171#T2', compared_quantity='six system fracture load groups plus experimentally observed origin', refutes_us=True))
    dump('raw/R1_FIT.json', out)
    state('R1_complete', 'FAIL_identifiability_and_mode: one study, one area, no connector origin', 'freeze R2 competing-mode lower bound and within-study response; build fresh 3-unit FE port')
    print('R1: full-rank gate fails; study transfer and absolute connector capacity UNKNOWN')
    return out
if __name__ == '__main__':
    run()
