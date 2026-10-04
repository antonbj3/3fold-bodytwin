from common import *
sys.path.insert(0, str(DATA / 'participant_code'))
from legacy.geometry import area_weights
from legacy.generators import submit
from legacy.checks import all_checks
from quality import roof_metrics, nondominated
from scipy.sparse import eye

def run():
    from integrity import verify
    verify()
    xy = np.concatenate([np.array([[0, 0], [1, 0], [1, 1], [0, 1]], float) + [2 * j, 0] for j in range(4)])
    faces = np.concatenate([np.array([[0, 1, 2], [0, 2, 3]]) + 4 * j for j in range(4)])
    gap = np.repeat([0.0625, 0.0625, 0.25, 0.25], 4)
    ref = 3 - gap
    t = dict(task_id='pareto_probe', frame='fixture_mm', family='molar_crown', xy=xy, faces=faces, preparation_z=np.zeros(16), ceiling=np.full(16, 3.0), weights=area_weights(xy, faces), A=eye(16, format='csr'), obstacle_b=np.full(16, 3.0), requirements=dict(film_min_mm=0.02, film_max_mm=0.15, wall_mm=0.25))
    inner = np.full(16, 0.0625)
    a = ref - 0.0625
    b = ref.copy()
    b[gap > 0.1] -= 0.5
    qa = roof_metrics(t, a, inner, ref)
    qb = roof_metrics(t, b, inner, ref)
    va = all_checks(t, submit(t, a, inner))['validity']
    vb = all_checks(t, submit(t, b, inner))['validity']
    front = nondominated([[qa['anatomy_rmse_mm'], qa['contact_symdiff_mm2']], [qb['anatomy_rmse_mm'], qb['contact_symdiff_mm2']]], [0.0, 0.0])
    series = np.array([[0.03125, 0.03125], [0.125, 0.125]])
    parallel = series.T.copy()
    mean_error = float(abs(series.mean() - parallel.mean()))
    hist_error = float(np.max(abs(np.sort(series.ravel()) - np.sort(parallel.ravel()))))

    def flow(h):
        return float(np.sum(1 / np.sum(1 / h ** 3, axis=0)))
    (qs, qp) = (flow(series), flow(parallel))
    checks = dict(both_L1=va == vb == 'PASS', anatomy_tradeoff=qa['anatomy_rmse_mm'] < qb['anatomy_rmse_mm'], contact_tradeoff=qa['contact_symdiff_mm2'] > qb['contact_symdiff_mm2'], both_nondominated=bool(front.all()), wrong_single_winner_rejected=front.tolist() != [True, False], cement_exact_mean=mean_error == 0, cement_exact_histogram=hist_error == 0, flow_discriminates=qp / qs > 2)
    out = dict(round='R6', claim_type='capability', decision='ACTUAL_METRIC_TRADEOFF_PASS' if all(checks.values()) else 'FAIL', checks=checks, design_A=dict(L1=va, quality=qa), design_B=dict(L1=vb, quality=qb), cement_fields_mm={'series': series, 'parallel': parallel}, mean_identity_error_mm=mean_error, histogram_identity_error_mm=hist_error, downstream_conductance_ratio=qp / qs, external_referent=dict(kind='closed_form', locator='https://doi.org/10.1098/rstl.1886.0005', compared_quantity='Series/parallel discrete Reynolds conductance; dental surfaces here are constructed probes', refutes_us=True), scope='An actual pair of L1-pass surfaces produces opposing anatomy/contact metrics. Pareto set only within these candidates, not a global optimality claim. Source-faithful ideal would dominate both; no method win.')
    dump(ROOT / 'rounds/R6.json', out)
    print(out['decision'])
if __name__ == '__main__':
    run()
