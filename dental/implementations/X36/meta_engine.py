"""Known-sampling-variance multilevel REML; no generic MixedLM residual rescale."""
import numpy as np
from scipy import linalg, optimize, stats

def design(rows, terms):
    cols = [np.ones(len(rows))]
    names = ['intercept']
    for (key, kind, reference) in terms:
        if kind == 'numeric':
            cols.append(np.array([r[key] for r in rows], float))
            names.append(key)
        else:
            values = sorted({r[key] for r in rows})
            assert reference in values, (key, reference, values)
            for v in values:
                if v != reference:
                    cols.append(np.array([r[key] == v for r in rows], float))
                    names.append(key + '=' + v)
    return (np.column_stack(cols), names)

def sampling(rows, rho=0.0, no_n=False):
    v = np.array([r['variance_no_n'] if no_n else r['variance'] for r in rows])
    S = np.diag(v)
    for (i, a) in enumerate(rows):
        for j in range(i):
            b = rows[j]
            if a['dataset'] == 'cement' and a['study'] == b['study'] and (a['arm'] == b['arm']):
                S[i, j] = S[j, i] = rho * np.sqrt(v[i] * v[j])
    return S

def fit(y, X, S, clusters, names, profile=True, dense_control=False, cell_residual=True):
    y = np.asarray(y, float)
    (n, p) = X.shape
    (_, sing, vt) = np.linalg.svd(X, full_matrices=False)
    rank = int(np.linalg.matrix_rank(X))
    labels = sorted(set(clusters))
    k = len(labels)
    support = {'n_rows': n, 'independent_clusters': k, 'p_fixed': p, 'rank': rank, 'confirmatory_support_pass': rank == p and k >= 10 * (p - 1), 'required_clusters_for_confirmatory': 10 * (p - 1), 'singular_values': sing.tolist()}
    if rank < p:
        aliases = [{name: float(v) for (name, v) in zip(names, z) if abs(v) > 1e-06} for z in vt[rank:]]
        return dict(status='UNKNOWN_RANK_DEFICIENT', support=support, nullspace=aliases)
    scale = max(float(np.std(y)), 0.001)
    ys = y / scale
    Ss = S / scale ** 2
    groups = np.array(clusters)
    B = (groups[:, None] == groups[None, :]).astype(float)
    I = np.eye(n)
    projector = I - X @ np.linalg.pinv(X)
    components = np.column_stack([(projector @ B @ projector).ravel(), projector.ravel()])
    component_rank = int(np.linalg.matrix_rank(components))
    variance_identification = dict(component_rank=component_rank, requested_components=2 if cell_residual else 1, study_component_annihilated=bool(np.linalg.norm(components[:, 0]) < 1e-09), note='Projected covariance component rank; low independent cluster count remains a separate limitation')
    upper = 100.0
    bounds = [(0, upper), (0, upper) if cell_residual else (0, 0)]

    def calc(q, dense=False):
        V = Ss + q[0] * B + q[1] * I
        if dense:
            (sign, lv) = np.linalg.slogdet(V)
            W = np.linalg.inv(V)
            WX = W @ X
            wy = W @ ys
            H = X.T @ WX
            (signh, lh) = np.linalg.slogdet(H)
            if sign <= 0 or signh <= 0:
                return (np.inf, None)
            C = np.linalg.inv(H)
            beta = C @ (X.T @ wy)
            r = ys - X @ beta
            wr = W @ r
        else:
            L = linalg.cho_factor(V, lower=True, check_finite=False)
            WX = linalg.cho_solve(L, X, check_finite=False)
            wy = linalg.cho_solve(L, ys, check_finite=False)
            H = X.T @ WX
            Hc = linalg.cho_factor(H, lower=True, check_finite=False)
            C = linalg.cho_solve(Hc, np.eye(p), check_finite=False)
            beta = linalg.cho_solve(Hc, X.T @ wy, check_finite=False)
            r = ys - X @ beta
            wr = linalg.cho_solve(L, r, check_finite=False)
            lv = 2 * np.log(np.diag(L[0])).sum()
            lh = 2 * np.log(np.diag(Hc[0])).sum()
        return (0.5 * (lv + lh + r @ wr), (beta, C, wr, V))

    def obj(q):
        try:
            return calc(q)[0]
        except (np.linalg.LinAlgError, ValueError):
            return 1e+100
    fits = [optimize.minimize(obj, start, method='L-BFGS-B', bounds=bounds, options={'ftol': 1e-12, 'gtol': 1e-07, 'maxiter': 500}) for start in [[0.1, 0.1], [1.0, 0.01], [0.01, 1.0]]]
    best = min(fits, key=lambda z: z.fun)
    refined = optimize.minimize(obj, best.x, method='Powell', bounds=bounds, options={'ftol': 1e-11, 'xtol': 1e-10, 'maxiter': 200})
    if refined.fun < best.fun:
        best = refined
    q = best.x
    for j in range(2):
        if q[j] < 1e-07:
            z = q.copy()
            z[j] = 0
            if obj(z) <= obj(q) + 1e-09:
                q = z
    (val, (b, C, wr, V)) = calc(q)
    b = b * scale
    C = C * scale ** 2
    (tau, omega) = q * scale ** 2
    scores = np.array([X[groups == g].T @ wr[groups == g] for g in labels]) * scale
    meat = scores.T @ scores
    bread = C / scale ** 2
    robust = bread @ meat @ bread * k / (k - 1) if k > 1 else np.full((p, p), np.nan)
    se = np.sqrt(np.maximum(np.diag(C), 0))
    rse = np.sqrt(np.maximum(np.diag(robust), 0))
    ct = stats.t.ppf(0.975, k - 1) if k > 1 else np.nan
    W = linalg.solve(S, np.eye(n), assume_a='pos')
    P = W - W @ X @ linalg.solve(X.T @ W @ X, X.T @ W, assume_a='pos')
    typical = (n - p) / np.trace(P)
    denom = tau + omega + typical
    coeff = [dict(term=name, estimate=float(bb), se_conditional=float(ss), ci95_normal_conditional=[float(bb - 1.95996398454 * ss), float(bb + 1.95996398454 * ss)], se_cluster_sandwich=float(rs), ci95_cluster_t=[float(bb - ct * rs), float(bb + ct * rs)], inference_status='ILLUSTRATIVE_ONLY' if not support['confirmatory_support_pass'] else 'MODEL_CONDITIONAL') for (name, bb, ss, rs) in zip(names, b, se, rse)]
    out = dict(status='EXPLORATORY_FIT' if not support['confirmatory_support_pass'] else 'MODEL_CONDITIONAL_FIT', support=support, coefficients=coeff, tau2_study=float(tau), omega2_cell=float(omega), I2_total_percent=float(100 * (tau + omega) / denom), I2_study_percent=float(100 * tau / denom), I2_cell_percent=float(100 * omega / denom), typical_sampling_variance=float(typical), conditional_row_QE=float(y @ P @ y), row_df=n - p, optimizer=dict(success=bool(best.success), objective=float(val), message=str(best.message), boundary_upper=bool(np.any(q > 0.999 * upper))), variance_component_identification=variance_identification, cell_residual_enabled=cell_residual, coefficient_covariance=C.tolist(), cluster_covariance=robust.tolist(), beta=b.tolist(), fitted=(X @ b).tolist(), residual=(y - X @ b).tolist())
    if profile:

        def prof(t):
            if not cell_residual:
                return obj([t, 0]) - val
            z = optimize.minimize_scalar(lambda w: obj([t, w]), bounds=(0, upper), method='bounded', options={'xatol': 1e-09})
            return min(z.fun, obj([t, 0])) - val
        cutoff = stats.chi2.ppf(0.95, 1) / 2
        lo = 0.0 if prof(0) <= cutoff else optimize.brentq(lambda a: prof(a) - cutoff, 0, q[0], xtol=1e-08)
        hi = None if prof(upper) <= cutoff else optimize.brentq(lambda a: prof(a) - cutoff, q[0], upper, xtol=1e-08)
        out['tau2_profile_ci95'] = [float(lo * scale ** 2), None if hi is None else float(hi * scale ** 2)]
        out['profile_note'] = 'Gaussian likelihood profile, asymptotic cutoff; small-k coverage unverified; null upper means exceeds fixed search bound'
    if dense_control:

        def dobj(z):
            try:
                return calc(z, dense=True)[0]
            except (np.linalg.LinAlgError, ValueError):
                return 1e+100
        ref = optimize.minimize(dobj, q + [0.02, 0.02], method='L-BFGS-B', bounds=bounds, options={'ftol': 1e-12, 'gtol': 1e-07, 'maxiter': 500})
        db = calc(ref.x, True)[1][0] * scale
        out['dense_control'] = dict(objective_at_same_theta_abs_error=float(abs(calc(q, True)[0] - val)), refit_objective_abs_error=float(abs(ref.fun - val)), beta_max_abs_error=float(np.max(abs(db - b))), variance_relative_max_error=float(np.max(abs(ref.x - q) / (1 + q))), success=bool(ref.success))
    return out

def egger(y, se):
    y = np.asarray(y)
    se = np.asarray(se)
    k = len(y)
    if k < 3:
        return dict(status='NOT_ESTIMABLE', k=k)
    X = np.column_stack([np.ones(k), 1 / se])
    z = y / se
    b = np.linalg.lstsq(X, z, rcond=None)[0]
    r = z - X @ b
    cov = np.linalg.inv(X.T @ X) * (r @ r) / (k - 2)
    s = np.sqrt(cov[0, 0])
    t = b[0] / s
    return dict(k=k, intercept=float(b[0]), se=float(s), p_two_sided=float(2 * stats.t.sf(abs(t), k - 2)), status='DIAGNOSTIC_ONLY_NOT_PUBLICATION_BIAS_TEST', formal_gate_pass=k >= 10, limitation='Unlike estimands and source selection can create asymmetry; unknown within-study covariance; no publication-bias conclusion')
