"""Case-cluster uncertainty; tasks never count as independent patients."""
import numpy as np, math

def cluster_summary(matrix, seed=6102, replicates=2000):
    a = np.asarray(matrix, float)
    if a.ndim != 2 or not np.isfinite(a).all() or np.any((a < 0) | (a > 1)):
        raise ValueError('bounded cluster outcome matrix required')
    (k, m) = a.shape
    if k < 2:
        raise ValueError('at least two case clusters required')
    y = a.mean(1)
    mu = float(y.mean())
    rng = np.random.default_rng(seed)
    samples = np.empty(replicates)
    for i in range(replicates):
        samples[i] = y[rng.integers(0, k, k)].mean()
    ci = np.quantile(samples, [0.025, 0.975])
    msb = m * float(np.var(y, ddof=1))
    msw = float(np.sum((a - y[:, None]) ** 2) / (k * (m - 1))) if m > 1 else 0.0
    den = msb + (m - 1) * msw
    raw = (msb - msw) / den if den > 0 else None
    rho = max(0.0, raw) if raw is not None else None
    n = k * m
    neff = n / (1 + (m - 1) * rho) if rho is not None else None
    radius = math.sqrt(math.log(40) / (2 * k))
    return dict(mean_case_pass_fraction=mu, case_cluster_bootstrap_95=ci.tolist(), bootstrap_replicates=replicates, bootstrap_degenerate=bool(ci[0] == ci[1]), bounded_case_Hoeffding_95=[max(0.0, mu - radius), min(1.0, mu + radius)], case_clusters=k, independent_patient_case_records=k, task_rows=n, tasks_per_case=m, icc_anova_raw=raw, icc_nonnegative=rho, n_eff_tasks_estimated=neff, n_eff_tasks_sensitivity={str(r): n / (1 + (m - 1) * r) for r in [0.0, 0.1, 0.5, 1.0]}, resolution='POPULATION', uncertainty_scope='Marginal within-dataset case sampling only; scanner, label, cross-dataset identity and domain bias not included. K case units, never n_eff_task patients.')

def paired_case_difference(a, b, seed=6102):
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    if a.shape != b.shape:
        raise ValueError('paired case shapes differ')
    d = (a - b).mean(1)
    rng = np.random.default_rng(seed)
    v = [d[rng.integers(0, len(d), len(d))].mean() for _ in range(2000)]
    return dict(mean=float(d.mean()), case_bootstrap_95=np.quantile(v, [0.025, 0.975]).tolist(), case_clusters=len(d), interpretation='descriptive paired capability contrast, no algorithm-superiority gate')
