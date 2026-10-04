from patient_geometry import P, write
from fractions import Fraction as Q
import json, numpy as np, time, resource

def main():
    tick = time.perf_counter()
    r = json.load(open(P / 'rounds/R1/results.json'))
    pairs = r['pairs']
    g = [Q(p['minimum_axial_gap_mm']) for p in pairs]
    d = max(g) + 1
    it = [i for (i, p) in enumerate(pairs) if p['lower_fdi'] == 36]
    other = next((i for (i, p) in enumerate(pairs) if p['lower_fdi'] != 36))

    def world(ia, extra):
        f = [Q(1)] * len(g)
        f[ia] += extra
        c = [(d - x) / y for (x, y) in zip(g, f)]
        res = [x - d + z * y for (x, z, y) in zip(g, c, f)]
        return dict(total=sum(f), tooth36=sum((f[i] for i in it)), forces=f, compliances=c, max_residual=max((abs(x) for x in res)), injected_c_error_refuted=g[ia] - d + (c[ia] + Q(1, 1000)) * f[ia] != 0)
    A = world(it[0], 100 - len(g))
    B = world(other, 100 - len(g))
    C = world(it[0], 100)
    D = world(it[1], 100)
    pts = [np.array(pairs[i]['lower_point_cbct_mm']) for i in it]
    mdiff = np.cross(pts[1] - pts[0], [0, 0, 100.0])
    out = dict(round='R4', claim_type='capability', closure='Hypothetical independent positive linear normal-contact compliances; not calibrated tooth/PDL model', same_geometry_and_virtual_closure=True, virtual_closure_mm=float(d), force_summary_test=dict(total_A_N=float(A['total']), total_B_N=float(B['total']), exact_summary_identity_error_N=float(abs(A['total'] - B['total'])), float64_identity_error_N=abs(float(A['total']) - float(B['total'])), tooth36_A_N=float(A['tooth36']), tooth36_B_N=float(B['tooth36']), downstream_difference_N=float(abs(A['tooth36'] - B['tooth36'])), minimal_extension='Measured FDI36 force channel for this scalar force decision'), wrench_summary_test=dict(total_C_N=float(C['total']), total_D_N=float(D['total']), tooth36_C_N=float(C['tooth36']), tooth36_D_N=float(D['tooth36']), exact_summary_identity_error_N=float(abs(C['total'] - D['total']) + abs(C['tooth36'] - D['tooth36'])), float64_identity_error_N=abs(float(C['total']) - float(D['total'])) + abs(float(C['tooth36']) - float(D['tooth36'])), moment_difference_Nmm=mdiff.tolist(), moment_difference_norm_Nmm=float(np.linalg.norm(mdiff)), minimal_extension='Two transverse coordinates of force-weighted centroid for normal resultant moment; pressure field for local stress'), all_compliances_positive=all((x > 0 for w in [A, B, C, D] for x in w['compliances'])), exact_equilibrium_error=max((float(w['max_residual']) for w in [A, B, C, D])), all_injected_c_errors_rejected=all((w['injected_c_error_refuted'] for w in [A, B, C, D])), worlds=[dict(forces_N=[str(x) for x in w['forces']], compliance_mm_per_N=[str(x) for x in w['compliances']]) for w in [A, B, C, D]], external_referent=dict(kind='closed_form', locator='g-d+c*lambda=0 with c>0, lambda>0; rational arithmetic on stored digital gaps', compared_quantity='Conditional exact contact equilibrium; no physical patient force', refutes_us=True), cost=dict(wall_s=time.perf_counter() - tick, maxrss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024))
    out['all_pass'] = out['all_compliances_positive'] and out['exact_equilibrium_error'] == 0 and out['all_injected_c_errors_rejected'] and (out['force_summary_test']['float64_identity_error_N'] == 0) and (out['wrench_summary_test']['float64_identity_error_N'] == 0)
    write(P / 'rounds/R4/results.json', out)
    print(json.dumps({k: out[k] for k in ['all_pass', 'force_summary_test', 'wrench_summary_test']}, indent=2))
if __name__ == '__main__':
    main()
