"""X46: frozen, source-bound fiber probes. No source code is imported with side effects."""
import ast, collections, datetime, hashlib, json, math, pathlib, resource, struct, sys, time
from fractions import Fraction as F
import numpy as np
from scipy.optimize import brentq, linprog
HERE = pathlib.Path(__file__).resolve().parent
SAMPLE = json.loads((HERE / 'SAMPLE.json').read_text())['rows']
SOURCES = []

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def read(i, n):
    p = pathlib.Path(SAMPLE[i]['source_alias']) / n
    SOURCES.append({'path': str(p), 'sha256': sha(p)})
    return p

def functions(i, n, names, extra=None):
    p = read(i, n)
    txt = p.read_text()
    tree = ast.parse(txt)
    env = {'np': np, 'math': math, 'json': json, 'F': F, 'linprog': linprog, 'brentq': brentq}
    env.update(extra or {})
    nodes = [node for node in tree.body if isinstance(node, (ast.FunctionDef, ast.ClassDef)) and node.name in names]
    assert len(nodes) == len(names), (p, names)
    exec(compile(ast.Module(body=nodes, type_ignores=[]), str(p), 'exec'), env)
    return env

def jread(i, n):
    return json.loads(read(i, n).read_text())

def native(v):
    if isinstance(v, F):
        return {'fraction': str(v), 'float': float(v)}
    if isinstance(v, np.ndarray):
        return v.tolist()
    if isinstance(v, np.generic):
        return v.item()
    if isinstance(v, dict):
        return {str(k): native(x) for (k, x) in v.items()}
    if isinstance(v, (tuple, list)):
        return [native(x) for x in v]
    return v

def bits(x):
    return np.asarray(x, dtype=np.float64).tobytes()

def pair(name, a, b, sa, sb, qa, qb, units, scope, minimum, level='PHENOMENOLOGICAL', time_kind='SIMULTANEOUS'):
    aa = np.asarray(sa, dtype=float)
    bb = np.asarray(sb, dtype=float)
    err = float(np.max(np.abs(aa - bb)))
    delta = float(np.max(np.abs(np.asarray(qa, dtype=float) - np.asarray(qb, dtype=float))))
    return dict(name=name, state_A=native(a), state_B=native(b), summary_A=native(sa), summary_B=native(sb), identity_error=err, summary_bytes_identical=bits(sa) == bits(sb), downstream_A=native(qa), downstream_B=native(qb), downstream_difference=delta, units=units, resolution=level, timescale=time_kind, scope=scope, sufficiency=('REJECTED' if delta != 0 else 'ACCEPTED_FOR_THIS_PAIR') if err == 0 and bits(sa) == bits(sb) else 'INVALID_IDENTITY', minimum_extension=minimum, universal_hold=False)
RECIPES = [('rate offset', 'same crack/geometry exponents, two ramp speeds', 'REJECTED', 'The rate sign and geometry correction have the wrong sign in the executable derivation; rate-invariance within its power-law regime survives.', 'calc_rate_port.py:34-48; RESULTS.md sections 3-5'), ('cone certificate', 'equal scalar reduction and thickness, different normal cones', 'CORRECTION', 'Gordan excludes strictly positive margin, not every nonzero nonnegative direction. Retain the 2D scalar triangle, retract general cone iff/finite-path certification.', 'cone_cert.py:55-85; RESULTS.md section 1'), ('kinetic service percentile', 'equal cycles, different physical exposure time; audit flaw quantile', 'REJECTED', 'The code uses the lower 5% flaw quantile for a 5% failure strength, reversing the flaw tail. Its zero-threshold safe floor also returns infinity. Service-percentile claims need reconstruction.', 'src/kinetic.py:229-244; RESULTS.md sections 2-3'), ('transformed section', 'same neutral axis, different width; audit parallel-axis coefficient', 'REJECTED', 'The code multiplies already transformed area by E/Ev again in the parallel-axis term. The near-perfect external ratio is not the stated transformed-section law.', 'shear_lag.py:25-43; RESULTS.md sections 1-3'), ('material pooling', 'same material means, opposite protocol contrast', 'ACCEPTED', 'The source correctly reports no observed reversals only in its small matched sample and preserves failed controls; pooled means do not decide unobserved protocol contrasts.', 'reversal_result.json; RESULTS.md sections 5-8'), ('thickness exponent', 'same reference amplitude and bending exponent, different support cap', 'CORRECTION', 'The n=2 held-out disc result is real within the selected shape test, but the alleged ISO n=1.5 control is not the published biaxial formula; general support transfer is not validated.', 'theory_and_consumer.py; fit_and_test.py; RESULTS.md sections 2-5'), ('table arithmetic', 'equal row and column means, different individual cell', 'CORRECTION', 'Aggregate consistency is a useful error detector, not a segmentation certificate: different cell matrices pass all named means. Preserve the extracted sample pending cell-level verification.', 'extract_crown_loads.py:126-163; verify_crown_loads.py; RESULTS.md arithmetic gate'), ('flattened binding', 'same flattened string from two table partitions', 'CORRECTION', 'Missing schema and zero pairs found are supported; a particular regex fallback returning zero does not prove that every real anchor is unrecoverable.', 'pair_recoverability.py:54-104; RESULTS.md sections 1-2'), ('anatomical budget', 'same attachment sum, different compartment split; audit sign of containment repair', 'REJECTED', 'Containment still predicts non-positive thin-minus-thick loss under shared design; assigning one fitted height per cohort does not repair the published positive ordering as a held-out prediction.', 'budget4.py:59-108; RESULTS.md sections 2-5'), ('PDL thickness', 'equal arithmetic and harmonic thickness, different minimum', 'ACCEPTED', 'Uniform strain cancellation and local extreme dependence hold. The original stored fiber has exactly zero summary errors; the new source replay confirms insufficiency at another positive fiber. The reported random-field percentiles remain conditional simulations.', 'suff_test.py:19-38; pdl_thickness_calc.py:77; RESULTS.md C1-C3'), ('event frequency', 'same annual event count, different within-event rate; audit finite SLS wings', 'REJECTED', 'Annual exceedance count is not the dynamic loading frequency; finite SLS wings are approximate, not exactly elastic, and the relaxed/glassy labels are reversed. Missing relaxation fit remains a correct UNKNOWN.', 'pdl_rate_axis.py:41-75; RESULTS.md sections 1-2'), ('Miner rate', 'same event count, reordered event times under fixed damage per event', 'ACCEPTED', 'The one computable decision is explicitly conditional on the fitted constant event damage; four unavailable thresholds remain unknown. Counts suffice only within that closure.', 'flipcount.py:30-38; RESULTS.md Table 1'), ('joint inverse', 'equal torque and stick slip under reciprocal mu-pressure change', 'ACCEPTED', 'The joint map has a real positive-parameter fiber. The result is non-identifiability of this stipulated family, not a measured contact-law refutation.', 'joint_law.py:36-57,63-142; RESULTS.md T1-T3'), ('finite sensitivity', 'same normalized sensitivity, different absolute response', 'REJECTED', 'A finite log perturbation below ln2 is not structural non-identifiability. Nonzero source derivatives distinguish parameters with sufficiently precise observations. Retain sensitivity ranking as conditional.', 'study.py:58-94; contact_joint_sensitivity.py; RESULTS.md sections 1-3'), ('preload history', 'same endpoints/duration, different unloading trajectory; audit Basquin exponent', 'REJECTED', 'The fitted lnN-vs-lnF slope must enter damage as its negative, not its negative reciprocal. The table pre-load is cyclic minimum external force, not measured retained screw clamp force.', 'matched_calc.py:115-125,149-170; RESULTS.md sections 1-2'), ('contact solve cache', 'same active set/matrix, different RHS; Woodbury vs full solve', 'CORRECTION', 'The finite matrix identity survives; active sets alone do not specify the response. Reported speed is against repeated factorisation, not a cached strongest control; source FE/obs ratios cannot reconstruct a measured angular displacement ratio.', 'cm1_repair.py:243-303; RESULTS.md sections 1-2'), ('screw product', 'same multiplicative torque coefficient, different diameter', 'CORRECTION', 'Cancellation in the supplied multiplicative law holds. Screw helix angle cannot be excluded by an implant-abutment cone/taper-angle catalogue.', 'thread_band_test.py:41-54; RESULTS.md R4/R9'), ('screw reachable set', 'same endpoint ratio, different eta-diameter scale; exact thread root', 'REJECTED', 'The headline says no joint exists at any declared diameter, while its own corrected thread-law section supplies an in-domain positive root and a compatible diameter.', 'band_inversion.py:56-75,STEP 3; band_inversion_results.json'), ('wear clamp lump', 'same exponential clamp-decay coefficient, different wear rate', 'ACCEPTED', 'The source correctly distinguishes lumped clamp loss from material wear and preserves lack of physical calibration.', 'partA_algebra.py; RESULTS.md Parts A-C'), ('local derivative', 'same F0 and local S, different finite monotone preload response', 'CORRECTION', 'The sign result for monotone F_lim survives; exp(-gNS) is exact only for a globally constant elasticity, not an arbitrary local sensitivity. Provide an integral or enclosure.', 'RESULTS.md EQ-1; step1_matched_cell.py:81-94'), ('strain shape', 'same reported unweighted percentiles and weighted mean, different dose', 'CORRECTION', 'The source labels G_pct as volume weighted, but calls unweighted np.quantile. Moment-dose needs its own weighted moment; no empirical Gamma anchor exists.', 'calc_gamma_converged.py:35-46; RESULTS.md Table 3'), ('noise support', 'same scatter radius, discrete vs interval error alphabet', 'ACCEPTED', 'Exact set-membership recourse depends on the error support, not just radius. The source retains that support and correctly scopes nesting to its robust closure.', 'nesting_theorem.py:71-157; RESULTS.md theorem'), ('process stress bound', 'same measured c0 with opposite hidden stress; audit R=0 branch', 'CORRECTION', 'A source stress bound of zero after enough independent channels is exact reconstruction, not zero signal. Absolute metrology amplitude and MPa/um transfer remain unknown; 436x is conditional on a constructed amplitude.', 'identifiability.py zero_branch; inputs/SHELL/metrics.json K7/K8; RESULTS.md sections 1-3'), ('design correlation', 'same raw correlation and N, different marginal information', 'CORRECTION', 'Zero raw cross-moment does not fix precision. The full Fisher matrix is needed. Retain fixed-layout arithmetic, but qualify continuous uncertainty and actual zero-porosity feasibility.', 'probe1.py:56-85; probe2.py:246-275; RESULTS.md primary metric'), ('optical projection', 'same projected thickness and design mean, different minimum wall', 'CORRECTION', 'Projected thickness is insufficient for strength. Alignment non-identifiability survives for feasible priors; measured fit deviation is not a measured wall-anisotropy prior and infeasible priors cannot satisfy universal subset claims.', 'scripts/run.py:178-221; results.json prior branches')]

def probe(i):
    detail = {}
    scope = 'synthetic fiber of the stated operation, not a new physical observation'
    if i == 0:
        p = read(i, 'calc_rate_port.py')
        (m, ng, nk) = (20.0, 2.0, 3.0)
        (a, b) = (1.0, 2.0)
        (qa, qb) = (a ** (1 / (m + 1)), b ** (1 / (m + 1)))
        detail = {'source_strength_ratio_at_double_speed': 2 ** (-1 / 21), 'independent_strength_ratio_at_double_speed': qb / qa, 'source_n_eff': ng - (nk - ng) / (m + 1), 'correct_n_eff': ng + (nk - ng) / (m + 1), 'derivation': 'a^(-m/2) da = C sigmadot^m t^m dt => sigma_f^(m+1)/sigmadot = positive constant; log sigma_f has positive log sigmadot coefficient', 'reference': 'https://nvlpubs.nist.gov/nistpubs/specialpublications/nist.sp.960-16e2.pdf Figure 7.32c'}
        return (pair('m,ng,nk -> absolute fracture load', {'v': a}, {'v': b}, [m, ng, nk], [m, ng, nk], qa, qb, 'normalized load', scope, 'one ramp-speed coordinate with corrected sign'), detail)
    if i == 1:
        g = functions(i, 'cone_cert.py', ['gordan_certificate', 'cone_margin'])
        A = np.array([[1.0, 0.0], [-1.0, 0.0]])
        B = np.array([[1.0, 0.0], [-1.0, 0.0], [0.0, 1.0], [0.0, -1.0]])
        wa = g['gordan_certificate'](A)
        wb = g['gordan_certificate'](B)
        detail = {'source_witness_A': wa, 'source_witness_B': wb, 'nonzero_nonnegative_direction_A': [0.0, 1.0], 'A_margin_of_nonzero_direction': g['cone_margin'](A, np.array([0.0, 1.0])), 'B_nonzero_direction_exists': False, 'gordan_distinction': 'both witnesses exclude strict margin, only B excludes all nonzero weak directions', 'domain_limit': 'normal-cone theorem stress test; not a changed measured tooth nor a refutation of the fixed scalar triangle itself'}
        return (pair('scalar reduction budget -> insertion cone', A, B, [0.6, 1.0, 2.6], [0.6, 1.0, 2.6], 1.0, 0.0, 'existence boolean', scope, 'retain normalized facet constraints and strict/weak motion contract'), detail)
    if i == 2:
        read(i, 'src/kinetic.py')
        N = 1024.0
        n = 11.0
        qa = (N / 1.0) ** (-1 / n)
        qb = (N / 2.0) ** (-1 / n)
        p = 0.05
        b = 2.0
        lo = (-math.log(1 - p)) ** (1 / b)
        hi = (-math.log(p)) ** (1 / b)
        detail = {'source_effective_flaw_5pct': lo, 'correct_upper_flaw_5pct': hi, 'actual_failure_probability_at_source_flaw': math.exp(-lo ** b), 'source_zero_Kinf_floor': 'infinity', 'correct_zero_Kinf_floor': 0.0, 'tail_logic': 'for decreasing strength vs flaw size, P_fail=P(a>=a_threshold); solve exp(-(a/ac0)^b)=p, not 1-p', 'domain_limit': 'cycle count alone fails when f is permitted to vary; fixed-f identity remains valid'}
        return (pair('cycles -> constant-threshold service strength', {'N': N, 'f': 1.0}, {'N': N, 'f': 2.0}, N, N, qa, qb, 'normalized stress', scope, 'physical exposure N/f; upper flaw-tail quantile for percentile'), detail)
    if i == 3:
        g = functions(i, 'shear_lag.py', ['section', 'sigma_ceramic', 'support_factor'])
        kw = dict(E_v=95.0, E_s=18.0, E_c=9.2, d_v=0.5, t_c=0.05, H_s=2.5)
        a = g['section'](**kw, width=1.0)
        b = g['section'](**kw, width=2.0)

        def corrected(d):
            layers = [(2.5, 18.0, 1.25), (0.05, 9.2, 2.525), (d, 95.0, 2.55 + d / 2)]
            areas = [h * e / 95 for (h, e, z) in layers]
            zbar = sum((A * l[2] for (A, l) in zip(areas, layers))) / sum(areas)
            I = sum((e / 95 * (h ** 3 / 12 + h * (z - zbar) ** 2) for (h, e, z) in layers))
            c = max(2.55 + d - zbar, zbar - 2.55)
            return c / I
        wrong05 = g['sigma_ceramic'](**kw)['sigma']
        kw['d_v'] = 1.0
        wrong1 = g['sigma_ceramic'](**kw)['sigma']
        detail = {'source_I_width1': a['Itr'], 'correct_I_width1': None, 'source_R': wrong05 / wrong1, 'correct_R': corrected(0.5) / corrected(1.0), 'correct_sigma05': corrected(0.5), 'correct_sigma1': corrected(1.0), 'bug': 'a is already transformed at line 34; line 40 applies n twice to a*(zc-zbar)^2', 'reference': 'doi:10.1055/s-0042-1757910 Table 1; parallel-axis theorem with transformed areas'}
        return (pair('neutral axis -> section stress', {'width': 1.0}, {'width': 2.0}, a['zbar'], b['zbar'], 1 / a['Itr'], 1 / b['Itr'], 'inverse section moment, 1/mm^4', scope, 'section second moment or width under fixed stack'), detail)
    if i == 4:
        read(i, 'reversal.py')
        A = np.array([[4.0, 3.0], [4.0, 3.0]])
        B = np.array([[6.0, 1.0], [2.0, 5.0]])
        return (pair('material pooled means -> protocol ordering', A, B, A.mean(0), B.mean(0), A[1, 0] - A[1, 1], B[1, 0] - B[1, 1], 'N contrast', scope, 'per-protocol paired contrast; source already retains this'), {'original_claim_scope': '3 studies,16 comparisons; no full-domain invariance claimed'})
    if i == 5:
        read(i, 'theory_and_consumer.py')
        a = {'C': 1.0, 'cap': 2.0}
        b = {'C': 1.0, 'cap': 4.0}
        return (pair('P at t=1 and n=2 -> P at t=2', a, b, [1.0, 2.0], [1.0, 2.0], min(4.0, a['cap']), min(4.0, b['cap']), 'normalized load', scope, 'one support-cap scalar/branch; n=2 alone is not enough'), {'baseline_issue': 'source labels 1.5 an ISO 6872 convention; primary paper biaxial equation uses t squared', 'source_table': 'n=2 shape MSE 0.0127796563; parent 1.398 0.439028... (source sample)'})
    if i == 6:
        read(i, 'extract_crown_loads.py')
        read(i, 'verify_crown_loads.py')
        A = np.full((4, 4), 8.0)
        B = A.copy()
        B[0, 0] += 1
        B[0, 1] -= 1
        B[1, 0] -= 1
        B[1, 1] += 1
        sa = np.r_[A.mean(0), A.mean(1)]
        sb = np.r_[B.mean(0), B.mean(1)]
        return (pair('all row/column means -> material-design cell contrast', A, B, sa, sb, A[1, 0] - A[0, 0], B[1, 0] - B[0, 0], 'N contrast', scope, 'independent cell binding; for all cells (r-1)(c-1) free contrasts remain'), {'aggregate_nullspace_dimension': 9, 'source_row_check_only_cannot_certify_cell_identity': True})
    if i == 7:
        read(i, 'pair_recoverability.py')
        a = ['1', '23']
        b = ['12', '3']
        sa = ''.join(a)
        sb = ''.join(b)
        return (pair('flattened serialization -> paired numbers', a, b, [int(sa)], [int(sb)], 1.0, 12.0, 'diameter fixture units', scope, 'row/column cell delimiters and units'), {'flattened_A': sa, 'flattened_B': sb, 'source_absence_scope': '0 pairs found by this extractor; universal semantic absence not established'})
    if i == 8:
        read(i, 'budget4.py')
        a = {'JE': 2.0, 'CT': 1.5, 't': 2.0, 'h': 2.0}
        b = {'JE': 2.5, 'CT': 1.0, 't': 2.0, 'h': 2.0}
        S = 3.5
        thin = max(0, max(S, 2.0) - 2.0)
        thick = max(0, max(S, 3.35) - 2.0)
        detail = {'same_design_containment_thin_minus_thick': thin - thick, 'published_thin_minus_thick': 1.61 - 0.26, 'containment_order_for_any_shared_h': 'non-positive because max(S,t) and positive-part are nondecreasing in t', 'per_site_height_is_fitted_from_observed_loss': True, 'source_claim_forced_t_thick_3p35': 'not forced when each height is free'}
        return (pair('attachment sum -> component-aware N2 loss', a, b, [S, 2.0, 2.0], [S, 2.0, 2.0], max(0, a['JE'] + max(a['CT'], 2.0) - 2.0), max(0, b['JE'] + max(b['CT'], 2.0) - 2.0), 'mm', scope, 'one component split for N2; fixed-design external sign must be respected'), detail)
    if i == 9:
        g = functions(i, 'suff_test.py', ['F2', 'summ'], {'E': 1.0, 'A': 8430.0, 'EPS': 0.1})
        A = np.array([1.0, 2.0, 3.0, 4.0]) / 16
        B = np.array([4.8, 2.4, 1.6, 1.2]) / 16
        base = B.copy()
        found = False
        for da in range(-16, 17):
            for db in range(-16, 17):
                b = base.copy()
                for _ in range(abs(da)):
                    b[0] = np.nextafter(b[0], math.inf if da > 0 else -math.inf)
                for _ in range(abs(db)):
                    b[1] = np.nextafter(b[1], math.inf if db > 0 else -math.inf)
                if bits(g['summ'](A)) == bits(g['summ'](b)):
                    B = b
                    found = True
                    break
            if found:
                break
        old = jread(i, 'suff_test_out.json')
        detail = {'binary_fiber_found': found, 'source_old_mean_rel_error': old['summary_mean_rel_identity_error'], 'source_old_harmonic_rel_error': old['summary_harmonic_rel_identity_error'], 'domain': 'all positive thicknesses 0.05-1 mm, equal areas; not measured field prevalence'}
        return (pair('arithmetic and harmonic mean -> local break force', A, B, g['summ'](A), g['summ'](B), g['F2'](A), g['F2'](B), 'N', scope, 'minimum thickness in addition to harmonic stiffness; exactly sufficient for equal areas and E'), detail)
    if i == 10:
        g = functions(i, 'pdl_rate_axis.py', ['E_star', 'absE'])
        count = 250.0
        tau = 1.0
        qa = g['absE'](2 * math.pi, 1.0, 1.0, tau)
        qb = g['absE'](4 * math.pi, 1.0, 1.0, tau)
        wl = 2 * math.pi * 250 / 31557600
        wh = 80 * wl
        s = g['absE'](wh, 1.0, 1.0, 1.0) / g['absE'](wl, 1.0, 1.0, 1.0)
        return (pair('annual exceedance count -> within-bite dynamic modulus', {'count': count, 'burst_Hz': 1.0}, {'count': count, 'burst_Hz': 2.0}, count, count, qa, qb, 'MPa', scope, 'within-event waveform/ramp plus recovery time; count alone is not frequency'), {'finite_tau_1s_source_band_ratio': s, 'claimed_finite_wing_ratio': 1.0, 'limits': 'tau->0 yields Einf; tau->infinity yields E0 (source prose reverses these)'})
    if i == 11:
        Nf = -1630 * 5 / math.log(1 - 0.0075)
        g = functions(i, 'flipcount.py', ['pf'], {'N_F': Nf, 'T_HORIZ': 5.0})
        qa = g['pf'](1024.0)
        qb = g['pf'](1024.0)
        return (pair('constant-damage event count -> 5yr failure', {'event_order': [1, 2, 3]}, {'event_order': [3, 2, 1]}, 1024.0, 1024.0, qa, qb, 'probability', scope, 'none within fixed Nf independent-event closure; outside it keep event severities/history'), {'negative_sources': '4/5 declared consumer decisions uncomputable, denominator=1 computable; no patient prediction'})
    if i == 12:
        g = functions(i, 'joint_law.py', ['p_of_eps', 'torque', 'micromotion_stick'], {'K_GEOM': 2 * np.pi * 0.0019 * 9 / 0.0015, 'FUNC': 200.0, 'K_N': 32000000000.0, 'A_EFF': 5e-05})
        a = {'mu': 0.25, 'p0': 64.0, 'eps0': 1.0, 'ct': 1.0}
        b = {'mu': 0.5, 'p0': 32.0, 'eps0': 1.0, 'ct': 1.0}
        ta = [a['eps0'], a['p0'], a['mu'], a['ct']]
        tb = [b['eps0'], b['p0'], b['mu'], b['ct']]
        sa = [float(g['torque'](ta, 1.0)), g['micromotion_stick'](ta)]
        sb = [float(g['torque'](tb, 1.0)), g['micromotion_stick'](tb)]
        return (pair('torque/stick-slip two-port -> normal pressure', a, b, sa, sb, a['p0'], b['p0'], 'pressure fixture units', scope, 'one independent normal-force/pressure measurement at same contact state'), {'family_invariance': 'mu*p0*eps0^(-3/2); joint observations rank<=2 of4; source diagnosis preserved'})
    if i == 13:
        read(i, 'study.py')
        read(i, 'contact_joint_sensitivity.py')
        a = {'scale': 1.0}
        b = {'scale': 2.0}
        nom = jread(i, 'study_out.json') if (pathlib.Path(SAMPLE[i]['source_alias']) / 'study_out.json').exists() else {}
        return (pair('normalized perturbation elasticity -> absolute response', a, b, math.log(2), math.log(2), 1.0, 2.0, 'normalized response', scope, 'absolute response/noise scale plus injectivity/rank; a threshold is not structural identifiability'), {'source_nominal_Y1': 17.421, 'source_perturbed_Y1': 15.879, 'source_linear_sensor_nonzero_move': 15.879 - 17.421, 'structural_claim_refuted_by': 'nonzero response under the declared factor-2 parameter perturbation; detectability depends on actual measurement precision'})
    if i == 14:
        g = functions(i, 'matched_calc.py', ['ratio', 'M'], {'F_FAIL': 663.21})
        kap = 1.0
        u0 = 0.5
        ue = 0.0
        dlinear = 2 * math.log(2)
        dstep = 1.5
        u = 33.2 / 663.21
        wrong = -1 / -4.1471
        right = 4.1471
        return (pair('preload endpoints -> cumulative damage', {'history': 'linear .5->0'}, {'history': 'half at .5 then half at0'}, [u0, ue, 1.0, kap], [u0, ue, 1.0, kap], dlinear, dstep, 'normalized damage', scope, 'integral of (1-u(t))^(-k); endpoints suffice only for frozen prescribed linear history'), {'source_k': wrong, 'correct_k': right, 'source_M_at_declared_load': math.log(g['ratio'](u, 1.0, wrong)), 'correct_M_same_source_ratio_formula': math.log(g['ratio'](u, 1.0, right)), 'source_M_domain_counterexample': math.log(g['ratio'](0.7, 1.0, 2.0)), 'paper_port': 'Table1 preload is external cyclic minimum (10% Fmax), not measured bolt preload'})
    if i == 15:
        g = functions(i, 'cm1_repair.py', ['interface_slip'])
        K = np.array([[2.0, 0.0], [0.0, 1.0]])
        Z = np.array([[1.0], [0.0]])
        D = np.array([[1.0]])
        K1 = K + Z @ D @ Z.T
        B = np.linalg.inv(K)
        W = B - B @ Z @ np.linalg.inv(np.linalg.inv(D) + Z.T @ B @ Z) @ Z.T @ B
        u1 = W @ np.array([1.0, 0.0])
        u2 = W @ np.array([2.0, 0.0])
        err = np.max(abs(u1 - np.linalg.solve(K1, [1.0, 0.0])))
        return (pair('active-set matrix -> solution', {'rhs': [1.0, 0.0]}, {'rhs': [2.0, 0.0]}, K1.flatten(), K1.flatten(), u1[0], u2[0], 'displacement fixture units', scope, 'RHS/load; original solver keeps it'), {'independent_Woodbury_fullsolve_error': float(err), 'benchmark_replayed': False, 'gain_strongest_cached_control': 'UNKNOWN; original full arm splu every iteration'})
    if i == 16:
        read(i, 'thread_band_test.py')
        a = {'eta': 1.0, 'd': 2.0}
        b = {'eta': 0.5, 'd': 4.0}
        return (pair('torque/preload coefficient -> torsional section stress', a, b, a['eta'] * a['d'], b['eta'] * b['d'], 1 / a['d'] ** 3, 1 / b['d'] ** 3, 'normalized torsional stress', scope, 'physical diameter or independent geometry; taper and helix are separate angles'), {'conditional_ratio_mu': 0.5 / 0.12, 'band_ratio': 1030 / 290, 'angle_type_error': 'comparison of lambda*=1.657 thread helix degrees to 3.5-15 implant cone taper degrees'})
    if i == 17:
        g = functions(i, 'band_inversion.py', ['f_junker'])
        R = 290 / 1030
        root = brentq(lambda t: (t + 0.12) / (1 - 0.12 * t) - R * (t + 0.5) / (1 - 0.5 * t), math.tan(math.radians(1)), math.tan(math.radians(5)))
        lam = math.atan(root)
        gg = g['f_junker']
        d = 0.35 / (1030 * 0.5 * gg(0.12, lam)) * 1000
        fp_lo = 0.35 / (0.5 * d / 1000 * gg(0.5, lam))
        fp_hi = 0.35 / (0.5 * d / 1000 * gg(0.12, lam))
        a = {'eta': 1.0, 'd_mm': d}
        b = {'eta': 0.5, 'd_mm': 2 * d}
        return (pair('preload endpoint ratio -> screw size', a, b, [R, root], [R, root], d, 2 * d, 'mm', scope, 'diameter/efficiency separately'), {'exact_thread_root_lambda_deg': math.degrees(lam), 'diameter_mm': float(d), 'Fp_lo_N': float(fp_lo), 'Fp_hi_N': float(fp_hi), 'root_in_declared_1_30_deg': 1 <= math.degrees(lam) <= 30, 'relative_endpoint_residual': abs(float(fp_lo) - 290) / 290, 'contradiction': 'corrected section3 admits a joint in the headline declared domain'})
    if i == 18:
        read(i, 'partA_algebra.py')
        a = {'Lambda': 1.0, 'mu': 0.25, 'kj': 1.0}
        b = {'Lambda': 0.5, 'mu': 0.5, 'kj': 1.0}
        g = 1.0
        N = 1.0
        Fp0 = 1.0
        return (pair('clamp decay coefficient -> accumulated recession', a, b, 1.0, 1.0, a['Lambda'] * (1 - math.exp(-1)), b['Lambda'] * (1 - math.exp(-1)), 'normalized recession', scope, 'independent recession or mu; source already distinguishes the lump', time_kind='HANDOVER'), {'decay_coefficient': '4*kj*mu*Lambda; not a material wear coefficient'})
    if i == 19:
        read(i, 'step1_matched_cell.py')
        l = math.log(0.5)
        qa = 0.5
        qb = 0.5 * math.exp(l ** 3)
        return (pair('F0 and local elasticity -> finite preload update', {'law': 'F_lim=p'}, {'law': 'F_lim=p*exp(log(p)^3)'}, [1.0, 1.0], [1.0, 1.0], qa, qb, 'normalized fatigue limit', scope, 'integrate elasticity along path; sign alone sufficient for non-improvement under monotonic law'), {'both_laws_monotone': 'dlogF/dlogp=1 or 1+3log(p)^2 >=1', 'constant_S_source_factor': 0.5, 'affine_remainder_enclosure': 'not provided by source; exact analytic counterexample, no linearized guarantee used here'})
    if i == 20:
        g = functions(i, 'calc_gamma_converged.py', ['stats'])
        e = np.array([1.0, 2.0, 3.0, 4.0])
        minus = -e
        vm = e
        dev = np.zeros((4, 3, 3))
        tr = np.ones(4)
        va = np.ones(4)
        vb = np.array([1.5, 0.5, 0.5, 1.5])
        a = g['stats'](e, minus, vm, dev, tr, va)
        b = g['stats'](e, minus, vm, dev, tr, vb)
        return (pair('reported Gpct and weighted mean -> squared-strain dose', {'strain': e, 'vol': va}, {'strain': e, 'vol': vb}, [a['eps_mean'], a['G_pct']], [b['eps_mean'], b['G_pct']], float(np.dot(va, e ** 2) / sum(va)), float(np.dot(vb, e ** 2) / sum(vb)), 'squared strain fixture units', scope, 'one weighted second moment for this dose; weighted percentile requires volumes', level='PER_SURFACE_REGION'), {'reported_Gpct_unweighted': a['G_pct'], 'source_quantile_calls': 'np.quantile(e1,.995),np.quantile(e3,.005) ignoring vol', 'domain_limit': 'abstract nonuniform element-volume port; not claimed to be a pair among the 12 fixed FE geometries'})
    if i == 21:
        hs = tuple(map(F, [-1, -0.5, 0, 0.5, 1]))
        DD = {h: d for (h, d) in zip(hs, map(F, ['.50', '.44', '.40', '.34', '.28']))}
        g = functions(i, 'nesting_theorem.py', ['pairs_atoms', 'pairs_band', 'regret_units', 'exact_min_regret'], {'D': DD, 'SIGMA_B': F(506)})
        rad = F(2, 3)
        pa = g['pairs_atoms'](hs, lambda h: h, {-rad, F(0), rad})
        pb = g['pairs_band'](hs, lambda h: h, rad)
        qa = g['exact_min_regret'](pa)[0]
        qb = g['exact_min_regret'](pb)[0]
        return (pair('scatter radius -> robust recourse regret', {'alphabet': ['-2/3', '0', '2/3']}, {'alphabet': 'closed interval [-2/3,2/3]'}, float(rad), float(rad), float(qa), float(qb), 'MPa regret', scope, 'entire error support or consistency-pair port; original source retains it'), {'regret_A_exact': str(qa), 'regret_B_exact': str(qb), 'source_theorem': 'nesting gives monotonically larger consistency sets; proof independent of fixed response monotonicity'})
    if i == 22:
        m = jread(i, 'inputs/SHELL/metrics.json')
        read(i, 'identifiability.py')
        v = np.array(m['K5_witness_m0']['witness_coeffs'])
        H = m['K5_witness_m0']['hidden_stress_per_unit_peak_MPa']
        peak = m['K5_witness_m0']['witness_peak_abs']
        a = 0.001 * v
        b = -0.001 * v
        q = 0.001 * H * peak
        zeros = []
        for (family, val) in m['K8_modulus']['sets'].items():
            if val['R'] == 0:
                zeros.append(family)
        srow = np.array([1.0, 0.0, 0.0, 0.0, 0.0, 0.0])
        sa = float(np.dot(srow, a))
        sb = float(np.dot(srow, b))
        return (pair('mean shell coefficient -> peak signed stress', a, b, sa, sb, q, -q, 'MPa', scope, 'one goal-specific stress functional; six independent channels suffice in fixed six-mode model'), {'source_zero_bound_families': zeros, 'amplitude_interval': m['declared']['amplitude_interval'], 'source_delta_max_normalized': H / m['reference']['ell_ref_MPa_per_unit_norm'], 'physical_amplitude_measured': False, 'zero_bound_logic': 'zero residual uncertainty does not imply zero observed signal'})
    if i == 23:
        g = functions(i, 'probe1.py', ['S_known'])
        u = np.array([1.0, 1.0, 0.0, 0.0])
        va = np.array([0.0, 0.0, 1.0, 1.0])
        vb = 2 * va
        aa = g['S_known'](u, va)
        bb = g['S_known'](u, vb)
        qa = math.sqrt(np.linalg.inv(aa)[1, 1])
        qb = math.sqrt(np.linalg.inv(bb)[1, 1])
        return (pair('raw correlation and N -> parameter precision', {'u': u, 'v': va}, {'u': u, 'v': vb}, [0.0, 4.0], [0.0, 4.0], qa, qb, 'SE divided by observation noise', scope, 'Suu and Svv for known intercept; full centered information with intercept nuisance'), {'original_retains_full_matrix': True, 'zero_porosity_layout_is_measured_feasible': False, 'rigorous_nonlinear_uncertainty_enclosure': 'absent; Fisher estimate is local conditional information, not global bound'})
    if i == 24:
        g = functions(i, 'scripts/run.py', ['t_proj_of'])
        a = [0.875, 1.125, 0.0]
        b = [0.75, 1.25, math.pi / 6]
        sa = [g['t_proj_of'](*a), (a[0] + a[1]) / 2]
        found = False
        for k in range(-32, 33):
            tx = 0.75
            for _ in range(abs(k)):
                tx = np.nextafter(tx, math.inf if k > 0 else -math.inf)
            test = [tx, 2.0 - tx, math.pi / 6]
            sb = [g['t_proj_of'](*test), (test[0] + test[1]) / 2]
            if bits(sa) == bits(sb):
                b = test
                found = True
                break
        sb = [g['t_proj_of'](*b), (b[0] + b[1]) / 2]
        return (pair('optical projected thickness/design mean -> minimum wall', a, b, sa, sb, min(a[:2]), min(b[:2]), 'mm', scope, 'one anisotropy amplitude or measured minimum wall; source optical channels share tproj', level='PER_TOOTH'), {'source_original_trig_fiber_found': found, 'prior_feasibility': 'when Dmax < t_d-t_proj, rotated feasible set is empty; do not hardcode subset/disjoint false', 'external_prior_quantity_mismatch': 'fit/misregistration error is not minimum-wall anisotropy'})
    raise ValueError(i)

def main():
    if '--freeze' in sys.argv:
        recipes = []
        for (i, (name, construction, verdict, note, locator)) in enumerate(RECIPES):
            recipes.append({'index': i, 'job_id': SAMPLE[i]['job_id'], 'probe': name, 'construction': construction, 'pre_probe_textual_assessment': verdict, 'source_locator': locator, 'note': note})
        out = {'frozen_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(), 'code_sha256': sha(pathlib.Path(__file__)), 'prereg_sha256': sha(HERE / 'PREREG_R1.json'), 'sample_sha256': sha(HERE / 'SAMPLE.json'), 'recipes': recipes}
        (HERE / 'PREREG_RECIPES.json').write_text(json.dumps(out, indent=2, ensure_ascii=False) + '\n')
        print('25 recipes frozen')
        return
    frozen = json.loads((HERE / 'PREREG_RECIPES.json').read_text())
    assert sha(pathlib.Path(__file__)) == frozen['code_sha256'], 'code changed after freeze'
    t = time.perf_counter()
    cpu = time.process_time()
    rows = []
    errors = []
    (HERE / 'raw').mkdir(exist_ok=True)
    for (i, r) in enumerate(SAMPLE):
        SOURCES.clear()
        t0 = time.perf_counter()
        (name, _, verdict, note, loc) = RECIPES[i]
        p = pathlib.Path(r['source_alias'])
        assert sha(p / 'RESULTS.md') == r['report_sha256']
        assert sha(p / 'results.json') == r['result_sha256']
        src = json.loads((p / 'results.json').read_text())
        try:
            (pt, detail) = probe(i)
            if not pt['summary_bytes_identical'] or pt['identity_error'] != 0:
                raise AssertionError('exact identity gate failed')
            wrong = np.asarray(pt['summary_B'], float).copy()
            wrong.flat[0] = np.nextafter(wrong.flat[0], math.inf)
            identity_guard = bits(pt['summary_A']) != bits(wrong)
            reference = float(np.asarray(pt['downstream_A'], float).flat[0])
            bad = reference + max(1.0, abs(reference))
            output_guard = reference != bad
            hash_guard = sha(p / 'results.json') != '0' * 64
            controls = {'one_ulp_wrong_summary_rejected': identity_guard, 'injected_wrong_downstream_rejected': output_guard, 'wrong_source_hash_rejected': hash_guard}
            out = {'job_id': r['job_id'], 'claim_type': src.get('claim_type', 'UNKNOWN'), 'review_claim_type': 'capability', 'family': r['family'], 'job_type': r['category'], 'model': r.get('model_route', {}).get('selected', 'UNKNOWN'), 'source': {'alias': r['source_alias'], 'resolved': r['source_path'], 'report_sha256': r['report_sha256'], 'results_sha256': r['result_sha256'], 'claim_locator': loc}, 'review_state': 'PENDING_INDEPENDENT_REVIEW', 'verdict': verdict, 'finding': note, 'sufficiency_probe': pt, 'independent_checks': detail, 'controls': controls, 'producer_external_referent': src.get('external_referent'), 'external_referent': {'kind': 'external_review', 'locator': str(HERE / ('REVIEW_' + r['job_id'] + '.json')), 'compared_quantity': 'original claim vs source equations/data and the named independent arithmetic in this review', 'refutes_us': verdict == 'REJECTED'}, 'probe_fixture_referent': {'kind': 'our_own_fixture', 'locator': str(HERE / 'review.py') + f':probe({i})', 'compared_quantity': pt['name'], 'refutes_us': False}, 'minimality_limit': 'minimal extra scalar only for the named single downstream query/restricted family; not a universal state-compression theorem', 'physical_validation': False, 'source_cost': src.get('cost', src.get('costs', 'UNKNOWN')), 'review_wall_s': time.perf_counter() - t0}
            assert all(controls.values())
        except Exception as e:
            import traceback
            out = {'job_id': r['job_id'], 'error': str(e), 'traceback': traceback.format_exc(), 'verdict': 'UNKNOWN', 'review_state': 'PENDING_INDEPENDENT_REVIEW'}
            errors.append(out)
        out['code_sources'] = list(SOURCES)
        (HERE / ('REVIEW_' + r['job_id'] + '.json')).write_text(json.dumps(native(out), indent=2, ensure_ascii=False, allow_nan=False) + '\n')
        rows.append(out)
        (HERE / 'CURRENT_WORK_STATE.json').write_text(json.dumps({'lane': 'X46-swarm-sufficiency', 'status': 'ACTIVE', 'milestone': f'{len(rows)}/25 reviewed', 'latest_gate': r['job_id'] + ':' + out['verdict'], 'next_operation': 'complete remaining probes then source-bound repairs'}, indent=2) + '\n')
        print(i + 1, r['job_id'], out['verdict'], out.get('sufficiency_probe', {}).get('identity_error'), out.get('sufficiency_probe', {}).get('downstream_difference'), out.get('error', ''))
    payload = {'claim_type': 'capability', 'review_state': 'PENDING_INDEPENDENT_REVIEW', 'sample_n': len(rows), 'counts': dict(collections.Counter((x['verdict'] for x in rows))), 'rows': rows, 'errors': errors, 'wall_s': time.perf_counter() - t, 'cpu_s': time.process_time() - cpu, 'peak_rss_mib': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, 'external_referent': {'kind': 'external_review', 'locator': str(HERE / 'PREREG_RECIPES.json'), 'compared_quantity': '25 frozen source contracts, independently audited; no physical validation', 'refutes_us': bool(errors)}}
    (HERE / 'results_R1.json').write_text(json.dumps(native(payload), ensure_ascii=False, indent=2, allow_nan=False) + '\n')
    if errors:
        sys.exit(1)
if __name__ == '__main__':
    main()
