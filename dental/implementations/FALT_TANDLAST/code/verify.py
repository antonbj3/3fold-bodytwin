"""Independent quantity/certificate replay and adversarial model probes."""
from copy import deepcopy
from itertools import product
from collections import Counter
import json, time, resource
from model import *
from motion_engine.ncp.gap_box import GapProblem, enclose_gap_box, verify_gap_enclosure
from motion_engine.ncp.latent_contact import make_family, certify_family, verify_family

def independent_case(row):
    """Worst-corner global minima test; does not call classify."""
    (e, p) = edges_from_x21(row['case'])
    if row['source_gap_sha256'] != sha(p):
        return False
    entries = [(int(t), b, {int(t)}) for (t, b) in row['teeth'].items()]
    for (jaw, regions) in row['regions'].items():
        entries += [(r, b, {t for t in FDI if (t < 30) == (jaw == 'upper') and region(t) == r}) for (r, b) in regions.items()]
    for (key, b, S) in entries:
        inc = [i for (i, x) in enumerate(e) if x[0] in S or x[1] in S]
        corner = [x[3] if i in inc else x[2] for (i, x) in enumerate(e)]
        minimum = min(corner)
        necessary = any((i in inc for (i, x) in enumerate(corner) if x == minimum))
        cls = 'NEVER' if not inc else 'MUST' if necessary else 'CAN'
        (lo, hi) = map(Q, b['share_hull_pp'])
        (f0, f1) = map(Q, b['force_hull_N'])
        outside = len(e) - len(inc)
        expected = (Q(0), Q(0)) if not inc else (Q(100), Q(100)) if not outside else (Q(0), Q(100))
        expectedN = (Q(0), Q(0)) if not inc else TOTAL if not outside else (Q(0), TOTAL[1])
        margin = None if not inc or not outside else min((x[2] for (i, x) in enumerate(e) if i not in inc)) - min((e[i][3] for i in inc))
        if b['classification'] != cls or (lo, hi) != expected or (f0, f1) != expectedN or (b['license'] != 'AXIAL_POSITIVE_SUPPORT_V1') or (b['physical_classification'] != 'UNCERTAIN') or (b['can_bear'] != bool(inc)) or (b['must_bear'] != (cls == 'MUST')) or (b['never_bears'] != (cls == 'NEVER')) or (b['incident_count'] != len(inc)) or (b['outside_count'] != outside) or ((None if b['necessity_margin_mm'] is None else Q(b['necessity_margin_mm'])) != margin) or (b['infimum_attained'] != (cls != 'MUST' or not outside)):
            return False
    return True

def main():
    stopped()
    start = time.perf_counter()
    checks = Counter()
    points = []
    rows = json.load(open(ROOT / 'raw/CASE_CERTIFICATES.json'))
    assert len(rows) == 993
    for r in rows:
        assert independent_case(r)
        checks['full_case_replays'] += 1
    manifest = json.load(open(ROOT / 'raw/SOURCE_MANIFEST.json'))['sha256']
    for (path, value) in manifest.items():
        assert sha(path) == value, path
        checks['source_hash_replays'] += 1
    panel = json.load(open(X21 / 'raw/INPUT_MANIFEST.json'))['numerical_panel']
    for case in panel:
        stopped()
        (e, p) = edges_from_x21(case)
        gap = [(x[2] + x[3]) / 2 for x in e]
        for T in (Q(205), Q(879), Q(1963)):
            ans = exact_point(e, gap, T, {t: Q(1, 750) for t in FDI})
            assert replay_point(e, ans)
            for t in FDI:
                cls = classify(e, [t])['classification']
                value = ans['tooth_force_N'].get(t, Q(0))
                assert value > 0 if cls == 'MUST' else value >= 0
                if cls == 'NEVER':
                    assert value == 0
            points.append(dict(case=case, point=ans))
            checks['real_exact_QP_queries'] += 1
    e = [(11, 41, Q(0), Q(0)), (12, 42, Q(1, 5), Q(1, 5))]
    approach = []
    for M in (Q(2), Q(20), Q(200), Q(2000)):
        first = Q(6, 5) / (M + 1)
        second = 1 - first
        d = M * first
        pb = GapProblem.make([[1, 1, 0, 0], [0, 0, 1, 1]], [M / 2000, M / 2000, Q(1, 2000), Q(1, 2000)], [-d / 2000] * 4, [[0, 0], [Q(1, 5000)] * 2], 1)
        answer = enclose_gap_box(pb, max_nodes=32)
        assert verify_gap_enclosure(pb, answer)
        assert answer['bounds']['normal_force:0']['lo'] == first
        assert answer['bounds']['normal_force:1']['lo'] == second
        assert first > 0 and classify(e, [11])['classification'] == 'MUST'
        approach.append(dict(M=M, force1_N=first, force2_N=second, archive_replay=True))
        checks['strict_positive_zero_infimum_probes'] += 1
    tie = [(11, 41, Q(0), Q(0)), (12, 42, Q(0), Q(0))]
    ans = exact_point(tie, [Q(0), Q(0)], Q(205), {t: Q(1, 750) for t in (11, 12, 41, 42)})
    assert replay_point(tie, ans) and all((f > 0 for f in ans['tooth_force_N'].values()))
    assert all((classify(tie, [t])['classification'] == 'MUST' for t in (11, 12, 41, 42)))
    assert classify(tie, [18])['classification'] == 'NEVER'
    checks['tie_and_isolated_probes'] += 2
    fixtures = [[(11, 41, Q(0), Q(1)), (12, 42, Q(1), Q(2))], [(11, 41, Q(-1), Q(2)), (11, 42, Q(0), Q(1)), (12, 41, Q(-2), Q(0)), (12, 42, Q(1), Q(3))]]
    for edges in fixtures:
        for t in (11, 12, 41, 42):
            all_min = []
            for corner in product(*[(e[2], e[3]) for e in edges]):
                m = min(corner)
                all_min.append(any((t in edges[i][:2] for (i, v) in enumerate(corner) if v == m)))
            assert classify(edges, [t])['must_bear'] == all(all_min)
            checks['licensed_corner_exhaustions'] += 1
    pb = make_family([[2, 1], [1, 2]], [-5, -6], [7, -7], [-1, 1], step='1/100', source_id='published-Siconos-matrix:declared-extension', source_sha256='0' * 64, unit_system='EXTERNAL_UNSPECIFIED')
    family = certify_family(pb)
    assert verify_family(pb, family)
    assert family['bounds']['total_normal_force']['lo'] == Q(1100, 3)
    assert family['bounds']['total_normal_force']['lo'] < 600
    bad = deepcopy(family)
    bad['bounds']['total_normal_force']['lo'] = Q(600)
    assert not verify_family(pb, bad)
    checks['unlicensed_force_corner_rejections'] += 1
    base = next((r for r in rows if any((b['classification'] == 'MUST' for b in r['teeth'].values()))))
    for kind in ('class', 'bound', 'source', 'license', 'physical', 'margin', 'attainment'):
        bad = deepcopy(base)
        t = next((t for (t, b) in bad['teeth'].items() if b['classification'] == 'MUST'))
        b = bad['teeth'][t]
        if kind == 'class':
            b['classification'] = 'CAN'
        elif kind == 'bound':
            b['share_hull_pp'][1] = '99'
        elif kind == 'source':
            bad['source_gap_sha256'] = '0' * 64
        elif kind == 'license':
            b['license'] = '3D_COULOMB_CORNER'
        elif kind == 'physical':
            b['physical_classification'] = 'MUST'
        elif kind == 'margin':
            b['necessity_margin_mm'] = str(Q(b['necessity_margin_mm']) + 1)
        elif kind == 'attainment':
            b['infimum_attained'] = not b['infimum_attained']
        assert not independent_case(bad)
        checks['case_mutations_rejected'] += 1
    (e, _) = edges_from_x21(points[0]['case'])
    a = deepcopy(points[0]['point'])
    a['edge_force_N'][0] += 1
    assert not replay_point(e, a)
    checks['QP_mutations_rejected'] += 1
    write(ROOT / 'raw/POINT_CONTROLS.json', points)
    write(ROOT / 'raw/STRICT_POSITIVITY.json', dict(points=approach, force_infimum=0, minimum_attained=False, explanation='All finite M yield strictly positive load; M->infinity yields zero infimum.', archive_unit_transport='mm -> m, algebraic quasi-static compliance represented by GapProblem with h=1; not physical tooth inertia'))
    write(ROOT / 'raw/VERIFICATION.json', dict(status='PASS_SCOPED', checks=dict(checks), wall_s=time.perf_counter() - start, peak_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, independent_review='PENDING; producer replay only', proof='QUANTITY_LICENSE.md', source_drift=False, physical_validation='NOT_RUN'))
    print(dict(status='PASS_SCOPED', checks=dict(checks), wall_s=time.perf_counter() - start))
if __name__ == '__main__':
    main()
