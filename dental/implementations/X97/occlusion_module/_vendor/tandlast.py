from fractions import Fraction as Q
TOTAL = (Q(205), Q(1963))

def classify(edges, teeth):
    """Exact ALL-gap/ALL-positive-support quantifiers for one arch subset."""
    teeth = frozenset(teeth)
    incident = [i for (i, e) in enumerate(edges) if e[0] in teeth or e[1] in teeth]
    outside = [i for i in range(len(edges)) if i not in incident]
    if not incident:
        (cls, share, force, margin) = ('NEVER', (Q(0), Q(0)), (Q(0), Q(0)), None)
    elif not outside:
        (cls, share, force, margin) = ('MUST', (Q(100), Q(100)), TOTAL, None)
    else:
        a = min((edges[i][3] for i in incident))
        b = min((edges[i][2] for i in outside))
        margin = b - a
        cls = 'MUST' if margin >= 0 else 'CAN'
        (share, force) = ((Q(0), Q(100)), (Q(0), TOTAL[1]))
    return dict(classification=cls, can_bear=bool(incident), must_bear=cls == 'MUST', never_bears=cls == 'NEVER', force_hull_N=force, share_hull_pp=share, infimum_attained=cls != 'MUST' or not outside, necessity_margin_mm=margin, incident_count=len(incident), outside_count=len(outside), license='AXIAL_POSITIVE_SUPPORT_V1', physical_classification='UNCERTAIN', graph_coverage='DECLARED_PAIRS_ONLY', zero_force_physical_interpretation='UNKNOWN' if not incident else None)
