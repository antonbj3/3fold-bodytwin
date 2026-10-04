from common import *

def collapse(rows):
    groups = {}
    for r in rows:
        groups.setdefault((r['image_group'], r['fdi']), []).append(r)
    out = []
    for ((g, f), rr) in sorted(groups.items()):
        accepted = [r for r in rr if r['nominal_tf2_lower_mm'] >= 2]
        out.append(dict(image_group=g, fdi=f, accepted=bool(accepted), below=any((r['union_upper_mm'] < 2 for r in accepted)), tf2_below=any((r['tf2_upper_mm'] < 2 for r in accepted)), uncertain=any((r['union_lower_mm'] < 2 <= r['union_upper_mm'] for r in accepted)), paired=any((r['paired_revisions'] for r in rr)), rows=rr))
    return out

def proportion_ci(rows, value='below', accepted=True, seed=9601):
    r = [x for x in rows if x['accepted']] if accepted else rows
    if not r:
        return dict(n=0, count=0, fraction=None, ci95=None)
    groups = {}
    for x in r:
        groups.setdefault(x['image_group'], []).append(float(x[value]))
    nums = np.array([sum(v) for v in groups.values()])
    dens = np.array([len(v) for v in groups.values()])
    rng = np.random.default_rng(seed)
    ix = rng.integers(0, len(nums), size=(2000, len(nums)))
    rates = nums[ix].sum(1) / dens[ix].sum(1)
    patient = nums / dens
    return dict(n=len(r), count=int(nums.sum()), fraction=float(nums.sum() / dens.sum()), ci95=np.quantile(rates, [0.025, 0.975]).tolist(), image_groups=len(groups), image_group_mean_site_fraction=float(patient.mean()), image_group_mean_ci95=np.quantile(patient[ix].mean(1), [0.025, 0.975]).tolist(), image_groups_any_breach=int((nums > 0).sum()), resolution='POPULATION', ci_scope='image-group empirical bootstrap; convenience sample, fixed digital model; not guide population uncertainty')

def summarize(depth, rows):
    cc = collapse(rows)
    return dict(extra_mm=depth, pose_queries=len(rows), unique_sites=len(cc), all=proportion_ci(cc), paired_revisions=proportion_ci([r for r in cc if r['paired']]), unpaired=proportion_ci([r for r in cc if not r['paired']]), tf2_only=proportion_ci(cc, 'tf2_below'), unknown_numeric=sum((r['uncertain'] for r in cc)), nominal_below2=sum((not r['accepted'] for r in cc)), alias_pose_duplicates=len(rows) - len(cc))
