from common import *
import collections

def run():
    rows = read(ROOT / 'raw/FUNCTIONAL_ROWS.json')
    groups = collections.defaultdict(list)
    for r in rows:
        if r['status'] == 'SCORED':
            groups[r['family'], r['participant']].append(r)
    table = []
    for ((family, participant), rs) in sorted(groups.items()):
        counts = [r['function']['0']['contact_regions'] for r in rs if r['function'].get('0', {}).get('contact_regions') is not None] if rs and '0' in rs[0]['function'] else [r['function']['0.0']['contact_regions'] for r in rs if r['function']['0.0']['contact_regions'] is not None]
        table.append(dict(family=family, participant=participant, requested=len(rs), supported=len(counts), zero_contact_fraction=float(np.mean(np.asarray(counts) == 0)) if counts else None, connected_region_count_mean=float(np.mean(counts)) if counts else None, connected_region_count_range=[min(counts), max(counts)] if counts else None, resolution='POPULATION', input_resolution='PER_SURFACE_REGION', individual_result_resolution='PER_TOOTH', population_norm_admission='NOT_IDENTIFIED_METHOD_AND_TOOTH_MATCH'))
    refs = [dict(locator='https://pubmed.ncbi.nlm.nih.gov/37301414/', quantity='prevalence of no static habitual occlusal contact by natural posterior tooth', method='silicone registration, GEDAS II', resolution='POPULATION', reported_cells=[dict(tooth_type='maxillary_first_molar', FDI=[16, 26], zero_fraction=0.097), dict(tooth_type='mandibular_second_molar', FDI=[37, 47], zero_fraction=0.2)], scope='Observed distribution, not an ideal contact requirement; full protocol mismatch prevents numerical admission'), dict(locator='https://pmc.ncbi.nlm.nih.gov/articles/PMC8156897/', quantity='digital total contact area versus clinical articulating-foil contact count', resolution='POPULATION', method='digital casts versus8um foil', finding='Digital area did not significantly correlate with clinical contact count in this sample; numerical values not transferred to generated teeth'), dict(locator='https://doi.org/10.1046/j.1365-2842.1997.00596.x', quantity='occlusal number/distribution in normal dentitions', resolution='POPULATION', finding='Wide variation; classical theoretical fixed-number proposals not supported'), dict(locator='https://pmc.ncbi.nlm.nih.gov/articles/PMC13222086/', quantity='contacts per tooth in posterior teeth pooled over types', method='8um paper / IOS', resolution='POPULATION', reported_mean=2.89, reported_SD=1.32, range=[0, 8], rejected_reason='pooled posterior count is not a tooth-type norm; paper-versus-distance regime mismatch; cannot assign 2.89 to each crown')]
    out = dict(comparison_rows=table, published_referents=refs, dropout={'norm_candidate_sources': 4, 'sources_retained_as_context': 4, 'sources_admitted_as_individual_pass_threshold': 0, 'individual_threshold_rejection_fraction': 1.0, 'reasons': ['method/pose mismatch', 'population versus individual target', 'pooled tooth types', 'no unique normative count']}, conclusion='Connected geometric contact regions per tooth type are now quantified. Universal normative per-tooth pass counts are not supported by these primary sources; matched clinical contact acquisition needed.')
    dump(ROOT / 'raw/CONTACT_NORMS.json', out)
    return out
if __name__ == '__main__':
    run()
