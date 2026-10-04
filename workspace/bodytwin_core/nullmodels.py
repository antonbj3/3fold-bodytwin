"""B24 reference plus N1/N1g nulls; N1 formulas/parameters remain source scoped."""
from .config import paths
import json

def b24_lifting(bodyweight_N, config=None):
    p=paths(config)
    rows=json.loads((p['b24_results']/'results.json').read_text())['facets']
    return {r['facet_id']:float(bodyweight_N*r['population_null']['coefficient_pctBW']/100)
            for r in rows if 'Lifting' in r['facet_id']}

def n1_knee(grf_peak_N, coefficient=2.2):
    if coefficient not in (2.2,2.4):raise ValueError('N1 coefficients recorded for N12 are 2.2 or 2.4')
    return coefficient*grf_peak_N

def n1g_knee(grf_peak_N, activity_coefficient):
    """N1g = activity-group coefficient × |GRF|; coefficient must be externally trained."""
    return float(activity_coefficient)*grf_peak_N
