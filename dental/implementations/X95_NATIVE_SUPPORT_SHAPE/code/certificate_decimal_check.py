from common import *
from native_port import model
from decimal import Decimal, localcontext

def run():
    (a, g, m) = model()
    r = read(R / 'raw/R4_PORT_CERTIFICATE.json')
    point = [Decimal.from_float(x) for x in r['point_mm']]
    tri = m.Vd[m.Fd_master[:, 2:5]].astype(float)
    lower = []
    upper = []
    with localcontext() as ctx:
        ctx.prec = 100
        for t in tri:
            axes = [[Decimal.from_float(float(v)) for v in t[:, k]] for k in range(3)]
            delta = [max(min(axes[k]) - point[k], point[k] - max(axes[k]), Decimal(0)) for k in range(3)]
            lower.append(sum((x * x for x in delta)).sqrt())
            for v in t:
                dd = [Decimal.from_float(float(v[k])) - point[k] for k in range(3)]
                upper.append(sum((x * x for x in dd)).sqrt())
        exact_lo = min(lower)
        exact_hi = min(upper)
        lo_down = exact_lo.next_minus()
        hi_up = exact_hi.next_plus()
        bound = r['rigorous_distance_enclosure_mm']
        ok = Decimal.from_float(bound[0]) <= lo_down and Decimal.from_float(bound[1]) >= hi_up
    out = {'reference': '100-digit scalar Decimal calculation on exact float32-derived binary coordinates, sqrt bracketing via next_minus/next_plus', 'box_lower_reference': str(exact_lo), 'vertex_upper_reference': str(exact_hi), 'binary64_bounds': bound, 'pass': ok, 'scope': 'Independent scalar arithmetic check of the fixed geometry enclosure; not physical validation'}
    dump(R / 'raw/R4_DECIMAL_CHECK.json', out)
    if not ok:
        raise ValueError('Enclosure misses independent scalar bounds')
    print(json.dumps(out, indent=2))
if __name__ == '__main__':
    run()
