"""Frozen third-band direct transport audit with explicit seed and backend groups."""
import numpy as np
from atlas_oxygen_bounds import FRACTION, BACKGROUND
from atlas_spectral_transport_probe import assess, PACKET


def absorption(saturation):
    blood=np.log(10.)*15/64500.*(816*saturation+761.72*(1-saturation))
    return np.asarray([BACKGROUND+blood*FRACTION],dtype=np.float32)


def evaluate(rows):
    groups={(backend,seed):[r for r in rows if r['backend']==backend and r['seed']==seed]
            for backend in ('H100','L4') for seed in (20260912,20260913)}
    if len(rows)!=8 or any(len(g)!=2 for g in groups.values()):
        raise ValueError('Exactly two complete captures per declared backend and seed required')
    if any(r['wavelength']!=800 or r['saturation']!=.8 or r['photons']!=10000 for r in rows):
        raise ValueError('Frozen arm parameters changed')
    gates=dict(full_repeat=all(g[0]==g[1] for g in groups.values()),
               energy=all(r['absorbed']+r['escaped']+r['residual']==r['photons']*PACKET for r in rows),
               zero_failures=all(r['leaks']==r['caps']==r['residual']==0 for r in rows),
               direct_path_agreement=all(r['intensity_difference']<=r['bound'] for r in rows),
               cross_backend_exact=all(groups['H100',seed][0]['hashes']==groups['L4',seed][0]['hashes'] for seed in (20260912,20260913)))
    return dict(rows=rows,gates=gates,scope='Direct third-band transport under fixed assumed optical properties; no clinical, anatomy, scattering-spectrum or universal hardware acceptance.')


if __name__=='__main__':
    import argparse,json
    from pathlib import Path
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('capture',type=Path)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();report=evaluate(json.loads(args.capture.read_text())['rows'])
    args.output.write_text(json.dumps(report,indent=2)+'\n')
    raise SystemExit(0 if all(report['gates'].values()) else 2)
