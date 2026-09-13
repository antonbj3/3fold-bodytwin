"""Measure the external cube60 reference before a new photon solver design."""
import hashlib
import importlib.metadata
import json
from pathlib import Path
import numpy as np
import pmcx

ROOT=Path(__file__).resolve().parents[2]


def main():
    rows=[];fields=[];configs=[]
    for leg in range(2):
        cfg=pmcx.mcxcreate('cube60')
        result=pmcx.run(cfg)
        flux=np.asarray(result['flux']);fields.append(flux.copy())
        stats=result.get('stat',{})
        stable_stats={k:float(stats[k]) for k in ('energytot','energyabs','normalizer') if k in stats}
        rows.append({'shape':list(flux.shape),'dtype':str(flux.dtype),
                     'finite':bool(np.isfinite(flux).all()),'minimum':float(flux.min()),'maximum':float(flux.max()),
                     'fluence_integrated_absorption':float(np.sum(flux,dtype=np.float64)*cfg['tstep']*cfg['prop'][1][0]),
                     'statistics':stable_stats,'flux_sha256':hashlib.sha256(flux.tobytes()).hexdigest()})
        configs.append({k:v for k,v in cfg.items() if k!='vol'})
        print(json.dumps(rows[-1]),flush=True)
    gates={'finite_nonnegative':all(r['finite'] and r['minimum']>=0 for r in rows),
           'configuration_identical':configs[0]==configs[1],
           'two_fields_bit_identical':fields[0].tobytes()==fields[1].tobytes()}
    report={'legs':rows,'gates':gates,'config':configs[0],
            'maximum_repeat_difference':float(np.max(np.abs(fields[0]-fields[1]))),
            'pmcx_version':importlib.metadata.version('pmcx'),
            'scope':'external cube60 reference measurement; no new solver or throughput claim'}
    (ROOT/'reports').mkdir(exist_ok=True)
    (ROOT/'reports/mcx_reference.json').write_text(json.dumps(report,indent=2)+'\n')
    np.savez_compressed(ROOT/'reports/mcx_reference_fields.npz',first=fields[0],second=fields[1])
    print(json.dumps(report),flush=True)
    return 0 if all(gates.values()) else 1


if __name__=='__main__':raise SystemExit(main())
