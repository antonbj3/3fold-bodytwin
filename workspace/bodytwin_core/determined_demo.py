"""Render the frozen four-person knee-load comparison."""
import json
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parents[1]/'results/CX-DETERMINED'

def run():
    data=json.loads((HERE/'results.json').read_text())
    new=json.loads((HERE/'new_activity_results.json').read_text())
    print('trial                    person activity       frames law band set   law/N1g RMSE BW   LP status')
    for r in data['trials']:
        frac=f"{r['set_fraction_all']:.1%}" if r['set_fraction_all'] is not None else '—'
        print(f"{r['key'][:24]:24} {r['person']:6} {r['activity'][:13]:13} {r['n_frames']:6} {r['law_fraction']:7.1%} {frac:>5} {r['law_rmse_BW']:.2f}/{r['n1g_rmse_BW']:.2f} {r['set_status']}")
    for r in new['rows']:
        print(f"{('JW4__'+r['trial'])[:24]:24} JW     {r['activity'][:13]:13} {r['frames']:6} {r['law_fraction']:7.1%}     — {r['law_rmse_BW']:.2f}/{r['n1g_rmse_BW']:.2f} {r['set_status']}")
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    x=json.loads((HERE/'figure_data.json').read_text())
    t=np.asarray(x['time_s']);band=np.asarray(x['law_band_N']);setband=np.asarray(x['set_N'],float)
    fig,ax=plt.subplots(figsize=(10,5))
    ax.fill_between(t,band[:,0],band[:,1],color='#258a9f',alpha=.30,label='Determined law: lever-arm scenarios')
    ax.fill_between(t,setband[:,0],setband[:,1],color='#b17c2a',alpha=.23,label='Undetermined feasible set')
    ax.plot(t,x['law_N'],color='#08677b',linewidth=1.5,label='F_law nominal')
    ax.plot(t,x['measured_N'],color='black',linewidth=1.4,label='Implant measured')
    ax.plot(t,x['n1g_N'],color='#962f75',linestyle='--',linewidth=1.2,label='N1g reference')
    ax.plot(t,x['selected_N'],color='#6b5120',linestyle=':',linewidth=1.2,label='L1 S1b selected point')
    ax.set(xlabel='Time (s)',ylabel='Total knee compression (N)',title=x['key'])
    ax.legend(loc='upper left',fontsize=8)
    ax.set_ylim(bottom=0)
    fig.tight_layout();path=HERE/'determined.png';fig.savefig(path,dpi=170);plt.close(fig)
    print(f'Figure: {path}')
