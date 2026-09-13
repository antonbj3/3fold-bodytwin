"""Instrument a frozen photon kernel without editing its transport baseline."""
import hashlib
from pathlib import Path
PIN = 'f78dfc8eb15289e6dd0972ed228e4bab7fa032e394833ab234e28f5f9bb0a271'


def instrument(source, max_events=256):
    if type(max_events) is not int or max_events not in (256,50000):
        raise ValueError("Explicit prefix or complete event capacity required")
    if hashlib.sha256(source.encode()).hexdigest() != PIN:
        raise ValueError('Frozen photon source changed')
    edits = {
        'exits:wp.array2d(dtype=wp.float64)):' : 'exits:wp.array2d(dtype=wp.float64),query:wp.array3d(dtype=wp.float64),selected:wp.array2d(dtype=int)):',
        'face=_candidate_hit(mesh,start,ray,distance)': '''face=_candidate_hit(mesh,start,ray,distance)
        if events<=256:
            for axis in range(3):
                query[tid,events-1,axis]=start[axis]
                query[tid,events-1,axis+3]=ray[axis]
            query[tid,events-1,6]=distance
            selected[tid,events-1]=face''',
        '    wp.launch(_trace,photons,inputs=': '''    query=wp.zeros((photons,256,7),dtype=wp.float64,device=device)
    selected=wp.full((photons,256),-2,dtype=int,device=device)
    wp.launch(_trace,photons,inputs=''',
        'absorption,terminal,counters,paths,exits],device=device)': 'absorption,terminal,counters,paths,exits,query,selected],device=device)',
        "'exits':exits.numpy()}": "'exits':exits.numpy(),'query':query.numpy(),'selected':selected.numpy()}",
    }
    for old,new in edits.items():
        if source.count(old) != 1:
            raise ValueError('Instrumentation anchor is not unique')
        source=source.replace(old,new)
    if max_events == 50000:
        source=source.replace("events<=256", "events<=50000").replace("(photons,256,7)", "(photons,50000,7)").replace("(photons,256)", "(photons,50000)")
    return source
