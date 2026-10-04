"""Runnable independent-method port example: visible contralateral reflection.
Usage: python3 example_method.py MASKED_CONTEXT.npz PREDICTED_SURFACE.npz
"""
import sys
from crownbench import *
with np.load(sys.argv[1]) as d:
    v = d['vertices']
    f = d['faces']
    labels = d['labels']
    k = int(d['target_fdi'])
assert k not in np.unique(labels)
t = tooth_dict(v, f, labels)
q = 1 if k // 10 < 3 else 3
(n, c, _) = reflection(t, q)
(tv, tf) = t[swap(k)]
write_surface(sys.argv[2], reflect(tv, n, c), tf[:, ::-1])
