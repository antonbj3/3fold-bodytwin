"""Retain tiny triangles as conservative longest-edge capsule covers.

Every original triangle is retained. Tiny areas excluded only from random
normal sampling, and their collision geometry is covered, never deleted.
"""
import numpy as np
from codesign import Scene, EPS
from collision import segment_triangle, segment_segment

class SliverScene(Scene):

    def __init__(self, vertices, faces):
        self.v = np.asarray(vertices, float)
        self.f = np.asarray(faces, int)
        if self.v.ndim != 2 or self.v.shape[1] != 3 or (not np.isfinite(self.v).all()):
            raise ValueError('bad vertices')
        if self.f.ndim != 2 or self.f.shape[1] != 3 or (not len(self.f)) or (self.f.min() < 0) or (self.f.max() >= len(self.v)):
            raise ValueError('bad faces')
        self.tri = self.v[self.f]
        self.lo = self.tri.min(1)
        self.hi = self.tri.max(1)
        self.area = np.linalg.norm(np.cross(self.tri[:, 1] - self.tri[:, 0], self.tri[:, 2] - self.tri[:, 0]), axis=1) / 2
        self.bounds = np.array([self.v.min(0), self.v.max(0)])
        self.sliver = self.area < 1e-14
        self.valid = ~self.sliver
        t = self.tri[self.sliver]
        if len(t):
            lengths = np.array([np.sum((t[:, (k + 1) % 3] - t[:, k]) ** 2, axis=1) for k in range(3)]).T
            k = lengths.argmax(1)
            ix = np.arange(len(t))
            self.a = t[ix, k]
            self.b = t[ix, (k + 1) % 3]
            c = t[ix, (k + 2) % 3]
            u = self.b - self.a
            uu = np.sum(u * u, axis=1)
            alpha = np.clip(np.divide(np.sum((c - self.a) * u, axis=1), uu, out=np.zeros(len(t)), where=uu > 0), 0, 1)
            self.width = np.linalg.norm(c - self.a - alpha[:, None] * u, axis=1)
        else:
            self.a = self.b = np.empty((0, 3))
            self.width = np.empty(0)

    def clearance(self, a, b, r):
        a = np.asarray(a)
        b = np.asarray(b)
        lo = np.minimum(a, b) - r - EPS
        hi = np.maximum(a, b) + r + EPS
        keep = self.valid & np.all(self.hi >= lo, axis=1) & np.all(self.lo <= hi, axis=1)
        val = float(segment_triangle(a, b, self.tri[keep]).min() - r) if keep.any() else np.inf
        if len(self.a):
            val = min(val, float((segment_segment(a, b, self.a, self.b) - self.width - r).min()))
        return val
