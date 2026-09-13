"""Shared-origin batches beside the frozen scalar tetrahedral ray walker.

Complete cell/interval records in input order; no implicit file writes.
"""
import numpy as np
from bodytwin.geometry.tetra_ray_walk_v1 import TetraRayWalk


class TetraRayBatch(TetraRayWalk):
    def trace_batch(self, origin, directions):
        origin=np.asarray(origin,dtype=np.float64)
        directions=np.asarray(directions,dtype=np.float64)
        if origin.shape!=(3,) or directions.ndim!=2 or directions.shape[1]!=3:
            raise ValueError('Three-component origin and direction rows required')
        if not np.isfinite(origin).all() or not np.isfinite(directions).all():
            raise ValueError('Finite three-component ray required')
        if not len(directions):return []
        local=np.einsum('nij,nj->ni',self.inverse,origin-self.origins)
        bary=np.column_stack((1-local.sum(axis=1),local));inside=np.flatnonzero(np.min(bary,axis=1)>0)
        if len(inside)!=1:raise ValueError('Source must be strictly inside exactly one cell')
        return [self._walk(bary,int(inside[0]),direction) for direction in directions]

    def _walk(self, bary, initial_cell, direction):
        length=np.linalg.norm(direction)
        if not np.isfinite(length) or length==0:raise ValueError('Nonzero finite ray direction required')
        direction=direction/length
        cell=initial_cell;entry=-1;parameter=0.;records=[];visited=set()
        for _ in range(len(self.cells)+1):
            if cell in visited:raise ValueError('Traversal revisited a convex cell')
            visited.add(cell)
            q=self.inverse[cell]@direction;delta=np.r_[-q.sum(),q]
            candidates=[]
            for face in range(4):
                if face!=entry and delta[face]<0:
                    value=-bary[cell,face]/delta[face]
                    if value<parameter:raise ValueError('Negative cell interval')
                    candidates.append((float(value),face))
            if not candidates:raise ValueError('No forward exit from bounded cell')
            candidates.sort();exit_parameter,face=candidates[0]
            if len(candidates)>1 and candidates[1][0]==exit_parameter:raise ValueError('Ambiguous edge or vertex exit')
            if exit_parameter<=parameter:raise ValueError('Zero-length cell interval')
            records.append((cell,parameter,exit_parameter))
            neighbor=int(self.neighbors[cell,face])
            if neighbor<0:return np.asarray(records,dtype=np.float64)
            entry=int(self.entry_faces[cell,face]);cell=neighbor;parameter=exit_parameter
        raise ValueError('Cell traversal bound exceeded')

if __name__=='__main__':
    import runpy
    from pathlib import Path
    runpy.run_path(str(Path(__file__).resolve().parents[3]/'probes/optics/tetra_ray_batch_probe.py'),run_name='__main__')
