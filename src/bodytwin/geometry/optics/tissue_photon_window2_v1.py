"""Exact prepared photon transport with measured two-length search windows."""
import numpy as np
import warp as wp
from bodytwin.geometry.optics.tissue_photon_segmented_v1 import _window_hit, PACKET
from bodytwin.geometry.optics.tissue_photon_prepared_v1 import PreparedPhotonScene as PreparedBaseline

@wp.func
def _candidate_hit(mesh:wp.uint64,start:wp.vec3d,ray:wp.vec3d,limit:wp.float64):
    if limit<=wp.float64(2.0):
        return _window_hit(mesh,start,ray,wp.float64(0.0),limit)
    lower=wp.float64(0.0)
    for step in range(128):
        upper=wp.min(limit,wp.float64(step+1)*wp.float64(2.0))
        face=_window_hit(mesh,start,ray,lower,upper)
        if face>=0:
            return face
        if upper>=limit:
            return int(-1)
        lower=upper
    return _window_hit(mesh,start,ray,wp.float64(0.0),limit)


@wp.kernel
def _trace(mesh:wp.uint64,front:wp.array(dtype=int),back:wp.array(dtype=int),
           normals:wp.array(dtype=wp.vec3),plane_points:wp.array(dtype=wp.vec3d),plane_normals:wp.array(dtype=wp.vec3d),props:wp.array2d(dtype=float),
           source:wp.vec3,direction:wp.vec3,initial_region:int,seed:int,reflect:int,
           origin:wp.vec3,pitch:float,nx:int,ny:int,nz:int,window_ns:float,
           absorption:wp.array(dtype=wp.int64),terminal:wp.array2d(dtype=wp.int64),
           counters:wp.array(dtype=wp.int64),paths:wp.array2d(dtype=wp.float64),exits:wp.array2d(dtype=wp.float64)):
    tid=wp.tid();rng=wp.rand_init(seed,tid)
    p=wp.vec3d(source);d=direction;region=initial_region
    weight=wp.int64(1073741824);absorbed=wp.int64(0);escaped=wp.int64(0)
    tau=-wp.log(wp.float64(1.0)-(wp.float64(wp.randf(rng))+wp.float64(2.98023223876953125e-8)))
    continuous=wp.float64(1073741824.0)
    clock=wp.float64(0.0);events=int(0);hits=int(0);leak=int(0)
    while weight>wp.int64(0) and events<50000 and leak==0:
        events+=1
        mua=props[region,0];mus=props[region,1];g=props[region,2];index=props[region,3]
        distance=wp.float64(1.0e20)
        if mus>0.0:distance=tau/wp.float64(mus)
        start=wp.vec3d(p);ray=wp.vec3d(d)
        face=_candidate_hit(mesh,start,ray,distance)
        hit=int(0);normal=wp.vec3(0.0);next_region=int(0)
        if face>=0:
            hit=1;hits+=1;normal=normals[face]
            distance=wp.dot(plane_points[face]-start,plane_normals[face])/wp.dot(ray,plane_normals[face])
            if wp.dot(d,normal)>0.0:
                next_region=front[face]
                if back[face]!=region:leak=1
            else:
                next_region=back[face];normal=-normal
                if front[face]!=region:leak=1
        if distance<=wp.float64(0.0) or distance>wp.float64(1.0e19):leak=1
        if leak!=0:break
        traveled=wp.float64(0.0);segments=int(0)
        coordinate=(start-wp.vec3d(origin))/wp.float64(pitch)
        cell=wp.vec3i(int(wp.floor(coordinate[0])),int(wp.floor(coordinate[1])),int(wp.floor(coordinate[2])))
        for axis in range(3):
            if ray[axis]<wp.float64(0.0) and coordinate[axis]==wp.floor(coordinate[axis]):cell[axis]-=1
        while traveled<distance and segments<1024 and weight>wp.int64(0):
            segments+=1
            ix=cell[0];iy=cell[1];iz=cell[2]
            if ix<0 or iy<0 or iz<0 or ix>=nx or iy>=ny or iz>=nz:
                leak=1;break
            crossing=wp.vec3d(wp.float64(1.0e30))
            for axis in range(3):
                if ray[axis]!=wp.float64(0.0):
                    boundary=cell[axis]
                    if ray[axis]>wp.float64(0.0):boundary+=1
                    plane=wp.float64(origin[axis])+wp.float64(boundary)*wp.float64(pitch)
                    crossing[axis]=(plane-start[axis])/ray[axis]
            next_t=wp.min(distance,wp.min(crossing[0],wp.min(crossing[1],crossing[2])))
            length=next_t-traveled
            if length<wp.float64(0.0):leak=1;break
            remaining=continuous*wp.exp(-wp.float64(mua)*length)
            loss=weight-wp.int64(wp.floor(remaining+wp.float64(.5)))
            loss=wp.min(weight,wp.max(wp.int64(0),loss))
            within=wp.clamp((wp.float64(window_ns)-clock)*wp.float64(299.792458)/wp.float64(index),wp.float64(0.0),length)
            score=loss
            if within<length:
                score=weight-wp.int64(wp.floor(continuous*wp.exp(-wp.float64(mua)*within)+wp.float64(.5)))
                score=wp.min(loss,wp.max(wp.int64(0),score))
            wp.atomic_add(absorption,(ix*ny+iy)*nz+iz,score)
            continuous=remaining
            weight-=loss;absorbed+=loss
            paths[tid,region-1]+=length
            clock+=length*wp.float64(index)/wp.float64(299.792458)
            traveled=next_t
            if traveled<distance:
                for axis in range(3):
                    if crossing[axis]==next_t:
                        if ray[axis]>wp.float64(0.0):cell[axis]+=1
                        else:cell[axis]-=1
        if traveled<distance and weight>wp.int64(0) and leak==0:leak=1
        if leak!=0 or weight==wp.int64(0):break
        p=start+ray*distance
        if hit==1:
            tau=wp.max(wp.float64(0.0),tau-wp.float64(mus)*distance)
            transmit=int(1)
            if reflect!=0:
                outgoing=props[next_region,3]
                cosine=wp.clamp(wp.dot(d,normal),0.0,1.0)
                eta=index/outgoing;sine2=eta*eta*(1.0-cosine*cosine)
                if sine2>=1.0:transmit=0
                else:
                    ct=wp.sqrt(wp.max(0.0,1.0-sine2))
                    rs=(index*cosine-outgoing*ct)/(index*cosine+outgoing*ct)
                    rp=(outgoing*cosine-index*ct)/(outgoing*cosine+index*ct)
                    if wp.randf(rng)<.5*(rs*rs+rp*rp):transmit=0
                    else:d=wp.normalize(eta*d+(ct-eta*cosine)*normal)
            if transmit==0:d=wp.normalize(d-2.0*wp.dot(d,normal)*normal)
            else:region=next_region
            if region==0:
                for axis in range(3):
                    exits[tid,axis]=p[axis];exits[tid,axis+3]=wp.float64(d[axis])
                exits[tid,6]=clock
                escaped=weight;weight=wp.int64(0);break
            # Explicit microscopic offset at a handled interface; no implicit miss escape.
            if transmit==0:p-=wp.vec3d(normal)*wp.float64(1.0e-5)
            else:p+=wp.vec3d(normal)*wp.float64(1.0e-5)
        else:
            u=wp.randf(rng);cosine=1.0-2.0*u
            if wp.abs(g)>1.0e-6:
                ratio=(1.0-g*g)/(1.0-g+2.0*g*u)
                cosine=wp.clamp((1.0+g*g-ratio*ratio)/(2.0*g),-1.0,1.0)
            sine=wp.sqrt(wp.max(0.0,1.0-cosine*cosine));phi=6.283185307179586*wp.randf(rng)
            scatter_axis=wp.vec3(0.0,0.0,1.0)
            if wp.abs(d[2])>.9:scatter_axis=wp.vec3(1.0,0.0,0.0)
            t1=wp.normalize(wp.cross(scatter_axis,d));t2=wp.cross(d,t1)
            d=wp.normalize(cosine*d+sine*(wp.cos(phi)*t1+wp.sin(phi)*t2))
            tau=-wp.log(wp.float64(1.0)-(wp.float64(wp.randf(rng))+wp.float64(2.98023223876953125e-8)))
    terminal[tid,0]=absorbed;terminal[tid,1]=escaped;terminal[tid,2]=weight
    terminal[tid,3]=wp.int64(leak)
    wp.atomic_add(counters,0,wp.int64(events));wp.atomic_add(counters,1,wp.int64(hits))
    if events>=50000 and weight>wp.int64(0):wp.atomic_add(counters,2,wp.int64(1))



class PreparedPhotonScene(PreparedBaseline):
    def simulate(self, *, source, direction, initial_region, origin, shape, pitch=1.,
                 photons=1000000, seed=20260912, reflect=False, window_ns=5.):
        if not 0<initial_region<self._regions or photons<1 or photons*PACKET>np.iinfo(np.int64).max:raise ValueError('Invalid launch/energy count')
        shape=tuple(int(x) for x in shape)
        if len(shape)!=3 or min(shape)<=0 or pitch<=0 or window_ns<=0:raise ValueError('Positive grid/window required')
        direction=np.array(direction,dtype=float,copy=True);direction/=np.linalg.norm(direction)
        device=self._device
        absorption=wp.zeros(int(np.prod(shape)),dtype=wp.int64,device=device)
        terminal=wp.zeros((photons,4),dtype=wp.int64,device=device);counters=wp.zeros(3,dtype=wp.int64,device=device)
        paths=wp.zeros((photons,self._regions-1),dtype=wp.float64,device=device)
        exits=wp.zeros((photons,7),dtype=wp.float64,device=device)
        wp.launch(_trace,photons,inputs=[self._mesh.id,self._front,self._back,
            self._normals,self._plane_points,self._plane_normals,self._properties,wp.vec3(*source),wp.vec3(*direction),
            initial_region,seed,int(reflect),wp.vec3(*origin),float(pitch),*shape,float(window_ns),absorption,terminal,counters,paths,exits],device=device,block_dim=64)
        return {'absorption':absorption.numpy().reshape(shape),'terminal':terminal.numpy(),'counters':counters.numpy(),'paths':paths.numpy(),'exits':exits.numpy()}


if __name__=='__main__':
    import subprocess
    import sys
    from pathlib import Path
    root=Path(__file__).resolve().parents[4]
    raise SystemExit(subprocess.call([sys.executable,str(root/'probes/optics/photon_segmented_transport_probe.py'),'--window2',*sys.argv[1:]]))
