"""Observe frozen v3 boundary failures without changing photon paths.
Separate photon variant with strictly interior optical-depth RNG samples.

This experimental module requires the documented physical gates before use.
Source position and initial region are explicit; launch energy is inside tissue.
"""
import numpy as np
import warp as wp

PACKET=1 << 30


@wp.kernel
def _trace(mesh:wp.uint64,front:wp.array(dtype=int),back:wp.array(dtype=int),
           normals:wp.array(dtype=wp.vec3),plane_points:wp.array(dtype=wp.vec3d),plane_normals:wp.array(dtype=wp.vec3d),props:wp.array2d(dtype=float),
           source:wp.vec3,direction:wp.vec3,initial_region:int,seed:int,reflect:int,
           origin:wp.vec3,pitch:float,nx:int,ny:int,nz:int,window_ns:float,
           absorption:wp.array(dtype=wp.int64),terminal:wp.array2d(dtype=wp.int64),
           counters:wp.array(dtype=wp.int64)):
    tid=wp.tid();rng=wp.rand_init(seed,tid)
    p=source;d=direction;region=initial_region
    weight=wp.int64(1073741824);absorbed=wp.int64(0);escaped=wp.int64(0)
    tau=float(-wp.log(wp.float64(1.0)-(wp.float64(wp.randf(rng))+wp.float64(2.98023223876953125e-8))))
    clock=wp.float64(0.0);events=int(0);hits=int(0);leak=int(0);zero_losses=int(0);tail_zero=int(0)
    while weight>wp.int64(0) and events<50000 and leak==0:
        events+=1
        mua=props[region,0];mus=props[region,1];g=props[region,2];index=props[region,3]
        distance=wp.float64(1.0e20)
        if mus>0.0:distance=wp.float64(tau/mus)
        start=wp.vec3d(p);ray=wp.vec3d(d)
        query=wp.mesh_query_ray(mesh,p,d,float(distance))
        hit=int(0);normal=wp.vec3(0.0);next_region=int(0)
        if query.result:
            hit=1;hits+=1;normal=normals[query.face]
            distance=wp.dot(plane_points[query.face]-start,plane_normals[query.face])/wp.dot(ray,plane_normals[query.face])
            if wp.dot(d,normal)>0.0:
                next_region=front[query.face]
                if back[query.face]!=region:leak=1
            else:
                next_region=back[query.face];normal=-normal
                if front[query.face]!=region:leak=1
        if distance<=wp.float64(0.0) or distance>wp.float64(1.0e19):leak=2
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
                leak=3;break
            crossing=wp.vec3d(wp.float64(1.0e30))
            for axis in range(3):
                if ray[axis]!=wp.float64(0.0):
                    boundary=cell[axis]
                    if ray[axis]>wp.float64(0.0):boundary+=1
                    plane=wp.float64(origin[axis])+wp.float64(boundary)*wp.float64(pitch)
                    crossing[axis]=(plane-start[axis])/ray[axis]
            next_t=wp.min(distance,wp.min(crossing[0],wp.min(crossing[1],crossing[2])))
            length=next_t-traveled
            if length<wp.float64(0.0):leak=4;break
            loss=wp.int64(wp.floor(wp.float64(weight)*(wp.float64(1.0)-wp.exp(-wp.float64(mua)*length))+wp.float64(.5)))
            loss=wp.min(weight,wp.max(wp.int64(0),loss))
            within=wp.clamp((wp.float64(window_ns)-clock)*wp.float64(299.792458)/wp.float64(index),wp.float64(0.0),length)
            score=loss
            if within<length:
                score=wp.int64(wp.floor(wp.float64(weight)*(wp.float64(1.0)-wp.exp(-wp.float64(mua)*within))+wp.float64(.5)))
                score=wp.min(loss,wp.max(wp.int64(0),score))
            wp.atomic_add(absorption,(ix*ny+iy)*nz+iz,score)
            if loss==wp.int64(0):
                zero_losses+=1;tail_zero+=1
            else:tail_zero=0
            weight-=loss;absorbed+=loss
            clock+=length*wp.float64(index)/wp.float64(299.792458)
            traveled=next_t
            if traveled<distance:
                for axis in range(3):
                    if crossing[axis]==next_t:
                        if ray[axis]>wp.float64(0.0):cell[axis]+=1
                        else:cell[axis]-=1
        if traveled<distance and weight>wp.int64(0) and leak==0:leak=5
        if leak!=0 or weight==wp.int64(0):break
        p=wp.vec3(start+ray*distance)
        if hit==1:
            tau=wp.max(0.0,tau-mus*float(distance))
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
                escaped=weight;weight=wp.int64(0);break
            # Explicit microscopic offset at a handled interface; no implicit miss escape.
            p+=d*1.0e-5
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
            tau=float(-wp.log(wp.float64(1.0)-(wp.float64(wp.randf(rng))+wp.float64(2.98023223876953125e-8))))
    terminal[tid,0]=absorbed;terminal[tid,1]=escaped;terminal[tid,2]=weight
    terminal[tid,3]=wp.int64(leak)
    if leak!=0 or weight>wp.int64(0):
        for axis in range(3):
            terminal[tid,4+axis]=wp.int64(wp.float64(p[axis])*wp.float64(1.0e9))
            terminal[tid,7+axis]=wp.int64(wp.float64(d[axis])*wp.float64(1.0e9))
        terminal[tid,10]=wp.int64(distance*wp.float64(1.0e9));terminal[tid,11]=wp.int64(hit)
        terminal[tid,12]=wp.int64(events);terminal[tid,13]=wp.int64(hits)
        terminal[tid,14]=wp.int64(zero_losses);terminal[tid,15]=wp.int64(tail_zero)
    wp.atomic_add(counters,0,wp.int64(events));wp.atomic_add(counters,1,wp.int64(hits))
    if events>=50000 and weight>wp.int64(0):wp.atomic_add(counters,2,wp.int64(1))


def simulate(vertices,faces,properties,face_front,face_back,*,source,direction,initial_region,
             origin,shape,pitch=1.,photons=1000000,seed=20260912,reflect=False,window_ns=5.,device='cuda:0'):
    v=np.asarray(vertices,dtype=np.float32);f=np.asarray(faces,dtype=np.int32)
    prop=np.asarray(properties,dtype=np.float32);front=np.asarray(face_front,dtype=np.int32);back=np.asarray(face_back,dtype=np.int32)
    if v.ndim!=2 or v.shape[1]!=3 or f.ndim!=2 or f.shape[1]!=3:raise ValueError('Triangle arrays required')
    if not np.isfinite(v).all() or f.min()<0 or f.max()>=len(v):raise ValueError('Invalid geometry')
    if prop.ndim!=2 or prop.shape[1]!=4 or not np.isfinite(prop).all():raise ValueError('Optical rows must be mua,mus,g,n')
    if np.any(prop[:,:2]<0) or np.any(np.abs(prop[:,2])>1) or np.any(prop[:,3]<=0):raise ValueError('Invalid optical properties')
    if len(front)!=len(f) or len(back)!=len(f) or min(front.min(),back.min())<0 or max(front.max(),back.max())>=len(prop):raise ValueError('Invalid face region labels')
    if not 0<initial_region<len(prop) or photons<1 or photons*PACKET>np.iinfo(np.int64).max:raise ValueError('Invalid launch/energy count')
    shape=tuple(int(x) for x in shape)
    if len(shape)!=3 or min(shape)<=0 or pitch<=0 or window_ns<=0:raise ValueError('Positive grid/window required')
    normal=np.cross(v[f[:,1]]-v[f[:,0]],v[f[:,2]]-v[f[:,0]])
    lengths=np.linalg.norm(normal,axis=1)
    if np.any(lengths==0):raise ValueError('Degenerate triangle')
    normal/=lengths[:,None]
    plane_point=v[f[:,0]].astype(np.float64)
    plane_normal=np.cross(v[f[:,1]].astype(np.float64)-plane_point,v[f[:,2]].astype(np.float64)-plane_point)
    plane_normal/=np.linalg.norm(plane_normal,axis=1)[:,None]
    direction=np.asarray(direction,dtype=float);direction/=np.linalg.norm(direction)
    mesh=wp.Mesh(points=wp.array(v,dtype=wp.vec3,device=device),indices=wp.array(f.ravel(),dtype=int,device=device))
    absorption=wp.zeros(int(np.prod(shape)),dtype=wp.int64,device=device)
    terminal=wp.zeros((photons,16),dtype=wp.int64,device=device);counters=wp.zeros(3,dtype=wp.int64,device=device)
    wp.launch(_trace,photons,inputs=[mesh.id,wp.array(front,dtype=int,device=device),wp.array(back,dtype=int,device=device),
        wp.array(normal,dtype=wp.vec3,device=device),wp.array(plane_point,dtype=wp.vec3d,device=device),wp.array(plane_normal,dtype=wp.vec3d,device=device),wp.array(prop,dtype=float,device=device),wp.vec3(*source),wp.vec3(*direction),
        initial_region,seed,int(reflect),wp.vec3(*origin),float(pitch),*shape,float(window_ns),absorption,terminal,counters],device=device)
    return {'absorption':absorption.numpy().reshape(shape),'terminal':terminal.numpy(),'counters':counters.numpy()}



if __name__ == "__main__":
    import runpy
    from pathlib import Path
    runpy.run_path(str(Path(__file__).resolve().parents[2]/"probes/optics/photon_boundary_failure_probe.py"),run_name="__main__")
