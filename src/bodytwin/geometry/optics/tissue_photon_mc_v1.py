"""Mesh photon transport with explicit region labels and integer energy tally.

This experimental module requires the documented physical gates before use.
Source position and initial region are explicit; launch energy is inside tissue.
"""
import numpy as np
import warp as wp

PACKET=1 << 30


@wp.kernel
def _trace(mesh:wp.uint64,front:wp.array(dtype=int),back:wp.array(dtype=int),
           normals:wp.array(dtype=wp.vec3),props:wp.array2d(dtype=float),
           source:wp.vec3,direction:wp.vec3,initial_region:int,seed:int,reflect:int,
           origin:wp.vec3,pitch:float,nx:int,ny:int,nz:int,window_ns:float,
           absorption:wp.array(dtype=wp.int64),terminal:wp.array2d(dtype=wp.int64),
           counters:wp.array(dtype=wp.int64)):
    tid=wp.tid();rng=wp.rand_init(seed,tid)
    p=source;d=direction;region=initial_region
    weight=wp.int64(1073741824);absorbed=wp.int64(0);escaped=wp.int64(0)
    tau=-wp.log(wp.max(1.0-wp.randf(rng),1.0e-30))
    clock=float(0.0);events=int(0);hits=int(0);leak=int(0)
    while weight>wp.int64(0) and events<50000 and leak==0:
        events+=1
        mua=props[region,0];mus=props[region,1];g=props[region,2];index=props[region,3]
        distance=float(1.0e20)
        if mus>0.0:distance=tau/mus
        query=wp.mesh_query_ray(mesh,p,d,distance)
        hit=int(0);normal=wp.vec3(0.0);next_region=int(0)
        if query.result:
            hit=1;hits+=1;distance=query.t;normal=normals[query.face]
            if wp.dot(d,normal)>0.0:
                next_region=front[query.face]
                if back[query.face]!=region:leak=1
            else:
                next_region=back[query.face];normal=-normal
                if front[query.face]!=region:leak=1
        if distance<=0.0 or distance>1.0e19:leak=1
        if leak!=0:break
        remaining=distance;segments=int(0)
        while remaining>0.0 and segments<1024 and weight>wp.int64(0):
            segments+=1
            sample=(p-origin+d*1.0e-5)/pitch
            ix=int(wp.floor(sample[0]));iy=int(wp.floor(sample[1]));iz=int(wp.floor(sample[2]))
            if ix<0 or iy<0 or iz<0 or ix>=nx or iy>=ny or iz>=nz:
                leak=1;break
            length=remaining
            for axis in range(3):
                voxel=int(wp.floor(sample[axis]))
                face=origin[axis]+float(voxel)*pitch
                if d[axis]>0.0:face+=pitch
                if wp.abs(d[axis])>1.0e-12:
                    crossing=(face-p[axis])/d[axis]
                    length=wp.min(length,wp.max(crossing,1.0e-5))
            loss=wp.int64(wp.floor(wp.float64(weight)*(wp.float64(1.0)-wp.exp(-wp.float64(mua)*wp.float64(length)))+wp.float64(.5)))
            loss=wp.min(weight,wp.max(wp.int64(0),loss))
            within=wp.clamp((window_ns-clock)*299.792458/index,0.0,length)
            score=loss
            if within<length:
                score=wp.int64(wp.floor(wp.float64(weight)*(wp.float64(1.0)-wp.exp(-wp.float64(mua)*wp.float64(within)))+wp.float64(.5)))
                score=wp.min(loss,wp.max(wp.int64(0),score))
            wp.atomic_add(absorption,(ix*ny+iy)*nz+iz,score)
            weight-=loss;absorbed+=loss
            clock+=length*index/299.792458;p+=d*length
            remaining=wp.max(0.0,remaining-length)
        if remaining>0.0 and weight>wp.int64(0) and leak==0:leak=1
        if leak!=0 or weight==wp.int64(0):break
        if hit==1:
            tau=wp.max(0.0,tau-mus*distance)
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
            tau=-wp.log(wp.max(1.0-wp.randf(rng),1.0e-30))
    terminal[tid,0]=absorbed;terminal[tid,1]=escaped;terminal[tid,2]=weight
    terminal[tid,3]=wp.int64(leak)
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
    direction=np.asarray(direction,dtype=float);direction/=np.linalg.norm(direction)
    mesh=wp.Mesh(points=wp.array(v,dtype=wp.vec3,device=device),indices=wp.array(f.ravel(),dtype=int,device=device))
    absorption=wp.zeros(int(np.prod(shape)),dtype=wp.int64,device=device)
    terminal=wp.zeros((photons,4),dtype=wp.int64,device=device);counters=wp.zeros(3,dtype=wp.int64,device=device)
    wp.launch(_trace,photons,inputs=[mesh.id,wp.array(front,dtype=int,device=device),wp.array(back,dtype=int,device=device),
        wp.array(normal,dtype=wp.vec3,device=device),wp.array(prop,dtype=float,device=device),wp.vec3(*source),wp.vec3(*direction),
        initial_region,seed,int(reflect),wp.vec3(*origin),float(pitch),*shape,float(window_ns),absorption,terminal,counters],device=device)
    return {'absorption':absorption.numpy().reshape(shape),'terminal':terminal.numpy(),'counters':counters.numpy()}
