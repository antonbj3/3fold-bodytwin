"""Analytic synthetic ocular boundaries beside the frozen tissue transport.

The fixed spherical/planar prescription requires its own documented gates.
Launch energy starts in explicit air; the retinal plane terminates packets.
"""
import numpy as np
import warp as wp

PACKET=1 << 30


@wp.func
def _analytic_hit(start:wp.vec3d,ray:wp.vec3d,limit:wp.float64,z:wp.array(dtype=wp.float64),curvature:wp.array(dtype=wp.float64),edge_z:wp.array(dtype=wp.float64),radius:wp.float64):
    best=limit;normal=wp.vec3(0.0);front=int(-1);back=int(-1);face=int(-1)
    a=wp.dot(ray,ray)
    for i in range(6):
        c=curvature[i]
        for branch in range(2):
            t=wp.float64(-1.0)
            if c==wp.float64(0.0):
                if branch==0 and ray[2]!=wp.float64(0.0):t=(z[i]-start[2])/ray[2]
            else:
                radius_s=wp.float64(1.0)/c;center=wp.vec3d(wp.float64(0.0),wp.float64(0.0),z[i]+radius_s)
                oc=start-center;b=wp.dot(ray,oc);disc=b*b-a*(wp.dot(oc,oc)-radius_s*radius_s)
                if disc>=wp.float64(0.0):
                    root=wp.sqrt(disc)
                    if branch==0:t=(-b-root)/a
                    else:t=(-b+root)/a
            if t>=wp.float64(0.0) and t<best:
                p=start+t*ray
                valid=p[0]*p[0]+p[1]*p[1]<=radius*radius
                if c!=wp.float64(0.0):valid=valid and (p[2]-z[i]-wp.float64(1.0)/c)*c<=wp.float64(0.0)
                if valid:
                    best=t;face=i;back=i;front=i+1
                    if i==5:front=0
                    normal=wp.vec3(0.0,0.0,1.0)
                    if c!=wp.float64(0.0):normal=wp.normalize(wp.vec3(wp.vec3d(-p[0]*c,-p[1]*c,wp.float64(1.0)-(p[2]-z[i])*c)))
    axy=ray[0]*ray[0]+ray[1]*ray[1]
    if axy>wp.float64(0.0):
        bxy=start[0]*ray[0]+start[1]*ray[1];disc=bxy*bxy-axy*(start[0]*start[0]+start[1]*start[1]-radius*radius)
        if disc>=wp.float64(0.0):
            for branch in range(2):
                t=(-bxy-wp.sqrt(disc))/axy
                if branch==1:t=(-bxy+wp.sqrt(disc))/axy
                if t>=wp.float64(0.0) and t<best:
                    p=start+t*ray
                    for region in range(1,6):
                        if p[2]>=edge_z[region-1] and p[2]<=edge_z[region]:
                            best=t;face=5+region;back=region;front=0
                            normal=wp.normalize(wp.vec3(wp.vec3d(p[0],p[1],wp.float64(0.0))))
    return face,best,normal,front,back


@wp.kernel
def _trace(z:wp.array(dtype=wp.float64),curvature:wp.array(dtype=wp.float64),edge_z:wp.array(dtype=wp.float64),radius:wp.float64,props:wp.array2d(dtype=float),
           sources:wp.array2d(dtype=float),directions:wp.array2d(dtype=float),launches:int,initial_region:int,seed:int,reflect:int,
           origin:wp.vec3,pitch:float,nx:int,ny:int,nz:int,window_ns:float,
           absorption:wp.array(dtype=wp.int64),terminal:wp.array2d(dtype=wp.int64),
           counters:wp.array(dtype=wp.int64),paths:wp.array2d(dtype=wp.float64),exits:wp.array2d(dtype=wp.float64),history:wp.array2d(dtype=int)):
    tid=wp.tid();rng=wp.rand_init(seed,tid)
    launch=tid%launches
    p=wp.vec3d(wp.float64(sources[launch,0]),wp.float64(sources[launch,1]),wp.float64(sources[launch,2]));d=wp.vec3(directions[launch,0],directions[launch,1],directions[launch,2]);region=initial_region
    reflections=int(0);exit_face=int(-1)
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
        face,boundary_distance,normal,next_front,next_back=_analytic_hit(start,ray,distance,z,curvature,edge_z,radius)
        hit=int(0);next_region=int(0)
        if face>=0:
            hit=1;hits+=1;distance=boundary_distance
            if wp.dot(d,normal)>0.0:
                next_region=next_front
                if next_back!=region:leak=1
            else:
                next_region=next_back;normal=-normal
                if next_front!=region:leak=1
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
            if reflect!=0 and face!=5:
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
            if transmit==0:
                reflections+=1
                d=wp.normalize(d-2.0*wp.dot(d,normal)*normal)
            else:region=next_region
            if region==0:
                exit_face=face
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
    history[tid,0]=reflections;history[tid,1]=exit_face
    terminal[tid,0]=absorbed;terminal[tid,1]=escaped;terminal[tid,2]=weight
    terminal[tid,3]=wp.int64(leak)
    wp.atomic_add(counters,0,wp.int64(events));wp.atomic_add(counters,1,wp.int64(hits))
    if events>=50000 and weight>wp.int64(0):wp.atomic_add(counters,2,wp.int64(1))


def simulate(properties,*,sources,directions,photons=1000000,seed=20260913,reflect=True,
             origin=(-4.,-4.,-1.5),shape=(32,32,104),pitch=.25,window_ns=5.,device='cuda:0'):
    z=np.array([-1.,0.,.6,3.6,7.6,23.6],np.float64)
    curvature=np.array([0.,1/8,1/7,1/10,-1/6,0.],np.float64);radius=3.
    edge_z=z+curvature*radius**2/(1+np.sqrt(1-curvature**2*radius**2))
    prop=np.asarray(properties,dtype=np.float32);source=np.asarray(sources,dtype=np.float32);direction=np.array(directions,dtype=np.float32,copy=True)
    if prop.shape!=(6,4) or not np.isfinite(prop).all() or np.any(prop[:,:2]<0) or np.any(np.abs(prop[:,2])>1) or np.any(prop[:,3]<=0):raise ValueError('Six finite optical rows required')
    if source.ndim!=2 or source.shape[1]!=3 or not len(source) or direction.shape!=source.shape or not np.isfinite(source).all() or not np.isfinite(direction).all():raise ValueError('Finite matching launch vectors required')
    if np.any(np.sum(source[:,:2]**2,axis=1)>=radius**2) or np.any(source[:,2]<=-1) or np.any(source[:,2]>=0):raise ValueError('Sources must be strictly inside the launch air slab')
    lengths=np.linalg.norm(direction,axis=1)
    if np.any(lengths==0):raise ValueError('Nonzero directions required')
    direction/=lengths[:,None]
    if photons<1 or photons*PACKET>np.iinfo(np.int64).max:raise ValueError('Bounded positive photon count required')
    shape=tuple(int(x) for x in shape)
    if len(shape)!=3 or min(shape)<=0 or pitch<=0 or window_ns<=0:raise ValueError('Positive grid/window required')
    absorption=wp.zeros(int(np.prod(shape)),dtype=wp.int64,device=device)
    terminal=wp.zeros((photons,4),dtype=wp.int64,device=device);counters=wp.zeros(3,dtype=wp.int64,device=device)
    paths=wp.zeros((photons,5),dtype=wp.float64,device=device);exits=wp.zeros((photons,7),dtype=wp.float64,device=device);history=wp.zeros((photons,2),dtype=int,device=device)
    wp.launch(_trace,photons,inputs=[wp.array(z,dtype=wp.float64,device=device),wp.array(curvature,dtype=wp.float64,device=device),wp.array(edge_z,dtype=wp.float64,device=device),wp.float64(radius),wp.array(prop,dtype=float,device=device),
        wp.array(source,dtype=float,device=device),wp.array(direction,dtype=float,device=device),len(source),1,seed,int(reflect),wp.vec3(*origin),float(pitch),*shape,float(window_ns),absorption,terminal,counters,paths,exits,history],device=device)
    return {'absorption':absorption.numpy().reshape(shape),'terminal':terminal.numpy(),'counters':counters.numpy(),'paths':paths.numpy(),'exits':exits.numpy(),'history':history.numpy()}


if __name__=='__main__':
    import runpy
    from pathlib import Path
    runpy.run_path(str(Path(__file__).resolve().parents[4]/'probes/optics/ocular_transport_probe.py'),run_name='__main__')
