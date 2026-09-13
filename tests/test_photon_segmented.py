"""Boundary, tie and fallback controls for segmented photon candidate queries."""
import numpy as np
import pytest
wp = pytest.importorskip('warp')
from bodytwin.geometry.optics import tissue_photon_paths_v1 as baseline
from bodytwin.geometry.optics import tissue_photon_segmented_v1 as candidate
from bodytwin.geometry.optics import tissue_photon_window2_v1 as window_candidate
from bodytwin.geometry.optics import tissue_photon_cell_v1 as cell_candidate

@wp.kernel
def compare(mesh:wp.uint64,q:wp.array2d(dtype=wp.float64),out:wp.array2d(dtype=int)):
    i=wp.tid();p=wp.vec3d(q[i,0],q[i,1],q[i,2]);d=wp.vec3d(q[i,3],q[i,4],q[i,5])
    out[i,0]=baseline._candidate_hit(mesh,p,d,q[i,6])
    out[i,1]=candidate._candidate_hit(mesh,p,d,q[i,6])


@wp.kernel
def compare_window2(mesh:wp.uint64,q:wp.array2d(dtype=wp.float64),out:wp.array2d(dtype=int)):
    i=wp.tid();p=wp.vec3d(q[i,0],q[i,1],q[i,2]);d=wp.vec3d(q[i,3],q[i,4],q[i,5])
    out[i,0]=baseline._candidate_hit(mesh,p,d,q[i,6])
    out[i,1]=window_candidate._candidate_hit(mesh,p,d,q[i,6])


@wp.kernel
def compare_cell(mesh:wp.uint64,q:wp.array2d(dtype=wp.float64),out:wp.array2d(dtype=int)):
    i=wp.tid();p=wp.vec3d(q[i,0],q[i,1],q[i,2]);d=wp.vec3d(q[i,3],q[i,4],q[i,5])
    out[i,0]=baseline._candidate_hit(mesh,p,d,q[i,6])
    out[i,1]=cell_candidate._candidate_hit(mesh,p,d,q[i,6])


def check_boundaries(device='cpu', window2=False, cell_certified=False):
    vertices=np.array([[x,y,z] for x in (8,16) for y,z in ((-1,-1),(1,-1),(1,1),(-1,1))],dtype=np.float32)
    faces=np.array([[0,1,2],[0,2,3],[4,5,6],[4,6,7],[0,1,2]],dtype=np.int32)
    queries=np.array([[0,0,0,1,0,0,8-1e-9],[0,0,0,1,0,0,8],[0,0,0,1,0,0,8+1e-9],
                      [8,0,0,1,0,0,20],[8+1e-6,0,0,1,0,0,20],[16,0,0,-1,0,0,20],
                      [0,0,0,0,1,0,1e20],[0,2,0,1,0,0,2000],[0,0,0,1e-8,0,0,1e20],
                      [-1000,0,0,1,0,0,2000]],dtype=np.float64)
    expected=np.array([-1,0,0,0,2,2,-1,-1,0,0])
    mesh=wp.Mesh(points=wp.array(vertices,dtype=wp.vec3,device=device),indices=wp.array(faces.ravel(),dtype=int,device=device))
    q=wp.array(queries,dtype=wp.float64,device=device);out=wp.zeros((len(q),2),dtype=int,device=device)
    arrays=[]
    for _ in range(2):
        wp.launch(compare_cell if cell_certified else compare_window2 if window2 else compare,len(q),inputs=[mesh.id,q,out],device=device);arrays.append(out.numpy())
    assert np.array_equal(arrays[0][:,0],expected)
    assert np.array_equal(arrays[0][:,1],expected)
    assert arrays[0].tobytes()==arrays[1].tobytes()
    return arrays[0]


def test_segmented_boundary_tie_and_fallback():
    check_boundaries('cpu')


def check_prepared(device='cpu'):
    from bodytwin.geometry.label_interfaces_v1 import label_interfaces
    from bodytwin.geometry.optics.tissue_photon_prepared_v1 import PreparedPhotonScene
    mesh=label_interfaces(np.ones((3,3,3),dtype=np.uint8))
    v,f,front,back=(mesh[k] for k in ('vertices','faces','front','back'))
    properties=np.array([[0,0,1,1],[.01,.5,.8,1.37]])
    kwargs=dict(source=[1.5,1.5,1.5],direction=[.2,.3,.7],initial_region=1,origin=[0,0,0],shape=(3,3,3),photons=32,seed=20260913,reflect=True)
    reference=candidate.simulate(v,f,properties,front,back,device=device,**kwargs)
    scene=PreparedPhotonScene(v,f,properties,front,back,device=device)
    first=scene.simulate(**kwargs);before={k:a.tobytes() for k,a in first.items()}
    # The scene owns a snapshot and retained outputs are independent of later calls.
    v.fill(99);f.fill(0);properties.fill(0);front.fill(0);back.fill(0)
    second=scene.simulate(**kwargs)
    for k in reference:
        assert reference[k].tobytes()==before[k]==first[k].tobytes()==second[k].tobytes()
    scene.simulate(**dict(kwargs,photons=16,seed=20260914))
    for k in first:assert first[k].tobytes()==before[k]
    for change in ({'initial_region':0},{'photons':0},{'shape':(0,3,3)},{'pitch':0},{'window_ns':0}):
        with pytest.raises(ValueError):scene.simulate(**dict(kwargs,**change))
    return before


def test_prepared_snapshot_and_outputs():
    check_prepared('cpu')


def test_window2_boundary_tie_and_fallback():
    check_boundaries('cpu',window2=True)


def test_cell_certificate_refuses_unaligned():
    vertices=np.array([[.25,0,0],[.25,1,0],[.25,0,1]])
    with pytest.raises(ValueError,match='integer axis plane'):
        cell_candidate.PreparedPhotonScene(vertices,[[0,1,2]],[[0,0,0,1],[.1,1,0,1]],[0],[1],device='cpu')


def test_cell_boundary_tie_and_fallback():
    check_boundaries('cpu',cell_certified=True)
