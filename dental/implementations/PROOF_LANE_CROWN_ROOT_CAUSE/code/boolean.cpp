#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/corefinement.h>
#include <CGAL/Polygon_mesh_processing/self_intersections.h>
#include <fstream>
#include <iostream>
#include <vector>
using K=CGAL::Exact_predicates_exact_constructions_kernel;
using M=CGAL::Surface_mesh<K::Point_3>;
namespace PMP=CGAL::Polygon_mesh_processing;
bool read(const char* p,M& m){std::ifstream in(p);size_t n,f;in>>n>>f;std::vector<M::Vertex_index> v;for(size_t i=0;i<n;i++){double x,y,z;in>>x>>y>>z;v.push_back(m.add_vertex(K::Point_3(x,y,z)));}for(size_t i=0;i<f;i++){size_t a,b,c;in>>a>>b>>c;if(m.add_face(v[a],v[b],v[c])==M::null_face())return false;}return true;}
int main(int argc,char**argv){M a,b,c;if(!read(argv[1],a)||!read(argv[2],b))return 4;
 bool ok=PMP::corefine_and_compute_difference(a,b,c);
 if(!ok){std::cout<<"{\"success\":false}";return 2;}
 c.collect_garbage();std::ofstream out(argv[3]);out.precision(17);out<<c.number_of_vertices()<<" "<<c.number_of_faces()<<"\n";
 for(auto v:c.vertices()){auto p=c.point(v);out<<CGAL::to_double(p.x())<<" "<<CGAL::to_double(p.y())<<" "<<CGAL::to_double(p.z())<<"\n";}
 for(auto f:c.faces()){auto h=c.halfedge(f);out<<c.target(h).idx()<<" "<<c.target(c.next(h)).idx()<<" "<<c.target(c.next(c.next(h))).idx()<<"\n";}
 std::vector<std::pair<M::Face_index,M::Face_index>> pp;PMP::self_intersections(c,std::back_inserter(pp));
 std::cout<<"{\"success\":true,\"exact_intersection_count\":"<<pp.size()<<",\"closed\":"<<(CGAL::is_closed(c)?"true":"false")<<"}";
}
