#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/self_intersections.h>
#include <CGAL/Polygon_mesh_processing/repair_self_intersections.h>
#include <fstream>
#include <iostream>
#include <vector>
using K=CGAL::Exact_predicates_inexact_constructions_kernel;
using M=CGAL::Surface_mesh<K::Point_3>;
int main(int argc,char**argv){
 std::ifstream in(argv[1]); size_t n,f;in>>n>>f;M m;std::vector<M::Vertex_index> vv;
 for(size_t i=0;i<n;i++){double x,y,z;in>>x>>y>>z;vv.push_back(m.add_vertex(K::Point_3(x,y,z)));}
 for(size_t i=0;i<f;i++){size_t a,b,c;in>>a>>b>>c;if(m.add_face(vv[a],vv[b],vv[c])==M::null_face()){std::cerr<<"Invalid face "<<i;return 4;}}
 bool ok=CGAL::Polygon_mesh_processing::experimental::remove_self_intersections(m,CGAL::parameters::number_of_iterations(3).preserve_genus(true));
 m.collect_garbage();std::ofstream out(argv[2]);out.precision(17);out<<m.number_of_vertices()<<" "<<m.number_of_faces()<<"\n";
 for(auto v:m.vertices())out<<m.point(v)<<"\n";
 for(auto f:m.faces()){auto h=m.halfedge(f);out<<m.target(h).idx()<<" "<<m.target(m.next(h)).idx()<<" "<<m.target(m.next(m.next(h))).idx()<<"\n";}
 std::cout<<"{\"success\":"<<(ok?"true":"false")<<"}";
}
