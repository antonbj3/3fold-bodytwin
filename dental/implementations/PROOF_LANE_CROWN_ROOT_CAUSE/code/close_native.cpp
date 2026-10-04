#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/triangulate_hole.h>
#include <CGAL/Polygon_mesh_processing/border.h>
#include <CGAL/Polygon_mesh_processing/self_intersections.h>
#include <fstream>
#include <iostream>
#include <vector>
using K=CGAL::Exact_predicates_inexact_constructions_kernel;
using M=CGAL::Surface_mesh<K::Point_3>;
namespace PMP=CGAL::Polygon_mesh_processing;
int main(int argc,char**argv){
 std::ifstream in(argv[1]);size_t nv,nf;in>>nv>>nf;M m;std::vector<M::Vertex_index> vv;
 for(size_t i=0;i<nv;i++){double x,y,z;in>>x>>y>>z;vv.push_back(m.add_vertex(K::Point_3(x,y,z)));}
 for(size_t i=0;i<nf;i++){size_t a,b,c;in>>a>>b>>c;if(m.add_face(vv[a],vv[b],vv[c])==M::null_face())return 4;}
 std::vector<M::Halfedge_index> borders;PMP::extract_boundary_cycles(m,std::back_inserter(borders));
 for(auto h:borders){std::vector<M::Face_index> faces;PMP::triangulate_hole(m,h,std::back_inserter(faces));}
 std::vector<std::pair<M::Face_index,M::Face_index>> pairs;PMP::self_intersections(m,std::back_inserter(pairs));
 std::ofstream out(argv[2]);out.precision(17);out<<m.number_of_vertices()<<" "<<m.number_of_faces()<<"\n";
 for(auto v:m.vertices())out<<m.point(v)<<"\n";
 for(auto f:m.faces()){auto h=m.halfedge(f);out<<m.target(h).idx()<<" "<<m.target(m.next(h)).idx()<<" "<<m.target(m.next(m.next(h))).idx()<<"\n";}
 std::cout<<"{\"border_count\":"<<borders.size()<<",\"intersection_count\":"<<pairs.size()<<",\"original_facets\":"<<nf<<",\"final_facets\":"<<m.number_of_faces()<<",\"closed\":"<<(CGAL::is_closed(m)?"true":"false")<<"}";
}
