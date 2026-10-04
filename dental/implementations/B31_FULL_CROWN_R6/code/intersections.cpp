#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/self_intersections.h>
#include <fstream>
#include <iostream>
#include <vector>
using K=CGAL::Exact_predicates_inexact_constructions_kernel;
using M=CGAL::Surface_mesh<K::Point_3>;
int main(int argc,char**argv){
 std::ifstream in(argv[1]); size_t n,f;in>>n>>f;M m;std::vector<M::Vertex_index> vv;
 for(size_t i=0;i<n;i++){double x,y,z;in>>x>>y>>z;vv.push_back(m.add_vertex(K::Point_3(x,y,z)));}
 for(size_t i=0;i<f;i++){size_t a,b,c;in>>a>>b>>c;if(m.add_face(vv[a],vv[b],vv[c])==M::null_face()){std::cerr<<"Invalid face "<<i;return 4;}}
 std::vector<std::pair<M::Face_index,M::Face_index>> pairs;
 CGAL::Polygon_mesh_processing::self_intersections(m,std::back_inserter(pairs));
 std::cout<<"{\"count\":"<<pairs.size()<<",\"pairs\":[";bool first=true;
 for(auto p:pairs){if(!first)std::cout<<",";first=false;std::cout<<"["<<p.first.idx()<<","<<p.second.idx()<<"]";}
 std::cout<<"]}";
}
