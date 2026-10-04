#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_triangle_primitive.h>
#include <CGAL/squared_distance_3.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Side_of_triangle_mesh.h>
#include <fstream>
#include <iostream>
#include <vector>
using K=CGAL::Exact_predicates_exact_constructions_kernel;
using Tri=K::Triangle_3;using It=std::vector<Tri>::iterator;
using Tree=CGAL::AABB_tree<CGAL::AABB_traits<K,CGAL::AABB_triangle_primitive<K,It>>>;
int main(int argc,char**argv){
 if(argc<3)return 2;bool surface=argc>3;std::ifstream in(argv[1]);size_t nv,nf;in>>nv>>nf;std::vector<K::Point_3> v;
 for(size_t i=0;i<nv;i++){double x,y,z;in>>x>>y>>z;v.emplace_back(x,y,z);}std::vector<Tri> tris;
 for(size_t i=0;i<nf;i++){size_t a,b,c;in>>a>>b>>c;Tri t(v[a],v[b],v[c]);if(t.is_degenerate()){std::cerr<<"degenerate triangle";return 4;}tris.push_back(t);}
 CGAL::Surface_mesh<K::Point_3> mesh;std::vector<CGAL::Surface_mesh<K::Point_3>::Vertex_index> ids;for(auto p:v)ids.push_back(mesh.add_vertex(p));std::ifstream again(argv[1]);size_t an,af;again>>an>>af;double tmp;for(size_t i=0;i<an*3;i++)again>>tmp;for(size_t i=0;i<af;i++){size_t a,b,c;again>>a>>b>>c;if(mesh.add_face(ids[a],ids[b],ids[c])==mesh.null_face())return 6;}std::unique_ptr<CGAL::Side_of_triangle_mesh<decltype(mesh),K>> side;if(!surface)side.reset(new CGAL::Side_of_triangle_mesh<decltype(mesh),K>(mesh));Tree tree(tris.begin(),tris.end());tree.accelerate_distance_queries();std::ifstream q(argv[2]);size_t n;q>>n;
 std::cout.precision(17);std::cout<<"{\"backend\":\"CGAL EPECK exact binary64-coordinate separation\",\"rows\":[";
 for(size_t i=0;i<n;i++){double x,y,z,r;q>>x>>y>>z>>r;K::Point_3 p(x,y,z);auto cp=tree.closest_point_and_primitive(p);auto d2=CGAL::squared_distance(p,cp.first);K::FT rr=K::FT(r)*K::FT(r);if(i)std::cout<<",";
 std::cout<<"{\"pass\":"<<(d2>=rr?"true":"false")<<",\"inside\":"<<(surface?"null":((*side)(p)==CGAL::ON_BOUNDED_SIDE?"true":"false"))<<",\"facet\":"<<std::distance(tris.begin(),cp.second)<<",\"distance_squared_exact\":\""<<CGAL::exact(d2)<<"\",\"radius_squared_exact\":\""<<CGAL::exact(rr)<<"\"}";}
 std::cout<<"]}\n";
}
