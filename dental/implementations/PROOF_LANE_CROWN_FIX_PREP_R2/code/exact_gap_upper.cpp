#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_triangle_primitive.h>
#include <CGAL/squared_distance_3.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <algorithm>
using K=CGAL::Exact_predicates_exact_constructions_kernel;using P=K::Point_3;using Tri=K::Triangle_3;using It=std::vector<Tri>::iterator;using Tree=CGAL::AABB_tree<CGAL::AABB_traits<K,CGAL::AABB_triangle_primitive<K,It>>>;
std::vector<Tri>read(const char*p){std::ifstream in(p);size_t n,f;in>>n>>f;std::vector<P>v;for(size_t i=0;i<n;i++){double x,y,z;in>>x>>y>>z;v.emplace_back(x,y,z);}std::vector<Tri>out;for(size_t i=0;i<f;i++){size_t a,b,c;in>>a>>b>>c;out.emplace_back(v[a],v[b],v[c]);}return out;}
size_t cells=0,unknown=0,witnesses=0;int deepest=0;K::FT bad_d=0;size_t badface=0;Tree*tree;K::FT bound;
bool cover(const Tri&t,int dep,size_t idx){cells++;deepest=std::max(dep,deepest);std::vector<It>cand;P center=CGAL::centroid(t);std::vector<P>pts={t[0],t[1],t[2],center};for(auto p:pts){auto cp=tree->closest_point_and_primitive(p);K::FT d=CGAL::squared_distance(p,cp.first);if(d>bound){if(witnesses++==0){bad_d=d;badface=idx;}return false;}cand.push_back(cp.second);}for(auto it:cand){bool ok=true;for(int k=0;k<3;k++)if(CGAL::squared_distance(t[k],*it)>bound){ok=false;break;}if(ok)return true;}if(dep>=6){unknown++;return false;}P a=CGAL::midpoint(t[0],t[1]),b=CGAL::midpoint(t[1],t[2]),c=CGAL::midpoint(t[2],t[0]);bool ok=true;ok=cover(Tri(t[0],a,c),dep+1,idx)&&ok;ok=cover(Tri(a,t[1],b),dep+1,idx)&&ok;ok=cover(Tri(c,b,t[2]),dep+1,idx)&&ok;ok=cover(Tri(a,b,c),dep+1,idx)&&ok;return ok;}
int main(int argc,char**argv){if(argc<3)return 2;auto source=read(argv[1]);auto query=read(argv[2]);Tree tr(source.begin(),source.end());tr.accelerate_distance_queries();tree=&tr;bound=K::FT(36)/K::FT(10000);size_t bad=0;for(size_t i=0;i<query.size();i++){if(query[i].is_degenerate())return 3;if(!cover(query[i],0,i))bad++;}std::cout<<"{\"backend\":\"CGAL EPECK continuous convex-target facet cover\",\"all_pass\":"<<(bad==0?"true":"false")<<",\"query_triangles\":"<<query.size()<<",\"failed_triangles\":"<<bad<<",\"cells\":"<<cells<<",\"maximum_depth\":"<<deepest<<",\"unresolved_cells\":"<<unknown<<",\"exact_refuting_points\":"<<witnesses<<",\"first_bad_query_face\":"<<badface<<",\"first_bad_point_squared_distance_exact\":\""<<CGAL::exact(bad_d)<<"\",\"certified_upper_squared_exact\":\""<<CGAL::exact(bound)<<"\"}\n";}
