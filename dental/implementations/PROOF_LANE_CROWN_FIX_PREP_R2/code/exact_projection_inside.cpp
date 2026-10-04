#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/intersections.h>
#include <fstream>
#include <iostream>
#include <vector>
using K=CGAL::Exact_predicates_exact_constructions_kernel;
int main(int argc,char**argv){
 std::ifstream in(argv[1]);size_t n;double dx,dy,dz;in>>n>>dx>>dy>>dz;std::vector<K::FT>d={dx,dy,dz};int j=0;for(int k=1;k<3;k++)if(CGAL::abs(d[k])>CGAL::abs(d[j]))j=k;int a=(j+1)%3,b=(j+2)%3;
 auto point=[&](std::istream& s){double x,y,z;s>>x>>y>>z;std::vector<K::FT>q={x,y,z};return K::Point_2(q[a]-d[a]/d[j]*q[j],q[b]-d[b]/d[j]*q[j]);};
 CGAL::Polygon_2<K> poly;for(size_t i=0;i<n;i++)poly.push_back(point(in));std::ifstream q(argv[2]);size_t nv,nf;q>>nv>>nf;std::vector<K::Point_2>p;for(size_t i=0;i<nv;i++)p.push_back(point(q));bool ok=poly.is_simple();size_t bad_points=0,bad_segments=0;
 if(ok){for(auto &x:p)if(poly.bounded_side(x)!=CGAL::ON_BOUNDED_SIDE)bad_points++;for(size_t i=0;i<nv;i++)for(size_t k=i+1;k<nv;k++)if(p[i]!=p[k]){K::Segment_2 e(p[i],p[k]);for(auto it=poly.edges_begin();it!=poly.edges_end();it++)if(CGAL::do_intersect(e,*it))bad_segments++;}}
 std::cout<<"{\"all_pass\":"<<(ok&&!bad_points&&!bad_segments?"true":"false")<<",\"simple\":"<<(ok?"true":"false")<<",\"bad_points\":"<<bad_points<<",\"bad_segments\":"<<bad_segments<<",\"query_vertices\":"<<nv<<"}\n";
}
