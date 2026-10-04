#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include <fstream>
#include <iostream>
#include <vector>
using K=CGAL::Exact_predicates_exact_constructions_kernel;
int main(int argc,char**argv){std::ifstream in(argv[1]);size_t n;double dx,dy,dz;in>>n>>dx>>dy>>dz;std::vector<K::FT>d={dx,dy,dz};int j=0;for(int k=1;k<3;k++)if(CGAL::abs(d[k])>CGAL::abs(d[j]))j=k;int a=(j+1)%3,b=(j+2)%3;CGAL::Polygon_2<K>poly;for(size_t i=0;i<n;i++){double x,y,z;in>>x>>y>>z;std::vector<K::FT>q={x,y,z};poly.push_back(K::Point_2(q[a]-d[a]/d[j]*q[j],q[b]-d[b]/d[j]*q[j]));}std::cout<<"{\"simple\":"<<(poly.is_simple()?"true":"false")<<",\"vertices\":"<<n<<",\"area_exact\":\""<<CGAL::exact(poly.area())<<"\"}\n";}
