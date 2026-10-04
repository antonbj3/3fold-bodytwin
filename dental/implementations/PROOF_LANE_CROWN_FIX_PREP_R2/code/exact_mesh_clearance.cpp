#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_triangle_primitive.h>
#include <CGAL/squared_distance_3.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Side_of_triangle_mesh.h>
#include <CGAL/Polygon_mesh_processing/connected_components.h>
#include <fstream>
#include <iostream>
#include <vector>
using K=CGAL::Exact_predicates_exact_constructions_kernel;using P=K::Point_3;using Tri=K::Triangle_3;using It=std::vector<Tri>::iterator;using M=CGAL::Surface_mesh<P>;
using Tree=CGAL::AABB_tree<CGAL::AABB_traits<K,CGAL::AABB_triangle_primitive<K,It>>>;
bool read(const char* path,std::vector<Tri>& tt,M& m){std::ifstream in(path);size_t nv,nf;if(!(in>>nv>>nf))return false;std::vector<P>v;std::vector<M::Vertex_index>ids;for(size_t i=0;i<nv;i++){double x,y,z;in>>x>>y>>z;v.emplace_back(x,y,z);ids.push_back(m.add_vertex(v.back()));}for(size_t i=0;i<nf;i++){size_t a,b,c;in>>a>>b>>c;Tri t(v[a],v[b],v[c]);if(t.is_degenerate())return false;tt.push_back(t);if(m.add_face(ids[a],ids[b],ids[c])==M::null_face())return false;}return true;}
K::FT decimal(std::string s){bool neg=s[0]=='-';if(neg)s=s.substr(1);auto j=s.find('.');long den=1;if(j!=std::string::npos){for(size_t i=j+1;i<s.size();i++)den*=10;s.erase(j,1);}return K::FT(std::stol(s)*(neg?-1:1))/K::FT(den);}
int main(int argc,char**argv){if(argc<4)return 2;std::vector<Tri>T,Q;M tm,qm;if(!read(argv[1],T,tm)||!read(argv[2],Q,qm))return 3;bool surface=argc>4;K::FT w=decimal(argv[3]),w2=w*w;Tree tree(T.begin(),T.end());size_t count=0,bad=0,firstq=0,firstt=0;K::FT firstd=0;std::vector<It>near;for(size_t i=0;i<Q.size();i++){auto t=Q[i];K::FT lo[3],hi[3];for(int k=0;k<3;k++){lo[k]=hi[k]=t[0][k];for(int j=1;j<3;j++){lo[k]=CGAL::min(lo[k],t[j][k]);hi[k]=CGAL::max(hi[k],t[j][k]);}lo[k]-=w;hi[k]+=w;}K::Iso_cuboid_3 box(P(lo[0],lo[1],lo[2]),P(hi[0],hi[1],hi[2]));near.clear();tree.all_intersected_primitives(box,std::back_inserter(near));for(auto it:near){count++;K::FT d=CGAL::squared_distance(t,*it);if(d<w2){if(bad==0){firstq=i;firstt=it-T.begin();firstd=d;}bad++;break;}}}
 bool inside=true;size_t tc=0,qc=0;if(!surface){if(!CGAL::is_closed(tm)||!CGAL::is_closed(qm))return 4;CGAL::Side_of_triangle_mesh<M,K> side(tm);inside=side(Q[0][0])==CGAL::ON_BOUNDED_SIDE;auto tf=tm.add_property_map<M::Face_index,size_t>("f:cc",0).first;auto qf=qm.add_property_map<M::Face_index,size_t>("f:cc",0).first;tc=CGAL::Polygon_mesh_processing::connected_components(tm,tf);qc=CGAL::Polygon_mesh_processing::connected_components(qm,qf);}
 std::cout.precision(17);std::cout<<"{\"backend\":\"CGAL EPECK exact triangle-triangle; exact expanded-box pruning\",\"all_pass\":"<<((bad==0&&inside&&(surface||(tc==1&&qc==1)))?"true":"false")<<",\"boundary_clearance_pass\":"<<(bad==0?"true":"false")<<",\"inside_witness\":"<<(surface?"null":inside?"true":"false")<<",\"support_boundary_components\":"<<tc<<",\"query_boundary_components\":"<<qc<<",\"query_triangles\":"<<Q.size()<<",\"tested_pairs\":"<<count<<",\"bad_query_triangles\":"<<bad<<",\"required_squared_exact\":\""<<CGAL::exact(w2)<<"\",\"first_bad_query\":"<<firstq<<",\"first_bad_support\":"<<firstt<<",\"first_bad_squared_exact\":\""<<CGAL::exact(firstd)<<"\"}\n";
}
