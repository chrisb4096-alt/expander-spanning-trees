// Enumerate all spanning trees of a small graph (m <= 64) as edge bitmasks.
// stdin: "n m" then m lines "u v" (0-based). argv[1]: output file (raw uint64 little-endian).
// Backtracking over edges with union-find rollback; prunes by remaining-edge count and by
// a cheap connectivity test of (chosen ∪ undecided) edges.
#include <cstdio>
#include <cstdint>
#include <vector>
#include <numeric>
using namespace std;
int n, m;
vector<int> U, V;
vector<int> par, sz;
vector<pair<int,int>> hist;
vector<uint64_t> out;
int find(int x){ while(par[x]!=x) x=par[x]; return x; }
bool connected_with(uint64_t excluded){
  // is graph minus excluded edges connected?
  vector<int> p(n); iota(p.begin(),p.end(),0);
  auto f=[&](int x){ while(p[x]!=x){ p[x]=p[p[x]]; x=p[x];} return x; };
  int comps=n;
  for(int i=0;i<m;i++) if(!((excluded>>i)&1ULL)){ int a=f(U[i]), b=f(V[i]); if(a!=b){p[a]=b; comps--;} }
  return comps==1;
}
void rec(int i, int cnt, uint64_t mask, uint64_t excl){
  if(cnt==n-1){ out.push_back(mask); return; }
  if(m-i < n-1-cnt) return;
  int a=find(U[i]), b=find(V[i]);
  if(a!=b){
    if(sz[a]>sz[b]) swap(a,b);
    par[a]=b; sz[b]+=sz[a];
    rec(i+1,cnt+1,mask|(1ULL<<i),excl);
    par[a]=a; sz[b]-=sz[a];
  }
  uint64_t ex2 = excl|(1ULL<<i);
  if(connected_with(ex2)) rec(i+1,cnt,mask,ex2);
}
int main(int argc, char** argv){
  if(scanf("%d %d",&n,&m)!=2) return 1;
  U.resize(m); V.resize(m);
  for(int i=0;i<m;i++) if(scanf("%d %d",&U[i],&V[i])!=2) return 1;
  par.resize(n); sz.assign(n,1); iota(par.begin(),par.end(),0);
  rec(0,0,0ULL,0ULL);
  FILE* fo=fopen(argv[1],"wb"); fwrite(out.data(),8,out.size(),fo); fclose(fo);
  printf("%zu\n", out.size());
  return 0;
}
