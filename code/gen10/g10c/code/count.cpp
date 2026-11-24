// Plain exact model counter (DPLL + unit propagation, no shortcuts). Reads DIMACS on stdin. G10c.
#include <bits/stdc++.h>
using namespace std;
int nv, nc; vector<vector<int>> cl; vector<vector<int>> occ; vector<signed char> val; vector<int> trail;
typedef unsigned __int128 u128;
bool setlit(int lit, vector<int>& q){ int v=abs(lit); signed char s= lit>0?1:-1; if(val[v]==s) return true; if(val[v]==-s) return false; val[v]=s; trail.push_back(v); q.push_back(lit); return true; }
bool propagate(vector<int>& q){ size_t h=0; while(h<q.size()){ int lit=q[h++]; int v=abs(lit); for(int ci: occ[v]){ int unas=0, last=0; bool sat=false; for(int l: cl[ci]){ int a=abs(l); if(val[a]==0){unas++; last=l;} else if((val[a]>0)==(l>0)){sat=true;break;} } if(sat) continue; if(unas==0) return false; if(unas==1){ if(!setlit(last,q)) return false; } } } return true; }
u128 rec(){ int v=0; for(int i=1;i<=nv;i++) if(val[i]==0){v=i;break;} if(!v) return 1; u128 tot=0; for(int s=-1;s<=1;s+=2){ size_t mark=trail.size(); vector<int> q; bool ok=setlit(s*v,q)&&propagate(q); if(ok) tot+=rec(); while(trail.size()>mark){ val[trail.back()]=0; trail.pop_back(); } } return tot; }
int main(){ string s; while(cin>>s){ if(s=="c"){ getline(cin,s); continue;} if(s=="p"){ cin>>s>>nv>>nc; break; } } cl.resize(nc); occ.assign(nv+1,{}); val.assign(nv+1,0); bool ok=true; vector<int> q;
 for(int i=0;i<nc;i++){ int x; while(cin>>x && x) cl[i].push_back(x); for(int l: cl[i]) occ[abs(l)].push_back(i); }
 for(int i=0;i<nc;i++){ if(cl[i].empty()) ok=false; if(cl[i].size()==1) ok = ok && setlit(cl[i][0],q); }
 ok = ok && propagate(q); u128 r = ok? rec():0; string out; if(r==0) out="0"; while(r){ out.push_back('0'+(int)(r%10)); r/=10;} reverse(out.begin(),out.end()); cout<<out<<"\n"; }
