/* Geometrie de l'espace des solutions : 3-SAT vs 3-XORSAT.
   Ilots = composantes connexes du graphe (solutions, voisines si distance de Hamming 1). */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
typedef uint64_t u64;
#define CAP 400000
int n,m; int cl[400][3]; u64 sols[CAP]; int S; int par[CAP];
int val[64];
static int find(int x){while(par[x]!=x){par[x]=par[par[x]];x=par[x];}return x;}
static int cmp(const void*a,const void*b){u64 x=*(u64*)a,y=*(u64*)b;return x<y?-1:x>y;}
/* ---- 3-SAT : enumeration par retour arriere avec verification partielle ---- */
static int clause_dead(int c){for(int j=0;j<3;j++){int l=cl[c][j],v=val[abs(l)]; if(v==0) return 0; if((l>0)==(v>0)) return 0;} return 1;}
int **cv; int *ncv;
static void enum_sat(int v,u64 mask){
  if(S>=CAP) return;
  if(v>n){sols[S++]=mask;return;}
  for(int b=1;b>=-1;b-=2){val[v]=b; int ok=1;
    for(int k=0;k<ncv[v];k++) if(clause_dead(cv[v][k])){ok=0;break;}
    if(ok) enum_sat(v+1, b>0? mask|(1ULL<<(v-1)) : mask);}
  val[v]=0;}
/* ---- 3-XORSAT : x_a ^ x_b ^ x_c = rhs, enumeration de l'espace affine ---- */
u64 row[400]; int rhs[400];
static void enum_xor(void){
  int r=0, pivcol[64]; u64 R[400]; int B[400];
  for(int i=0;i<m;i++){R[i]=row[i];B[i]=rhs[i];}
  for(int c=0;c<n && r<m;c++){int p=-1; for(int i=r;i<m;i++) if(R[i]>>c&1){p=i;break;}
    if(p<0) continue; u64 t=R[p];R[p]=R[r];R[r]=t; int tb=B[p];B[p]=B[r];B[r]=tb;
    for(int i=0;i<m;i++) if(i!=r && (R[i]>>c&1)){R[i]^=R[r];B[i]^=B[r];}
    pivcol[r++]=c;}
  for(int i=r;i<m;i++) if(B[i]) {S=0;return;}
  int isp[64]={0}; for(int i=0;i<r;i++) isp[pivcol[i]]=1;
  int fr[64],nf=0; for(int c=0;c<n;c++) if(!isp[c]) fr[nf++]=c;
  if(nf>18){nf=18;} /* borne */
  S=0;
  for(u64 f=0; f<(1ULL<<nf) && S<CAP; f++){u64 x=0;
    for(int j=0;j<nf;j++) if(f>>j&1) x|=1ULL<<fr[j];
    for(int i=0;i<r;i++){int b=B[i]; u64 o=R[i]&~(1ULL<<pivcol[i])&x; b^=__builtin_popcountll(o)&1; if(b) x|=1ULL<<pivcol[i];}
    sols[S++]=x;}}
/* ---- analyse : ilots + histogramme des distances ---- */
static void analyse(int *ncomp,double *big,double *gapfrac){
  qsort(sols,S,sizeof(u64),cmp);
  for(int i=0;i<S;i++) par[i]=i;
  for(int i=0;i<S;i++) for(int b=0;b<n;b++){u64 y=sols[i]^(1ULL<<b);
    u64*f=bsearch(&y,sols,S,sizeof(u64),cmp); if(f){int j=f-sols; int a=find(i),c=find(j); if(a!=c) par[a]=c;}}
  static int cnt[CAP]; memset(cnt,0,sizeof(int)*S); int nc=0,mx=0;
  for(int i=0;i<S;i++){int r=find(i); if(cnt[r]++==0) nc++; if(cnt[r]>mx) mx=cnt[r];}
  *ncomp=nc; *big=(double)mx/S;
  /* fraction de paires a distance "intermediaire" (entre 15% et 35% de n) */
  long tot=0,mid=0; for(int t=0;t<20000;t++){int i=rand()%S,j=rand()%S; int d=__builtin_popcountll(sols[i]^sols[j]);
    tot++; if(d>0.15*n && d<0.35*n) mid++;}
  *gapfrac=(double)mid/tot;}
int main(int argc,char**argv){
  int type=atoi(argv[1]); n=atoi(argv[2]); int T=atoi(argv[4]); double a=atof(argv[3]); m=(int)(a*n+0.5);
  srand(777+(int)(a*1000)+type);
  cv=malloc(sizeof(int*)*(n+1)); ncv=malloc(sizeof(int)*(n+1)); for(int v=0;v<=n;v++) cv[v]=malloc(sizeof(int)*m);
  double sc=0,sb=0,sg=0,ss=0; int got=0;
  for(int t=0;t<T*20 && got<T;t++){
    for(int c=0;c<m;c++){int x[3]; for(int j=0;j<3;j++){int ok;do{ok=1;x[j]=1+rand()%n;for(int i=0;i<j;i++) if(x[i]==x[j]) ok=0;}while(!ok);}
      for(int j=0;j<3;j++) cl[c][j]=(rand()&1)?x[j]:-x[j];
      row[c]=(1ULL<<(x[0]-1))|(1ULL<<(x[1]-1))|(1ULL<<(x[2]-1)); rhs[c]=rand()&1;}
    S=0;
    if(type==0){ memset(ncv,0,sizeof(int)*(n+1)); memset(val,0,sizeof val);
      for(int c=0;c<m;c++){int mxv=0; for(int j=0;j<3;j++) if(abs(cl[c][j])>mxv) mxv=abs(cl[c][j]); cv[mxv][ncv[mxv]++]=c;}
      enum_sat(1,0);}
    else enum_xor();
    if(S<2||S>=CAP) continue;
    int nc; double big,g; analyse(&nc,&big,&g); sc+=nc; sb+=big; sg+=g; ss+=S; got++;}
  if(!got){printf("%s n=%d alpha=%.2f : aucune instance utilisable\n",type?"XORSAT":"3-SAT ",n,a);return 0;}
  printf("%s n=%d alpha=%.2f | solutions moy=%8.0f | ilots moy=%7.1f | plus grand ilot=%5.1f%% | paires a distance moyenne=%5.1f%%\n",
    type?"XORSAT":"3-SAT ",n,a,ss/got,sc/got,100*sb/got,100*sg/got);
  return 0;}
