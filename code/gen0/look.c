#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#define MAXN 400
#define MAXM 2000
int n,m, cl[MAXM][3];
int occ[2*MAXN+2][MAXM], nocc[2*MAXN+2];
int val[MAXN+1], trail[MAXN+1], tlen;
long long props, nodes;
static int idx(int l){return l>0?2*l:2*(-l)+1;}
static int litval(int l){int v=val[abs(l)]; return l>0?v:-v;}
static int setlit(int l){int lv=litval(l); if(lv==1) return 1; if(lv==-1) return 0; val[abs(l)]=l>0?1:-1; trail[tlen++]=l; return 1;}
static void undo(int to){while(tlen>to){tlen--; val[abs(trail[tlen])]=0;}}
static int propagate(int q){
  props++;
  while(q<tlen){int l=trail[q++]; int ni=idx(-l);
    for(int k=0;k<nocc[ni];k++){int c=occ[ni][k], un=0, ul=0, sat=0;
      for(int j=0;j<3;j++){int x=litval(cl[c][j]); if(x==1){sat=1;break;} if(x==0){un++;ul=cl[c][j];}}
      if(sat) continue; if(un==0) return 0; if(un==1) setlit(ul);}}
  return 1;}
/* refuted by depth-d lookahead ? */
static int fails(int d){
  if(d==0||tlen==n) return 0;
  for(int v=1;v<=n;v++) if(!val[v]){
    int s=tlen; int f1=!(setlit(v)&&propagate(s))||fails(d-1); undo(s);
    if(f1){int f2=!(setlit(-v)&&propagate(s))||fails(d-1); undo(s); if(f2) return 1;}
  }
  return 0;}
/* polynomial candidate: 0=UNSAT 1=SAT 2=UNKNOWN */
static int candidate(int k){
  tlen=0; memset(val,0,sizeof val); int guessed=0;
  while(tlen<n){
    int forced=0;
    for(int v=1;v<=n;v++) if(!val[v]){
      int s=tlen;
      int tok=setlit(v)&&propagate(s)&&!fails(k-1); undo(s);
      int fok=setlit(-v)&&propagate(s)&&!fails(k-1); undo(s);
      if(!tok&&!fok) return guessed?2:0;
      if(!tok){setlit(-v);propagate(s);forced=1;}
      else if(!fok){setlit(v);propagate(s);forced=1;}
    }
    if(!forced&&tlen<n){
      int v=1; while(val[v]) v++;
      int sc=nocc[idx(v)]-nocc[idx(-v)];
      int s=tlen; guessed=1;
      if(!(setlit(sc>=0?v:-v)&&propagate(s))) return 2;
    }
  }
  return 1;}
static int dpll(void){
  nodes++; if(tlen==n) return 1;
  int v=1; while(val[v]) v++;
  int s=tlen;
  if(setlit(v)&&propagate(s)&&dpll()) return 1; undo(s);
  if(setlit(-v)&&propagate(s)&&dpll()) return 1; undo(s);
  return 0;}
static void gen(void){
  memset(nocc,0,sizeof nocc);
  for(int c=0;c<m;c++){int a[3];
    for(int j=0;j<3;j++){int ok; do{ok=1;a[j]=1+rand()%n; for(int i=0;i<j;i++) if(a[i]==a[j]) ok=0;}while(!ok);
      int l=(rand()&1)?a[j]:-a[j]; cl[c][j]=l; int ix=idx(l); occ[ix][nocc[ix]++]=c;}}}
int main(int argc,char**argv){
  int k=atoi(argv[1]); int T=atoi(argv[3]); n=atoi(argv[2]); m=(int)(4.26*n+0.5);
  srand(12345+n);
  int ok=0,sat=0,okS=0,okU=0; double tc=0; long long tn=0;
  for(int t=0;t<T;t++){
    gen(); tlen=0; memset(val,0,sizeof val); nodes=0; int truth=dpll(); tn+=nodes; sat+=truth;
    clock_t c0=clock(); int r=candidate(k); tc+=(double)(clock()-c0)/CLOCKS_PER_SEC;
    if(r==1&&truth){ok++;okS++;} if(r==0&&!truth){ok++;okU++;}
  }
  printf("k=%d n=%3d  formules=%d  SAT=%2d  candidat correct=%2d/%d (%3.0f%%)  temps candidat moyen=%.3fs  noeuds DPLL moyens=%lld | SAT trouves=%d/%d UNSAT prouves=%d/%d\n",
    k,n,T,sat,ok,T,100.0*ok/T,tc/T,tn/T,okS,sat,okU,T-sat);
  fflush(stdout); return 0;}
