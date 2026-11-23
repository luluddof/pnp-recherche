// Membership of e_target in the row module of a sparse matrix over Z_q, q in {2,3,4} (q=p^k).
// Input (stdin, binary-ish text): q nrows ncols target ; then nnz triples row col val.
// Howell-form style elimination for Z_4; plain Gauss for primes.
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
static int q;
static int val(int a){ // 2-adic valuation for q=4; for prime q: 0 if a!=0
  if(a==0) return 99; if(q==4) return (a&1)?0:1; return 0; }
static int inv(int a){ if(q==2) return 1; if(q==3) return a; /*q=4*/ return a; /*1->1,3->3*/ }
int main(){
  int nr,nc,tgt; long nnz;
  if(scanf("%d %d %d %d %ld",&q,&nr,&nc,&tgt,&nnz)!=5) return 1;
  int cap = nr + nc + 16; // extra rows from saturation (at most one per pivot)
  uint8_t **R = malloc(sizeof(uint8_t*)*cap);
  for(int i=0;i<nr;i++){ R[i]=calloc(nc,1); }
  for(long t=0;t<nnz;t++){ int r,c,v; scanf("%d %d %d",&r,&c,&v); R[r][c]=(uint8_t)(((R[r][c]+v)%q+q)%q); }
  int nrows=nr;
  // echelon: pivot row per column
  int *pivrow = malloc(sizeof(int)*nc); for(int j=0;j<nc;j++) pivrow[j]=-1;
  uint8_t **E = malloc(sizeof(uint8_t*)*nc);
  int top=0; // rows [top..nrows) are remaining
  for(int j=0;j<nc;j++){
    int best=-1,bv=99;
    for(int i=top;i<nrows;i++){ int v=val(R[i][j]); if(v<bv){bv=v;best=i; if(v==0)break;} }
    if(best<0) continue;
    uint8_t *t=R[best]; R[best]=R[top]; R[top]=t; uint8_t *P=R[top]; top++;
    // normalize: multiply by unit so leading = p^bv
    int lead=P[j]; int u = (q==4)? ((lead&1)? inv(lead): ((lead==2)?1:1)) : inv(lead);
    if(u!=1) for(int c=j;c<nc;c++) P[c]=(P[c]*u)%q;
    int pl=P[j]; // 1 or 2
    for(int i=top;i<nrows;i++){ int a=R[i][j]; if(!a) continue; int f=a/pl; // divisible
      uint8_t *X=R[i]; for(int c=j;c<nc;c++) X[c]=(uint8_t)((X[c]+ (q-f)*P[c])%q); }
    E[j]=P; pivrow[j]=1;
    if(q==4 && pl==2){ // saturation: 2*P has zero at col j -> add as new remaining row
      uint8_t *S=calloc(nc,1); int nz=0; for(int c=j+1;c<nc;c++){ S[c]=(2*P[c])%4; nz|=S[c]; }
      if(nz){ R[nrows++]=S; } else free(S);
    }
  }
  // reduce target
  uint8_t *b=calloc(nc,1); b[tgt]=1; int ok=1;
  for(int j=0;j<nc;j++){ if(!b[j]) continue; if(pivrow[j]<0){ok=0;break;}
    int pl=E[j][j]; if(b[j]%pl){ok=0;break;} int f=b[j]/pl;
    for(int c=j;c<nc;c++) b[c]=(uint8_t)((b[c]+(q-f)*E[j][c])%q); }
  printf("%d\n",ok); return 0;
}
