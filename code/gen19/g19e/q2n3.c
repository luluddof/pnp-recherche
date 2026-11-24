// Exhaustive over F_2, N=3: min rank(F - G1(x)G1' - G2(x)G2'), G1 in rank normal form.  (crossing, q = 2)
#include <stdio.h>
#include <stdlib.h>
static int kr[8][8];
static inline int rank9(unsigned short *r){ int rk=0; for(int c=0;c<9;c++){ int piv=-1; for(int i=rk;i<9;i++) if((r[i]>>c)&1){piv=i;break;} if(piv<0) continue; unsigned short t=r[piv]; r[piv]=r[rk]; r[rk]=t; for(int i=0;i<9;i++) if(i!=rk && ((r[i]>>c)&1)) r[i]^=r[rk]; rk++; } return rk; }
int main(){
  for(int a=0;a<8;a++) for(int b=0;b<8;b++){ int m=0; for(int u=0;u<3;u++) for(int v=0;v<3;v++) if(((a>>u)&1)&&((b>>v)&1)) m|=1<<(u*3+v); kr[a][b]=m; }
  unsigned short Frow[9]; for(int x=0;x<3;x++) for(int y=0;y<3;y++) Frow[x*3+y]=1<<(y*3+x);  // [x=v][y=u]: u=y, v=x
  long hist[10]={0}; int best=9;
  for(int rho=1;rho<=3;rho++){
    int G1[3]={0,0,0}; for(int i=0;i<rho;i++) G1[i]=1<<i;
    long h[10]={0};
    for(int g1p=0;g1p<512;g1p++) for(int g2=0;g2<512;g2++) for(int g2p=0;g2p<512;g2p++){
      unsigned short r[9];
      for(int x=0;x<3;x++) for(int y=0;y<3;y++)
        r[x*3+y]=Frow[x*3+y]^kr[G1[x]][(g1p>>(3*y))&7]^kr[(g2>>(3*x))&7][(g2p>>(3*y))&7];
      h[rank9(r)]++;
    }
    printf("rank(G1)=%d: histogram of rank(F-C):",rho); for(int i=0;i<10;i++){ printf(" %ld",h[i]); hist[i]+=h[i]; if(h[i]&&i<best) best=i;} printf("\n"); fflush(stdout);
  }
  printf("N=3, F_2, q<=2: exact min p = %d\n",best);
  return 0; }
