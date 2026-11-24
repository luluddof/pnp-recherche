import numpy as np
L3=np.log2(3)
def H(x): return 0 if x<=0 or x>=1 else -x*np.log2(x)-(1-x)*np.log2(1-x)
def KL(a,b): return a*np.log2(a/b)+(1-a)*np.log2((1-a)/(1-b))
print("one-bit Kronecker (window [0,D]):  kappa theta T theta+T circleT")
for k in np.arange(0.25,0.76,0.05):
    th=H(k)+(1-k)*L3; T=H(k)+k*L3
    print(f"{k:.2f} {th:.3f} {T:.3f} {th+T:.3f} {np.sqrt(max(0,4-th*th)):.3f}")
print("centred window (AW-type) degree D=kappa n, window [n/4-D/2, n/4+D/2]:")
for k in [0.01,0.02,0.05,0.1,0.15,0.2,0.3,0.4,0.5]:
    th=2-KL(0.25+k/2,0.25); T=H(k)+k*L3
    print(f"{k:.2f} {th:.3f} {T:.3f} {th+T:.3f} circleT={np.sqrt(max(0,4-th*th)):.3f}  lineT={2-th:.3f}")
