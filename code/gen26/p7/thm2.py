import math
z2=math.pi**2/6
def v2(n): return (n&-n).bit_length()-1
P=0; worst=0.0  # P = sum_{k<K} popcount(k)
for K in range(1,2**20+1):
    P+=bin(K-1).count('1')
    if K>=2: worst=max(worst,P/((K/2)*math.log2(K)))
print("max P(K)/((K/2)log2K), K<=2^20:",worst)
# (2) sum_{k<K} s_k - 2K^2 >= 0 and identity = K(m-1)-2P(K); (3) A(Y)>=0 min; (4) G(theta) <= z2/(4theta)
worstG=0
for m in range(1,19):
    D=2**m; s=[]; vk=0
    for k in range(D//2):
        if k: vk+=v2(k)
        s.append(m+1+2*k+2*vk)
    assert all(s[k]<s[k+1] for k in range(len(s)-1))
    S=0; P=0; mn=10**18
    for K in range(1,D//2+1):
        S+=s[K-1]; P+=bin(K-1).count('1')
        assert S-2*K*K==K*(m-1)-2*P
        mn=min(mn,S-2*K*K)
    # G(theta)
    for e in range(-60,41):
        th=2.0**(e/4)
        G=sum(-math.log1p(-math.exp(-th*x)) for x in s if th*x<700)
        worstG=max(worstG,G*4*th/z2)
    print(m,"min_K(sum s - 2K^2)=",mn,"end value",S-2*(D//2)**2)
print("max G*4theta/zeta2 over m<=18, theta in 2^[-15,10]:",worstG)
