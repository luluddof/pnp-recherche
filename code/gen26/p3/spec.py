from load import *
for f in FILES:
    N,A,Bm=load(f); S,K=mats(N,A,Bm); E=(S-K).astype(float)
    s=np.linalg.svd(E,compute_uv=False); s=s[s>1e-9]
    tr=np.trace(S.T@E)  # = N^2 - tr(Swap K)
    print(f.split('/')[-1],'p',len(s),'q',len(A),'trSE',tr,'nuc^2/F^2 %.2f'%(s.sum()**2/(s**2).sum()),'F^4/S4^4 %.2f'%((s**2).sum()**2/(s**4).sum()),
      'tr^2/F^2 %.2f'%(tr**2/(s**2).sum()), 'sv',np.round(s,2))
