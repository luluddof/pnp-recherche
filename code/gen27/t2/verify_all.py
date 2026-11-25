# run G26-P15's independent verifier (pure Python Fraction + Bareiss) on files; also ranks mod two primes (flint nmod).
import sys, json, importlib.util
spec=importlib.util.spec_from_file_location('mv','myverify_lib.py'); mv=importlib.util.module_from_spec(spec); spec.loader.exec_module(mv)
from lib import load, swapm, Kmat, osr, modrank
for fn in sys.argv[1:]:
    b=fn.split('/')[-1].replace('.json','').split('_'); N,p,q=int(b[-3]),int(b[-2]),int(b[-1])
    N2,r1,r2,nt,o=mv.check(fn,p,q)
    _,E,C=load(fn); M=swapm(N)-Kmat(E,C); R=osr(E,C)
    mods=[(modrank(M,pr),modrank(R,pr)) for pr in (1000003,2147483647)]
    ok=(N2==N and r1==r2==p and o==q)
    print(f"({p},{q})@{N}: rank_Q frac={r1} bareiss={r2} terms={nt} osr={o} mod-primes={mods} -> {'CONFIRMED' if ok else 'CHECK'}",flush=True)
