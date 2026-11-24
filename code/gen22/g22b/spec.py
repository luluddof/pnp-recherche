import numpy as np, sympy as sp
from anat import load, Rmat, files
def PiTheta(N, G, Gp):
    R, C = Rmat(N, G, Gp); R4 = R.reshape((N,)*4)
    # Pi[(a,b),(c,d)] = C(x=a,y=d,u=b,v=c);  Theta[(a,b),(c,d)] = R(x=a,y=b,u=d,v=c)
    Pi = C.transpose(0, 2, 3, 1).reshape(N*N, N*N); Th = R4.transpose(0, 1, 3, 2).reshape(N*N, N*N)
    return Pi, Th
if __name__ == "__main__":
    for f in files:
        N, G, Gp = load(f); Pi, Th = PiTheta(N, G, Gp)
        print("=====", f.split('/')[-1], " Pi symmetric:", (Pi == Pi.T).all(), " Theta symmetric:", (Th == Th.T).all())
        for nm, M in (("Pi", Pi), ("Theta", Th)):
            Ms = sp.Matrix(M.tolist()); cp = sp.factor(Ms.charpoly().as_expr())
            print("  charpoly", nm, ":", cp, "| rank", Ms.rank(), "| rank of M^2:", (Ms*Ms).rank(), "| trace", Ms.trace())
