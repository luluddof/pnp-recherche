# enumerate monomials (multisets of d index tuples (x,y,u,v)) with prescribed contents per leg
import itertools
def parts(d, n=3):
    out = []
    def rec(rem, mx, cur):
        if len(cur) == n:
            if rem == 0: out.append(tuple(cur))
            return
        for a in range(min(rem, mx), -1, -1): rec(rem - a, a, cur + [a])
    rec(d, d, []); return out
def monos(cont, N=3):
    """cont = 4 compositions (tuples of length N); returns list of sorted tuples of variable indices (var = x*27+y*9+u*3+v)"""
    d = sum(cont[0]); res = []
    vars_ = list(range(N**4))
    def dec(v): return (v//27, (v//9)%3, (v//3)%3, v%3)
    rem = [list(c) for c in cont]
    def rec(k, start, cur):
        if k == d: res.append(tuple(cur)); return
        for v in range(start, N**4):
            t = dec(v)
            if all(rem[l][t[l]] > 0 for l in range(4)):
                for l in range(4): rem[l][t[l]] -= 1
                cur.append(v); rec(k+1, v, cur); cur.pop()
                for l in range(4): rem[l][t[l]] += 1
    rec(0, 0, []); return res
if __name__ == "__main__":
    for d in range(1, 8):
        P = parts(d)
        for mu in P:
            for nu in P:
                if nu > mu: continue
                n = len(monos((mu, nu, nu, mu)))
                print(d, mu, nu, n, flush=True)
