# Exact-versus-border gap for the strengthening (c):  m = n = 2, Q = 2.
# Exact: osr(B) < 4  ==> rank(A) >= Q = 2.   Border: P = lim (A_eps + B_eps) with rank(A_eps) = 1 and osr(B_eps) <= 3.
import sympy as sp, itertools
e = sp.symbols('eps')
m = n = Q = 2
def idx_out(g, ee, w): return (g*m + ee)*Q + w      # G (x) E (x) W
def idx_in(ee, g, w): return (ee*n + g)*Q + w       # E (x) G (x) W
Dm = m*n*Q
P = sp.zeros(Dm, Dm)
for ee, g, w in itertools.product(range(m), range(n), range(Q)): P[idx_out(g, ee, w), idx_in(ee, g, w)] = 1
# c_eps = alpha_1 beta_1^T + eps alpha_2 beta_2^T  (functional on site-0 operators |g><e| : c(|g><e|) = c[g,e])
c = sp.Matrix([[1, 0], [0, e]])
# x = g_1 (x) u_1 + g_2 (x) u_2/eps,  u_w = e_1 (x) w in E (x) W ;  y = e_1 (x) v_1 + e_2 (x) v_2, v_w = g_1 (x) w in G (x) W
x = sp.zeros(Dm, 1); y = sp.zeros(Dm, 1)
x[idx_out(0, 0, 0)] = 1; x[idx_out(1, 0, 1)] = 1/e
y[idx_in(0, 0, 0)] = 1; y[idx_in(1, 0, 1)] = 1
A = x * y.T
def contract(T, c):     # (c (x) 1) T : rest_in (g',w') -> rest_out (e,w)
    R = sp.zeros(m*Q, n*Q)
    for g, ee in itertools.product(range(n), range(m)):          # site-0 operator |g><e'|: out g, in e'
        if c[g, ee] == 0: continue
        for eo, w, gi, wi in itertools.product(range(m), range(Q), range(n), range(Q)):
            R[eo*Q + w, gi*Q + wi] += c[g, ee] * T[idx_out(g, eo, w), idx_in(ee, gi, wi)]
    return R
Dl = contract(P - A, c).applyfunc(sp.simplify)
print("(c_eps (x) 1)(P - A_eps) =", Dl.tolist())
# N with c(N) = 1: N = |g=0><e=0|.  B = P - A - N (x) Dl
B = (P - A).copy()
for eo, w, gi, wi in itertools.product(range(m), range(Q), range(n), range(Q)):
    B[idx_out(0, eo, w), idx_in(0, gi, wi)] -= Dl[eo*Q + w, gi*Q + wi]
def realign(T):
    R = sp.zeros(n*m, m*Q*n*Q)
    for g, ee, eo, w, gi, wi in itertools.product(range(n), range(m), range(m), range(Q), range(n), range(Q)):
        R[g*m + ee, ((eo*Q + w)*n + gi)*Q + wi] = T[idx_out(g, eo, w), idx_in(ee, gi, wi)]
    return R
err = (A + B - P).applyfunc(sp.simplify)
print("rank A_eps =", A.rank(), "; osr(B_eps) =", realign(B).rank(), "(over Q(eps)); mn =", m*n, "; Q =", Q)
print("A_eps + B_eps - P has entries:", sorted(set(map(str, err))), " -> 0 as eps -> 0")
assert A.rank() == 1 and realign(B).rank() <= 3 and all(sp.limit(v, e, 0) == 0 for v in err)
print("BORDER decomposition with r = 1 < Q and b = 3 < mn EXISTS; exact ones do not (T7).  Border inequalities (a),(b) hold: ",
      all(k*Q*n <= 1 + k*k*Q*3 and k*n*Q <= 3*Q + k*k*1 for k in range(m+1)))
