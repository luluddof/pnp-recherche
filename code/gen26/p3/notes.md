# G26-P3 notes (monomial-C model MC, beyond-line attempts)
Tags [P] proof (unreviewed), [E] exact, [N] numeric, [H] heuristic.
02:05 start. load.py: all 5 exact optima on file ARE in MC (A_l, B_l signed partial perms; rank(Swap-K) re-verified exactly):
 (5,3)@4, (7,4)@5, (6,6)@5, (7,7)@6, (12,3)@6.  [E]
spec.py: singular values of E = Swap-K are nearly flat on all optima (nuc^2/F^2 within 10% of p).

## Key idea found (02:45): symmetrization + compression
Theta = K Swap = sum_l P_l (twisted signed partial maps), V = ker(1-Theta), m = dim V = N^2 - p.
Theta_s = (Theta + Theta^T)/2.  For v,w in V: <v,Theta_s w> = <v,w>  => compression of Theta_s to V is I
=> Theta_s has >= m eigenvalues >= 1 (Courant-Fischer).   Odd words in {P_l, P_l^T} are twisted => |tr| <= N
=> |tr Theta_s^j| <= q^j N for odd j.  ||Theta_s|| <= q.
Odd polynomial f ~ sign on [1,q]: (1-eps)m - (1+eps)(n-m) <= tr f(Theta_s) <= N * sum|f_j| q^j  => m <= n/2 + ...
 => THEOREM A: p >= (1-eps)N^2/2 - N L(q,eps), L <= 2(k+1)(1+sqrt2)^k, k = O(q log(1/eps)).
    => p >= (1/2 - o(1)) N^2 whenever q = o(log N).   NOT an exponent gain (needs q <~ log N).
 THEOREM B (monotone Fejer-type g): p >= c N^2 ln N / q uniformly (log gain, like MC2 but all q <= N^{1/2-}).
Obstruction [H]: any argument using only odd traces + norm needs degree k >~ a q for gain factor a, and coefficient
 blow-up e^{ck} <= N forces a <~ log N / q... i.e. log gains only.

## [E] thetas.py on the 5 exact optima (02:35)
 #eig(Theta_s) >= 1 equals m = N^2 - p exactly on 4 of 5 (25 vs 24 on (12,3)@6): the eigenvalue count is TIGHT there.
 ||Theta_s|| = 2 on all optima (q up to 7) -- far below q. Odd traces tiny vs q^j N.
 Note Theta_s = I - H, H = (E S + (E S)^T)/2, rank H <= 2p, so n - 2p eigenvalues are exactly 1.

## 02:35 results
[E] cert.py (exact rational certificates, Sturm root counts): odd g with |g|<=1 on [0,1], g>=a on [1/q,1]
   (03:00: first version checked only g<=1; constraint g>=-1 added to LP + verifier, ALL certificates re-verified, same L):
   q=3, a=0.90, deg 41: L=18.01  => p >= 0.4737 N^2 - 9.48 N      (MC, any N, q<=3)
   q=3, a=0.95, deg 61: L=42.68  => p >= 0.4872 N^2 - 21.89 N
   q=3, a=0.99, deg 81: L=420.8  => p >= 0.4975 N^2 - 211.4 N
   q=2, a=0.99, deg 61: L=17.68  => p >= 0.4975 N^2 - 8.89 N   (sanity: exact value N(N-1)/2 known, Thm 1 g24g)
[N] lpcert.py (LP, deg<=41, float): best odd-trace bound q=10: N=1e6 0.4246 N^2, N=1e9 0.472 N^2;
   q=30: N=1e3 2.95 N^2/q, N=1e6 6.28 N^2/q   (grows like ln N / q: log gain, as predicted).
[E] grpthm.py / grpex.py: Theorem G (group-regular A_l=+-L(g), B_l=+-L(h)): p >= (N^2 - sum d^3)/2, never violated
   (Z6,S3,Q8,D4,Z2^3,A4,D6, q<=6 random; S3,Z6 exhaustive q<=3: min p = 15 = N(N-1)/2 >= bound 13/15).
[E] dsum.py: direct sums of the optima: p ~ N^2/2 - O(N) (e.g. (12,3)@6 x3: N=18, p=156 > 153).
[E] more certificates: q=4,a=.9: L=91.6 (p>=.4737N^2-48.2N); q=5: L=471 (-248N); q=6: L=2603 (-1370N);
    q=6,a=.5,deg41: L=10.16 (p >= N^2/3 - 6.8N); q=10,a=.5,deg61: L=93.4 (p >= N^2/3 - 62.2N).
[N] bn_check.py: the odd-trace bound |tr Theta_s^j| <= q^j N holds (and is attained, degenerate G'=G^{-1}) for
    ORTHOGONAL factors -> Theorems A/B hold in the bounded-norm model BN (see below).

# ===================== WRITE-UP (02:55) =====================
Setting. Swap = R + C, C = sum_{l<=q} G_l (x) G'_l, K[x,y,u,v] = sum_l G_l[x,u] G'_l[y,v], rank R = p (rows (x,y)).
Theta := K S (S = Swap as N^2 x N^2 permutation). S^2 = 1 so p = rank(S - K) = rank(1 - Theta). n := N^2.
MC: G_l, G'_l signed partial permutation matrices.  BN ("bounded norm"): ||G_l||_op ||G'_l||_op summed: Q := sum_l ||G_l|| ||G'_l||.
MC is contained in BN with Q <= q.  Field: R (or C, with transposes -> adjoints).

LEMMA 1 (compression). V := ker(1 - Theta), m := dim V = n - p, Theta_s := (Theta + Theta^T)/2. Then <v,Theta_s w> = <v,w> for
 v,w in V, hence (Courant-Fischer) Theta_s has at least m eigenvalues >= 1.
 Proof: <v,Theta_s w> = (<v,Theta w> + <Theta v,w>)/2 = <v,w>.  lambda_m(Theta_s) = max_{dim U=m} min_{u in U,|u|=1} <u,Theta_s u> >= 1.
LEMMA 2 (odd traces). ||Theta_s|| <= Q and |tr Theta_s^j| <= Q^j N for every odd j.
 Proof: S (X(x)Y) S = Y(x)X, so Theta = sum_l (G_l(x)G'_l) S and Theta^T = sum_l (G'_l^T (x) G_l^T) S. A product of an odd number of
 factors (X_i (x) Y_i) S equals (M (x) M') S with M = X_1 Y_2 X_3 ..., M' = Y_1 X_2 Y_3 ... (alternating), ||M|| ||M'|| <=
 prod ||X_i|| ||Y_i||.  tr((M(x)M')S) = tr(M M') so |tr| <= N ||M|| ||M'||.  Expanding Theta_s^j = 2^{-j}(sum of 2q terms)^j gives
 |tr Theta_s^j| <= 2^{-j} (2Q)^j N.  (In MC, simply: an odd word in the twisted cell maps P_l^{+-} is a twisted signed partial map
 (a,b) -> +-(mu b, nu a), which fixes <= N cells.)  Norm: ||Theta_s|| <= ||K|| <= Q.
LEMMA 3 (counting). If g is an odd real polynomial, |g(y)| <= 1 on [-1,1] and g(y) >= a > 0 on [1/Q,1], and L = sum_j |g_j|
 (monomial coefficients), then  p >= (a n - L N)/(1 + a).
 Proof: eigenvalues of Theta_s/Q lie in [-1,1]; at least m of them are >= 1/Q (Lemma 1) where g >= a; at all others g >= -1. So tr g(Theta_s/Q) >= a m - (n - m). And tr g(Theta_s/Q) = sum_j g_j Q^{-j} tr Theta_s^j <= L N
 (Lemma 2). Hence m (1+a) <= n + L N, p = n - m >= (a n - L N)/(1+a).  QED
THEOREM A [P]. For every Q and a < 1 there is an explicit L(Q,a) with  p >= (a N^2 - L(Q,a) N)/(1+a)  in BN (hence MC with q = Q).
 In particular p >= (1/2 - eps) N^2 - O_{q,eps}(N) for fixed q, and p >= (1/2 - o(1)) N^2 whenever q = o(log N).
 Explicit (cert.py, exact rational Sturm verification):
   q=2: a=.99 L=17.68 | q=3: a=.90 L=18.01, a=.95 L=42.68, a=.99 L=420.8 | q=4: a=.9 L=91.6 | q=5: a=.9 L=471 | q=6: a=.9 L=2603,
   a=.5 L=10.16 | q=10: a=.5 L=93.4.     E.g. MC, q <= 3:  p >= 0.4737 N^2 - 9.48 N  and  p >= 0.4975 N^2 - 211.5 N.
 Existence for all (Q,a) and the o(log N) statement: take g = Gilyen-Su-Low-Wiebe sign approximant (odd, |g|<=1 on [-1,1],
 |g - sgn| <= eps on |y| >= 1/Q, degree k = O(Q log(1/eps))) [literature, Lemma 25 of GSLW STOC'19]; or self-contained:
 g(sin phi) = G(phi)/G(pi/2), G(phi) = int_0^phi cos(psi) F_K(2psi)^r dpsi (F_K Fejer kernel): odd-frequency, monotone, |g|<=1,
 1 - g(y) <= u^{2r-1} pi/(2r-1) for y >= delta with u = pi^2/(4(K+1)delta); degree 2Kr+1; Chebyshev coefficients sum <=
 (8/pi)(1+ln k)(pi^2/4)^r, and sum of monomial coefficients of T_t is <= (1+sqrt2)^t.  So L <= exp(0.89 k + O(r + log k)),
 k = O(Q log(1/eps)); L N = o(N^2) iff Q log(1/eps) <= c log N.
THEOREM B [P]. In BN with ln N <~ Q:  p >= c N^2 ln N / Q  (c > 0 absolute).  Proof: Lemma 3 with r = 1 (Fejer), K+1 ~ 0.28 ln N
 (so (1+sqrt2)^{2K+1} <= sqrt N), a = g(1/Q) >= 16(K+1)/(pi^4 Q) (valid when K+1 <= pi Q/2), L <= 2pi(1+ln(2K+2))(1+sqrt2)^{2K+1}.
 Gives p >= 0.04 N^2 ln N/Q - O(N^{3/2} loglog N).  LP (lpcert.py, deg<=41) gives much better constants: q=30: 2.95 N^2/q at N=1e3,
 6.3 N^2/q at N=1e6.  Still only a LOG gain (pq >= c N^2 log N), uniformly for q <= N^{1/2-o(1)}.
REMARK (obstruction for Lemmas 1-3, [H]): the certificate must have g(1/Q) >= a with |g| <= 1 near 0, which forces degree >~ aQ and
 coefficient sum L >= e^{c a Q} (Chebyshev/Bernstein: e.g. sin(aQy) has a and L = sinh(aQ)); L N < a N^2 forces aQ <~ log N.
 So any argument using only {eigenvalue count of Theta_s, norm, odd traces} gives at most p ~ N^2 log N / Q: no exponent.
 Since Lemmas 1-3 hold in BN, an exponent gain in MC MUST use monomiality beyond the norm bound.

THEOREM G [P] (term-wise group symmetry => no gain). Let alpha be a representation of a finite group Gamma on C^N,
 alpha = sum_rho rho^{(+) m_rho} (d_rho = dim rho). If every G_l and G'_l commutes with alpha(Gamma) (after the symmetry
 (G_l,G'_l) -> (X G_l Y^{-1}, Y G'_l X^{-1}), X,Y invertible, which preserves rank R and q), then for ANY q:
        p >= ( N^2 - sum_rho d_rho m_rho^2 ) / 2.
 Cases: Gamma acting regularly (G_l, G'_l in the left-regular algebra, e.g. +-left multiplications in a group of order N):
 m = d, p >= (N^2 - sum d^3)/2 >= (N^2 - N^{3/2})/2.  Free action with orbits of size k (abelian, d=1, m=N/k): (N^2 - N^2/k)/2.
 Commutative families (simultaneously diagonalisable; alpha = the torus/abelian group they generate, if multiplicity-free): (N^2-N)/2.
 Proof. Put beta(x,y) = alpha(x) (x) alpha(y) on C^N (x) C^N. K commutes with beta(x,y); S beta(x,y) = beta(y,x) S. So
 Theta = K S satisfies Theta beta(x,y) = beta(y,x) Theta, i.e. Theta maps the beta-isotypic component I(rho,pi) (dim d_rho d_pi m_rho m_pi)
 into I(pi,rho).  (i) rho != pi: on I(rho,pi) (+) I(pi,rho), 1 - Theta = [[1, -T12],[-T21, 1]], rank >= dim I(rho,pi).
 (ii) rho = pi: I(rho,rho) = (V (x) V) (x) (M (x) M), V = rep space, M = C^{m_rho}; let F swap the two V factors; F^{-1} Theta commutes with
 beta on I(rho,rho), V(x)V is irreducible for Gamma x Gamma, so by Schur F^{-1}Theta = 1_{V(x)V} (x) Theta', i.e. Theta = F (x) Theta'.
 F = +1 on Sym^2 V (dim d(d+1)/2), -1 on Alt^2 V (dim d(d-1)/2), so rank(1-Theta)|_{I(rho,rho)} = d(d+1)/2 rank(1-Theta') +
 d(d-1)/2 rank(1+Theta') >= d(d-1)/2 * m^2 (since rank(1-T)+rank(1+T) >= dim).  Sum: p >= sum_{rho<pi} d_rho d_pi m_rho m_pi +
 sum_rho d_rho(d_rho-1) m_rho^2/2 = [(sum d m)^2 - sum d^2 m^2]/2 + [sum d^2 m^2 - sum d m^2]/2 = (N^2 - sum d m^2)/2.  QED
 Checks [E]: grpthm.py/grpex.py (regular reps of Z6,S3,Q8,D4,Z2^3,A4,D6; random q<=6; S3,Z6 exhaustive q<=3), grp2.py
 (free Z_k actions, N<=16): never violated.  NB the Q8-invariant ST (8,8)@8 array (p = N) is NOT term-wise symmetric (Q8
 permutes the terms), so it does not contradict G: symmetric-but-term-permuting constructions escape G.

THEOREM C [P] (Schatten-4 / relation energy). With X = 1 - Theta (rank p), tr X = n - tr Theta >= n - qN, Hoelder:
 n - qN <= ||X||_* <= p^{3/4} ||X||_{S4}, and
   ||X||_{S4}^4 = n + 2 tr Theta^2 + 4||Theta||_F^2 + ||Theta^T Theta||_F^2 - 4 tr Theta - 4 tr(Theta Theta^T Theta).
 In MC: ||Theta||_F^2 = sum_{l,m} <A_l,A_m><B_l,B_m>, tr Theta^2 = sum_{l,m} tr(A_l B_m) tr(B_l A_m),
 ||Theta^T Theta||_F^2 = sum_{l,m,l',m'} tr(A_l^T A_m A_m'^T A_l') tr(B_l^T B_m B_m'^T B_l') =: n * Energy; odd terms <= 4(q+q^3)N.
 The 2q^2 'trivial' quadruples contribute <= 2q^2 n. Hence if all pair-coincidences and the NONTRIVIAL quadruple energy are
 <= q^{3-eta} n, then  p >= c n / q^{1 - eta/3}:  exponent gain eta/3 (e.g. eta=1: p >= c N^2/q^{2/3}).
 (S_{2k} version: if tr((X^T X)^k) <= C_k q^k n, p >= n/(C_k^{1/(2k-1)} q^{k/(2k-1)}) -> ~ n/(2 sqrt q): ST's exponent.)
DICHOTOMY PICTURE (not a theorem): low relation energy -> Thm C exponent gain; exact group structure -> Thm G (~N^2/2);
 the missing piece is a ROBUST inverse theorem: Energy >= q^{3-eta} n  =>  (after X,Y) the families are close, up to rank
 o(N/q) per term, to families commuting with a group action with small sum d m^2.  Abelian case = Lemma T/Fourier; nonabelian
 approximate groups (BGT) are where it is open.  The perturbation step is fine: changing G_l on r rows changes rank R by <= rN.
q=30 a=.2667 deg41 L=1468 => p>=0.2105N^2-1159N ; q=100 a=.08 deg41 L=1116 => p>=0.0741N^2-1034N (verified)

## 03:05 checks
[N] fejer.py: Theorem-B kernel (K=5..20, Q=20..100): |g|<=1, monotone, a = g(1/Q) ~ 0.64(K+1)/Q (claimed >= 0.164(K+1)/Q),
    L ~ e^{0.71(2K+1)} (claimed <= 2pi(1+ln(2K+2))(1+sqrt2)^{2K+1}). Claims hold.
[N] S4 expansion of Thm C verified numerically (731 = 731); tr Theta^2 = sum tr(A_l B_m) tr(B_l A_m) verified.
[E] kronI.py: A (x) I_m families of all optima: N^2/2 - p stays O(N) (q fixed), as Theorem A predicts.
Data fit (not a theorem): all 5 optima satisfy p >= N^2/2 - q N/3 (c = 1/4..1/3).  CONJECTURE (MC strong lemma):
    p >= N^2/2 - O(q N) in MC.  Would give p >= N^2/4 for q <= N/C: a huge exponent gain. Thm A proves the o(N^2) version
    for q = o(log N) only.
RIGOROUS LIMIT OF LEMMA 3 [P, sketch]: if |g| <= 1 on [-1,1], g(0)=0, g(1/Q) >= a, then M := max_{|z|=1}|g| >= e^{c a Q}
 (Cauchy estimate for g' on a circle of radius rho ~ 1/log M around [0,1/Q] + two-constant theorem in the slit disk, where the
 harmonic measure of the outer circle is O(rho)), and L = sum|g_j| >= M. So Lemma 3 gives p >~ a n only if e^{caQ} < aN, i.e.
 a <~ log N / Q:  the method is capped at p ~ N^2 log N / Q.  (Fake spectra by LP duality: [H].)
[N] fejer_r.py: self-contained Theorem-A kernel G' = cos(psi) F_K(2psi)^r: |g|<=1, 1-g(1/Q) far below the proved bound
    pi u^{2r-1}/(1.4(2r-1)), u = pi^2/(4(K+1)delta); L below 1.82(1+ln k)(pi^2/4)^r (1+sqrt2)^k, k = 2Kr+1.
    => with K+1 = ceil(e pi^2 Q/4), r ~ ln(1/eps)/2: k ~ 6.7 Q ln(1/eps), L <= eps^{-5.9Q - O(1)} poly(Q): so
       p >= (1/2 - eps) N^2  as soon as  N >= eps^{-6Q - O(1)}  (no external citation needed).

## FINAL SUMMARY (03:00)
PROVED (all in char 0; MC subset of BN with Q <= q):
 A. p >= (a N^2 - L(Q,a) N)/(1+a) for every a<1; hence p >= (1/2-eps)N^2 - O_{q,eps}(N) (fixed q) and (1/2-o(1))N^2 for
    q = o(log N). Exact certificates: q<=3: p >= 0.4737N^2 - 9.48N, p >= 0.4975N^2 - 211.5N; q<=6: N^2/3 - 6.8N;
    q<=10: N^2/3 - 62.3N; q<=30: 0.2105N^2 - 1159N; q<=100: 0.0741N^2 - 1034N.   Holds in BN: Q = sum ||G_l|| ||G'_l||
    (indeed for any C with projective (op (x)_pi op) norm <= Q: number of terms irrelevant).
 B. p >= c N^2 ln N / q for ln N <~ q (uniform; MC2 needed q <= N^{1/3}/3).  LOG gain only.
 G. term-wise group symmetry: p >= (N^2 - sum d m^2)/2 (regular: >= (N^2 - N^{3/2})/2), any q.
 C. conditional exponent: p >= c n/q^{1-eta/3} if relation energy <= q^{3-eta} n (S4 Hoelder).
NOT PROVED: any exponent gain p >= N^{2-beta+eps} at q = N^beta.  Obstruction: Lemmas 1-3 hold in BN and are provably capped at
 N^2 log N/Q; acyclic/independent-set methods capped at N^2 log q/q on random-like instances [H].  Exponent needs monomial
 input beyond norms: a robust inverse theorem (high relation energy => near term-wise group symmetry), open.
Files: load.py spec.py thetas.py cert.py lpcert.py grpthm.py grpex.py grp2.py dsum.py kronI.py thmC.py bn_check.py fejer.py fejer_r.py
COROLLARY of G [P]: families G_l (x) I_m, G'_l (x) I_m (base size N0, N = N0 m, ANY base, any q, any field of char 0):
 they commute with I_{N0} (x) rho(Gamma) for a finite group with an irreducible m-dim rep (m >= 1: e.g. Gamma = a group with
 an m-dim irrep; for the bound only d = m, multiplicity N0 matter): p >= (N^2 - m N0^2)/2 = N^2/2 - N0 N/2.
 So the A (x) I_m route of leaf S2.3.a2 can give at most a LINEAR (in N) gain below N^2/2, for every base. (Checked kronI.py.)
 [Needs a finite Gamma with an m-dim irrep: exists for every m (e.g. dihedral D_{2m'}... take Gamma = S_{m+1} standard rep).]
Remark: BN at Q = N already allows p = 0 (Swap = N^{-1} sum_P +-P (x) P over real Paulis, N = 2^k), and Q = 1 gives N(N-1)/2;
 Pauli-diagonal partial sums are inefficient (tau(sigma) = transpose sign is bent). Whether BN admits p ~ N^2 polylog/Q is open,
 so "no exponent from BN-only arguments" remains [H] beyond the Lemma-3 cap.
[E] q=3 more: a=.98 deg81 L=143.6 (p>=0.4949N^2-72.5N); a=.995 deg101 L=1435.7 (p>=0.4987N^2-719.6N).
    (1-a, L): (.1,18) (.05,42.7) (.02,143.6) (.01,420.8) (.005,1436): at N=1e4 best gain bound N^2/2-p <= ~1.2e6 = 0.012 N^2.
    So in MC at q=3 the gain below N^2/2 is o(N^2) (theory: <= N^{2-c/q}); bounded-vs-linear (leaf a2) still open.
All my processes finished (other python jobs on the machine belong to other agents).
