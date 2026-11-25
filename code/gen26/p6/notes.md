# G26-P6 notes (agent, unreviewed)
02:25 UTC. Read g18d notes, smith_cauchy.py, flat_packet.py, paper (D.7) lines 3361-3400.
Kernel (D.7): W_xy(T) = C(T)/(C'(0)(T-(y-x))), S={-L..L-1}, C(T)=prod_{q in S}(T-q), L=2^a, alpha=2^{a+1}.
Crossing block at node t (v2(t)=i<a), T=t+alpha*zeta: E = (W_xy(T)) for x,y in B_b={x in [0,L): bit_i(x)=b}.
PROOF SKETCH (found 02:25): E = c * Cauchy(1/(T+x-y)), c=C(T)/C'(0), v2(c)=a+1+v2(zeta) (Lemma D.4 counting).
Minor(R,C) valuation = k v2(c) + V(R)+V(C) - sum_{R x C} min(i, v2(x-y)); decompose by residue rho mod 2^i,
cross terms form PSD quadratic form Q(e), e_rho=r_rho-c_rho; within-class: superfactorial bound.
=> d_k = k(v2(c)-i) + 2*(sum of k smallest of {j+v2(j!)}, mult 2^i). QED (to be written up).
02:30 verify_smith.py: 856 cases (a=1..7, all/many t with v2(t)<a, zeta in {1,2,3,-5,12,6}, b=0,1) all match formula;
  brute force over all minors (a=3,4) confirms min minor valuation = prefix formula. Output verify_smith_out.txt.
Paper D.4.5 read: crossing = packets with differing symbols, rows x/cols y with equal bit i; tiles y0 mod 2^{i+1};
  paper weight w_Z=b+1 (true a-i+1+b). Row H modes = Smith modes of exactly my E (submatrix for real ports: weights only grow).
CAVEAT found: Smith sum G(th)=sum_k f(th s_k) = zeta(2)/(4 th) + O(log^2(1/th)) (s_k=4k+c0-2popcount(k)), so the
  "n -> n/4" exponent halving is asymptotic; need ratio rho(D,y) at P3's parameters. Computing (rho.py).
02:40 rho.py (rho_out.txt): exact Smith GF exponent / (D.15) exponent < 0.5 at all tested (D,y).
  PROVED: G(th)=sum_{k<D/2} f(th s_k) <= zeta(2)/(4 th) for all th>0 (D=2^m). Proof: N(y)=#{k: s_k<=y};
  A(Y)=int_0^Y (y/4-N) dy >= 0 since min over Y = sum_{k<K} s_k - 2K^2 = K(m-1) - 2P(K) >= 0, P(K)=sum_{k<K} popcount(k)
  <= (K/2)log2 K <= K(m-1)/2 (K<=2^{m-1}); then zeta2/(4th)-G = int A(y) w''(y) dy >= 0 (w=f(th y) convex, boundary 0).
  => (D.15) with n -> n/4 EXACTLY (no o(1)): row H count bound rigorous.
C1 analysis: upper bound R_p(flat) <= Smith-monomial count: easy [P]. Lower bound for p <= 2 s_0 via submatrix
  [[1,0,0],[0,E,0],[0,0,E'^T]] after (1,z) basis change on Y: R_p >= 1+R_p(E)+R_p(E') [P] => C1 exact for p<=2(a-i+1).
  C1 as stated FAILS far from saturation: a=3,i=1,p=10: 35 vs 55 (rank cap 256); a=3,i=2,p=6: 37 vs 45.
  Observation: for i=a-1 (one column per tile) data = SQUARE-FREE Smith-monomial count exactly (37,93,163,219,255).
  Testing: F ~ M (x) M' with M[S,X'] = prod_{x in S}(E X')_x ? (c1_test.py)
02:30 c1_test.py (c1_test_out.txt): R_p(F) = R_p(M1 (x) M2) EXACTLY for all (a,i) in {(2,0),(2,1),(3,0),(3,1),(3,2)}, all p,
  where M[S,X'] = prod_{x in S}(E X')_x is the one-direction "Boolean product matrix" (2^{L/2} x 2^{L/2}).
  Smith(M_{3,1}) = {0,3,5,7}^{(x)2} (= tile profile of a'=2), Smith(M_{3,2}) = {0,2}^{(x)4}: tensor over tiles.
  Smith(M_{2,0}) = {0,3,5,7}; Smith(M_{3,0}) = {0,4,6,9,10,10,12,13,13,14,15,17,19,20,20,22}
  vs full monomials in modes (4,6,10,12): 0,4,6,8,10,10,12,12,12,...; Boolean squares shifted (+1).

=== THEOREM 1 (Smith profile of the crossing block) — PROVED 02:25-02:30 ===
Setting: a>=1, L=2^a, S={-L..L-1}, C(T)=prod_{q in S}(T-q), alpha=2^{a+1}; O = ring of integers of an unramified
ext. of Q_2 (or Z_(2)), zeta in O\{0}; t in S, i=v2(t) <= a-1; T = t+alpha*zeta; b in {0,1};
B_b = {x in [0,L): bit_i(x)=b}; E = (W_xy(T))_{x,y in B_b}, W_xy(T) = C(T)/(C'(0)(T-(y-x)))  [(D.7)], an (L/2)x(L/2)
matrix over O. Claim: Smith form of E over O = diag(2^{s}) with multiset
  { s_k := (a-i+1) + v2(zeta) + 2k + 2 v2(k!) : 0 <= k < 2^{a-i-1} }, each value repeated 2^i times.
Proof.
(1) v2(c)=a+1+v2(zeta), c := C(T)/C'(0): C(T) = alpha*zeta*prod_{q!=t}(t-q+alpha*zeta), and v2(t-q)<=a<a+1 for q!=t in S
    (complete residue system mod 2^{a+1}), so v2 C(T) = a+1+v2(zeta)+v2 C'(t); v2 C'(t)=v2 C'(0) (Lemma D.4's count).
(2) For x,y in B_b: v2(T+x-y) = min(i, v2(x-y)) (v2(x-y)=inf -> i).  Write x = rho_x + b2^i + 2^{i+1}x', rho in [0,2^i).
    If rho_x != rho_y: v2(x-y)=v2(rho_x-rho_y)<i and v2(t+x-y)=v2(x-y); alpha*zeta has valuation >= a+1. If rho_x=rho_y:
    t+x-y = 2^i(odd + 2(x'-y')) so valuation i (and then v2(x-y) >= i+1).  In particular T+x-y != 0.
(3) E = c * (1/(X_x - Y_y)) with X_x = T+x, Y_y = y: a Cauchy matrix. For R,C subsets of B_b, |R|=|C|=k,
    v2 det E[R,C] = k v2(c) + V(R) + V(C) - sum_{x in R, y in C} min(i, v2(x-y)),  V(R) := sum_{pairs in R} v2(x-x~).
(4) Group by residue rho: r_rho=|R_rho|, c_rho=|C_rho|, e_rho=r_rho-c_rho. Pairs in different classes: valuation
    m(rho,rho')=v2(rho-rho'); same class: i+1+v2(x'-x~'). Algebra gives
      f(R,C) := V(R)+V(C)-sum min(...) = Q(e) - i k + sum_rho [ g(R_rho) + g(C_rho) ],
      Q(e) = (1/2) sum_{rho,rho'} m(rho,rho') e_rho e_rho'  with m(rho,rho)=i,
      g(X) = C(|X|,2) + sum_{pairs in X} v2(x'-x~').
    m(rho,rho') = sum_{l=1..i} [rho == rho' mod 2^l], so Q(e) = (1/2) sum_l sum_{classes mod 2^l} (sum e)^2 >= 0.
    Superfactorial lemma: for n distinct integers z, prod_{j<l}(z_l-z_j) / prod_{j<n} j! = det[binom(z_l, j)] in Z,
    so g(X) >= G(n) := sum_{j<n} h(j), h(j) = j + v2(j!), with equality for consecutive x'.
(5) h strictly increasing => min over (r_rho) with sum r_rho = k of sum_rho G(r_rho) = H_k := sum of the k smallest
    elements of the multiset {h(j): 0<=j<2^{a-i-1}} x 2^i. Hence min f = 2H_k - ik, attained at R=C= union of
    initial segments (e=0). So d_k(E) = min_{R,C} v2 det = k(v2(c)-i) + 2H_k.
(6) Over the DVR O, d_k = s_(0)+...+s_(k-1) (sorted invariants); so the invariants are v2(c)-i+2h(j), j<2^{a-i-1},
    each 2^i times, = (a-i+1)+v2(zeta)+2j+2v2(j!).  QED.  (zeta=0: E=0.)
Numerics: verify_smith.py, 856 cases a<=7 incl. negative t, t=3*2^i etc., zeta in {1,2,3,-5,6,12}; brute-force minors.

=== THEOREM 2 (count) — PROVED ===
For D=2^m, s_k = m+1+2k+2v2(k!) (k<D/2) [this is s_k for zeta unit, m=a-i], f(x) = -ln(1-e^{-x}):
  sum_{k<D/2} f(theta s_k) <= zeta(2)/(4 theta)  for every theta>0.
Proof in 02:40 entry. Consequence: the Smith-mode system at node t has <= 6M/D modes with each index k (as (D.13)),
real-port submatrices have termwise larger invariants (R_w(sub) <= R_w), zeta in {0} u O^x, so
#{beta: sum s beta < p_t} <= exp(theta p_t + (6M/D) zeta(2)/(4 theta)) = exp(pi sqrt(2 (6M/(4D)) p_t/3)):
exactly (D.15) with n_t -> n_t/4, i.e. the exponent is halved, no o(1).
02:40 PROPOSITION 3 [P] (any kernel W): F ~ M (x) M'^T over O. Proof: change Y-basis (l0,l1)->(1,z) (unimodular), entry
  = prod_{x in S}(A_x(X_g)+u_x(X_g')) prod_{x' in S'}(A'_x'(X_g')+u'_x'(X_g)); expand both products over S1 c S, S1' c S':
  F' = P_g N P_h^T, P_g[(S,X_g),(S1,X~)] = [X~=X_g][S1 c S] prod_{S\S1} A_x(X_g) (unitriangular in inclusion order),
  N[(S1,X~g),(S1',X~g')] = M[S1,X~g'] M'[S1',X~g] = Kronecker M (x) M'^T. Hence Smith(F) = {sigma+sigma'}.
  Explains c1_test exact equality. C1 reduces to Smith profile of M[S,X'] = prod_{x in S}(E X')_x.

=== C1 STATUS (02:45) ===
Definitions. One packet (two symbols differ in coordinate i), node T as in Thm 1. Groups g = {rows bit_i=0, cols bit_i=1},
g' = {rows bit_i=1, cols bit_i=0}. Packet tensor K(Y,X) = prod_x l_{Y_x}((W(T)X)_x). Flattening F: rows (Y_g,X_g), cols (Y_g',X_g')
(2^L x 2^L). R_p(F) = #{Smith invariants of F with valuation < p}.  E (g'->g crossing, rows bit 0, cols bit 0) and E'
(rows bit 1, cols bit 1) are the two crossing blocks of Thm 1.  M(E)[S,X'] = prod_{x in S}(E X')_x, S c rows, X' in {0,1}^cols.
(a) [P] Prop. 3: F ~ M(E) (x) M(E')^T over O, for ANY kernel W (only affineness of l_Y used). Hence
    Smith(F) = {sigma+sigma' : sigma in Smith M(E), sigma' in Smith M(E')}. Confirmed numerically (c1_test, 5 cases, all p).
(b) [P] Upper bound, all p: R_p(F) <= #{(beta,beta'): sum s beta + sum s beta' < p} (Smith-mode monomials). Proof: expand
    prod l(local + P Delta z) — every Z^beta coefficient divisible by 2^{w(beta)}, each term rank one across the cut.
    (This is the direction the ALGORITHM uses; it is Thm 1 + expansion, no C1 needed.)
(c) [P] Lower bound, linear range: rows (S=0,X_g=0),(S={x},X_g=0),(S=0,X_g=e_y); cols analogous give the submatrix
    [[1,0,0],[0,E,0],[0,0,E'^T]] (after one column reduction) => R_p(F) >= 1 + R_p(E) + R_p(E'). With (b): C1 EXACT for p <= 2 s_0
    = 2(a-i+1+v2 zeta).
(d) [P] Level i = a-1 (D=2): E = 2G, G in GL(O) (Thm 1: all invariants 2). M(E) = diag(2^{|S|}) N, N[S,X'] = prod_{S}(GX')_x;
    mod 2, N is block-triangular in the Moebius basis with diagonal blocks = compounds of G-bar (per = det in char 2),
    so N is unimodular: Smith M(E) = {2|S|}. Hence R_p(F) = #{(S,S'): 2(|S|+|S'|) < p} = SQUARE-FREE Smith-monomial count.
(e) [N, REFUTATION] C1 as literally stated (= full Smith-monomial count until saturation) is FALSE, already at L=8 (a=3,i=0,
    zeta=1, the case cited as support): R_9(F) = 6 < 8 (Z0^2, Z0'^2 sit at weight 2s_0+1 = 9, not 8; g18d only tested even p);
    a=3,i=1,p=10: 35 < 55; a=3,i=2,p=6: 37 < 45, all far below the rank cap 256 (c1_odd_out.txt, flat_out.txt).
    Cause: Boolean columns (X^2 = X): z_k^2 == z_k mod 2, so repeated modes lose weight structure.
(f) [N] Corrected conjecture C1' (sandwich): termwise (sorted)
       full-monomial weights  <=  Smith(M(E))  <=  square-free subset sums of Smith(E).
    Verified for (a,i) in {(1,0),(2,0),(3,0),(4,0),(3,1),(4,1),(4,2),(4,3)} (m_profile_out.txt; a=4 matrices 256x256).
    Left inequality = (b) [P]. Right inequality proved for i=a-1 (d) and in range p<=2s_0 (c); open in general.
    Tile-tensor guess Smith M_{a,i} = (Smith M_{a-i,0})^{(x)2^i}: true for (3,1),(4,2),(4,3), FALSE for (4,1).
(g) Consequence if C1'-right holds: by Prop 3 + Lemma T(iii) the KF lower bound is >= #{0/1 tuples of Smith modes with weight < p_t}
    = exp(pi sqrt(n p/12)(1-o(1))) against the algorithm's exp(pi sqrt(n p/6)) (Thm 2): mode enumeration optimal in KF up to
    a factor sqrt 2 in the exponent constant (not up to 1 as C1 claimed).
02:42 generalE.py (generalE_out.txt): for RANDOM E = P diag(2^s) Q (P,Q in GL_n(Z_2), n<=5, 300 trials) the left half
  (Smith M >= full monomials) never fails (it is proved), but the RIGHT half fails in 93/300 trials, e.g. n=2, s=(0,5):
  Smith M = (0,0,1,6) vs square-free (0,0,5,5): R_6 = 3 < 4. So C1'-right is NOT a general Smith-form fact; it must use the
  Cauchy structure (Lambda = span of square-free products is not invariant under row operations on E, which is why the
  Smith factor P matters). Plausible mechanism for the Cauchy block: Mahler modes z_k = sum_y binom(y,k) X_y; products
  binom(y,j)binom(y,k) expand in binom(y,l) with l >= max(j,k) (weights go UP under the Boolean collapse). Not proved.

=== FINAL STATUS / ROW H ===
* Thm 1 [P]: exact Smith form of every crossing block of (D.7) (all a, all t in S with v2(t)<a, all zeta in O, b=0,1):
  invariants (a-i+1)+v2(zeta)+2k+2v2(k!), k<2^{a-i-1}, multiplicity 2^i.  (G18d's formula, now a theorem; 856 extra checks.)
* Thm 2 [P]: the Smith-mode tuple count obeys (D.15) with n_t replaced by n_t/4 exactly (G(th) <= zeta(2)/(4th), all th).
* => P3 row H: the running-time bound of the Smith-mode variant is rigorous (relative to App. D + P3's rows A-G, which are
  P3's own unreviewed re-proofs): c = 37.5 in E <= c k m / sqrt(log2 s) (P3's closed form 75 with n -> n/4). Thresholds
  2.25e5 / 1.5e6 are P3's numbers (not recomputed here). Row H no longer depends on C1: C1 is only the optimality claim.
* C1: literal form refuted (L=8, p=9); reduced exactly to the one-direction Boolean product matrix (Prop 3); proved for
  p <= 2 s_0 and for the top level i=a-1; corrected sandwich conjecture C1' verified a<=4; its lower-bound half is false for
  generic E, so any proof must use the Cauchy/Mahler structure.  If C1'-right holds: KF optimality up to sqrt 2 in the constant.
Row I (kappa=2) untouched.  No processes left running.  Nothing written to /home/user/pnp-recherche.
