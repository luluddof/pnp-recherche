# G9a — can the gain rho of Appendix B go below 1/sqrt(2) with a slower-growing graded space?

Verdict: NO in the Appendix-B ("graded window") framework. Two independent reasons, one a
3-line rank count (proved), one structural (locality breaks; proved for the simplest case,
numerically confirmed). alpha = 1/2 (Hardy-Ramanujan) is exactly the optimum of the framework
among power-law growths. Status marks: [P] proved here, [N] numerical, [H] heuristic.

## 1. What Appendix B uses (abstract axioms)

Data, for each lambda (N = 2^lambda):
(G) a graded space W = (+)_{E>=0} W_E with fixed bases, d(E) = dim W_E, L(E) = log2 max_{w<=E} d(w);
(S) weights E0 (for y,u), h_a, h_b, and N basis vectors of W_{E0}: d(E0) >= N;
(O) for x, v in [N] and 0<=i<=H operators a_x^(i)[r] : W_E -> W_{E+r}, b_v[r] : W_E -> W_{E+r};
(L9a) "both orders": Q_s(x,v) := sum_{i=0..H} (-1)^i a_x^(i)[s-h+i] b_v[h-s-i]
                              = sum_{i=0..H} (-1)^i b_v[h-s-i] a_x^(i)[s-h+i]   on W_{E0};
(L9bc+L10) Q_s = 0 outside an interval I of shifts, and sum_{s in I} Q_s(x,v) = d [x=v] Id, d != 0;
(F) |I|, H, E0 = poly(lambda); entries have poly(lambda) bits; tables computable in 2^{O(lambda)}.
In the paper: W = Q[t_1,t_2,...], E0 = h_a = h_b = h, H = 2h, I = [-2h,2h], d = h!, and
(V) "state-field": the operators come from N independent states of weight <= h, with
    H = wt a + wt b (pole order of the full contraction), hence H >= 2 L^{-1}(lambda).
Consequences [P] (this is (54)-(56) of the paper): A_s(x,y,u,v) = <u|Q_s(x,v)|y> factors
row-wise through (+)_i W_{E0+s-h+i}  and column-wise through (+)_i W_{E0-s+h-i}; with x := max
row weight, the max column weight is M - x, where M = 2 E0 + H (= 4h in the paper). So
   c_row(s) + c_col(s) <= L(x) + L(M-x) + 2 log(H+1),   max over s is 2 L(M/2) for concave L.
With L(E) = c E^alpha and M = 4h, L(h) = lambda:  2 L(2h) = 2^{alpha+1} lambda, i.e.
rho = 2^{alpha-1}. The paper's sentence in Appendix H is correct AS A CALCULATION.
Theorem A then gives alpha' = log2(1/rho) = 1 - alpha, time 2^{O(m/(log m)^{1-alpha})}
(K is squared per level, e_l ~ rho^l, 2^l ~ log K_l / (3x0+6)); h ~ lambda^{1/alpha} would only
affect the constant c of Definition 7.

## 2. Obstruction A: rank count  [P]

The crossing [x=v][y=u], as a matrix (x,y)|(u,v), is a permutation matrix of rank N^2
(same for (x,u)|(y,v)). Rank is subadditive, so for ANY pieces
      sum_s rank_row(A_s) >= N^2   and   sum_s rank_col(A_s) >= N^2.
In the window framework rank_row(A_s) <= (H+1) max_{w <= M} d(w). Hence
      L(M) >= 2 lambda - log(|I|(H+1)) = 2 lambda - O(log lambda).            (*)
(In the paper this is tight: the piece s = 2h passes row-wise through weight 4h and
log p(4h) = 2 log p(h) + O(log h): Hardy-Ramanujan growth is exactly borderline.)

Corollary 1 (Appendix-B geometry, M = 4h, L(h) ~ lambda). (*) reads L(4h) >= 2 L(h) - O(log).
For L = c E^alpha this is 4^alpha >= 2, i.e. alpha >= 1/2. So NO graded space with
log dim ~ E^alpha, alpha < 1/2 (squares: alpha = 1/3, k-th powers: 1/(k+1), powers of 2:
log^2, bounded parts: log) can satisfy Lemma 9 + Lemma 10 in this geometry, whatever the
operators are. The value rho = 2^{alpha-1} < 1/sqrt2 is never attained.

Corollary 2 (free geometry: E0, H, M arbitrary, I an interval of shifts, L concave).
Some piece has row weight x1 with L(x1) >= 2 lambda - O(log), some piece has column weight
M - x2 with L(M - x2) >= 2 lambda - O(log); the shifts in between are pieces too, and
x -> L(x) + L(M-x) is concave and symmetric, so some piece costs >= min(4 lambda, 2 L(M/2)) - O(log).
 (a) power law: L(M) >= 2 lambda forces M >= 2^{1/alpha} h1 (h1 = L^{-1}(lambda)), so
     2 L(M/2) >= 2^{2-alpha} lambda, i.e. rho >= 2^{-alpha}; with (V) (M >= 4 h1) also
     rho >= 2^{alpha-1}. Together   rho >= max(2^{alpha-1}, 2^{-alpha}) >= 1/sqrt2,
     with equality only at alpha = 1/2.  Partitions are optimal among power laws.
 (b) more generally, if log L is concave in log E (local growth exponent non-increasing; true
     for H2-type and n-species gas profiles I checked) and (V) holds, then rho >= 1/sqrt2:
     g(u) = log L(e^u) concave, g(u_M) - g(u_1) >= log 2, u_M - u_1 >= 2 log 2  ==>
     g(u_M - log 2) - g(u_1) >= log 2 (1 - log2/(u_M-u_1)) >= (log 2)/2.
 (c) merely concave L with (V): rho >= 2/3 (interpolate (h1,lambda), (M,2 lambda) at M/2, M>=4h1).
     A profile reaching 2/3 would be L(E) = min(lambda E/h, (2/3)lambda + lambda E/(3h)): slope
     drops by 3 exactly at h and the local exponent then INCREASES from 1/3 to 2/3 on [h,4h].
     Non-concave "staircase" profiles (d = N on [h,3h], N^2 above) are not excluded by counting.
     I know no algebraic structure with such a profile and a Lemma-9 identity. [open, speculative]

## 3. Obstruction B: locality breaks on a sub-Fock space  [P for the base case, N in general]

The modes nabla_{+-k}, k in S, do form a Heisenberg Lie algebra acting on W_S = Q[t_k: k in S].
But Lemma 9 is not a statement about the Lie algebra; it is locality of FIELDS. For a = t_n,
a[r] = binom(r-1,n-1) nabla_{-r} for ALL r: the fields contain every mode, W_S is not invariant
(t_n[r] 1 = binom(r-1,n-1) t_r), and W_S is not a vertex subalgebra (the subalgebra generated by
t_n contains all t_r, r >= n, because T = translation is injective in positive weight).
Compressing (dropping the modes outside S): for a = b = t_1, H = 2,
   Q_s(orig) - Q_s(exch) = sum_i (-1)^i binom(2,i) [a[s-1+i], b[1-s-i]] = -(Delta^2 delta)(s),
   delta_r = r * 1[|r| in S].
A second difference of delta vanishes identically iff delta is affine on Z iff S = all of N.
In general the proof of Lemma 9 needs kappa = sum_{k in S} k w1^{-k-1} w2^{k-1} to be the two
expansions of ONE rational function, i.e. (w1-w2)^H must kill the two-sided sequence
k 1_S(|k|) p(k); a two-sided linear-recurrence sequence with H consecutive zeros is zero, and by
Skolem-Mahler-Lech its support has positive density (S = symmetric union of residue classes:
rescaled/twisted Heisenberg, growth exp(c sqrt E) again, alpha = 1/2). [P modulo SML; sketch]
Numerics (heisS.py comp): for S = {1,2}, {1,2,3}, {1,2,4,8}, {1,4}, {1,3}, odd, h = 4..8:
both-orders equality FALSE and sum_s A_s = h! crossing FALSE in every case.
Tried repair [P, negative]: states of the full Heisenberg VOA whose fields skip the modes in a gap
(n, R] exist (a[r] = p(r) nabla_{-r}, p vanishing on the gap) but then wt a >= 2(R-n), H >= 4(R-n),
the windows reach weight > R, and each letter can create one quantum of a mode > R: the
intermediate space gains ~ 2m log R >= 2 lambda bits. No gain.
Other candidates [H]: rank-n Heisenberg (n colours), fermions / Clifford (distinct parts),
Virasoro / W-algebra / lattice modules: all have log dim ~ c sqrt(E) (Cardy growth, c_eff > 0)
or faster (local exponent in [1/2, 1]); T-covariance forces the number of modes at weight k to be
non-decreasing in k, so a vertex algebra with a nondegenerate pairing cannot grow slower than
partitions-like. Not proved in general (would need: invariant form nondegenerate on V_h ==>
growth >= exp(c sqrt E)).

## 4. Numerics

heisS.py full (x,y,u,v restricted to Par_S(h), operators of the full VOA; identity holds by
restriction; exact ranks over Q):

 S          h   N    identity  sum_s min(rank)  N^2   log_N
 all        4   5    True      203              25    3.301   (= G4a)
 all        6   11   True      1139             121   2.935   (= G4a)
 {1,2}      8   5    True      587              25    3.961
 {1,2}      12  7    True      1749             49    3.837
 {1,2,4}    6   6    True      530              36    3.501
 {1,2,4,8}  8   10   True      1798             100   3.255
 (restriction keeps the number of pieces and nearly full ranks: strictly worse than S = all.)
heisS.py comp: see section 3 (all False).
profile.py (out_profile.txt): for S = squares, rho_calc = L(2h)/2L(h) = 0.676 (h=1000),
0.654 (h=10^4), 0.642 (h=10^5) -> 2^{-2/3} = 0.63, but L(4h)/L(h) = 1.81, 1.70, 1.65 -> 1.587 < 2
and the exact inequality (4h+1)(2h+1) max_{w<=4h} p_S(w) >= p_S(h)^2 FAILS for the first time at h = 15269
(slack -9.3 bits at h = 3*10^4, -34 bits at 10^5). For S = all the slack stays positive
(L(4h)/L(h) -> 2 from above).

## 5. Adversarial self-review
- Obstruction A assumes every factorization of every piece goes through the graded components
  of ONE space W with weights bounded by M, and counts ranks only. It does not exclude pieces of a
  different shape (the general floor remains rho >= 1/2).
- Corollary 2 assumes the shift set is an interval and L concave; (V) is an assumption taken
  from vertex algebras, not a theorem about all conceivable operator families.
- Obstruction B is proved only for linear fields (one-oscillator states); the multi-particle case
  rests on the same contraction function kappa and on the numerics.
- The "no slow-growing vertex algebra" statement is heuristic.
- Nothing here improves the published bound; nothing here is a lower bound for #3-SAT.

## 6. Extra profile tables (out_profile2.txt)
 S            first h violating the rank count   rho_calc at largest h   L(4h)/L(h) there
 all          none (slack grows)                 0.715 (h=10^4)          2.040
 odd parts    none (slack grows)                 0.716 (h=10^4)          2.044
 squares      15269                              0.642 (h=10^5)          1.645
 cubes        > 30000 (slack 18.5 -> 7.5, falling) 0.624 (h=3*10^4)      1.548
 powers of 2  3694                               0.577 (h=3*10^4)        1.319
 parts <= 8   440                                0.553 (h=10^4)          1.214
 parts <= 30  603                                0.581 (h=10^4)          1.326
(profile.py as first run used an O(h^2) scan and was stopped after 'squares'; the file now has the
prefix-max fix; out_profile.txt = all + squares, out_profile2.txt = the rest.)
