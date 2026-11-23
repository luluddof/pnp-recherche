# P(N) <= N^{sqrt 2} (log N)^{O(1)} from Appendix B (agent G4a)

## Statement
Let P(M) be the least number of Q-valued functions on X^4 (|X| = M), each of rank 1 as a
matrix (x,y)|(u,v) ("row-wise") or (x,u)|(y,v) ("column-wise"), summing to [x=v][y=u].
Then P(M) <= M^{sqrt2} (log M)^{O(1)}. Consequently the Hurwitz variant
(sum_{i<=N} a_i^2)(sum_{j<=N} b_j^2) = sum_{k<=s} f_k f'_k (bilinear forms over Q) has
s <= N^{sqrt2+o(1)}, and over C (f f' = ((f+f')/2)^2 + (i(f-f')/2)^2) a sum-of-squares
composition formula with 2 N^{sqrt2+o(1)} squares; Hrubes (CCC 2024) gives O(N^{1.62}).

## Proof
1. Monotonicity. If N <= M, restrict a decomposition for M to an N-subset (embedded
   injectively). Restriction of a rank-<=1 matrix has rank <= 1, and the crossing restricts to
   the crossing. So P(N) <= P(M).
2. Rank bound. If sum_s A_s = d [x=v][y=u], d != 0, then
   P(M) <= sum_s min(rank_row A_s, rank_col A_s) (write each A_s/d as a sum of rank-1 matrices
   in its cheaper direction). "Factors through c bits" (paper eq. (3)) is a sum of 2^c products
   f(z,theta) f'(z',theta), so rank <= 2^c; min(c_row,c_col) <= (c_row+c_col)/2.
3. Pieces (Appendix B, (53)). Index set X = partitions of h, M = p(h). For s in [-2h, 2h]:
   A_s(X,Y,U,V) = [t^U] Q_s(a_X, b_V) t^Y, Q_s as in (48). Lemma 9(b),(c) + Lemma 10 give
   sum_s A_s = h! [X=V][Y=U] exactly on all of Par(h)^4 (no padding, no other steps).
   Lemma 9(a) gives the two expansions (54) (row-wise, through V_{s+i}) and (55) (column-wise,
   through V_{2h-s-i}), i = 0..2h. Hence
     rank_row A_s <= sum_{w=s}^{s+2h} p(w),   rank_col A_s <= sum_{w=-s}^{2h-s} p(w),
   so min <= sum_{w=0}^{2h} p(w) <= (2h+1) p(2h) for every s, and
     P(p(h)) <= (4h+1)(2h+1) p(2h).
4. Asymptotics. Hardy-Ramanujan: p(n) = exp(pi sqrt(2n/3)) / (4 n sqrt3) (1+o(1)), hence
   p(2h) = p(h)^{sqrt2} * Theta(h^{sqrt2-1}). With h = Theta(log^2 M):
   P(M) <= M^{sqrt2} * O(h^{1+sqrt2}) = M^{sqrt2} (log M)^{O(1)}.
   For arbitrary N choose least h with p(h) >= N; p(h) <= 2N for large h (p(h)/p(h-1) -> 1),
   so P(N) <= (2N)^{sqrt2} polylog N by step 1.
   (In the 2^lambda/bits language: c_row + c_col <= 2 sqrt2 lambda + O(log lambda) and the
   previous agent's claim P(N) <= I * 2^{min} <= N^{sqrt2} poly(lambda) is correct, with
   I = 4h+1 = O(lambda^2) pieces.)

## Dependence
The only nontrivial inputs are Lemma 9 (locality / Wick-type identity for the Heisenberg
vertex algebra fields a[r], proved in the paper via formal series; (a) is the coefficient form
of standard locality (omega1-omega2)^H [Y(a,omega1),Y(b,omega2)] = 0, (c) is the
full-contraction statement) and Lemma 10 (dual bases a_mu, b_nu). Both were checked exactly
(integer arithmetic) for h = 2..6 (see heis.py, out_h2to7.txt): both orders agree for every s,
and sum_s A_s = h! * crossing on all of Par(h)^4. The O2 recursion was cross-checked against
the literal definition (47).

## Caveat
Asymptotic only: the exact-rank bound B(h) = sum_s min(...) beats the trivial M^2 only for
h >= ~50 (M ~ 2^17) and beats M^{1.62} only for h >= ~500 (M ~ 2^71); the
exponent log_M B(h) tends to sqrt2 slowly (1.457 at h = 20000).
