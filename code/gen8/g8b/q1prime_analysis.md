# G8b: Q1' (counting independent sets at average degree d ~ sqrt(log n)): statements and proofs

Labels: [P] proved here, [H] heuristic, [L] literature, [Q] open. e := Theta(1/sqrt(log n)); logs base 2 unless ln.

## 1. What is actually needed (correcting the form of Q1')
Costs on a d-regular G viewed as a monotone 2-CNF (m = nd/2, D = ceil(d/2)):
* Claude, Prop. R (G6c): exponent C e n log(2d).
* Enumeration-type methods (max-degree branching, frontier/tree-decomposition DP whose states are independent
  sets of a bag): work >= number of leaves/states, about i(G). For random d-regular graphs
  ln i(G) ~ n (ln d)^2/(2d) [H, first moment, H(a) - d a^2/2 at a ~ ln d/d]; the branching recurrence in the
  G(n,p) model (branching keeps p, so degree ~ pk at k vertices) gives the same (n/2d)(ln d)^2 [H].
  Measured: slope ~ 0.22 n log^2 d / d (Section 4).
* min of the two is maximal where d/ln d ~ sqrt(log n): exponent Theta(n log log n/sqrt(log n)).
[P] **Q1' as stated is necessary at d ~ sqrt(log n) and sufficient for d >= sqrt(log n), but NOT sufficient**:
for d = (log n)^c with 0 < c < 1/2, 2^{O(n/d)} is far above the target and Claude pays c n e log log n.
Also min(n/d, e n log d) is maximal at d log d ~ sqrt(log n), value ~ n log log n/(2 sqrt(log n)).
The correct target is **Q1***: #IS of every n-vertex graph with n polylog(n) edges in time 2^{O(n/sqrt(log n))},
i.e. Prop. R with log D replaced by O(1) for monotone 2-CNFs (uniformly in d).
[P] For *average* degree d, a 2^{O(n/d)} bound for all d is far too strong: pad a cubic graph H on N vertices
with a disjoint random 2d-regular graph on N vertices (counted separately, i(G) multiplies); this would count
IS of cubic graphs in 2^{O(N/d)} for every d, and with d = N^c, via the #ETH-style reductions of
Dell-Husfeldt-Marx-Taslaman-Wahlen (#3-SAT -> #IS bounded degree, O(N) vertices), #3-SAT in 2^{O(N^{1-c})},
which is below the 2^{n/log n} "wall" of G4b/G4c. So only regular / min-degree versions are plausible.
Blow-ups G[K_t], G[bar K_t] (values Z_G(t), Z_G(2^t-1), interpolation gives Z_G) keep n/d fixed: n/d is
scale-invariant, Claude's cost is not; so blow-ups cannot move an instance into a cheaper regime [P].

## 2. Route (ii): Claude's network directly on the IS 2-CSP  [P]
For a network whose nodes hold vertex bits (clusters A_a), the edge (a,b) needs c_ab >= log rank M_ab,
M_ab[S,T] = [no edge between S and T]; restricted to an induced matching of H_ab (others 0) M_ab contains the
tensor power of [[1,1],[1,0]] (rank 2), so c_ab >= im(H_ab) (largest induced matching); a vertex cover gives
c_ab <= tau(H_ab). With g = 3K^8 clusters and K ~ n^{1/17}, cluster size s ~ n^{9/17}; in a random d-regular
graph and a uniformly random clustering, whp every H_ab has max degree <= 4 and almost all edges cross,
so im(H_ab) >= e(H_ab)/7 and U = Theta(nd) (exponent e n d) [P]; for an arbitrary clustering the same is
expected since small sets of a random regular graph span O(size) edges [H], worse than Prop. R by d/log d. Small constant g gives U >= Omega(g n) with e = e(K) constant.
So the direct network is worse; Appendix D's marks/parity checks are exactly what turns nd into n log d.

## 3. Route (c) (G6c's Q2) is FALSE on the benchmark family  [P]
**Lemma CA.** G on n vertices, m edges, max degree Delta; F any edge set; 2 <= t <= n;
q = t(t-1)/(n(n-1)). There is X subset V with E(X) ∩ F = ∅ and e(X) >= q (m - |F| (1 + (Delta-1)(t-2)/(n-2))).
*Proof.* X0 uniform t-subset; E e(X0) = qm. For each f in F inside X0 delete a fixed endpoint u_f; this
kills every F-edge and removes <= sum deg_{X0}(u_f) edges; E[deg_{X0}(u_f) 1[f in X0]] =
q (1 + (deg u_f - 1)(t-2)/(n-2)). Linearity. ∎ (Monte-Carlo check: check_avoid.py, n = 3000, d = 8,16,32.)
**Corollary.** Let G be d-regular, Phi_G its IS 2-CNF, n̄ >= 1 the cut-off, t = ceil(2 sqrt(n̄ n/d)) + 2 <= n.
Any list of checks (sets of clauses = edges) such that every V(xi) = E(xi) with |V(xi)| > n̄ meets a check
oddly must cover |F| >= m / (2(1 + d t/n)) edges. With Appendix D's n̄ = n/log(2d):
|F| >= (1 - o(1)) n sqrt(d log 2d)/8; for any n̄ = n/L: |F| >= (nd/4)/(1 + 2 sqrt(d/L) + o(1)).
Hence total incidence I >= |F| and, with r checks and M >= max(r, n) marks, the kind-1 cost
M log(3 s̄) >= min over r of (r+n) log(3|F|/r) >= n (1/2 log d - O(1)), for every choice of n̄ and r.
So on d-regular IS instances Appendix D's route can save at most a factor 2 on log D; Q2 (I = O(n 2^{O(k)}))
is false. Caveat: only checks made of edge clauses; adding redundant clauses (e.g. 3-clauses on triples
containing an edge) changes V(xi) and is not covered (and raises m, cf. Lemma S of G7d).

## 4. Routes (i), (iii), (iv)
(i) Separators: random d-regular graphs have balanced separators of size Theta(n); DP over IS of the separator
    = enumeration type, exponent ~ n log^2 d/d. Treewidth of random d-regular is Theta(n) [L].
(iii) Branch on a set S (|S| = s n/d, avg degree s inside), then Prop. R on the residual: typical leaf leaves
    ~ n/s vertices of degree ~ d/s, but the leaf I_S = ∅ keeps G - S (degree ~ d, cost e n log d); in the
    G(n,p) model (branching preserves p) any adaptive branch-then-Prop.R strategy stopping at size k costs
    (n/2d)[(ln d)^2 - (ln D)^2] + e' D (n/d) ln D, D = dk/n; minimized at D = min(d, ~1/e'), value
    Theta(n e log min(d, 1/e)) for d <= 1/e  [H]. No gain over min of the two.
(iv) Inclusion-exclusion over colour classes (chi ~ d/(2 ln d) classes of size 2n ln d/d) gives a chi-node
    network with U ~ chi n = nd/ln d: worse. Edge inclusion-exclusion: 2^m terms.
Best bound I can PROVE for d-regular IS: 2^{O(min(e n log d, n log d · H))}, i.e. nothing beyond
min(Prop. R, enumeration) (the enumeration side is proved only as O*(i(G)) with i(G) <= 2^{n(d+1)/(2d)}
[Kahn-Zhao] in the worst case, and ~ n log^2 d/d only for random graphs [H]).

## 5. Status
* Q1' and Q1* are not in the literature [L, searched]. ETH-based lower bounds (Traxler: d^{Omega(n)} for
  domain-d binary CSP, giving 2^{Omega(n log d/d)} for MIS at degree d via the microstructure graph; ETH-hardness
  of IS on d-regular graphs, arXiv 2008.09008; #ETH-hardness of #IS at bounded degree, DHMTW) all relied on
  (#)ETH, refuted by Claude's theorem; their reductions transfer quantitatively but only give
  "IS in 2^{o(n log d/d)} => 3-SAT in 2^{o(N)}", which Claude's bound already satisfies: no barrier to Q1'.
* Q1* is open; obstructions proved: direct network (Section 2), parity-check route (Section 3, factor <= 2).
