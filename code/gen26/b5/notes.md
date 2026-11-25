# G26-B5 — decompose into sub-problems, attack the smallest tractable leaves (working notes)
Tags: [P] proof here (unreviewed), [E] exact computation, [N] numerics, [H] heuristic, [O] open.
Read: gen25 BRIEF, P9, P10, g24g, g22b, PISTES 2-5.

## 0. Key reduction used to prune the tree
* 4 rho = sup_theta (theta + T(theta)) (sup over the frontier). Hence constructions only matter if they push the
  frontier down AT THE ARGMAX (balanced region); endpoint improvements (q = O(1) or q = N^beta, beta -> 0) are worthless
  for step 2 on the construction side.  LOWER bounds at ANY single theta beyond the line give rho > 1/2 (route closed).
  So: construction leaves -> balanced region only; lower-bound leaves -> wherever easiest (endpoints!).
* P9: composition (any recursion tree) moves along Phi-additivity; finite-N interior points (e.g. (8,8)@8, theta = T = 1)
  do NOT propagate.  So every construction leaf must be an ASYMPTOTIC family.

## 1. TREE (status: A = attempted (file), N = not attempted, X = closed)
ROOT S2: 3-SAT in 2^{O(n/log n)}
 S2.1 Claude route, per-level rho -> 1/2 (P9: the only lever besides leaves)
  S2.1.1 asymptotic piece family strictly inside circle at balanced point  [A: ~15 families, g20c,g22c,g23a,g24a,b1]
   a. monomial ST model: exponent in [4/3, 2]; leaves: ST (18,18)@12, (28,28)@16 [A g22b,g24d: undecided at 600 s];
      group ansatz [X g23a]
   b. "monomial-C, arbitrary R" model (E_l, C_l signed partial perms, R anything; contains ST and G24g's (12,3)@6)
      [N] -> attacked here (sect. 3)
   c. coset window [A G24f, P10: reduces to (a point strictly inside the circle at one slope)]
   d. border families ((5,5)@5, (8,3)@5) [A g24e,g19b: border only, no family]
   e. composition of finite points [X P9 Lemma U / Prop P]
  S2.1.2 leaf solvers with e0^2 log K0 -> 0  [N; = beating sqrt(log) on clique CSPs themselves; algebraic leaf FMM
         gives mu = 9.1 < 16, worse, P9]  hard
  S2.1.3 D-direction pieces rho(D) < D^{-1/2}  [A g20a numerics N <= 4; no sign]
 S2.2 OpenAI route: kernel with sublinear 2-adic rank profile [A g18d: nothing]; C1 (barrier side) [A g18d numerics L<=8]
 S2.3 CLOSE the route: a bound p q >= N^{2+eps} at ONE theta (=> rho > 1/2)
  a. endpoint, q = N^beta: need p >= N^{2-c beta}, c < 1. Known: line c = 1 (Thm C, Lean); ST c = 1/2 (Lean).
     a1. q = 1, 2: exact N(N-1)/2 [P g24g Thm 1, Lean]
     a2. q = 3: is the gain bounded or linear?  criterion: osr-3 A with (Swap A)^2 = 1, tr(Swap A) > k
         [O in g24g sect 8] -> attacked here (sect. 2)
     a3. strong lemma f(N,r) >= N(N-1)/2 - c r N for CONSTANT r [O g24g]; FALSE for r = N^beta (App. B asymptotically)
     a4. "pencil + low rank" [P10 Prop 4: empty for generic spans]
  b. spreading conjecture 2 theta + T >= 4 [A B2, P7: open]
  c. restricted-model bounds: ST [done], block model [A B2 partial], monomial-C [N, sect. 3]
  d. border version of Thm 1 (q = 2) [A P10 sect 7: failed]
ROOT S1: constant of step 1
 S1.1 Claude: c1 = 4.82 (LV) / 6.81 (det) [P9];
   a. optimal-order lemma zs = O(1/g) -> 4.17/5.89 [C in P9, N as proof]
   b. remove grouping factor t_l (gap 4.8 vs ideal 0.7 is all balancing overhead) [N]
   c. thresholds m* with the P9 levers [N: P9 says not recomputed]   easy, numeric
 S1.2 OpenAI: 75 -> 37.5 under Smith modes (C1) [A P3, g18d]; kappa 6 -> 2 [O]

## 2. Leaf S2.3.a2 (q = 3 linear gain criterion)  -- in progress
## 3. Leaf S2.1.1b / S2.3.c (monomial-C model) -- in progress

### 2.1 progress (01:50 UTC, after container restart)
invol.py (real LM) / involc.py (complex trf): k=3,r=3,t=5: best cost 1.52 (12 real trials); k=3,r=3,t=7: 5.8; k=4,r=3,t=6: 1.63.
Sanity: k=2,r=4,t=4 found (8e-19); k=3,r=4,t=5 -> 8e-8 (converging). So the solver finds r=4 solutions; r=3 none so far.
Reformulation [P]: with L = sum E_l (x) C_l in M_k (x) M_k (product (a(x)b)(c(x)d) = ac(x)bd) and sigma = tensor flip,
 (Swap A)^2 = 1  <=>  L sigma(L) = 1, and tr(Swap A) = sum_l tr(E_l C_l).  r = 1: E C = +-I so tr = +-k <= k [P].
### 3.1 monoC.py (SA over signed partial perms, objective = sum of singular values beyond target) reproduces (5,3)@4, (4,4)@4
 exactly (ranks mod 10007/3/2 equal).
 01:55: k=3,r=3,t=5: 60 more real trials best 1.52319; complex trf (7 trials) best 1.52323 -> the SAME positive minimum
 1.5232 real and complex: strong [N] evidence that no osr-3 involutive Theta with tr > k exists at k = 3 (over R or C).
 monoC (9,3)@5: not found in 240 s (stuck at rank 10, F = sqrt2).

## 3. Monomial-C model (C-part monomial with cancellations allowed, R arbitrary)
Model MC(N,q): K = sum_{l<=q} E_l (x) C_l, E_l, C_l signed partial permutation matrices; p = rank(Swap - K). Contains ST
(G22b), the non-strict optima (7,4)@5, (6,6)@5 (G22b sect 1: G, G' monomial, G' != G^T), and G24g's (12,3)@6.
In cell language Theta = K Swap = sum_l P_l, P_l: cell (a,b) -> +- cell(sigma_l(b), tau_l(a)) (twisted partial perms).
LEMMA MC1 [P]. In MC(N,q): p >= (N^2 - qN)/(q+1).  (Exact at q = 1: N(N-1)/2.)
 Proof. rank(1 - Theta) >= |S| for every cell set S with Theta[S,S] nilpotent (g24g principle). Digraph on cells: c -> P_l c.
 Out-degree <= q. Cells with a self-loop (P_l c = c for some l) number <= qN (for each l, P_l c = c forces tau_l sigma_l b = b,
 <= N cells). Random order of the cells; S = non-self-looped cells all of whose out-neighbours come later. Edges inside S go
 forward => acyclic => Theta[S,S] nilpotent. E|S| >= (N^2 - qN)/(q+1). QED
 Remarks: constant 1 instead of Theorem C's 2/9 but still the LINE (pq ~ N^2). Tight for general digraphs (disjoint
 bidirected (q+1)-cliques), so going beyond the line in MC needs the twisted structure: can the cells be covered by
 "bidirected (q+1)-cliques" of twisted perms?  q = 1: yes (pairs {(a,b),(b,a)}) -> exact. (12,3)@6 has p = N^2/q exactly.
LEMMA MC2 [P] (beyond the line by a log, in MC). Let H = undirected cell graph (c ~ c' iff c' = P_l c or c = P_l c'),
 max degree <= 2q.  (i) every triangle c, c'=Xc, c''=Yc, c''=Zc' (X,Y,Z in {P_l, P_l^{-1}}, twisted partial bijections
 (a,b) -> (mu b, nu a)) gives ZX c = Y c with ZX UNTWISTED ((a,b) -> (mu_Z nu_X a, nu_Z mu_X b)) and Y twisted: then
 b = mu_Y^{-1} mu_Z nu_X a is determined by a, so <= N base cells per (X,Y,Z): <= 8 q^3 N cells lie on triangles.
 (ii) delete them and the <= qN self-looped cells: H' triangle-free, max deg 2q, >= N^2 - 9 q^3 N cells.
 (iii) Shearer: alpha(H') >= |H'| f(2q), f(d) = (d ln d - d + 1)/(d-1)^2 ~ ln d / d. An independent S of non-looped cells
 has (1 - Theta)[S,S] = I, so
        p >= (N^2 - 9 q^3 N) f(2q)  ~  N^2 ln(2q)/(2q)      (q <= N^{1/3}/3, say).
 So in MC, p q >= (1/2 - o(1)) N^2 ln q: strictly beyond the line (by ln q, NOT by N^eps). The general fact behind it:
 odd closed walks are twisted maps and fix <= N cells, so H is locally bipartite-free of short odd cycles; even
 cycles can be everywhere (untwisted words may be the identity). Limit of the method: random-regular-like graphs have
 independence ratio ~ 2 ln d/d, so local arguments cannot give an exponent gain; ST's sqrt(q) used the PALETTE structure.

### 2.2 results (02:20 UTC)  [N]
 k=3: best cost 1.52319 (72 real trials, 7 complex) ; k=4,t=6: 1.63396 (52 real trials, same value repeatedly);
 k=5,t=7: 1.589 (15 trials). r=4 control: solutions found (k=2 exact Swap_2; k=3,t=5 converging).
 VERDICT (numerical, not a proof): no osr-3 A with (Swap A)^2 = 1 and tr(Swap A) > k for k = 3,4,5, over R (and C at k=3).
 So the A (x) I_m family cannot give a linear gain at r = 3 for k <= 5; consistent with "bounded gain at r = 3" (g24g data).
 k = 2 is trivial [P]: t = 4 forces Theta = 1, A = Swap_2, osr 4 > 3.
 Open: a proof for all k; idea: t > k forces L = -1 on a nonzero subspace of skew matrices and L = +1 on a
 (t+k)/2-dim space of symmetric ones, L(W) = sum_{l<=3} E_l W C_l^T.
### 3.2 monoC SA calibration: FAILED to find the known exact MC point (6,6)@5 in 240 s (stuck at rank 7); (9,3)@5,
 (8,3)@5 not found either. => SA negative results are uninformative beyond N = 4. Batch killed.
 tri_check.py [E]: triangle cells (random / cyclic / involutive shift MC instances, N <= 40, q <= 5): 16..290, linear in N,
 far below 8 q^3 N.  Constants of MC2: q f(2q) = 0.69 (q=3), 0.79 (4), 1.13 (10), 2.17 (100): beyond pq = N^2 only for
 q >= ~7; asymptotically pq >= (1/2 - o(1)) N^2 ln q.  So MC2 does NOT exclude pq = N^2 at q = 3 (the (12,3)@6 type).

## 4. Outcome of the attacked leaves
 L1 = S2.3.a2 (linear gain at r = 3 via A (x) I_m): FAILED to find; numerical no-go for k = 3,4,5 (robust repeated
      positive minima 1.523 / 1.634 / 1.589); trivial proof for k = 2. Not a proof.
 L2 = S2.3.c / S2.1.1b (monomial-C model, new): PARTIAL. Lemma MC1 (line with constant 1, exact at q = 1) and Lemma MC2
      (beyond the line by ln q, for q <= N^{1/3}/3) [P]. No exponent gain; shows the line is not tight in MC and that
      MC (which contains every known exact optimum for N <= 8) cannot keep pq = O(N^2) as q grows. SA search tool works at
      N = 4 only (fails to rediscover (6,6)@5): computational side inconclusive.
 (Third leaf: none attacked beyond these; see ranking.)

## 5. Ranked next leaves (value x tractability), with first step
 1. MC exponent: does MC satisfy p >= N^{2-c beta} (c < 1) at q = N^beta, i.e. pq >= N^{2+eps}?  First step: replace
    independent sets by acyclic sets + use EVEN cycles: an even closed walk is an untwisted map (alpha, beta) fixing
    Fix(alpha) x Fix(beta); show that many such words must be the identity on few cells unless the P_l "commute",
    then a Ramsey/girth argument. A pq >= N^{2+eps} bound in MC would be the first beyond-line exponent in a model that
    contains all known exact optima (ST already has it, but excludes them).
 2. Better MC search: exact SAT/ILP over F_3 or GF(2)-sign-free supports instead of SA (G24d's CaDiCaL encodings
    adapted: unknown supports of E_l, C_l + rank condition via an explicit rank-p factorisation R = U V^T). First step:
    reproduce (6,6)@5, then decide (5,5)@5 and (11,3)@6 inside MC.
 3. r = 3 involution no-go proof: use L = -1 on skew U2 != 0 and L = +1 on symmetric U1 of dim (t+k)/2 with
    L(W) = sum_{l<=3} E_l W C_l^T; first step: k = 3 by Groebner basis / exact elimination on the 54 unknowns with the
    gauge (E1,E2) -> (I, D) (g24g sect 3), which leaves ~30 unknowns.
 4. Strong lemma for constant r: f(N,3) >= N(N-1)/2 - c N. First step: prove it in MC (cell-graph argument with q = 3:
    show the non-independent part is O(N) cells using MC2(i)-type counting of 2-cycles).
 5. Constants of step 1 (S1.1c): recompute thresholds m* with P9's levers (row sweep + unbalanced levels) using p2's
    thmA_explicit.py; cheap, gives the article numbers.
 6. ST (18,18)@12 / (28,28)@16: still the best test of whether monomial arrays beat N^2/8; needs > 600 s SAT runs or
    symmetry breaking (G24d encodings).
Files: notes.md, invol.py, involc.py, monoC.py, tri_check.py, logs log_*.txt, mc_*.txt. All my processes killed.
