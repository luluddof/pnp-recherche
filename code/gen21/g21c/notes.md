# G21c — Other routes to "step 2" (3-SAT in 2^{O(n/log n)}): equivalences, conditional theorems, one route pushed

Tags: [K] known result (reference given; R = I read the text, A = abstract/snippet only, M = from memory, not re-checked),
[P] proof written here (short, unreviewed), [S] sketch, [C] conditional on the two Oct-2026 preprints being correct
(Claude 2610.00577v1, OpenAI "ETH is false"; neither re-verified here), [O] open.
Downloads (untrusted data): dl/ (arXiv abstract pages + PDF of 2410.02606). theoryofcomputing.org and cs.bme.hu are not
reachable from this session, so Marx's paper itself was NOT read; its statements are quoted through 2410.02606 (R) and
2311.05913 (A).

## 0. Two scales, and a correction to the task statement
* The preprints give 2^{O(m/sqrt(log m))} in the number m of clauses (Claude Thm A(c) / OpenAI Thm 1.1 at m <= Dn), but
  only 2^{O(n loglog n/sqrt(log n))} in n alone (both papers [R]).  "Step 1 = 2^{n/sqrt log n}" is exact only at linear
  density.  Passing from the m-scale to the n-scale goes through sparsification with eps = 1/f(n), whose clause blow-up
  c(eps) = (k/eps)^{O(k)} [CIP06, M] costs polylog factors.  So "step 2 in m" (2^{O(m/log m)}) does NOT give step 2 in n
  without a sparsification lemma with c(eps) = O(1) at eps ~ 1/log n (not known; I do not know a lower bound either). [O]
  Everything below is in the m-scale (linear-size instances) unless said otherwise.
* ETH <=> no f(k)N^{o(k)} for k-Clique is NOT a known equivalence.  [K, M] Chen-Huang-Kanj-Xia (JCSS 72, 2006) prove only
  ETH => no f(k)N^{o(k)}; the converse is open, and the 2026 algorithms do not give it (G9e: "gain lost for k-Clique").

## 1. Equivalences / implications with exact parameter translation

### 1a. k-Clique (multicoloured: k parts of N vertices)
Reduction [K, M; standard, e.g. CHKX06, Cygan et al. book Thm 14.21; described in 2410.02606 sec. 1, R]: split the m
clauses of a 3-CNF into k groups; a vertex of part i = an assignment of the <= 3m/k variables of group i satisfying it
(N <= 7^{ceil(m/k)}); edges = consistency.  k-cliques <-> satisfying assignments (bijection; counting preserved).
 => Prop 1a [P]: if MCC(k,N) is solved in time N^{k/h(k)} poly(kN) (uniformly in k), then 3-SAT with m clauses is solved in
    2^{2.81 m/h(k)} poly for every k <= m; take k = m/(C log m), N = m^{2.81C} = k^{O(1)}:
       step 1 <= MCC with N = k^{O(1)} in 2^{O(k sqrt(log k))}       (i.e. N^{O(k/sqrt log k)})
       step 2 <= MCC with N = k^{O(1)} in 2^{O(k)}                   (i.e. N^{O(k/log k)})
    against trivial N^k = 2^{Theta(k log k)} and best known N^{omega k/3} (Nesetril-Poljak 1985, [M]) = constant factor.
    Only for k(n) growing (k ~ m/log m, N polynomial in k) can a clique algorithm give step 2; fixed k never does.
 Converse: none.  MCC(k,N) as SAT needs k log N variables but Theta(k^2 N^2) clauses of width 2 log N (or a Tseitin blow-up
 of the same size): dense, wide; neither preprint's bound (M-bound needs few wire bits; N-bound needs constant width, R2 of
 G9e) gives anything.
 What the 2026 algorithms give DIRECTLY on MCC [C,P]: MCC is the binary CSP on K_k with domain N, i.e. (via Def. 2 of the
 Claude paper, each edge factoring through log N bits) a network with U = k(k-1) log N wire bits.  Running Lemma 1 +
 Thm A(b): 2^{O(k^2 log N/sqrt(log k))}: worse than N^k as soon as k > C sqrt(log k), i.e. always for growing k.
 Structural reason: the cost is paid per WIRE bit, and a clique instance carries each node's value on k-1 wires (U/n = k-1);
 it is the WORST density for these methods, not the best.  The replication [x=y]^{k-1} is the "reused equality" case of G9d
 sec. 2(c): break-even at best.  Dense/symmetric/small bundles do not help: the bundle per edge is log N bits, and symmetry
 is not used by either algorithm.

### 1b. CSP and treewidth (Marx)
Facts [K]:
 * Marx, "Can you beat treewidth?", ToC 6 (2010) 85-112 [A via 2410.02606 (R) and 2311.05913 (A)]: assuming ETH, no class of
   graphs with unbounded treewidth admits binary CSP in f(G)||I||^{o(tw/log tw)}; ColSub(H) for fixed H of treewidth t has no
   O(n^{alpha t/log t}) algorithm.  Upper bound n^{tw+1} (dynamic programming).
 * Curticapean-Doring-Neuen-Wang, "Can you link up with treewidth?", arXiv 2410.02606 [R, sec. 1 and proof of Thm 1.3]:
   linkage capacity gamma(H); Thm 1.3: ETH => no O(n^{alpha gamma(H)}) for ColSub(H).  The reduction: 3-Colouring on n vertices
   of max degree 4  ->  ColSub(H) with parts of size s <= 3^t, t = 12 ceil(n/gamma).  gamma(H) = Omega(t/log t) for tw t,
   gamma = O(tw); k-vertex max-degree-3 expanders have gamma = Theta(k/log k), tw = Theta(k); dense H (average degree
   k^beta) have gamma = Theta(k) (tight bounds there).  "Whether the logarithmic loss in Marx's lower bound can be avoided"
   is stated as a major open question (Conjecture 1.7, "You cannot beat treewidth").
Exact translation [P from the reduction as written]: ColSub(H) in time s^{gamma(H)/h} poly  =>  3-Col(deg<=4) and hence
 (linear reductions) 3-SAT(m) in 2^{O(m/h)}.  For expanders (gamma ~ k/log k):
       step 1 <= ColSub(expander_k) in s^{O(k/(log k)^{3/2})},   step 2 <= ColSub(expander_k) in s^{O(k/(log k)^2)}.
 So step 2 is NOT "beating treewidth by the log factor Marx leaves open".  Marx's open log is the congestion (tw/gamma),
 a routing quantity; s^{O(tw/log tw)} for expanders would give, through this reduction, only 3-SAT in 2^{O(m)}: nothing.
 The two logs are different logs.
What a fast SAT-type algorithm gives for CSP (other direction) [C,S]: see Theorem T1 below: the preprints' network form gives
 binary #CSP on max-degree-Delta graphs with v variables and domain s in s^{O(Delta v/sqrt(log v))}; on bounded-degree
 expanders this is s^{O(tw/sqrt(log tw))}: treewidth beaten by sqrt(log tw) (conditional).  A step-2 algorithm IN THE
 WIRE-BIT MODEL (networks with U wire bits in 2^{O(U/log U)}) would give s^{O(tw/log tw)} = s^{O(gamma)} on bounded-degree
 expanders: it would close Marx's gap from ABOVE (algorithmically) for that class, and by Thm 1.3 that is exactly the point
 where the reduction stops transferring anything.  Fate of Marx's theorem if ETH is false at the sqrt-log level: its proof
 still gives "if 3-SAT(m) needs 2^{Omega(m/a(m))} then ColSub(H) needs s^{Omega(gamma(H)/a)}"; with a = sqrt(log): no
 s^{o(k/(log k)^{3/2})} for expanders, against the (conditional) upper bound s^{O(k/sqrt(log k))}; the gap stays a factor
 log k (= tw/gamma) whatever a is.  Under SH / Weak ETH (Fan-?/..., arXiv 2611.02650 [A]: no 2^{n^{1-eps}}) Marx-type bounds
 become s^{gamma^{1-eps}}-type.  arXiv 2611.12321 [A] already uses "3-SAT has no 2^{o(N/log N)} algorithm" (= "step 2 is
 false") as a hypothesis.
Unconditional classes with CSP in N^{o(tw)}: none known to me for unbounded tw [M]; for fixed-domain/fixed-language CSPs
 the right parameter is c^{tw} (e.g. arXiv 2611.16729 [A]), a different regime (domain fixed, not N).  Planar etc.: tw is
 already sqrt(v), and N^{O(sqrt v)} is ETH-tight [M].
Reading Marx's embedding backwards as an algorithm: refuted below (Theorem T2/T3): any round trip SAT -> ColSub(bounded-degree
 H) -> network multiplies the wire bits by Omega(log k), which exactly eats a step-2 saving.

### 1c. Other reformulations (m-scale; "<=" = only this direction known)
| problem (size param) | equiv. to step 1 | equiv. to step 2 | relation | best known unconditional / [C] |
|---|---|---|---|---|
| 3-Col, max deg O(1) (n vertices) | 2^{O(n/sqrt log n)} | 2^{O(n/log n)} | linear reductions both ways [K,M] | [C] 2^{O(n/sqrt log n)}; dense: 2^{O(n loglog n/sqrt log n)} via N-bound (OpenAI Cor 1.3, R) |
| Indep. Set / Vertex Cover, max deg O(1) | same | same | both ways (IPZ01 SERF, [M]) | [C] same |
| Circuit-SAT, size s = O(n) | 2^{O(s/sqrt log s)} | 2^{O(s/log s)} | both ways (Tseitin) | size s = n g(n): gain only if g = o(sqrt log n) (G9e R1) |
| MCC (k parts, N = k^{O(1)}) | 2^{O(k sqrt log k)} | 2^{O(k)} | <= only (1a) | N^{omega k/3}; [C] nothing |
| ColSub(expander_k), parts s | s^{k/(log k)^{3/2}} | s^{k/(log k)^2} | <= only (1b) | s^{O(k)}; [C] s^{O(k/sqrt log k)} (T1) |
| 2-CSP, n vars, domain d, (d,2)-CSP | d^{n/sqrt log n} at d = n^{O(1)} | d^{O(n/log n)} = 2^{O(n)} at d = n^{O(1)} | <= only (group log n clauses per variable; Traxler IWPEC 2008 [M] for the ETH form) | d^n; [C] d^{O(Delta n/sqrt log n)} for max degree Delta (T1) |
| Subset Sum (n items, O(n)-bit numbers) | 2^{O(n/sqrt log n)} | 2^{O(n/log n)} | <= only [M, standard] | 2^{n/2} (Horowitz-Sahni/Schroeppel-Shamir) [M]; [C] nothing |
| CVP_p, dimension n | 2^{O(n/sqrt log n)} | 2^{O(n/log n)} | <= only (Bennett-Golovnev-Stephens-Davidowitz 2017, ETH form, [M]) | 2^{O(n)} [M] |
| k-OV: k sets of N vectors, dim d = ck log N | N^{k/sqrt log(kN)} | N^{O(k/log k)} at N = k^{O(1)} | <= only (split variables into k blocks; d = m) | k = 2, d = c log N: N^{2-1/O(log c)} (Abboud-Williams-Yu 2015, [M]) |
| k-SUM (N numbers) | N^{k/sqrt log k} at N = k^{O(1)} | N^{O(k/log k)} at N = k^{O(1)} | <= only (Patrascu-Williams SODA 2010: k-SUM in n^{o(k)} => ETH false, [M]) | N^{ceil(k/2)} [M] |
| chromatic number / Set Cover (n) | 2^{O(n/...)} | | <= only | 2^n poly; under ARC 1.99982^n (BCHKP, arXiv 2404.04987 [A]) |
Pattern: every problem with a LINEAR reduction in both directions (bounded-degree graph problems, linear circuits) is
equivalent at every step; every "parameterized" reformulation (Clique, k-SUM, k-OV, ColSub) is one-way, and the useful
regime is always k(n) ~ n/log n with polynomial N, where step 2 reads "2^{O(k)} instead of 2^{Theta(k log k)}".

## 2. Conditional routes through famous conjectures

### 2a. Asymptotic rank conjecture (ARC) / omega = 2
[K, A] Bjorklund-Kaski arXiv 2310.11926 (STOC 2024): ARC => set cover conjecture false (k-Set Cover in (2-eps)^n) and a
(2-delta)^n randomized directed Hamiltonicity.  Pratt arXiv 2311.02774 (STOC 2024): ARC => set cover (bounded sets) in
(3/2^{2/3}+eps)^n.  Bjorklund-Curticapean-Husfeldt-Kaski-Pratt arXiv 2404.04987: ARC (char 0) => chromatic number in
deterministic O(1.99982^n).  All abstract-only.  All are BASE-shaving (2^n -> c^n) for problems whose trivial bound is 2^n
over subsets of a universe; none mentions k-CNF-SAT; none gives 2^{o(n)} for anything.
Inside the 2026 recursion (G9d Lemma T, G18f sec. 4; made precise here):
 Theorem R [P, one line]: let T = T_{g,h,lambda} (EPR pairs on the diagonals of a g x h grid, S bits).  Any expansion
 T = sum_{r<=R} a_r^{(1)} (x) ... (x) a_r^{(h)} into terms of rank one over the h column parties has R >= max over column
 bipartitions J|J^c of the flattening rank = 2^{lambda g(g-1)|J||J^c|} >= 2^{(1/2 - O(1/h))S}  (for h = 3: R >= n^2, n = N^{g(g-1)};
 ARC/omega = 2 makes R = n^{2+o(1)}, i.e. 2^{(2/3+o(1))S}).
 Proof: flattening rank is subadditive and equals 1 on each term.  QED
 Consequences: (i) such a level must ENUMERATE >= 2^{S/2} terms; the recursion's compounding needs 2^{o(S)} terms per level
 (E1 asks exactly that).  So rank can be paid once, at one level, and multiplies the brute-force cost of that level by a
 constant, never the recursive cost.  Best use (G9d, [S]): at the leaves: 3-SAT in 2^{(2 omega/9) c m/sqrt(log m)}, i.e.
 under omega = 2 the constant drops by 4/9; the shape m/sqrt(log m) is unchanged.  (ii) "Exact missing step" asked in the
 task: a rank-R decomposition would be usable recursively only if each factor a_r^{(j)} were itself a cheap K_g-network on
 column j AND R = 2^{o(S)}; the second condition is impossible by Theorem R.  So the obstruction is not (only) that the
 summands are not contractible networks: it is the term count, forced by flattening, independent of any conjecture.
 (iii) Cleanest conditional statement: ARC (even in its strongest K-party "asymptotic rank = max flattening" form) gives
 at most a constant factor in the exponent of any algorithm of the 2026 type; ARC => step 2 is NOT available by this route.
 Does ARC refute (Weak) ETH-type statements elsewhere?  Not known to me [M]; BK/Pratt stay at c^n.

### 2b. Other famous statements (one paragraph each; all "no known implication")
* Walsh-Hadamard non-rigidity (Alman-Williams FOCS 2017 [M]; G20c): WH is not rigid enough for Valiant's programme; the
  mechanism (polynomial method on Kronecker powers) gives rank r = N^{1-Theta(eps^2)} with sparsity N^{O(eps log(1/eps))}; in
  G18f's reading this is an entangled term with alpha + beta > 1 [S].  Non-rigidity is about one cut; step 2 needs a joint
  (row,column) saving with 2^{o(S)} terms.  No implication known or visible.
* Log-rank conjecture: a statement about communication protocols (log of the number of rectangles vs log rank, polylog
  gap).  Our costs are already ranks; the conjecture concerns 0/1 matrices and polylog factors in log-scale, i.e. it would
  change exponents by polynomial factors in log, at the level of protocols, not of signed decompositions.  No implication.
* PFR (Gowers-Green-Manners-Tao 2023, arXiv 2311.05762 [M]): polynomial structure for sets of small doubling in F_2^n.  Used
  in fine-grained complexity via BSG-type arguments (3SUM, [M]); SAT instances have no additive structure to exploit and the
  losses are polynomial in K, i.e. constant factors in exponents at best.  No implication.
* Sunflowers (Alweiss-Lovett-Wu-Zhang 2019 [M]): the sparsification lemma (IPZ01) is sunflower-based, so better sunflower
  bounds can only improve c(eps) in sparsification, i.e. the translation m-scale -> n-scale of sec. 0.  That is the ONE place
  where a famous bound touches the n/log n scale: step 2 in m plus a sparsification lemma with c(eps) = O(1) at
  eps = 1/log n would give step 2 in n.  Whether such c(eps) is possible is open to me [O]; ALWZ bounds are for unbounded
  k-sets and do not obviously give it.
* Cap set / slice rank (Ellenberg-Gijswijt 2017, Tao [M]): tools for LOWER bounds on rank-like decompositions (G9d uses slice
  rank of MM >= subrank to show slice-type triangle pieces have rho = 1 - o(1)).  They are barrier tools, not algorithms.

### 2c. Circuit lower bounds (Williams' programme)
[K] Williams (JACM 2014 / Wil13 Thm 6.1, [A via both preprints]): ETH false => E^NP not in SIZE(O(n)).  OpenAI Cor 1.4 [R,
proof read]: from 3-SAT at density <= 5 in 2^{A V/sqrt(log V)}, via Jahanjou-Miles-Viola, E^NP needs > c n sqrt(log n)
gates infinitely often.  The proof (guess a circuit of s_n = eta n a(n) gates, test circuit -> 3-CNF of O(s_n) variables,
time 2^{3A eta n + o(n)} < 2^{n/2}) uses only a(.) with a(Theta(n a(n))) ~ a(n).  Hence [P, same proof]:
   3-SAT (linear density) in 2^{O(V/a(V))}  =>  E^NP has no circuits of size c n a(n) (i.o.), for any such a.
   step 1: c n sqrt(log n) (OpenAI);   step 2: c n log n.
 Neither touches superpolynomial bounds (those need SAT for poly-size circuits beating 2^n by n^{omega(1)}, e.g. NEXP not
 in ACC0 [M]).  An n log n lower bound for E^NP would still be far below S_2E's near-maximum bounds (Chen-Hirahara-Ren,
 Li 2024, cited by Claude [A]) and above the best E^NP bound 3.1n (Li-Yang 2022 [A]).  So step 2 is not "shocking" by its
 circuit consequences; its surprise is algorithmic (it is the m/log m wall of all product-piece schemes, G9d).

## 3. Route chosen: CSP / treewidth / linkage capacity (pushed to two theorems)

Why this route: it is the only one where the 2026 algorithms give something DIRECT and new-looking (T1), and where the
"embed and solve" idea can be decided (T2/T3).  ARC: closed by Theorem R.  Clique: closed by 1a.

Theorem T1 [C; S]  (binary #CSP via the network form).  Let H have v vertices, max degree Delta, domain s per vertex, and
 binary weight functions f_e : [s]^2 -> Z with |f_e| < 2^B computable in poly time.  Then sum_x prod_e f_e(x_u,x_v) is
 computed deterministically in time 2^{O(Delta v ceil(log s)/sqrt(log v))} poly(v, s, B).
 Proof sketch.  Follow the proof of the m-bound (Claude paper sec. 6, "Clusters ... Cost" [R] and Thm A(c)): put the v
 CSP variables into g clusters of <= ceil(v/g) variables; node setting = assignment of the cluster (P_a <= ceil(v/g) log s
 bits); an edge between clusters a,b factors through log s bits (theta = x_u; f_{a,b} = [theta = x_u],
 f_{b,a} = f_e(theta, x_v)), so c_{a,b} <= (#edges between a,b) ceil(log s), sum_b c_{a,b} <= Delta ceil(v/g) ceil(log s) =: P.
 Definition 2 holds; U <= Delta v ceil(log s), gP <= Delta(v+g) ceil(log s).  Lemma 1 + Thm A(b) give exponent
 e_l (U + gP/K) + 2gP/K with K = K_l <= v^{1/10}, e_l = Theta((log K_l)^{-1/2}): O(Delta v log s / sqrt(log v)).
 Gap: I did not re-run the paper's choice of l and the poly(lambda+B)^{C K^2} factor; it is the same calculation as for 3-CNF
 with 6m replaced by Delta v log s, valid when log s <= v^{O(1)}.  Not checked: that Thm A(c)'s proof uses nothing about
 the 3-CNF beyond (U, gP, B, tau).
 Corollary T1' [C; S]: ColSub(H) / binary CSP on bounded-degree expanders H (tw = Theta(v)) in s^{O(tw/sqrt(log tw))}:
 "you CAN beat treewidth by sqrt(log tw)" for this class, if the preprints are right.  Novelty: not established; neither
 preprint states it (they cite Mar10 only as an ETH consequence, [R]); G9e treated CSP only through binary encoding into
 k-CNF (rule R2: d^{o(n)} only for log d = o((log n)^{1/4})); T1 removes that restriction for bounded degree.
 For MCC (H = K_k, Delta = k-1) T1 is useless (1a).

Theorem T2 [P]  (linkage capacity of bounded-degree graphs).  If H has k vertices and max degree Delta >= 2 then
 gamma(H) * floor(log_Delta(gamma(H)/4)) <= 4k + o(1)   (t -> infinity);  hence gamma(H) = O(k log Delta / log k).
 Proof.  Let X be matching-linked in F = H [x] K_t (strong product), |X| = x = gamma t.  The F-vertices whose H-projection is
 within H-distance r of a given H-vertex number at most t(1 + Delta + ... + Delta^r) <= 2t Delta^r.  Greedily pair elements
 of X at H-distance > r: while more than 2t Delta^r + 1 elements are unpaired, an unpaired element has an unpaired partner
 at distance > r.  This gives a matching M on X with |M| >= (x - 2t Delta^r - 1)/2.  A path in F between two vertices whose
 projections are at H-distance > r has at least r internal vertices (each F-edge moves the projection by <= 1).  The |M|
 paths are vertex-disjoint and F has kt vertices: r(x - 2tDelta^r - 1)/2 <= kt.  Take r = floor(log_Delta(gamma/4)), so
 2Delta^r <= gamma/2: r gamma/4 <= k + o(1).  QED  (For expanders this matches gamma = Theta(k/log k) stated in 2410.02606;
 the general bounded-degree upper bound is probably folklore; I did not find it stated.)

Theorem T3 [P given T2 and the reduction of 2410.02606 Thm 1.3]  (the embedding round trip is neutral at step 2).
 Route: 3-Col(n, deg<=4) -> ColSub(H) (parts s <= 3^{12 ceil(n/gamma)}) -> network of T1.  Wire bits of the result:
 U = 2|E(H)| ceil(log s) >= (k-1) * 19 n/gamma(H) >= Omega(n log k / log Delta) for connected H of max degree Delta (T2).
 (Uses |E(H)| >= k-1 and T2.)  For dense H only |E|/gamma >= |E|/k = d_avg/2 is claimed.  So a network algorithm
 with saving factor a(U) (time 2^{O(U/a(U))}) yields for 3-Col, by this route, 2^{O(n log k/(a(U) log Delta))}: at step 2
 (a = log U, U <= poly(n)) this is 2^{O(n)} for every choice of bounded-degree H.  For large Delta (dense H) the loss is
 |E|/gamma >= Delta/O(1) instead.  So "embed the sparse CSP into a structure where a fast algorithm is known" cannot produce
 step 2 from any algorithm whose cost is linear in the wire bits; it would need a ColSub algorithm whose cost scales with
 gamma(H) (routing capacity) rather than |E(H)| (wires): that is exactly Marx's open question, read algorithmically.
 Caveat: T3 is about this reduction (linkage/topological-minor embeddings in blow-ups), not about all reductions.

First open step of the route [O]: find ANY graph family H_k with ColSub(H_k) (or binary CSP) in s^{o(|E(H_k)|/log k)} *
 polylog-free factors, i.e. an algorithm whose exponent is below (wire bits)/log; on families with gamma(H_k) =
 Theta(|E(H_k)|) (dense, Thm 1.4/1.5 of 2410.02606) any s^{o(gamma)} algorithm refutes ETH outright, so the only room is
 sparse families, where the wire count exceeds gamma by log k, and a step-2-type algorithm is exactly what is missing.

## 4. Verdict
* Equivalences: only linear-size problems (bounded-degree 3-Col/IS/VC, linear Circuit-SAT) are step-by-step EQUIVALENT to
  3-SAT; Clique, k-SUM, k-OV, (d,2)-CSP, ColSub, Subset Sum, CVP are one-way (fast X => fast SAT), with step 2 <=> "2^{O(k)}
  instead of 2^{k log k} at N = k^{O(1)}" or ColSub(expander) in s^{k/log^2 k}.
* Conditional theorems: T1 [C,S] (CSP on bounded degree in s^{O(Delta v/sqrt log v)}, beats treewidth by sqrt log on
  expanders); T2 [P]; T3 [P] (embedding round trip neutral at step 2); Theorem R [P] (any CP-rank level has >= 2^{S/2}
  terms; ARC/omega=2 buy a constant only; "ARC => step 2" unavailable by this route); circuit consequence [P, OpenAI's proof
  generalised]: step 2 => E^NP not in SIZE(c n log n) i.o.
* Most promising route: none reaches step 2.  The CSP route at least locates the missing object: an algorithm whose cost
  scales with routing capacity rather than wire count; first open step stated at the end of sec. 3.

## 5. Adversarial self-review
* T1 rests on the Claude paper's Lemma 1 / Thm A as written; I checked Definition 2 and the 3-CNF cluster construction [R]
  but not the recursion constants; if the paper is wrong, T1 is void.  "Beats treewidth" also needs tw(expander) = Theta(v)
  [M, standard].
* T2: the claim "each F-edge moves the projection by <= 1" uses the strong-product definition H [x] K_t (edges inside a
  blob or along H-edges) as in 2410.02606 [R].  Endpoints in X may sit in the same blob; handled by distance > r.  The
  constant is loose.  The definition of gamma there involves "for all large t" and floor(ct): fine for the inequality.
* Growing k: the general reduction of 2410.02606 runs in O(n^{f(k)}), f(k) = 3k^{k+2}, so for k = k(n) growing one must use
  their explicit Benes-type graphs (Thm 4.4) and I did NOT check that their routing is computable in poly(n,k); the step
  translations for ColSub(expander_k) assume it.
* T3 compares wire bits of one specific reduction; an algorithm whose cost is not linear in wire bits escapes it.
* The m- vs n-scale issue (sec. 0) means "step 2" must always be stated with the scale; I used m.
* Marx's paper not read (blocked host); its statements come from 2410.02606 (read) and abstracts.  Traxler, Patrascu-
  Williams, BGS17, AWY15, CIP06, Alman-Williams, ALWZ, GGMT, EG17: from memory, not re-checked here.
* Table entries "best known" for Subset Sum, CVP, k-SUM: from memory.
