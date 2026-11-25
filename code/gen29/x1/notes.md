# G29-X1 notes: is SWAP the only crossing operator? (start 12:50 UTC)

Tags: [P] proof written here, [N] numerics, [E] exact (integer construction, ranks mod 2 primes), [R] read in paper/notes.

## 0. Reading (12:50-13:00) [R]
- Paper Sec 2.5 / Prop 3: nodes of a (g^2,lam)-clique on a g x g grid. Edge {w,w'}: same row or same column = straight
  ("side"), otherwise a diagonal of exactly one rectangle (2 rows x 2 cols). Rectangle corners a=(r,j), a'=(r,j'),
  b=(r',j), b'=(r',j'); strings x=z_a(b'), y=z_a'(b), u=z_b(a'), v=z_b'(a); the two diagonals give [x=v][y=u].
- Row network: nodes = rows, legs (x,y) at top row, (u,v) at bottom row -> ROW flattening (x,y)|(u,v): X is the
  swap permutation matrix, rank N^2 (N = 2^lam). Column network: (x,u) at left column, (y,v) at right -> COLUMN
  flattening (x,u)|(y,v): rank of X there = operator-Schmidt rank (X = sum E_xv (x) E_yu over (x,u),(y,v)). So
  p = row rank of R (merge rows), q = column rank of C (merge columns). Term iota chooses rows or columns GLOBALLY.
- Def. 7: pieces sum to d*X, each factors through c_row (rows) and c_col (cols), poly(lam) many. Lemma 8 (balancing)
  needs pairwise network (Def. 2), c_{a,b} <= P/g.

## 1. Why SWAP (task 1) [P, design-level]
Two partitions (rows pi_R, columns pi_C) of the node set; WLOG a grid (common refinement; g x h allowed, Lemma U).
A strand (= one bit-string edge, an identity between its two legs) is
  type RC (diagonal): cut by both partitions; type R (column side edge): cut by rows only; type C (row side edge).
Requirements for a gadget (set of strands decomposed together into pieces):
 (G1) each piece must be a PAIRWISE function in the merged network of each direction (Def. 2), or be refactored as a
      pairwise network among the classes it touches (a k-party piece costs the total bond bits of such a network);
 (G2) the gadgets partition the strands, product of gadget identities = product of all edge tests.
A single RC strand has the SAME bipartition {x}|{v} in both directions: any R + C = I has p + q >= N, every piece costs
the same in both directions -> no gain. Gain needs two DIFFERENT bipartitions of the gadget's legs. Inside 2 rows x 2 cols
(the only 2-party-in-both-directions configuration) the RC edges are exactly the 2 diagonals => crossing
[x=v][y=u] = SWAP in the row flattening and also in the column flattening (realignment of SWAP = SWAP); the third
pairing {x,v}|{y,u} would cut nothing (partial transpose of SWAP = rank 1).
Classification of 2-party (both directions) gadgets: legs sit in the 4 cells TL,TR,BL,BR; a strand is
 main (TL-BR) or anti (TR-BL) [RC], vertical (TL-BL, TR-BR) [R], horizontal [C]. All main strands merge into ONE identity
 (dimension = product), same for anti => every 2-party gadget is  X_{N1,N2} (x) (side strands),
 X_{N1,N2}(x,y,u,v) = [x=v][y=u], x,v in [N1], y,u in [N2]  ("rectangular swap").
Gadgets spanning >= 3 rows or columns: multi-party in that direction (see 3-cycle below).

## 2. Alternatives (task 2) -- see sections below as they are filled

### 2.0 General reduction lemma (BIPARTITE REDUCTION) [P, design level]  (13:10)
Setting: 2-direction layout (rows/cols, any g x h). A gadget = set of RC strands on cells (row class, col class);
its ROW GRAPH has an edge r1-r2 per strand, its COLUMN GRAPH c1-c2. A piece's cost in a direction = total bond bits of a
pairwise network among the classes it touches (Def. 2 admits exactly this: bond index copied at both ends + equality test).
Lemma. If both graphs are bipartite, then h_G(mu) >= h_SWAP(mu) for all mu (support functions of normalized frontiers),
hence alpha*(G) <= alpha*(SWAP).
Proof. 2-colour row classes and column classes. Every strand joins different colours in both => under the 2-party grouping
every strand is RC: G = X_{N^a,N^b} (a = #main strands, b = #anti strands, all in one identity each). A network cost is
>= log of the rank across any grouping (bonds crossing <= all bonds), baselines are equal (each strand costs lam in each
direction in both pictures). So h_G >= h_{X_{a,b}}. Then X_{a,b} (x) X_{b,a} = SWAP_{N^{a+b}} (main strands of both form
one identity, anti likewise), tensor products of piece families average the normalized costs (h of product <= mean of
h's), and X_{b,a} is the left-right mirror of X_{a,b} (exchanges main/anti, keeps both bipartitions): same h. QED
Consequences: rectangular swaps, cyclic permutation in its 2-vs-1 all-crossing flattening (= X_{N,N^2}), all EVEN cycles,
any gadget on 2 rows whose column digraph has bipartite underlying graph: no better than SWAP.
Only gadgets with an ODD cycle in the row graph or column graph escape (smallest: 2 rows x 3 cols, directed 3-cycle).

### 2.1 Rectangular swap X_{N1,N2} (unequal bit groupings of the 2 diagonals) [P + E]
Valid layout (gadget i takes k_i bits of diagonal 1, l_i of diagonal 2). Operator: rows (x,y)|(u,v), cols (x,u)|(y,v),
rank N1N2 = OSR. Frontier symmetric (realignment exchanges y,u, both in [N2]).
Theorem-C type [P]: restrict strands to k1<=N1, k2<=N2: p >= k1 k2 - q min(k1,k2)^2. With N1=N^a <= N2=N^b:
  pq >= N2^2/4 for q >= N2/(2N1)  =>  normalized u+v >= 2b/(a+b) > 1 (STRICTLY above the line in the middle);
  q < N2/(2N1) => p >= N1N2/2: flat at the corner (u = 1-o(1) for v < (b-a)/(a+b)): corner exponent 0.
  (Consistent with the lemma: an anti strand of excess dimension is a lone identity, cost (1,1).)
3-register cyclic permutation (x,y,z)->(y,z,x), "2-vs-1": it is the swap of register x with register (y,z), i.e.
X_{N,N^2}: balanced u >= 2/3 proved, corner flat; upper (1+2u_SWAP)/3. Strictly worse than SWAP's possibilities.

### 2.2 Partial swap (swap of some bits only), side strands [P]
[x1=v1][y1=u1][x2=u2][y2=v2] = X_{N1} (x) (vertical strands TL-BL, TR-BR) = SWAP (x) column side edges (cut by rows only).
Geometrically: bits that do not cross are side edges; in a g x g grid they are a 2/(g+1) fraction (vanishing); in g x h
grids their load is exactly Lemma U's straight-edge terms 1/h, 1/g, already in the alpha* formula. Bundling them into the
gadget: restriction (fix side strands to one vector) maps any bundled (p,q) to a SWAP (p,q); the Theorem-C bound gets
multiplied by the side dimension M (p >= M (k1k2 - q min^2)), the trivial upper bound (p_s M, q_s) matches. "Tensor
power mix" frontier = w F_SWAP + (1-w){(1,0)}: no new trade-off, only a different load split (Lemma U).

### 2.3 Copy crossing [x=v][y=u][x=w] [P]
Copies never occur in the clique instance (each string is held by its two endpoints); a copy is internal to a node table.
If w sits in the cell of v (or of x): w = local injective map of a leg => T ~ SWAP exactly (same frontier).
If w sits in BL (or TR): summing w is a local map in BL => any T decomposition gives a SWAP one with the same (p,q), same
baselines (N^2, N^2) => T at least as hard as SWAP; and T = (copy-merge at TL)(SWAP (x) side strand) => T <= SWAP(x)side.
Creating copies adds strings (extra lam bits on extra edges): worse.

### 2.4 Half-crossings [P]
(i) Routing a diagonal through a corner node (a -> a' -> b'): replaces an RC strand by one C strand + one R strand: cost
(lam, lam) unchanged, but the crossing disappears => loses the only source of saving. Cut argument: any path a ~> b'
crosses both the row cut and the column cut, so a single strand costs >= lam in each direction whatever the routing.
(ii) sqrt(SWAP) = aI + bX, X = H.H: composition needs an intermediate string layer (more bits), and H - aI = bX with I
column-free (I = vertical strands): each half is as hard as X. (iii) swap of half the bits = X_{sqrtN} (x) X_{sqrtN}:
tensor product, same asymptotic frontier.

### 2.5 Hexagonal / triangular lattice (3 line families, nodes i+j+k = 0 mod g) [P]
Two nodes share at most one coordinate. A 2-party-in-every-direction gadget containing an all-cut edge w-w' must lie in
{i,i'}x{j,j'}x{k,k'}: no third lattice point (it would force j=j' or i=i'). So all-cut edges (fraction 1-3/(g+1)) are
isolated single strands for 2-party gadgets, cut identically by all directions: NO gain
(multi-party gadgets there would be D=3 analogues of the odd-cycle case 2.7: not analysed). (2-cut edges form tetrahedra only for even g with
offsets g/2: product of 3 swaps each free in one direction; vanishing fraction.) Useless layout.

### 2.6 Cube / D = 3 (G20a) [R + P]
Body-diagonal gadget X_3 (4 strands, 3 cuts). Each pair of cuts sees X_3 as SWAP_{N^2} (for cuts (0,1): strands {0,1} main,
{2,3} anti; computed from G20a's cut lists) => the D=3 frontier contains the D=2 SWAP_{N^2} frontier on its faces.
alpha* in D dims: same formula with ||mu||_{1/(1-alpha)} on the simplex; step 2 <=> corner at a vertex. Near a vertex,
Theorem G gives nothing (even the line is not proved there); no construction beats the face (SWAP) values (product point
(1,1/2,1/2) normalized). Not reducible by the bipartite lemma (3 cuts), but no evidence of a better corner.

### 2.7 ODD CYCLE: 2 rows x 3 columns, directed 3-cycle P_3 (the only family NOT covered by 2.0) (13:20)
Layout: between two rows a,a' the RC edges form the complete digraph on the g columns (top col j -> bottom col j').
Decompose it into directed 3-cycles (Mendelsohn triple system: exists for g = 0,1 mod 3, g != 6). Gadget = 3 strands
s1: a1(L,top)->b2(M,bot), s2: a2(M,top)->b3(R,bot), s3: a3(R,top)->b1(L,bot).
Operator P_3[a;b] = [b2=a1][b3=a2][b1=a3] = cyclic shift of 3 registers.
 ROWS: 2-party, top (a1,a2,a3) | bottom (b1,b2,b3): matrix rank, full = N^3.
 COLUMNS: 3-party ring L-M-R; a piece's cost = log(b12 b23 b31) for a triangle network (node i holds (a_i,b_i)).
 In column language P_3 IS the matrix-multiplication tensor <N,N,N> (triangle of identities, bonds N).
Validity [P, design level]: ring pieces become pairwise networks (bond index copied at both ends + equality test,
 f_i multiplied into W_col i); per column pair and row pair 2 gadgets (one per orientation) => c_{j,j'} <= g^2 lam (+ factor 2
 if bonds of pieces exceed N; absorbed by P' = 2P, lower order), Lemma 8(ii) OK; Sum of pieces, integrality, poly(lam)
 pieces: as for any piece family. Asymmetric (rows 2-party, cols 3-party): handled by the unbalanced level map (Lemma U/B2).
 The "2-vs-1 flattening" of the task = grouping columns {L,M}|{R}: there P_3 = SWAP_N (s2,s3) (x) side strand s1.
Lower bound [P] (compression): restrict s3 to t dims; L top, R bottom, M top/bottom:
 rank_rows(C) <= N t^2 b31 min(b12,b23), operator rank N^2 t  => (rotate)  p >= N^3 / (4 b_(1) b_(2))  (two smallest bonds)
 => p q^{2/3} >= N^3/4, normalized u >= 1 - (2/3) v  (STRONGER than SWAP's line u >= 1 - v); balanced u >= 3/5.
 p3_lb.py (LP/grid over all compressions and side assignments) gives exactly this min over bond patterns (no better).
 Tight at N=2: bonds (1,1,4): p = 2 = 8/4 [E].
Corner B (u->0, cheap rows) REDUCES to SWAP [P]: kill the strand on the largest bond => path => bipartite => SWAP_N with
 rank <= p, OSR <= b_(1)b_(2) <= q^{2/3}:  (u,v) in F_P3 => (1.5u, <= v) in F_SWAP. A linear corner B of slope c gives a
 linear SWAP corner of slope 1.5c.
Corner A (u->1, cheap ring) does NOT reduce (restriction gives u_s = 1.5u > 1; 2-vs-1 grouping forces ring cost >= N...).
 => "triangle-network rigidity of the cyclic shift": C with ring bonds N^eta, rank(P_3 - C) = N^{3-delta}; known bound
 delta <= 2 eta (q = N^{3 eta}) still allows a LINEAR corner. No construction: C = I gives p = (2/3)(N^3-N) (constant factor), site-tensor-power
 truncations give the usual Gaussian (quadratic) corner [H].
Small exact points P_3 (rows rank p ; ring bonds) [E, flint, mod 1000003 and 998244353]:
 N=2: p=8 (C=0); p=4 @ (1,1,1) [C = I]; p=2 @ (1,1,4) [= SWAP_2 (1,1) x side strand, on the bound]; p=0 @ (2,2,2).
 Numerics (vp.py/cyc_scan.py; ring numerics FRAGILE, failed to find p=0 at (2,2,2) with 8 starts): no p<4 at (1,1,2),(1,2,2)
 [tails 1e-5..1e-8 at p=3: border-like], p=3 @ (1,2,3) numerically.

## 3. Small exact / numerical frontiers (13:25)
Solver calibration (vp.py, variable projection + L-BFGS): square SWAP reproduces F(2,1)=1, F(3,1)=F(3,2)=3, F(4,2)=6,
(4,4)@4 exact, (3,4)@4 border-like (tail 2e-3) -- matches the project's table.
Rectangular X_{N1,N2} (p = row rank, q = OSR); [E] = integer {-1,0,1} solution, ranks over Z and mod 1000003, 998244353
(mc_search.py); symmetric points (q,p) follow by realignment (exchange of y,u), so they are exact too.
 X_{2,3}: (6,0) (3,1)E (2,2)E (1,3)E (0,6)                       bound p >= k1k2 - q min^2: q=1 -> 2
 X_{2,4}: (8,0) (5,1)E (2,2)E (1,5)[realign] ; (4,1) border-like (tail 4e-9), E-search stuck at 5.  [= 3-cycle 2-vs-1, N=2]
 X_{3,4}: (12,0) (6,1)E (5,2)E (4,3)E (3,4) (2,5) (1,6) [realign]
 X_{2,5}: (10,0) (7,1)E (4,2)E (3,3)E (2,4) (1,7) [realign]
 Pattern N1=2: F(1) = 2N2-3 (bound 2N2-4).
Finite-chain support functions (chains.py, B2's log-sum-exp formula): every small frontier (SWAP N=2,3,4,7; X_{2,3},
 X_{2,4}, X_{3,4}, X_{2,5}; P_3 N=2) has h(1/2) > 1 (sum-ranks) and >= 1 (optimistic): useless, as for SWAP (B2).
 Two-piece formula 1/2 + log_{D0}(sqrt p + sqrt q): X_{2,4}(2,2): 1.0; P_3 N=2 (2;1,1,4): 1.09; SWAP octonion (8,8)@8: 0.93.
Copy crossing (w in BL), numerics: N=2: q=1 -> p=2 (SWAP: 1), q=2 -> 1. N=3: q=1 -> 6 (SWAP 3), q=2 -> 5 (SWAP 3),
 q=3 -> 3 (SWAP 1). Strictly worse than SWAP at small N, as the sandwich predicts.
P_3 at N=3..5 [E] (p3_constructions.py): C = I (bonds 1,1,1): p = (2/3)(N^3-N) = 16, 40, 80; C = register transposition
 T13 = [a1=b3][a2=b2][a3=b1] (bonds 1,1,N^2; = SWAP_N (F(N,1)) x side strand): p = (N^3-N^2)/2 = 9, 24, 50 (bound N^3/4).
 Numerics N=3 (cyc_pts.py): bonds (1,1,9): p=9, p=8 plateau 3.0 (optimal); (1,1,2),(1,2,2): 16 (p=15 tails 1e-4..1e-5,
 undecided); larger bonds: solver fails to reach even known values -> ring numerics unreliable beyond N=2.
 mc_search ({-1,0,1} annealing) useless for ring formats at N=3 (does not find known points).
Known ASYMPTOTIC frontier of P_3: only trivial points and constant-factor points (T(u) = 1 for u < 1) => alpha* = 0 from
 known constructions; SWAP-derived constructions (via {L,M}|R grouping) cost ring >= N^2 q_s^2, v >= 2/3: useless.

## 4. VERDICT (13:35)
(b)/(c) reduce to SWAP or worse:
 - rectangular swaps X_{N1,N2} and the 3-register cyclic permutation in its all-crossing 2-vs-1 flattening (= X_{N,N^2}):
   h >= h_SWAP (tensor with mirror = square swap); proven bounds strictly worse (u+v >= 2b/(a+b), flat corner).
 - partial swaps / side strands: SWAP (x) side strands; side load = Lemma U's straight-edge terms; no new trade-off.
 - all even cycles, any gadget with bipartite row and column strand graphs: h >= h_SWAP (lemma 2.0).
 - copy crossing: same as SWAP (copy in v's/x's cell) or at least as hard (other cells); numerically strictly worse.
 - half-crossings: rerouting removes the crossing (cost unchanged, saving lost); sqrt(SWAP) halves contain X.
 - hexagonal/triangular lattice: all-cut edges cannot be grouped into 2-party gadgets => no gain (multi-party not analysed).
Open but no evidence of advantage:
 - odd cycles (P_3, 2 rows x 3 cols; valid layout): proven u >= 1 - (2/3)v (stronger than SWAP's line; balanced >= 3/5),
   corner B => SWAP corner (slope x1.5); corner A = new "triangle-network rigidity of the cyclic shift" (<N,N,N>),
   linear corner not excluded, no construction (known asymptotic frontier trivial: alpha* = 0 vs 1/2 for SWAP).
 - D = 3 cube (G20a X_3): faces = SWAP_{N^2} (so at least as good as SWAP); vertex corner not reducible, lower bounds weaker
   (Theorem G), no construction beyond face/product points.
(a) none: no alternative has both a corner allowed by all bounds AND a construction.
alpha* (B2/R1 formula): for every reducing alternative alpha*(alt) <= alpha*(SWAP) (h larger); for P_3 the known frontier
gives alpha* = 0; no alternative frontier is known that would raise alpha*.
Processes: all mine killed at the end (see report).

## 5. Proof check of the P_3 bound (13:40) [P]
Restrict a3,b1 to [t] (submatrix: local, keeps the ring structure). Operator rank N^2 t. Put L,M on the row side, R on the
column side: crossing wires b1 (t), b2 (N), a3 (t), bonds b31, b23 => rank_rows(C) <= N t^2 b31 b23; M on the column side
=> N t^2 b12 b31. t = N/(2 b31 m), m = min(b12,b23) => p >= N^3/(4 b31 m). Rotating, compress the strand on the smallest
bond: p >= N^3/(4 b_(1) b_(2)) >= N^3/(4 q^{2/3}). Corner-B reduction: restrict the strand on the largest bond to one value
-> path s_i, s_{i+1} -> grouping {ends}|{middle node} -> SWAP_N, OSR <= q/b_(3) <= q^{2/3}, rank <= p.
Session ended ~13:42 UTC; no processes left running.
