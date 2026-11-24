# G18a — Can the 2-adic precision of App. D be cut by more than a constant using the actual kernels?

Agent notes (AI; not independently reviewed). Status tags: [PROVED] = short proof written here,
[NUM] = exact computation on toys, [HEUR] = heuristic, [OPEN].
Source read: paper.txt lines 3257-3660 (App. D.4). src_paper.tex in the dl folder is a 404 stub
(14 bytes), so the TeX source was NOT available; I used paper.txt.

## 1. Definition (task 1)

Setting of one call: grid S={-L..L-1}, L=2^a, alpha=2^{a+1}, w_d(T)=C(T)/(C'(0)(T-d)).
Packets P_1..P_k, packet j has real slot set Sigma_j in {0..L-1}; M = sum |Sigma_j|, N=M(2L-1).
Kernel K_T(Y,X) = prod_packets prod_{x in Sigma_j} l_{Y_x}( sum_{y in Sigma_j} X_y w_{y-x}(T) ).
F(T) = sum_{X,Y} A(Y)B(X) K_T(Y,X), F(0) = sum_X A(X)B(X).

 * Kernel lattice  Lam_ST := Z_(2)-span of { K_T(Y,X) } inside Q[T]_{<=N}.
   (With A,B arbitrary integral tensors F ranges exactly over Lam_ST tensor the coefficient ring;
   Smith invariants do not change under the unramified base change, so Z_(2) suffices.)
 * Black-box lattice Lam_BB := { F : deg<=N, F(t+alpha z) in Z_(2)[z] for all t in S } (Lemma D.5's class).
   Lam_ST is contained in Lam_BB (Lemma D.4).
 * Data actually available to the algorithm: App. D interpolates the ENTIRE polynomial F(t+alpha z)
   mod 2^{p_t} at each active node (v2(t)<u), so data_t(F) = all N+1 coefficients of F(t+alpha z).
   (The multiplicities q_t only enter formula (D.9); more coefficients can only lower the need.)
 * A precision vector (p_i)_{i<u} is FEASIBLE for a lattice Lam and target modulus 2^p iff
      for all F in Lam: [ data_t(F) = 0 mod 2^{p_{v2(t)}} for all active t ]  =>  F(0) = 0 mod 2^p.
   This is exactly "F(0) mod 2^p is determined by the truncated node data" on a linear class.
   Equivalently (p=1): no WITNESS F in Lam with F(0) odd and data_t(F) divisible by 2^{p_{v2(t)}}.
   need(Lam; shape) = 1 + max depth of a witness.  Feasible set is up-closed; "minimal (p_t)" is its
   Pareto frontier.  Cost-relevant scalar (count = exp(c sqrt(p_t * 6M/D_t))): lambda* = min max_i p_i/D_i.
 * Remark (non-linear class): if A,B are restricted to products of bounded-degree vertex tensors the
   class is a subset of Lam_ST, and determinacy on it can only need less.  With one hub vertex per
   side the class is all of Lam_ST.  I compute the linear (lattice) version. [definition, not a claim]

Computation: Smith normal form over Z_(2) of the (basis x data-coordinates) matrix with the target
F(0) carried along: if pivots have valuations v_i and transformed targets w_i then
   need(p) = max{ v_i + p - v2(w_i) : v2(w_i) < p }.   [PROVED: elementary, see g18lib.py]

## 2. Structural facts worked out before computing

 (S1) w_d = (C'(d)/C'(0)) * (Lagrange basis poly of node d on S), C'(d)/C'(0) a 2-adic unit.
      On the disc t+alpha z (t != d):  w_d(t+alpha z) = 2^{a+1} z G_t(z) / (t-d+alpha z), G_t unit poly,
      so its content is exactly 2^{a+1-v2(t-d)}.  [PROVED, = Lemma D.4 + rev1 R1]
 (S2) Boolean reduction X_y^2=X_y: Lam_ST of a packet is spanned by
      P_{R,Y'} = sum_{phi: R ->> Y'} prod_{x in R} w_{phi(x)-x}  (R subset of rows, Y' subset of columns),
      NOT by all monomials in the w_d.  Only differences d in Sigma-Sigma occur. [PROVED]
 (S3) Product structure: Lam_ST(P_1 u P_2) = Lam_ST(P_1) * Lam_ST(P_2); hence witness depth is
      SUPER-additive over packets: depth(P1 u P2) >= depth(P1)+depth(P2). [PROVED]
 (S4) Extreme sparse case, one wire per packet: Lam_ST = span{ w_0^k, k<=M }, and
      w_0(t+alpha z) = 2^{a+1-i} z * unit(z); the map c -> coefficients is triangular in z, so the
      witness depth at a level-i node is exactly <= M(a+1-i), versus Theta(M D_t/u) = Theta(M 2^{a-i}/u)
      black-box.  So the gain from structure CAN be unbounded (factor 2^{a-i}/(u(a+1-i))) -- but in that
      instance all wires cross at every node and the tuple enumeration is 2^M anyway. [PROVED]

## 3. Computations (task 2)  -- code: g18lib.py, run1.py ... ; raw outputs *_out.txt

Sanity (sanity_out.txt): my independent code reproduces G17a: a=1 [01]: 10; two packets: 19;
a=2 [01]: BB 9 (G17a with q=(6,6): 9) vs ST 3; [012]: 13 vs 12; [013]: 13 = 13.
NEW already here: (i) Lam_ST is NOT saturated in Lam_BB (index 2^5 for a=1 [01], 2^74 for the full a=2
packet) although the needs coincide when the Q-rank is full; (ii) for the sparse packet the need of the
saturation of Lam_ST (5 or 7) is strictly between ST (3 or 4) and BB (9 or 25).  So neither "Q-span" nor
"black box" explains the actual need: the INTEGRAL structure matters.

### 3.1 Raw observations so far (before the zero-threshold fix; ST and BB needs unaffected, see 3.2)
Uniform shape p_i = P (all active levels), target mod 2, full jets.  "need" = minimal P.
 a=1 (u=1): [01]^k: ST = BB = 9k+1 (k=1..5).  [0]^k: ST = 2k+1 = M(a+1)+1 vs BB 5,10,14,19,23.
            [01]+[0]: 14=14, [01]+[0]+[0]: 19=19, [01]^2+[0]: 23=23  (one 2-wire packet already restores BB).
 a=2 single packet (u=1 | u=2 uniform | u=2 balanced m, p=(2m,m)):   BB  vs  ST
   [0]      12 | 5 | 4      vs   4 | 3 | 2          rank 2/8
   [01]     25 | 9 | 7      vs   4 | 3 | 2          rank 7/15   (G17a's "4 vs 14" case)
   [02]     25 | 9 | 7      vs  12 | 4 | 2          rank 7/15
   [03],[12] as [01]
   [012]    37 | 13 | 10    vs  28 | 12 | 9         rank 19/22
   [013],[023]  37|13|10    vs  37 | 13 | 10        rank 22/22
   [0123]   50 | 17 | 13    vs  50 | 17 | 13        rank 29/29
 a=2 several packets:
   [0123]^2: 99|35|26 vs 99|35|26.
   [01]^k, k=1..5:  u=1: ST 4,29,61,88,110 vs BB 25,50,74,99,123  (gap tends to an additive ~13)
                    u=2 uniform: ST 3,7,11,15,19 (=2M-1) vs BB 9,17,26,35,43 ; balanced 2,7,11,15,19 vs 7,13,19,26,32
   [02]^k, k=1..4:  u=1: ST 12,24,48,73,(97) vs BB 25,50,74,99,(123); u=2 unif 4,12,18,27 vs 9,17,26,35
   [0]^k: u=1: 3k+1 ; u=2 uniform: 2k+1  (= M(a+1-i)+1, as proved in S4) vs BB 12,25,37,50,61,74 | 5,9,13,17,22,26

### 3.2 After the fix (rank/sat columns now right; K=2048 vs K=4096 give identical needs: run3_kcheck4096_out.txt)
Bug found by self-review: triangular insertion mod 2^K left junk rows of valuation ~K (inflated ranks,
e.g. "21/29" for [01]^2 is really 19/29; sat(ST) for those rows was garbage).  Fixed with a zero
threshold 2^{K/2}; ST and BB needs were not affected (re-run identical).
 a=2 mixes (BB vs ST;  u=1 | u=2 uniform | u=2 balanced):
   [01]^2  rank 19/29: 50|17|13 vs 29|7|7       [01]^3 rank 37/43: 74|26|19 vs 61|11|11
   [02]^2  rank 19/29: 50|17|13 vs 24|12|8      [02]^3 rank 37/43: 74|26|19 vs 48|18|13
   [01]+[02] 21/29: vs 31|10|6    [01]+[03] 21/29: vs 36|7|7    [01]+[23] 19/29: vs 29|7|7
   [012]^2 rank 41/43: 74|26|19 vs 68|24|18
   [01]+[02]+[03] rank 43/43: 74|26|19 vs 74|26|19   <- three 2-wire packets with all differences: = BB
   [0123]+[01] rank 43/43: 74|26|19 = ;   [0123]+[0] rank 36/36: 61|22|16 = .
 a=3 single packets (u=1|u=2|u=3 uniform), BB vs ST:
   [0] 27|11|5 vs 5|4|3 ; [01] 56|22|9 vs 5|4|4 ; [04] 56|22|9 vs 16|11|4 ;
   [012] 84|32|13 vs 15|7|6 ; [013] vs 14|7|6 ; [025] vs 15|7|6          (rank 19..22 of 46)
   [0123] 113|43|17 vs 65|25|7 (rank 44/61) ; [0136] vs 79|29|7 (51/61) ; [0246] vs 67|27|. (44/61)
   [01234] rank 70/76: u=1: 140 vs 125 (run cut by the 10-minute limit before u=2,3)
   [0145] 113|43|17 vs 80|21|12 (rank 49/61)
 a=3 several sparse packets [01]^k (k=1..4, M=2k), uniform need, BB vs ST:
   u=1: 56,113,169,226  vs  5,13,22,97      (rank 7/31, 19/61, 37/91, 61/121: grows quadratically in k)
   u=2: 22,43,64,85     vs  4,10,16,22  = 3M-2     (level-1 nodes +-2: all differences odd)
   u=3: 9,17,26,35      vs  4,7,11,15   = 2M-1     (level-2 nodes +-4)
 a=2 [012]^3 rank 62/64: 110|39|29 vs 104|37|27.

## 4. The law (task 3)

Notation: Delta = union of difference sets Sigma_j - Sigma_j; e_d(t) = a+1-v2(t-d) = 2-adic content
exponent of w_d on the disc of node t (S1); N+1 = M(2L-1)+1 = dimension of the black-box class.

 (L1) [NUM, 17 instances incl. a=3 [012457] with u=2,3; no counterexample] FULL RANK => NO GAIN AT ALL.  Whenever Lam_ST has Q-rank N+1,
      need_ST = need_BB exactly, for uniform and for balanced shapes (a=1: [01]^k k<=5, [01]+[0],
      [01]^2+[0], [01]+[0]^2; a=2: [013], [023], [0123], [0123]^2, [01]+[02]+[03], [0123]+[01],
      [0123]+[0]).  This is NOT because Lam_ST = Lam_BB: the index is 2^5 .. 2^285.  The deep witnesses of
      the black-box class survive in the sublattice.  UNPROVED.
 (L2) [PROVED for one-wire packets (S4); NUM otherwise] CUT NODES.  If at an active node t every nonzero
      difference d in Delta has v2(d) < v2(t) (then W(t)=0 on all packets: every wire is cut at t, and
      w_d = 2^{>=1} * w_0 * unit on that disc), the witness depth at t is <= M(a+1-v2(t)), instead of
      Theta(M 2^{a-v2(t)}/u).  Observed exactly: [0]^k: M(a+1-i)+1;  a=2 [01]^k, u=2: 2M-1;
      a=3 [01]^k: u=2: 3M-2, u=3: 2M-1.  This is G17a's unexplained "4 vs 14" (one packet {0,1} in L=4:
      Delta={-1,0,1}, the level-1 nodes +-2 are cut nodes, bound M(a+1-1)=4).
      Proof attempt: F = sum_m w_0^m (c_m + 2 B_m), F(0) = sum c_m; triangular in the order of vanishing
      for m_0 <= 2; for m_0 >= 3 my argument needs an independence property of the powers of the unit
      series Gamma_t(z) that I did not prove.  So: UNPROVED beyond one-wire packets.
 (L3) [NUM] RANK-DEFICIENT, NO CUT NODE: additive gain only.  For k copies of a deficient type the rank is
      N+1-delta with delta bounded (a=2: [01]^k: 6, [02]^k: 6, [012]^k: 2) and BB-ST stays bounded:
      a=2,u=1: [01]^k: 21,21,13,11,13;  [02]^k: 13,26,26,26,26;  [012]^k: 9,6,6.  Ratio -> 1.
      While the rank is far from full (a=3 [01]^k, rank ~ k^2 vs 30k+1) the need is small (5,13,22) and it
      jumps (97 at k=4) as the rank fills: the gain is a FEW-PACKETS effect.
      Algebraic reason [PROVED part]: the graded algebra generated by one packet's span V_1 is the
      coordinate ring of a rational curve of degree n_1 = sigma(2L-1) (V_1 contains 1, is base-point free,
      and is birational as soon as it contains two different w_d), so dim V_1^k = k n_1 + 1 - delta for
      large k (Hilbert polynomial); only one-wire packets (V_1 = <1,w_0>, a degree-(2L-1) cover of a line)
      keep rank M+1.  [UNPROVED part: that the depth follows the rank.]
 Scaling of the black-box optimum itself (= actual optimum for full packets), uniform shape:
      a=1,u=1: 4.5 M;  a=2: u=1: 12.3 M, u=2: 4.4 M;  a=3: u=1: 28 M, u=2: 10.6 M, u=3: 4.4 M
      i.e. need / (M D_{u-1}/u) = 2.3 | 3.1, 4.4 | 3.5, 5.3, 6.6  -> Theta(M D_t/u), linear in M (exactly
      additive in the number of full packets: 9k+1; 50,99; 17,35; balanced 13,26), constant between the
      proved witness bound (~2) and the LP value (~5.9) asymptotically.
      Paper's stated increment 8hD_t+8M(a+1) is 10-15x the optimum on these toys (a=2,u=2,[0123]^2:
      (577,385) against the balanced optimum (52,26)); the exact need of formula (D.9) with the paper's
      taper is 1.3-2.7x the optimum ((142,35) against (52,26); a=1 [01]^5: 89 against 46).
 (L2') [NUM] ONE NON-CUT PACKET RESTORES THE BLACK BOX (run4_out.txt).  a=2, u=2 uniform:
      [01]^3: 11 vs BB 26;  [01]^2+[02]: 24 vs 26;  [01]^3+[02]: 33 vs 35;  [0123]+[01]: 26 = 26.
      So the cut-node gain needs EVERY packet to be cut at that node.

## 5. Running time (task 4)

Count at an active node (D.15): variables: n_t = 6M/D_t modes per weight w >= 1; tuples with
sum w beta < p_t: ln(count) <= min_theta [theta p_t + n_t zeta(2)/theta] = 2 sqrt(zeta(2) n_t p_t)
= 2 pi sqrt(M p_t/D_t).  With p_t = kappa M D_t/u:   count = exp(2 pi sqrt(kappa) M/sqrt(u)).
 * Paper: kappa = 24 (8hD_t, h=3M/u; the terms p and 8M(a+1) are lower order because D_t >= 2^u).
 * Black-box optimum (G17a + this note): kappa* in [~2 (product witness, proved), ~5.9 (LP, attained by
   an LP-optimal taper)]; toys: 2.3-6.6 (uniform shape, not the balanced one).
 * Actual kernels, App. D's packets (L slots, >= D_arr = 64^{2a} packets, average filling >= 1/3, all
   differences present): by (L1)/(L2')/(L3) the actual optimum equals the black-box optimum up to an
   additive O(1) per packet type.  [NUM/conjectural]
 => re-optimised exponent: 2^{c M/sqrt(u)} with c smaller by a factor sqrt(24/kappa*) in [2.0, 3.5];
    the kernel structure adds NOTHING to that factor.  The recursion sum_d M_d/sqrt(u_d), M_d = Cs 8^-d,
    u_d = eps*lambda/2^d, is geometric and stays Theta(s/sqrt(lambda)); the leaf size H = s/sqrt(lambda)
    is untouched.  Final bounds keep their form: 2^{O(n loglog n/sqrt(log n))} in general and
    2^{O(n/sqrt(log n))} at bounded density; only the hidden constant of the non-leaf part drops (<= 3.5x).
 * What a real gain would have bought (NOT available): if p_t = M poly(a) were enough with the same mode
   supply, the count would be exp(O(M sqrt(a/D_t))) <= exp(O(M sqrt(a) 2^{-u/2})), i.e. 2^{M/2^{Omega(u)}},
   a polynomial saving in s^{eps'}.  The only instances where p_t = M(a+1-i) suffices are the cut-node
   ones (L2), and there the mode structure is trivial: every tile holds one real column, each wire is its
   own weight-(a+1-i) mode, the budget allows about M of them, and the enumeration is all 2^{Theta(M)}
   subsets.  [HEUR: "precision x supply" does not drop; shown only in the two extreme regimes.]

## 6. Verdict, obstruction, refuting experiment (task 5)

VERDICT: no gain beyond a constant, and for App. D's own packets no gain at all over the black-box
optimum.  The only improvable thing is the constant kappa (24 -> between 2 and 5.9).
OBSTRUCTION (precise, conjectural):
 Conjecture A (full rank => black box).  If rank_Q Lam_ST = M(2L-1)+1 then for every level set and every
   precision vector, (p_i) is feasible for Lam_ST iff it is feasible for Lam_BB.
 Conjecture B (one non-cut packet suffices).  If some packet has a nonzero difference d with
   v2(d) >= v2(t) for every active t, and the number of packets is >= L, then
   need_ST >= need_BB - O_{a}(1), hence need_ST = Theta(M D_t/u).
 Together with G17a's witness bound these give: every recovery rule that is a function of the truncated
 node data of App. D needs p_t >= c M D_t/u at some active node, even knowing the kernels.
 What is PROVED: S1-S4, the Hilbert-polynomial rank statement, the formula need = max(v_i+p-v(w_i)).
 What is NUMERIC: L1 (17 instances, a<=3), L2 beyond one-wire packets,
   L2', L3.
SMALLEST EXPERIMENT THAT WOULD REFUTE THE OBSTRUCTION (i.e. reopen a super-constant gain):
 (E1) one full packet at a=3 (M=8, N=120), u=3: need_ST < need_BB would refute Conjecture A.
      Not finished here (see 3.3); the six-slot substitute [012457] did NOT refute it.
 (E2) a=4, u=2 (the paper's a=2u), two full packets (M=32, N=992): compare slopes of need_ST and need_BB
      in M.  Not feasible with this pure-Python code inside the 10-minute limit; needs flint/C.
 (E3) for Conjecture B: [01]^k+[02] at a=3 with k up to 8 (rank full), u=2,3.

## 7. Adversarial self-review
 - "Actual kernel span" is the LINEAR span (arbitrary A,B).  A recovery rule that exploits the product
   structure of A and B (non-linear class) is outside what I computed; with one hub vertex per side the
   class is the full lattice, but for bounded-degree networks this is a real gap.  UNADDRESSED.
 - Target modulus 2 only (p=1); need(p) = max(v_i + p - v(w_i)) can differ from need(1)+p-1 if some
   transformed target is even.  Not tabulated.
 - Uniform and balanced (p_i proportional to D_i) shapes only; the full Pareto frontier was not computed.
 - Toys have a <= 3, u <= 3, M <= 12: "Theta(M D_t/u)" for u large is extrapolated from G17a's witness, not
   measured here.  The paper's regime a = 2u (cumulative) starts at a=4.
 - L3's "additive gap" rests on k <= 5 copies.
 - All arithmetic is in Z/2^2048 with a zero threshold 2^1024; checked against Z/2^4096 on 4 instances
   and against G17a's independent Fraction code on 6 instances.  A junk-row bug was found and fixed.
 - The src_paper.tex given as input is a 404 stub; I read paper.txt.  App. D's correctness not assessed.
 - Section 5's "precision x supply does not drop" is heuristic; I did not exhibit a lower bound for the
   number of rank-one terms mod 2^{p_t}.

### 3.3 a=3 with nearly full packets (experiment E1)
 - One FULL packet [01234567] (M=8, N=120): a random sublattice of Lam_ST reaches rank 121/121 after
   2000 generators (run3_a3full_partial.log), so the full packet has full rank.  The need computation did
   not finish inside the 10-minute limit (twice; the machine gives ~1 effective core).  NOT DONE.
 - Replacement: six of eight slots, N=90, sampled sublattice (need of a sublattice is a rigorous LOWER
   bound for need_ST, need_BB a rigorous upper bound): run5_*_out.txt.  RESULT (uniform shape):
     [012457]  rank 91/91 (full):  u=3: 26 <= need_ST <= BB = 26;   u=2: 64 <= need_ST <= BB = 64.
     [012345]  rank >= 87/91:      u=3: 25 <= need_ST <= 26;        u=2: 58 <= need_ST <= 64.
   So (L1) holds exactly in the first instance with a genuine 3-level taper (a=3, u=3), and a
   rank-deficient 3/4-full packet is within 1 (u=3) resp. 6 (u=2) bits of the black box.
   E1 for the full 8-slot packet itself remains NOT DONE (only its rank, 121/121, is established).

## 8. Files
 g18lib.py (library), run1.py..run5.py (drivers), mktable.py -> table.md (127 rows: formula, BB, ST),
 *_out.txt raw outputs, run3_a3full_partial.log.
 Caveat on files: run2_*_out.txt were produced before the junk-row fix; for [01]^k (k>=2) their rank and
 sat(ST) columns are wrong (true ranks 19, 37, 51, 65), the BB and ST needs are right (k=2,3 re-run).
