# G23c — does a dip of L/sqrt(E) survive with the TRUE pairing rank?  (working notes, saved as I go)
Tags: [P] proved here (short, unreviewed), [E] exact computation, [N] numerical, [H] heuristic, [L] literature from memory, [O] open.
Nothing here is a new algorithm or a new bound for 3-SAT.

## 0. Setting (from G21a/G22e)
Lattice VOA V_L (Gram G, rank r), module M = V_{lambda+L}, cost grading g_xi = L0 - xi(0) (xi in L(x)Q), weights
w_xi(p,k) = |p|^2/2 - <xi,p> + k.  Cost K = max_x [L_M(x) + L_M(B-x)], B = 2E0 + H, all in g_xi.
xi in the closed Voronoi cell of 0 => V^{xi}_{<0} = 0 and V^{xi}_0 = span{e^g : g in S_xi},
S_xi = {g in L : |g|^2/2 = <xi,g>}  (lattice points on the sphere through 0 centred at xi).

## 1. True pairing rank: channels  [P + E]   (trank.py, test1-3.py)
K = 0 forces c_xv = a_x(H-1)b_v in V^xi_0 = span{e^g : g in S_xi}  (all a at g_xi-weight ha, all b at hb, ha+hb = H;
locality a(n)b = 0 for n >= H is then automatic in the Voronoi cell).  Write c_xv = sum_g B_g(a_x,b_v) e^g, O_g = o_xi(e^g).
(A3) reads   sum_{g in S} B_g(x,v) <u|O_g|y> = d [x=v][y=u]        ("channel decomposition").
 - B_g(a,b) = coefficient of e^g in a(H-1)b; exact Wick formula for lattice VOAs (trank.Bval), cocycles drop out of ranks.
 - O_g maps charge sector p -> p+g, preserves g_xi-weight; O_g (g != 0) raises t := <xi, charge> by |g|^2/2 > 0: NILPOTENT.
Bounds [P]: N^2 <= sum_g rk(B_g|) rk(O_g|),  N <= sum_g min(R_g, Lab_g),  N <= dim M_E0.
 Grading transfer [P]: R_g(H) <= 2 D_V^{(g/2)}(H/2) (G22e lemma in grading L0 - (g/2)(0), where e^g has weight 0);
   Lab_g(E0) <= D_M^{(xi-g/2)}(E0) (min(dim M^p, dim M^{p+g}) <= dim at the average weight, which is the g_{xi-g/2} weight).
   Channel g splits xi = g/2 (pairing) + (xi - g/2) (labels).
VALIDATION [E]: App. B (Heisenberg, H = 2h): R_0(h,h) = Lab_0(h) = p(h) for h = 1..6 (N = p(h), cost 2 log P(2h) as in the paper).
R1 = A1, xi = alpha/2, S = {0, alpha}; exact ranks, h, E0 <= 7 (test2.py):
   R_0(h,h)  = 1,4,5,9,13          = Sym(h)  = coeff of sum q^{n^2+|n|}/eta   (G22e's "true pairing")
   R_a(h,h)  = 2,6,8,14,20 (h=1..5) = Dip(h) = dim V'_h  -- the FULL loose count: channel alpha pairs sector m with 1-m,
                                                            which have EQUAL g_xi-weights (self-conjugate in the cost grading)
   Lab_0(E0) = Dip(E0) (o(1) = Id),   Lab_a(E0) = Sym(E0) (blocks n -> n+1 have full rank min(dims)).
   B_0 + k B_a (k = 1, 3/7) is NONDEGENERATE on V'_h x V'_h, h = 0..6 (test3.py).
 So the 3.030 of G22e is the channel-0 value; it is not forced by a selection rule: channel alpha has the full pairing rank.
 What happens instead is a TRADE-OFF: in channel 0 the pairs are Sym and the labels Dip; in channel alpha the reverse.

## 2. Dichotomy theorem (the actual mechanism)  [P]
Let Y = span(labels), U = span(dual labels), P = projection onto Y along U^perp, Om_g = P O_g |_Y (N x N), J = span{Om_g : g != 0}.
 (i) If Id not in J: a functional tau on End(Y) with tau(Id)=1, tau(J)=0 applied to (A3) gives B_0(x,v) = d [x=v]:
     N <= R_0(H) = vacuum-channel rank <= 2 D_V^{conf}(H/2)   (charge-neutral pairing, CONFORMAL grading).
 (ii) If Id = sum k_g Om_g:  A := sum k_g O_g is nilpotent (all O_g raise t), so (A-1) is invertible, (A-1)Y c U^perp and
     Y cap U^perp = 0 force Y cap (A-1)Y = 0:  2N <= dim M_E0;  also N = rk(P A|_Y) <= sum_{g: k_g != 0} Lab_g.
 Single copy: N <= max( min(R_0, Lab_0), min(floor(D/2), sum Lab_g) ).  Both branches are REALISED for A1 (sect. 3).
 General adversary form [P]: for every tau in K^perp = {(Tr(T O_g))_g : T in Y(x)U}, tau(1) = 1:  N <= rk(sum_g tau_g B_g);
 tau may vanish on a set G' of channels iff Id is not in span{Om_g : g in G'}.  So N <= R_0 + min over hitting sets.
 Branch (ii) costs only a factor 2 in total, NOT per copy: for tensor-power families this theorem alone does not reach
 the per-copy rate.  It kills R1 for every FIXED algebra, but the family question stays open (sect. 4).

## 3. R1 family values, per copy (fam.py, chan.py) [N, Legendre entropies, eps <= 3.5]
 (d) loose (labels Dip, pairs Dip)                                   2.8157   <- upper bound on N, all my rigorous bounds reduce to it
 (a) channel 0 only (labels Dip, pairs Sym) [achievable]             3.0300
 (a') channel alpha only (labels Sym, pairs Dip) [achievable]        3.0300
 (c) branch (ii) per copy, tensored (labels Dip/2 per copy)          3.1129
 single channel g with a fraction pi of alpha's, n copies, RIGOROUS (Conf bounds): min over pi 2.902 (pi = 1/2), edge of grid
   (same with Sym in place of Conf: 3.014).
 Channel sums: N <= sum_g min(R_g, Lab_g) has rate max_pi [H(pi) + min(r_pi, l_pi)] >= rate of Dip, because Dip <= 2 Sym
 coefficientwise; so counting bounds alone CANNOT separate R1 from the loose 2.8157.  [P]

## 4. Scan (scan.py -> scan_r1.out, scan_r2.out, scan_r3.out; twist.py -> twist.out)  [N, per-copy Legendre rates]
Columns: LOOSE = labels+pairs in the cost grading (upper bound on N for ANY construction); B0 = vacuum channel (pairs in
conformal grading; rigorous whenever S_xi = {0}, and branch (i) otherwise); SC = best single channel g (or 0/g mixture),
rigorous Conf-type bounds.  Minima "at the grid edge" (e ~ 3r) = approach to 2sqrt2 from above, not a dip.
 lattice/module          xi (basis coords)     |S|  LOOSE   B0      SC
 A1 V                    0, .125, .25          1    2.91-2.89 2.91-2.93  -
 A1 V                    .375, .45 (interior!) 1    2.819, 2.812   2.946, 2.945   -     <- dips killed by charge neutrality
 A1 V / V_{a/2+L}        .5 (= R1)             2    2.817   2.932   2.908
 Z (odd, VOSA)           .5                    2    2.838   2.859   2.838
 A2 V                    0, generic            1    2.89-2.88 2.90-2.91  -
 A2 V, V_{w1+L}          root/2 [.5,0],[.5,.5] 2    2.809   2.926   2.893
 A2 V, V_{w1+L}          deep hole w1, w2      3    2.815   2.935   2.898
 A1+A1                   (.5,0) / (.25,.5) / (.5,.5)  2/2/4  2.799/2.809/2.817  2.942/2.950/2.932  2.909/2.917/2.908
 A3                      [.5,0,0], [.5,.5,.5]  2    2.803   2.906   2.884
 A3                      [.5,1,.5] (|S|=6), [.25,.5,.75] (|S|=4)   2.819, 2.815   2.920, 2.916   2.884, 2.884
 twisted sectors (fixed-point-free g, S = {0}, true bound = B0):  A1 -1: 2.926; A2 -1: 2.953; A2 order 3: 2.908;
   D4 -1: 2.913; D4 order 3: 2.916; E8 -1: 2.853; E8 order 3: 2.852 (all at grid edge; "own anchor" values 2.83-2.89 are NOT bounds).
 Not run: D4/E8 regradings (rank-4 lattice sums too slow in the time box), non-abelian orbifolds.
 NO value below 2sqrt2 with a true-rank bound.  LOOSE < 2sqrt2 occurs iff xi is near/on the Voronoi boundary.

## 5. Barrier / escape  [P + N]
 (a) SELECTION RULE (proved): if xi is in the INTERIOR of the Voronoi cell (S_xi = {0}), c_xv must be a multiple of the
     vacuum, charge conservation pairs sector p with -p, and the vacuum pairing rank is grading-free: R_0(H) <= 2 D_V^conf(H/2).
     The regrading then only moves the label/cost profile: N <= min(2D_V^conf(H/2), D_M^xi(E0)) (B0 column, >= 2.90 here).
     This is what kills R1-type dips at generic xi (A1 at xi = .45 alpha: loose 2.812 -> true >= 2.945).
 (b) On the Voronoi BOUNDARY (S_xi != {0}) the rule fails: non-neutral channels e^g (g in S_xi) appear in V^xi_0 and channel g
     has the full loose pairing rank (exact: A1 R_alpha = Dip).  What remains is the dichotomy of sect. 2: single channel
     g -> pairing in grading g/2, labels in grading xi - g/2, never both in xi (SC >= 2.88);  multi-channel -> factor-2
     label loss (branch ii), harmless at the rate level.
 (c) So the precise open question (the only escape left for lattice regradings): a MULTI-CHANNEL, non-product construction in
     a tensor-power family (e.g. A1^n or (A1+A1)^n at xi on the Voronoi boundary) with N at the loose rate.  Requirement:
     sum_g B_g (x) Om_g = I (x) I with Om_g compressions of the nilpotent O_g, and every "hitting set" of channels heavy in R.
     Upper bound 2.7989 (A1+A1, xi = (alpha1/2, 0)), best realised 3.030 (A1, channel 0 or alpha).  [O]

## 6. Verdict
 Best TRUE value: nothing below 2sqrt2 = 2.8284.  Realised constructions: App. B (-> 2.8284 from above); R1 channel 0 / alpha: 3.030.
 True-rank lower bounds on the cost: generic regradings and twisted sectors >= 2.85-2.95 (grid); boundary regradings,
 single channel >= 2.88.  Proved: the selection rule (interior xi) and the dichotomy (boundary xi).  Open: multi-channel
 non-product constructions on the Voronoi boundary in tensor-power families (gap 2.7989 .. 3.030).
 Adversarial: per-copy values are Legendre rates on a finite grid (eps <= ~3r); SC uses avg (Conf) bounds which are rigorous
 but not tight (exact per-copy ranks use Sym, giving 3.01 for A1); twisted ground dims for D4/E8 order 3 taken from
 |L/(1-g)L|^{1/2} [L, from memory]; branch (ii) single-copy realisability argued [P], not verified by an exact tensor check.
 Next step: decide (c) — either a rate-level strengthening of the dichotomy (e.g. a t-filtration argument showing that the
 compressions Om_g of heavy channels cannot express Id with N at the loose rate), or an explicit 2-copy A1^2 search
 (labels at E0 <= 4, exact ranks with trank.py) for N > max(single-channel, product) values.
Files: trank.py (exact evaluator), test1-3.py, fam.py, chan.py, scan.py, twist.py, *.out, this notes.md.
