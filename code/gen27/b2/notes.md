# G27-B2 notes (start 07:38 UTC)

## Reading (07:40-08:00)
- Claude paper Prop 3: level = g x g grid, per rectangle crossing X=[x=v][y=u] expanded into pieces A_i, per term iota
  merge ALL rows or ALL columns (global choice), cost min(U_row,U_col) <= average -> rho = sup (c_row+c_col)/(4 lam).
- Appendix A: K_{l+1} ~ K_l^2, e_{l+1} = rho e_l -> alpha' = log2(1/rho).
- G25-P9 (Prop P, Lemma U): unbalanced g x h grids; Phi=1/e^2 additive in x=log K for App B (circle). Says power change needs
  frontier "strictly inside the quarter circle for some weighting". Did NOT compute the exact alpha for a general frontier.

## Q1 analysis (08:00-)
Level map (Lemma U, design level): e' = max_{(u,v) in H} min(e_a u, e_b v), x' = x_a + x_b, H = convex hull of normalized
piece costs (u,v) = (theta_s/2, T_s/2) (worst iota = average over rectangles).
Self-similar ansatz e = c x^{-alpha}, split w: consistent iff H misses open quadrant {u > w^a, v > (1-w)^a}.
Separation (H convex, quadrant open up-set) => iff exists mu in [0,1]: h(mu) := max_H (mu u + (1-mu) v) <= mu w^a + (1-mu)(1-w)^a.
Max over w (Hoelder) = ||(mu,1-mu)||_q, q = 1/(1-alpha).
**FORMULA**: alpha* = sup{ alpha : exists mu, h(mu) <= ||(mu,1-mu)||_{1/(1-alpha)} }.
Checks: circle h = l2 norm for all mu -> alpha* = 1/2 (= Prop P). Balanced mu=1/2: h(1/2)=rho, ||.||_q = 2^{-alpha} -> alpha=log2(1/rho).
Theorem C (pairs: max_A theta + max_B T >= 2 for every split of pieces) => h(mu) >= max(mu,1-mu) => alpha* <= 1.
**alpha* = 1 iff inf_mu [h(mu) - max(mu,1-mu)] = 0.** For mu = c/(1+c) != 1/2 this is the CORNER condition
   T(theta) <= c (2 - theta) for theta in [2-2/c, 2]  (finite slope at the corner (2,0)), NOT the line.

### Numerical checks (dp_alpha2.py = DP over unbalanced Lemma-U trees with exact hull via support function;
### alpha_formula.py = self-similar formula)
| family (normalized u=theta/2, v=T/2) | balanced rho | log2(1/rho) | DP fitted alpha (x<=1e22) | formula alpha* (mu<=0.95) |
| circle (App B) | 0.7071 | 0.500 | 0.500 | 0.500 |
| line | 0.5 | 1.0 | 1.000 | 1.0 |
| l_1.5 sphere | 0.630 | 0.667 | 0.667 | 0.667 |
| circle + linear corner slope 4 (v=4(1-u) near (1,0)) | 0.7071 | 0.5 | 0.722 (pre-asymptotic) | 1.0 |
| circle + linear corner slope 10 | 0.7071 | 0.5 | 0.551 (pre-asympt.) | 1.0 |
| circle + power corner v ~ (1-u)^0.75 | 0.7071 | 0.5 | 0.611 | 0.683 (-> 0.75 as mu->1) |
=> SAME balanced rho as App B, but a different corner raises the power of log. Approach to alpha* is slow (needs
   splits x_b/x_a ~ c^{-1/(1-alpha)}).

**Theorem (design level, model M of G25-P9 with general piece family)**:
 alpha* = sup_{mu in [1/2,1)} alpha(mu), alpha(mu) = sup{alpha : h(mu) <= ||(mu,1-mu)||_{1/(1-alpha)}};
 alpha(1/2) = log2(1/rho); lim_{mu->1} alpha(mu) = gamma_corner where v_max(1-u=t) ~ t^gamma near (1,0).
 Hence alpha* >= max(log2(1/rho), gamma_corner); alpha* <= 1 (Theorem C).
 Step 2 (alpha* = 1) <=> exists c<inf: every piece has T_s <= c (2 - theta_s) (+o(1)), i.e. finite slope at a corner
 (c = 1 is the line).  Intermediate alpha <=> T(theta) = O((2-theta)^alpha) near theta -> 2 suffices.
 Corner problem = "Kronecker-rank rigidity of SWAP": K of operator-Schmidt rank q = N^eta with
 rank(SWAP - K) <= N^{2 - eta/c}; App B gives only N^{2 - eta^2/4} (gamma = 1/2).
 Theorem B barrier: rank methods cannot exclude slope 1, so the corner is as open as the line.
 MC lower bound p >= c N^2 ln N / q allows slope 1 too. q = o(log N): factor 2 only (constant, irrelevant).
Caveats: Lemma U is [P] design-level (G25-P9, unreviewed), extreme unbalance needs x_b >> log x (fine), pieces must be a
 nested family with poly(lambda) pieces; a fixed finite family has gap g = O(log S / log N0) near the corner (rank
 subadditivity), so the corner must be filled by lambda-dependent families (S -> inf), like App B.

### Exactness (both directions, design level, model M)
Duality: e' >= max_H min(e_a u, e_b v) = min_mu h(mu)/(mu/e_a + (1-mu)/e_b)  (worst iota).
If h(mu) >= ||(mu,1-mu)||_{1/(1-a)} for all mu then Phi=1/e satisfies Phi' <= ||(Phi_a,Phi_b)||_{1/a} (Hoelder),
so Phi^{1/a} is subadditive in x = log K -> power <= a. Conversely if h(mu) < ||.|| for one mu, the self-similar
induction (e(x) <= C x^{-a}, split x_b/x_a = kappa fixed, x_b >> log x) gives power >= a. So
   alpha* = inf{a : h(mu) >= ||(mu,1-mu)||_{1/(1-a)} for all mu}   (generalises G25-P9's Pythagoras a=1/2).
Corollary: alpha* = 1 <=> gamma_corner = 1 (corner exponent of the frontier at (theta,T)=(2,0)).
Pieces need NO nestedness: any chain of realizable pairs (p_s,q_s) gives pieces A_s = R_s - R_{s-1} = C_{s-1} - C_s
with ranks <= (p_s+p_{s-1}, q_s+q_{s-1}); so alpha* is a function of the asymptotic pair frontier T(theta) only.
selfsim.py: corner slope 4: self-similar ratio <1 for alpha=0.6..0.9 (0.9: 1-1e-7 at x_b/x_a=1e-6). DP (dp_long.py)
stalls at 0.666 because np.interp clamps e(x_a) for x_a in the last grid cell (artifact; the induction is exact).

### Monomial-shaped frontier (mono_check.py)
Frontier T = 4 - 2 theta on [4/3,2] (+ mirror) = shape of the monomial-model bound p^2 q ~ N^4:
 balanced rho = 2/3 (alpha_bal = 0.585) but every mu >= 2/3 satisfies h(mu) = mu  => alpha* = 1 (STEP 2).
=> The "spreading conjecture" wall (2theta+T >= 4 => rho >= 2/3, G25-B2) is NOT a wall on alpha in model M.
=> HWY conjecture (balanced point >= N^{1+eps}, eps < 1/3) is compatible with step 2 in model M.

### Q4 (targets.py): what alpha=0.585 needs
 balanced (mu=1/2): whole frontier under theta+T = 8/3: beat circle on theta in [0.86,1.81], max 13% in T.
 mu=0.8: beat circle on theta in [1.82,2], up to 46% in T; mu=0.9: [1.96,2], factor 1.72; corner: gamma >= 0.585.
 finite data (finite_chain.py): a fixed finite chain tensor-powered has h(mu) = log_{N0^2} sum_s r_s^mu c_s^{1-mu}
 (entropy of mixing = log-sum-exp). Exact N=7 frontier / N=8 points: rho > 0.93 even optimistically -> useless.
 Two-piece (Lemma 2 style) families: rho = 1/2 + log_{N0^2}(sqrt p0 + sqrt q0) >= 3/4 for every N0 (Theorem C).

### MC spectral remark (corner)
 tr(Theta) = sum_l tr(G_l G'_l) >= N^2 - (1+Q)p  and |tr(G G')| <= N Q  =>  p >= N^2/(2(1+Q)) or qQ >= N/2.
 So a corner point p = N^{2-delta}, q = N^{c delta} needs gauge-minimal norms Q >= N^{min(delta, 1-c delta)}/2.
 Bell-diagonal (Weyl/abelian-symmetric) model: r = s - k, |r^| = N off supp(k^) => sum r^2 >= N^2 - q; with entries
 bounded by B: p >= (N^2-q)/B^2 (no corner). Large entries are necessary there too.

## Q2 (OpenAI side) - reading g18e, g18s, g19a, ARTICLES Art. 3
- Cost of App. D root identity = 2^{Theta(M/sqrt u)} in model DF (A'' lower 0.564 M/sqrt u for near-uniform product
  loads; App. D upper ~12.8 M/sqrt u). u <= a/2 and (D.6) r^{4a} <= M: u = O(log M) is forced by the UNIFORM PAIR-CLASS
  LOADS of Lemma D.3 (every class of the label cube must be populated), not by precision (precision Theta(M D_t/u)
  even decreases with u).
- Step 2 on this side needs u ~ (log M)^2 at cost 2^{O(M/sqrt u)}: lower bounds do NOT forbid it (Theorem E:
  M/(2u) - u; A'' on circular labellings <= M/u). The obstacle is constructive: no DF identity with 2^{o(M)} terms is
  known on any non-uniform labelling (G19a). Smallest decisive object: staircase packet N_DF(u, m) (u bundles of m
  wires on the 2u-cycle, diameter cuts): need 2^{O(m sqrt u)} for growing u; u = 2 is the HWY/Hurwitz instance.
- The unbalanced-children mechanism of Q1 does NOT transfer: App. D is root-dominated (G18d), the root cost is a
  single term count, there are no children with different exponents to weight.

## Q3 / route list (see report for ranking)
R1 [NEW] Corner of the crossing frontier (Kronecker-rank rigidity of SWAP).
   F(N,q) = min rank(SWAP - K), Schmidt-rank(K) <= q;  delta(eta) = 2 - log_N F(N, N^eta).
   App B: delta >= eta^2/4 (gamma = 1/2). Theorem C / MC: delta <= eta (+o). F(N,1)=F(N,2)=N(N-1)/2 (Lean).
   Model M: alpha* >= gamma where delta(eta) ~ eta^{1/gamma}; step 2 <=> gamma = 1.
   Not excluded by Theorem C, Theorem B (rank methods cannot pass the line), MC bounds (linear), monomial bound
   (slope 2). Excluded in: clique-block MC model (T >= (5-2theta)/3: no corner), commuting / Bell-diagonal / group-
   symmetric families (factor <= 2, or needs entries >= N^{delta/2}), fixed vertex algebras (Cardy: growth sqrt ->
   gamma = 1/2), tensor products (corner slope gets worse), finite chains (entropy).
   Construction heuristic: in App-B-type graded constructions gamma = growth exponent of the side that is CHEAP at the
   corner (log dim(weight w) ~ w^beta); an ASYMMETRIC construction (cheap side with beta > 1/2) would suffice.
R2 [NEW] Monomial-shaped frontier p^2 q = N^{4+o(1)} on theta in [4/3,2] => step 2 (mu = 2/3), although balanced
   rho = 2/3 only. Obstacle: best monomial/MC constructions have exponent 2; G26 spectral bounds need large norms.
R3 Balanced route (Naslund P12 / HWY): line => step 2; rho=2/3 needs whole curve below theta+T=8/3. Now shown NOT
   necessary. HWY conjecture (balanced >= N^{1+eps}, eps < 1/3) is compatible with step 2 in model M.
R4 OpenAI more directions (u ~ log^2 M, non-product labellings / staircase packet). Open; u=2 instance = HWY.
R5 D-direction pieces / non-binary splits: same support-function law in R^D; corner analysis identical. Closed as an
   independent route (G18c, G20a, G6b) - adds nothing beyond R1.
R6 Amortisation / structured children: children see equality-test cliques; a restricted crossing is a full swap
   (G20e Prop R); harmonic law caps at step 2 (G9d). Closed.
R7 Re-embedding between levels (G25-B3 G1-2): dead (G25-P9). Closed.
R8 ARC / omega = 2 / fixed tensor ranks: constants only (G21c Theorem R). Closed.
R9 Classical rigidity (Valiant; Alman-Williams, Dvir-Edelman non-rigidity): sparsity is not Kronecker rank; Fourier-
   sparse (Bell-diagonal) perturbations need large entries (computed above). No implication either way.
R10 Border decompositions: already usable by interpolation (G19b, G21b); border frontier is what counts; small-N
   border values b = 3,5,6 (N=4,5,6). Corner of the border frontier equally open. Not a separate route.
R11 Williams-type consequences: step 2 would only push E^NP lower bounds slightly beyond linear (via local
   reductions); no known-hard statement implied. Not a barrier.
R12 Sat-preserving sparsification / decision: removes loglog in the n-scale only (G25-P1). Not step 2.

### Norm lemma (any decomposition written as in G26 MC: Theta = K Swap, Q = sum ||G_l|| ||G'_l||) [P, sketch]
 eigenvalue 1 has multiplicity >= N^2 - p, other |eigenvalues| <= Q, tr Theta = sum_l tr(G_l G'_l), |.| <= N Q
 => p >= (N^2 - N Q)/(1 + Q). With Q = N^kappa (kappa<1): theta >= 2 - kappa, whatever q is.
 => a corner point (2-delta, eta) needs gauge-minimal norms Q >= N^{delta}. (Same spirit as G26-P3.)

## Verdict on Q1
NOT necessary. The equivalence "step 2 <=> curve near the line" is an artefact of balanced levels (g = h, equal
children). With the unbalanced levels already validated in G25-P9 (Lemma U), the power is a functional of the
frontier's support function h(mu); step 2 <=> corner exponent 1. Every route inside model M still goes through
the crossing frontier, but NOT through Naslund P12's balanced point. Walls stated for balanced rho (spreading
conjecture 2theta+T>=4, "rho >= 2/3 with concavity") are not walls on alpha unless they hold for every weighting.
Processes: none of mine left (the running sscan.py / driver.py belong to other agents). Ended ~08:10 UTC.
