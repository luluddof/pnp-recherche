# G29-X2 notes: P_3 (odd-cycle gadget) vs matrix-multiplication constructions (start 13:35 UTC)
Tags: [P] proof here, [E] exact computation, [N] numeric, [H] heuristic, [R] read / literature.

## 0. Reading done (13:35-13:50): BRIEF, x1/notes + scripts, gen27/b2, gen29/u1, paper Prop 3 / Def 2,7 / Lemma 8.

## 1. Gadget (own re-derivation) -- see below
### 1.1 Gadget re-derivation [R + P design level] (~13:38-13:42)
Paper: level = grid; per rectangle the diagonal factors [x=v][y=u] are expanded into pieces A_i (Def. 7: sum = d*X,
each factors in rows through c_row bits and in columns through c_col bits); per iota a GLOBAL choice rows/cols, cost
min(U_row,U_col) <= average; Lemma 8 balances. P_3 layout (X1): between rows a,a' the RC edges = complete digraph on
the h columns (top col j -> bottom col j'); decompose into directed 3-cycles (Mendelsohn triple system, h = 0,1 mod 3,
h != 6). Gadget = 3 strands a1->b2, a2->b3, a3->b1 (a_i = string held at (a, j_i), b_i at (a', j_i)).
 Operator: P_3 = [b2=a1][b3=a2][b1=a3] = cyclic shift sigma on (C^N)^{(x)3} (N = 2^lam), sigma = S12*S23.
 ROW cost: 2-party (top row a holds a1,a2,a3; bottom row a' holds b1,b2,b3): c_row = log2 rank_{top|bottom}(piece).
 COLUMN cost: 3 columns, column j_i holds (a_i,b_i); a piece must be a pairwise network = periodic MPO
   piece(a,b) = Tr(A1[a1,b1] A2[a2,b2] A3[a3,b3]) with bond dims b12,b23,b31; c_col = log2(b12 b23 b31)
   (bond index stored in node settings at both ends + equality test, which factors through log b bits; local
   tensors multiplied into W_col). So the column side is the triangle tensor network: P_3 = <N,N,N> with bonds N.
 Baselines: 3 lam bits each side per gadget (rank N^3; bonds N,N,N). Normalised piece cost
   (u,v) = (log rank, log prod bonds)/(3 lam); trivial pieces (1,1). Same level map as SWAP (Lemma U, G25-P9,
   confirmed G29-U1) with H = hull of normalised P_3 piece costs (rows 2-party, cols 3-party, nothing else changes).
 => alpha*(P_3) = inf{a : h_P3(mu) >= ||(mu,1-mu)||_{1/(1-a)} for all mu}, h_P3 = support function of the P_3
   frontier in (u,v); balanced criterion rho_P3 = h_P3(1/2), alpha_bal = log2(1/rho_P3); corner exponents at
   A = (1,0) ("cheap ring") and B = (0,1) ("cheap rows"). Chain trick (B2/U1) gives pieces from any realisable chain.
 Validity caveats checked: bond indices enlarge column settings P by a constant factor (entry cost 2^{2gP/K}: lower
   order); Lemma 8(ii) c_{j,j'} <= P/h holds up to factor 2 when bonds exceed N. Agree with X1: VALID at design level.
 X1 bound p >= N^3/(4 b(1) b(2)): proof re-checked (compress strand on smallest bond to t values; cut putting two
   nodes on the top side crosses N t^2 b_e min(other)); consistent with all exact points (N=2: (4;1,1,1) bound 2;
   (2;1,1,4) bound 2 TIGHT; (3;1,2,3)[N] bound 1). Normalised u >= 1 - (2/3) v.

### 1.2 NEW: composite construction  SWAP => P_3  [P + E] (~13:42)
 sigma = S12 S23. If SWAP = R + K (rank R = p, osr K = q) put K12 = K(x)I, K23 = I(x)K. Then
   sigma - K12 K23 = R12 S23 + K12 R23   (rank <= 2pN),  and K12 K23 = sum_{l,m} G_l (x) G'_l H_m (x) H'_m is a
   periodic MPO with bonds (q, q, 1).
 => every SWAP point (p,q)@N gives the P_3 point (2pN ; q,q,1)@N. Normalised: u = (1+theta)/3, v = 2T/3.
 Exact check (composite.py, composite.out; flint over Z): reproduces C=I for q=1 (N=2..5: 4,16,40,80); with gen26
   exact SWAP points: N=7 (11,7)->148 [ring 49 vs 343], (20,3)->238; N=8 (20,6)->294, (22,5)->311; always <= 2pN.
 Asymptotic consequence [P, given App B]: App B circle theta^2+T^2=4 gives P_3 points
   u = (1+theta)/3, v = 2T/3: NONTRIVIAL (both <1) for theta in (1.32, 2), e.g. balanced App B point -> (0.805, 0.943).
   Corner A: 1-u = delta/3, v = 2 eta/3, delta = eta^2/4  ->  1-u = (3/16) v^2: QUADRATIC corner, gamma_A = 1/2.
 => X1's "known asymptotic P_3 frontier is trivial (alpha* = 0)" is WRONG: it is at least the image of SWAP's.
 alpha (alpha_p3.py/.out): composite frontier h(mu) = mu/3 + ||(2mu/3, 4(1-mu)/3)||_2 near mu->1 equals
   mu + (4/3)(1-mu)^2/mu > ||(mu,1-mu)||_2 -> alpha(mu) < 1/2, alpha(mu) -> 1/2 as mu -> 1 (0.42 @0.95, 0.465 @0.999,
   0.479 @1-1e-5; values at 1-1e-6 are float artefacts, the circle itself shows them). alpha*(P_3 composite) = 1/2
   as a supremum (not attained), i.e. NOT better than SWAP's App B (1/2 at every mu). [N + P]
 If SWAP reached its line, composite P_3 would have alpha* = 1 (finite corner slope) -- p3_comp_line. X1's bound
   u >= 1-(2v/3) as a frontier would give balanced 0.585 and alpha* = 1.
 STRUCTURAL: SWAP corner (linear) => P_3 corner A linear (composite). X1: P_3 corner B => SWAP corner.
   So P_3 corner A is a RELAXATION of the SWAP corner problem (weakly easier target), corner B a strengthening.

## 2. Matrix-multiplication decompositions in P_3 language (~13:44)
### 2.1 Dictionary [P]
 <N,N,N> = sum x_ij y_jk z_ki; node L=(a1,b1) <-> A_ij (i=b1, j=a1); M=(a2,b2) <-> B_jk (j=b2,k=a2); R=(a3,b3) <-> C (k=b3,i=a3).
 A CP (rank-1) term u(x)v(x)w = LOCAL PRODUCT OPERATOR U1(x)U2(x)U3 (U_i : a_i -> b_i). Costs: ring bonds (1,1,1);
   row rank = rk U1 * rk U2 * rk U3.
 A sum of s CP terms as a ring: bonds (s,s,1) (hub index routed through two edges) -> ring cost <= s^2.
 Conversely ring bonds (b12,b23,b31) => CP rank <= b12 b23 b31. So ring cost in [CP rank, CP rank^2].
 => corner A of P_3  <=>(slope up to factor 2)  "3-factor Kronecker-rank rigidity of the cyclic shift":
    C = sum_{r <= N^eta} U_r (x) V_r (x) W_r with rank(sigma - C) <= N^{3-delta}, delta >= c eta.
 Full MM decompositions: C = sigma, rank r = R(<N,N,N>) >= N^2 (flattening), = N^{omega+o(1)}:
   star realisation ring = r^2 = N^{2 omega} >= N^4 > N^3 = trivial ring. Useful only if omega < 3/2 (impossible).
   Border rank: same (border ring rank of <N,N,N> on the triangle = (N,N,N): flattening ranks N^2 bound border too).
   omega = 2 would give ring N^{4+o(1)}: STILL worse than trivial. The trivial ring N^3 is already optimal for sigma.
### 2.2 Strassen exact (strassen.py, strassen.out) [E]
 In P_3 language Strassen = sigma_2 = I + T6: term M1 = (A11+A22)(B11+B22)(C11+C22) IS the identity operator
   (row rank 8, ring 1); the other 6 terms are local products of rank-1 operators (row rank 1, ring 1).
 Natural split C = I: (p = 4; 1,1,1) = X1's C=I point. Kronecker powers N = 2^k: sigma = (I + T6)^{(x)k} = sum_j G_j
   (G_j = j non-identity factors). Truncations C_m = sum_{j<=m} G_j (CP rank sum_{j<=m} C(k,j) 6^j):
   k=2 (N=4): rank(sigma - C_m) = 40, 16, 0 for m = 0,1,2 (CP terms 1, 13, 49)
   k=3 (N=8): 336, 256, 64, 0 for m = 0..3 (CP terms 1, 19, 127, 343)
   Exact law (m >= 1): rank(sigma - C_m) = 8^k Pr[Bin(k,1/2) > m] (joint eigenbasis: T6 = sigma_2 - I has 4 zero
   eigenvalues out of 8; Sum_{j>m} G_j vanishes on eigenvectors with <= m nonzero factors).
   Flattening ranks (node|rest) of C_m: sum_{j<=m} C(k,j) 3^j (k=3: 1, 10, 37, 64) [E].
 Asymptotics m = beta k, beta in (1/2,1):  1-u = (1-H(beta))/3;
   star ring: v = 2(H(beta) + beta log2 6)/3 >= 1.53  (OUTSIDE the trivial square);
   any ring (flattening bound, optimistic): v >= (H(beta) + beta log2 3)/2 >= 0.896 (beta=1/2, where u=1);
   e.g. beta = 0.6 -> (u,v) >= (0.990, 0.96). No approach to (1,0): no corner even optimistically. Strassen gives NOTHING.
 m < k/2: rank(sigma - C_m) = N^3 (1 - o(1)).
### 2.3 Relation to known quantities [R + P]
 - Column cost = tensor-network (triangle / periodic-MPO) bond dimension; <N,N,N> = triangle graph tensor (EPR pairs on
   the 3 edges). Its triangle TN rank AND border TN rank is (N,N,N) (flattening ranks N^2). Border bond dimension
   concept: Christandl-Lucia-Vrana-Werner, SciPost Phys. 9 (2020), arXiv:1809.08185.
 - Row cost = rank in the flattening that SPLITS every party (in-legs | out-legs): the cyclic-shift permutation matrix.
   So the P_3 frontier = rigidity of the cyclic shift against periodic MPOs (equivalently, up to factor 2 in exponent,
   against sums of N^eta local product operators U(x)V(x)W). Not a quantity studied in MM theory.
 - omega enters only via ring <= CP^2 for C = sigma: 2 omega > 3 always. omega < 2.371339 (arXiv:2404.16349, 2024);
   an Aug-2026 preprint arXiv:2608.16884 reports 2.371177 (search snippet only, unverified). Irrelevant either way.
 - Direction mismatch: MM upper bounds degenerate the triangle INTO a GHZ/unit tensor (<N,N,N> <= <r>), laser method /
   subrank results show what big tensors CONTAIN; P_3 needs triangle(N) ~ triangle(b<<N) modulo an additive error of low
   rank in a party-splitting flattening. Additive low-rank errors are not monotone under restriction/degeneration, so
   asymptotic-spectrum / laser machinery does not transfer. (Same as G21c "ARC / omega: constants only".)

## 3. VERDICT (13:48 UTC; timestamps above approximate)
 - No MM construction (Strassen and powers, any rank / border-rank decomposition, omega=2, CW/laser) gives a
   nontrivial P_3 point: star realisation ring >= N^{2 omega} > N^3; truncated Strassen powers stay near u = 1 with
   v >= 0.9 even with an (unrealised) optimal ring, and >= 1.53 with the star ring. Obstruction: <N,N,N> is already the
   CHEAPEST object for the ring side (it is the network's own structure tensor); MM speed-ups are about CP/GHZ cost,
   which the pairwise column network never charges.
 - Nontrivial asymptotic P_3 points DO exist, from SWAP (composite, Sec 1.2): App B -> quadratic corner A (gamma = 1/2),
   alpha*(P_3) = 1/2 as a supremum, NOT above SWAP's 1/2. X1's "known frontier trivial, alpha* = 0" corrected.
 - Structural gain: SWAP linear corner => P_3 linear corner A, so P_3 corner A is a relaxation of the SWAP corner problem
   (weakly easier target for step 2; any native P_3 corner-A construction with delta >= c eta would suffice).
 No processes of mine left running (milp_reps.py seen in ps belongs to another agent). Nothing written to /home/user/pnp-recherche.
