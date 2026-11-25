# G28-O1 notes (start 09:44 UTC)
## Reading done: BRIEF, B2 notes, R1 notes, R1 chain/int/dp scripts, P13 dp13/mythr13, P14 review.
## Task 2 plan
Hypothetical frontier F_c = min(circle, c(1-u)) in normalized (u,v)=(theta/2,T/2), one corner at (1,0).
 (A short corner only on [1-s0,1] has effective concave-majorant slope max(c, sqrt(2/s0)): the corner line
  must reach the circle, so F_c is the right model.) c=1 => F_1 = full line.
Implement dp_corner.py = dp13 (splitiii) with circle replaced by hull H_c in levels AND row-sweep bases,
 both grid orientations (one-sided corner breaks symmetry). Validate with c=inf (must give 3.166).
## 09:53 dp_corner.py validated: c=inf reproduces P13 splitiii exactly (3.158/3.164/3.174/3.166).
Model remark: effective corner slope = slope of the hull. A family with true slope c but extent s0 only gives
 c_eff = max(c, v_A/s0) ~ sqrt(2/s0) (hull contains chord arc-end -> (1,0)); so "slope c" = line reaching the circle.
c=1 => full line frontier (rho = 1/2): e*x/log2x ~ 0.36-0.67 decreasing. c=2: e*x/log2x ~ 1.27-1.40 (R1 chain bound 4c log x/x = 8).
## 09:58 compare.py (ratio e_bal/e_c of the DP tables, x = log2 K_top ~ log2 m / 2), ONE-sided corner:
 c1 (= full line): gain >=1.1 from x=14, x2 from x=67, x4 from x=294.
 c2: >=1.1 from x=33, x1.5 at 183, x2 at 420, x4 at 2592.   c4: >=1.01 at 133, 1.1 at 299, 1.5 at 963, 2 at 2353.
 c8: 1.01 at 766, 1.1 at 1448, 1.5 at 11570.
 NOTE: frontier is symmetric under (p,q)<->(q,p) (relabel y<->u maps row flattening to column flattening, SWAP fixed),
 so a corner at theta=2 implies the mirror corner: TWO=1 runs are the natural model (c<=1+sqrt2: polygon (1,0),(k,k),(0,1)).
## Task 3 (finite_h.py): every chain of known points (N<=12, F(N,1) up to N=64), tensor-powered with the entropy
 formula h = log_{N^2} sum_s r_s^mu c_s^(1-mu): h_fin(mu) > h_B(mu) for ALL mu in [0,1] (min gap at the trivial endpoint
 mu=1, gap 0.16 (N=2) ... 1e-4 (N=32)), both optimistic nested pieces (p_s,q_{s-1}) and safe (p_s+p_{s-1},q_s+q_{s-1}).
 => no level of any shape can gain: e' = min_mu h/(mu/e_a+(1-mu)/e_b) is monotone in h. c1 unchanged (3.17).
 finite_hmax.py (entropy ignored, max instead of sum): chains would beat the circle (h(1/2)=0.50..0.60 for N=2..8)
 -> the multiplicity (entropy) of tensor-powered pieces is the whole obstruction, not the points.
## Task 1 (OpenAI side) -- reading of App. D (paper.txt 3007-3660) + g19a notes + ARTICLES Art. 3
Recursion structure. A call with M wires: for each node t in S (2L nodes, direction i = v2(t) < u) and each mode
tuple beta (count exp(O(sqrt(M p_t/D_t))) = exp(O(M/sqrt u)), (D.15)), F(t+alpha zeta) = prod_{g=1}^r Z_{g,beta}
(D.16): r children, each <= 3M/32 wires, contracted INDEPENDENTLY. Cost T(M) = sum_t N_t * sum_g T(M_g):
   log T(M) = log max_t N_t + max_g log T(M_g) + O(log)       ("divide" recursion: children are disjoint parts).
Claude's level: child = the WHOLE row (column) network, exponent e_a multiplies its full size ("compress" recursion),
   e' = max_H min(e_a u, e_b v): child efficiencies are WEIGHTS on the two directions -> support function, corner.
(a) Unequal direction groups / per-direction term counts N_i (i = valuation level/direction): App-D pays sum_i N_i
   ~ max_i N_i with EQUAL weights (children of every direction are the same kind of object and are negligible:
   root-dominated, 3M/32 shrink vs u halving). So only the equal-weight value (the analogue of mu = 1/2, i.e. the
   balanced rho) matters; the corner of a per-direction trade-off is irrelevant inside App D. For u = 2 the staircase
   packet with two bundles IS the swap crossing: direction 0 = row split (x,y)|(u,v), direction 1 = column split
   (x,u)|(y,v); a DF identity = R + C with N = p + q. App D's cost functional is p + q (balanced/HWY point), Claude's
   is the support function h(mu) of the (p,q) frontier. The corner question lives only in the weighted functional.
(b) Unbalanced groups (one big child (1-delta)M): log T = f(delta) M + log T((1-delta)M) -> exponent = f(delta)/delta,
   a RATIO, no compounding. Deforming only the delta M crossing wires lowers the degree of F to delta M (2L-1), so the
   precision scales too: f(delta) ~ delta/sqrt u: neutral. On expanders cut(delta-set) >= phi delta M and Theorem E
   (log N >= M_i/(2u) - u) gives ratio >= phi/(2u): never below the 1/u law.
(c) Precision/modes trade-off: tuple count ~ exp(2 sqrt(zeta(2) K p)), K = 6M/D mode variables per weight, p ~ 32 M D/u;
   K p = Theta(M^2/u) independent of the level i (taper equalises). Product K*p is the only invariant; black-box
   precision Theta(M D/u) optimal (G17a, const ~12) -> no corner (it is a product law, not a two-sided frontier).
(d) Hamming-distance profile (g19a sec. 4): per-wire A'' cost E|S_d|/u is CONCAVE in d; expanders force mean d >= phi' u;
   the LP optimum mixes the extremes d=0 and d=u (lower convex hull = chord) -> 0.8 phi'/sqrt(u). This IS a
   corner-type (mixture) computation and it gives nothing for exchangeable loads. Non-exchangeable (circular) labellings:
   A'' <= M/u, Theorem E >= M/(2u) - u: open (staircase packet).
(e) Why u = O(log M): each packet = one pair class, L = 2^a >= 4^u slots, pair classes ~ (r^2/2)^a, (D.6) M >= L D_arr
   => u <= a/2 <= (1/2) log2 M even for r = 2; also L <= s^{1/1000} for packet-tensor evaluation exp(O(L log s)).
   B2's statement confirmed (forced by the packet/pair-class design, not by the instance or by precision).
Analogue of alpha*: OpenAI exponent = M * f(u) at the root with u <= u_max(M); power of log = beta*kappa if
 f = u^-beta, u_max = (log M)^kappa. App D: 1/2 * 1. Step 2 needs beta*kappa = 1:
 (i) beta = 1, kappa = 1: per-wire 1/u on packets: FORBIDDEN for near-uniform product loads by A'' (0.564/sqrt u, Lean);
 (ii) beta = 1/2, kappa = 2: needs labellings outside Lemma D.3 (no packets); Theorem E allows (M/(2u)), flattening
      cannot decide (circular: A'' = M/u); constructive object = staircase packet N_DF(u,m) (R4 of B2).
 (iii) a compress-type recursion using App-D multi-direction pieces = Claude's architecture with D-direction pieces
      (B2 R5: same support-function law in R^D) -> the OpenAI "corner question" is Claude's corner question
      (u = 2: literally the swap crossing frontier) or its R^u version; nothing new.
Theorem E is the OpenAI analogue of Theorem C (alpha* <= 1): one-level identities on expanders cost >= M/(2u) - u, so
 with u = Theta(log M) step 2 is the ceiling; step 2 would need Hamming-1-like efficiency with every direction cutting
 >= phi' M wires (cut-conditioning can't: it costs min M_i >= phi' M).
## 10:10 thresholds (mythr13, III=1, S2q, r=4.26; log2 m where exponent/n <= 1 (2^n), log2 1.6423, log2 1.30704):
 balanced (c=inf): 1352 / 2683 / 9269 (reproduces P13).
 ONE-sided corner: c=1 (line): 229 / 315 / 592 ; c=2: 508 / 746 / 1518 ; c=4: 948 / 1462 / 3072 ; c=8: 1352 / 2362 / 5854.
 r=1, c=1: <=50 / 52 / 103 (bal: <=50 / 80 / 256).
R1-chain formula (e <= ~4c log2 x/x, brute-force small child, no App-B mixing) vs balanced 3.17/(2 sqrt(2x)):
 crossover x = 1.3e3 (c=1), 8.4e3 (c=2), 4.8e4 (c=4), 2.5e5 (c=8)  [log2 m ~ 2x]: pessimistic, the DP mixes
 circle+corner levels and gains much earlier.
Piece widths at the top: log2 lambda ~ 27-34 (N = 2^lambda): far beyond App B's finite-size regime (N <= 10^12 ~ lambda<=40
 is about the ALPHABET size, here lambda itself is ~1e9), so "App-B-like finite size" of the pieces is negligible here
 (rb^depth factor kept as in P13); all finite-size effects are recursion overheads (dp13 model).
## 10:13 test_lev.py: levH matches brute force over 2e5 boundary samples (rel 3e-6), saddle Ew+Ow = v to 1e-14 (TWO=0/1, c=2,4,8).
## 10:15 TWO=1 (mirror corner) tables are IDENTICAL to one-sided (to 6e-17): the one-sided runs already try both grid
 orientations, and the minimising orientation always puts the crossing on the corner segment. So "one corner" = "two
 corners" in this DP. c=2 two-sided = polygon (1,0),(2/3,2/3),(0,1) = B2's p2q / monomial-shaped frontier (route R2):
 its numbers are those of c=2.
## 10:20 topratio.out: top exponent/n (r=4.26) balanced vs corner, ratio bal/corner at log2 m = 100..25600:
 c1: 1.60 2.30 3.23 4.51 6.35 9.06 12.6 17.9 25.2
 c2: 1.14 1.22 1.51 1.91 2.50 3.28 4.26 5.35 6.70      (= p2q/monomial-shaped frontier, B2 route R2)
 c4: 1.00 1.00 1.02 1.15 1.38 1.76 2.19 2.54 2.86
 c8: 1.00 1.00 1.00 1.00 1.00 1.12 1.27 1.40 1.51
 Gain >= 10%: c1/c2 from log2 m ~ 100, c4 from ~ 700, c8 from ~ 3000. Factor 2: c1 ~150, c2 ~900, c4 ~5000, c8 >25600.
 Thresholds vs 2^n (r=4.26): 1352 (bal) -> 229 / 508 / 948 / 1352 (c = 1/2/4/8); vs 1.30704^n: 9269 -> 592/1518/3072/5854.
## Verdicts
 T1: B2 confirmed and sharpened: App D is a DIVIDE recursion (children disjoint, costs add in log), Claude's a COMPRESS
   recursion (child = whole network, exponents weight the two directions). Corner needs weights; App D's functional
   is the equal-weight sum over directions (= balanced/HWY point for u=2). Analogue of the corner question exists only
   in a compress-type use of D-direction pieces = Claude's corner question (R^u version, B2 R5). Inside App D the
   power of log is beta*kappa (f = u^-beta per wire, u <= (log M)^kappa); step 2 needs (beta,kappa) = (1,1) [blocked by
   A'' on product-uniform packets] or (1/2,2) [needs non-packet labellings; Theorem E allows; staircase packet open].
 T2: see above. T3: no finite piece helps (entropy); c1 stays 3.17.
Files: dp_corner.py, test_lev.py, compare.py, topratio.py, finite_h.py, finite_hmax.py, chain2.sh, *.out, sol_*.json.
No process left running.
