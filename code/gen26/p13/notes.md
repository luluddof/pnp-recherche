# G26-P13 notes (S2 variant: per-edge widths, second-order caps, covariance of disjoint group pairs)
## log
- 04:01 start. Read p10, p5, b7, p2 notes, paper Defs 1-5, Lemma 1, Lemma 2, Prop 3, Lemma 8, p9 Lemma U.
- 04:15 key ideas found:
  (a) caps via mu-weights mu_e = lambda_e + |lambda_e - c| and a RANDOM placement of the child's nodes on its grid
      (concentration) -> P, P/g of the row network bounded by (1+eps)(T + D)/g, (T+D)/g^2: the max width never enters.
      So P10's padding model (q = 1) is fully accounted: c1 = 3.312 survives (ii).
  (b) SPLITTING instead of padding: [x=v][y=u] with |x|=l1>l2=|y|: x=(x',x''), pieces on [x'=v'][y=u] (length l2),
      times [x''=v''] which factors through l1-l2 bits BOTH row- and column-wise. Rectangle cost
      w ea c_row + (1-w) eb c_col <= E_w*4 l2 + S_w*2(l1-l2)  (ordered), vs padding E_w*4 l1.
  (c) (iii) proved: Var(Y_e - Y_e') <= 2(1+eps_K) sigma^2 for disjoint group pairs, eps_K = (4K-2)/((K-2)(K-3)).

## 1. The S2 variant with per-edge widths: statement and proof  [P = proved here, R = inherited]
Model conventions as P2/P5/P10 (asymptotic: lambda -> inf, App. B pieces with (c_row, c_col) inside the quarter circle
of radius 2 lambda (HR leading term, r(lambda) -> 1), kappa/c_ev terms dropped). Everything is per "calibrated total":
for widths lambda = (lambda_e) on the K(K-1)/2 edges of a K-clique,
   T(lambda) = (K/(K-1)) * sum over ORDERED pairs of lambda_e   (= K^2 lambda for uniform width, the paper's unit N).
Deviation measures (all ordered, times K/(K-1)):
   D(lambda; c)  = (K/(K-1)) * 2 sum_e |lambda_e - c|            (any constant c),
   Dbar(lambda)  = (K/(K-1)) * E_pi sum_q 2|lambda_{e1(q)} - lambda_{e2(q)}|, pi = uniform random placement of the K nodes
                   on the solver's g x h grid, e1(q), e2(q) the two diagonals of rectangle q. Dbar <= D(.;c) for every c.

### Def. 3' (edge-weighted clique). As Def. 3, but the string on edge {w,w'} has its own length lambda_{ww'} >= 0;
c_{w,w'} = lambda_{ww'}, Psi = [z_w(w') = z_w'(w)]. (Factorisation through lambda_e bits by obs. (b) after Def. 1.)
### Def. 4' (two-parameter achievability). A solver has parameters (e, q, kind) if for every width vector it outputs the
value within 2^{e (T + q Dev) (1 + o(1)) + lower order}, Dev = Dbar for kind 'rect', Dev = min_c D for kind 'plain'.

### Step A (grouping produces edge-weighted children; replaces Lemma 1 / Lemma 8 (ii))  [P]
Row network on g = tK nodes with c_ab in [0, P/g], c_a <= P, total U, N = gP. Stratified grouping of R2F (groups of
exactly t nodes). The child is the edge-weighted K-clique with lambda_{gg'} := Y_{gg'} (NO appended bits): the
construction of Lemma 1 / (26) goes through verbatim with L replaced by Y_{gg'} on edge {g,g'} (blocks of c_{ab}
consecutive positions fill the string exactly; the factor [all appended bits are 0] disappears). Value unchanged.
Child totals: T_child = (K/(K-1)) sum_{ordered gg'} Y <= (K/(K-1)) U <= U + N/(K-1)  (deterministic).
E D(Y; EY) <= (K/(K-1)) K(K-1) sigma = K^2 sigma, sigma := max_pair sd(Y) <= sigma_S2 (P5's Lemma S2 [R, reviewed P10]),
  K^2 sigma_S2 / N = sdS2(t,K) = (1/sqrt2 + sqrt(V)/t)/t, V = sum_{k=1}^t (s2(sqrt(k-1)/2 + 1/2) + 1/2)^2.
E Dbar(Y) <= (K/(K-1)) 2|R| sqrt(2(1+eps_K)) sigma  (Step E below)  <= f_K K^2 sigma,  f_K = sqrt((1+eps_K)/2)
  (2|R| <= K(K-1)/2 rectangles-diagonals... precisely: the 2|R| diagonals are a subset of the K(K-1)/2 edges).
The child cost bound e_c (T_child + q_c Dev_child) is a computable function of the grouping, LINEAR in Dev_child; the
algorithm (exhaustive search, as Lemma 8) takes the grouping minimising it, which is <= its expectation:
   cost(child) <= e_c ( U + N/(K-1) + q_c f sdS2(t,K) N ),   f = f_K (kind rect) or 1 (kind plain / D-route).
So z(t,K) = q_c f sdS2 + 1/(K-1). Only expectations are used: no union bound, no max over pairs.

### Step B (caps without the max width: point (ii))  [P]
A level (Lemma U) receives an edge-weighted K'-clique, K' = g h, and places its nodes on the g x h grid. For a placement,
rectangle q costs at most 2 max(l1,l2) = l1 + l2 + |l1 - l2| bits per direction in total (padding: c_row_i(l1) <= 2 l1;
splitting: c_row_i(l2) + (l1 - l2) <= 2 l2 + l1 - l2 <= 2 max). So the row network has
   c_{a,a'} <= M(a,a') := sum of lambda_e over the h^2 edges between rows a, a' + sum_{q with rows a,a'} |l1(q) - l2(q)|,
   c_a <= sum_{a'} M(a,a').
Under a uniformly random placement, E M(a,a') is the same for all row pairs, and g^2 E M(a,a') <= (T + Dbar)(1 + O(1/g)).
M(a,a') is a function of a random bijection with bounded differences (a transposition changes it by <= 4h lambda_cap,
lambda_cap := t^2 P_par/g_par = N_par/K'^2 = the deterministic bound Y <= t^2 (P/g) inherited from the parent grouping):
McDiarmid for random permutations + union over g^2 pairs gives max M <= E M + O(h^{3/2} sqrt(log g)) lambda_cap, i.e.
   N_cap := g * (g max M) <= (1 + eps_pl)(T + Dbar) + eps'_pl N_par,  eps_pl, eps'_pl = O(sqrt(log g / h)) -> 0
(h >= K_b = 2^{x_b}). The same placement keeps the straight totals (vertical <= (1+eps)T/h, horizontal <= (1+eps)T/g),
the group setting lengths (entries term) and the realised Delta(pi) = sum_q 2|l1-l2| <= (1+eps) Dbar within (1+eps) of
their means (sums of very many bounded terms, same concentration); a placement meeting all of them exists (union bound)
and is found by search (poly(K') per level, counted like the grouping search). Hence the row and column networks satisfy
the hypotheses of Def. 2 / Lemma 8 (ii) / R2F / S2 / R3'' with P = N_cap/g, and every overhead term of the level
(straight 1/h, 1/g; grouping excess z; entries 2^{-x}; R3'' terms cA, Gam, 2/g of a sweep base) is (its uniform value)
x N_cap. The max width never enters: the caps are sums over h^2 edges / h^2 rectangles, which concentrate.
[eps_pl dropped in the numerics like P2's other lower-order terms; eps'_pl N_par is third order (times z of the level).]
Plain clique-sweep children need no caps (Step D).
=> (ii) is settled: with padding (q = 1) the recursion closes with exponent relative to T + Dbar (<= T + D), i.e. P10's model is exact
   and c1 = 3.312 (R3'') is fully accounted (no hidden max-width factor).

### Step C (splitting instead of padding)  [P]
Rectangle q, diagonals of lengths l1 >= l2 (strings x,v of length l1 on one diagonal, y,u of length l2 on the other).
Write x = (x', x''), v = (v', v''), |x'| = |v'| = l2. Then [x = v][y = u] = [x' = v'][y = u] * [x'' = v''].
Apply the pieces (Def. 7, length l2) to the first factor: sum_i A_i(x', y, u, v') = d [x'=v'][y=u]. The second factor
[x''=v''] factors through l1 - l2 bits row-wise (x'' in row a, v'' in row a') AND column-wise (x'' in column j, v'' in
column j') by obs. (b). By obs. (a): c_row(q) = c_row_i(l2) + (l1 - l2), c_col(q) = c_col_i(l2) + (l1 - l2).
(Also valid: Lemma 2(a)'s identity per rectangle; Prop. 3's expansion (34) and the row/column networks (35) are
rectangle-local; Lemma 2 and the piece bounds are applied per rectangle with its own l2.) With the minimax weight w:
   w ea 2 c_row(q) + (1-w) eb 2 c_col(q) <= E_w A_q + S_w Delta_q,  A_q = 4 l2, Delta_q = 2(l1 - l2),
   E_w = sqrt((w ea)^2 + ((1-w) eb)^2),  S_w = w ea + (1-w) eb  (<= sqrt2 E_w).
Padding would give E_w (A_q + 2 Delta_q)... precisely E_w * 4 l1 = E_w (A_q + 2 Delta_q)/... (ordered 4 l1 = A_q + 2Delta_q? no:
   4 l1 = 4 l2 + 4(l1-l2) = A_q + 2 Delta_q), i.e. padding charges E_w per unit of T_q = A_q + Delta_q plus E_w Delta_q extra;
splitting charges E_w per unit of T_q plus (S_w - E_w) Delta_q extra: extra ratio (S_w - E_w)/E_w <= sqrt2 - 1 = 0.414.

### Step D (level recursion with two parameters)  [P]
Children (e_a, q_a, kind_a), (e_b, q_b, kind_b); z_a = q_a f_a sdS2(t_a, K_a) + 1/(K_a - 1) (Step A), alpha = 1/h + z_a,
beta = 1/g + z_b, ent = 2^{-min(x_a, x_b)}. For each iota choose rows or columns (min <= w-average), so with
delta = Dbar/T, delta' = min(delta, 1) (Delta <= min((1+eps)Dbar, T)), and N_cap <= (1+delta) T (Step B, eps_pl dropped):
   cost/T <= (1 - delta') E_w + delta' S_w + (1 + delta)(O_w + ent),   O_w = w ea alpha + (1-w) eb beta.
At delta = 0 this is p9's Lemma U bound; with the minimax w it equals P10's e_lev = max_th min(ea(cos th + alpha),
eb(sin th + beta)) + ent (minimax over the convex hull of the arc; the max of the increasing min is on the arc).
Every chord slope in delta is <= S_w - E_w + O_w + ent, hence the level has parameters
   e = e_lev,   q = (S_w - E_w + O_w + ent)/e_lev,   kind = rect.
Padding instead: q = 1 exactly (cost/T <= (1+delta) e_lev).
Row-sweep base (R3'' ordering of the row network, P2/P10's sc()): cost/T = G(delta) = (1+delta) max_th min(
sc_g(xr_a), sc_h(xr_b)), xr_a = ((1-delta')cos th + delta')/(1+delta) + 1/h (xr per N_cap); G is nondecreasing in delta
(each term is), so q = sup_d (G(d) - G(0))/(d G(0)) <= max_i (G(d_{i+1}) - G(0))/(d_i G(0)) on a geometric grid
(1e-5..1e3, 160 points) [N; the interval (0, 1e-5] is covered by the chord at 1e-5 assuming G is right-differentiable
at 0 with a derivative not exceeding nearby chords - piecewise smooth max-min; not proved].
Plain clique sweep (s(K) cut): cut_j(lambda) <= s(K) c' + beta sum|lambda - c| with c' the empirical mean,
beta = 2 s(K)/(K(K-1)), giving q = 1 exactly with the D-route (kind plain, f = 1).
Top level: as P5/P10 (top grouping excess replaced by z = q_top f sdS2 + 1/(K-1) of the chosen top solver; the top
cluster network itself is P2's Lemma 8' [R]; P5's caveat that the top clique stays uniform except for the grouping
term applies).

### Step E (covariance of two disjoint group pairs; task (iii))  [P]
Claim: for distinct groups g1, g2, g3, g4, e = {g1,g2}, e' = {g3,g4}:
   Var(Y_e - Y_e') <= (1 + eps_K)(Var Y_e + Var Y_e') <= 2(1 + eps_K) sigma^2,   eps_K = (4K - 2)/((K - 2)(K - 3)).
Proof. F_k = sigma(pi_1..pi_k). Y_e depends only on (a_i, b_i) = (pi_i^{-1}g1, pi_i^{-1}g2), and strata are independent,
so E[Y_e | F_k] = E[Y_e | G_k] (P5's filtration): the Doob increments D_k of Y_e are those of Lemma S2/R2F, and D'_k
(of Y_e') likewise; both are F_k-martingale differences, so Var(Y_e - Y_e') = sum_k E(D_k - D'_k)^2.
Given F_{k-1}: D_k = F(a_k, b_k), D'_k = H(c_k, d_k) with F, H F_{k-1}-measurable, each of mean 0 over a uniform ordered
pair of distinct nodes of S_k, and (a_k, b_k, c_k, d_k) = pi_k^{-1}(g1..g4) is a uniform injection into S_k (|S_k| = K).
Then E[F H | F_{k-1}] = -(1/((K-2)(K-3))) E_{(a,b)}[F(a,b) m(a,b)], with (inclusion-exclusion, sum_{c != d} H = 0)
   m(a,b) = sum_{d != a} H(a,d) + sum_{d != b} H(b,d) + sum_{c != a} H(c,a) + sum_{c != b} H(c,b) - H(a,b) - H(b,a).
Each row/column sum is (K-1) E[H | one coordinate], whose rms is <= (K-1) rms(H) (Jensen); so rms(m) <= (4K - 2) rms(H)
and |E[F H | F_{k-1}]| <= eps_K rms(F) rms(H). Hence E[(D_k - D'_k)^2 | F_{k-1}] <= (1 + eps_K)(rms F^2 + rms H^2),
sum over k and take expectations. QED.
Consequence: E|Y_e - Y_e'| <= sqrt(Var(Y_e - Y_e')) <= sqrt(2(1+eps_K)) sigma (equal means by symmetry), against
E|Y_e - c| + E|Y_e' - c| <= 2 sigma: factor f_K = sqrt((1+eps_K)/2) -> 1/sqrt2. Under a uniformly random placement the
two diagonals of every rectangle join four distinct child nodes = four distinct groups, so E_grouping Dbar(Y)
<= (K/(K-1)) 2|R| sqrt(2(1+eps_K)) sigma <= f_K K^2 sigma (2|R| <= K(K-1)/2... number of diagonals <= number of edges).
(Remark: the linear part alone has E(Lambda - Lambda')^2 = 4 s^2 (1 + 1/(K-1))/K^2 vs E Lambda^2 = 2 s^2(1 - 1/(K-1))/K^2,
consistent with the lemma.)
- 04:13 first split run: row-sweep base q_split ranged 0.735..1.705 (> 1 for skewed d = 4 bases: splitting adds
  (l1-l2) to BOTH directions, bad for the low-share side of a non-linear sweep cost). Fix: per solver q = min(q_split, 1)
  (pad when better; padding gives q = 1 exactly). For levels S_w - E_w <= (sqrt2-1)E_w so splitting always wins there.
  Relaunched split runs at 04:14.
- 04:16 pad (q=1, f=1) = 3.312 (reproduces P10 exactly; thresholds 1514/2911/1.002e4 reproduced with mythr13 top S2q).
  padiii (q=1, f_K for rect children) = 3.284 [thr 1452 / 2849 / 9916].
- 04:17 split run #1: base q was inflated by a too-large grid margin (1e-6/(1e-5 G0) ~ 0.4) -> bases had q = 1 in the
  optimal trees; c1 3.274 (split) / 3.212 (splitiii) -- valid but loose. Fixed: base max-min now by bisection on the
  crossing (sc monotone), margin 1e-12. Base q now 0.47-0.60 (g = 2^8..2^20). Relaunched 04:18.
- 04:20 split (fixed base q) = 3.220; splitiii = 3.166 (a=2, [8000,16000]; ranges 3.158/3.164/3.174/3.166, so
  limsup over [1000,16000] is 3.174). a=4: 4.554 / 4.477. Level q: 0.37-0.43 (median 0.41); base q 0.47-0.79.
- 04:21 test_cov13.py/.out: (1) exact injection sums: |E F H| <= eps_K rms F rms H holds with ratio <= 0.15 (K = 5..12,
  random/correlated F, H); (2) Monte Carlo stratified grouping (K = 8, 16; t = 8; unif/bern/rank2 networks):
  Var(Ye - Ye')/(Var Ye + Var Ye') = 0.95-0.98 (bound 1 + eps_K), E|Ye - Ye'|/(E|Ye-EY| + E|Ye'-EY|) = 0.69-0.73
  (the lemma's factor f_K -> 0.707 is essentially tight on these examples).

## 2. Results (DP = P10's dp10.py model + (e, q, kind) per solver; R3'' base; a = 2; bins 1.02; lt step 1; c1 = limsup over
log m/a in [8000,16000]; thresholds = P5/P10 mythr5 with top grouping z = q_top f sdS2 + 1/(K-1) ("S2q", mythr13.py))
Formulas: z(t,K) = q_child * f * sdS2(t,K) + 1/(K-1); sdS2 = (1/sqrt2 + sqrt(V)/t)/t, V = sum_{k=1}^t (s2(sqrt(k-1)/2+1/2)+1/2)^2,
s2 = sqrt(2(1+1/(K-1))); f = sqrt((1+eps_K)/2), eps_K = (4K-2)/((K-2)(K-3)) (rect children, (iii)) or 1;
level q = min(1, (S_w - E_w + O_w + 2^{-min x})/e_lev); row-sweep base q = min(1, sup chord of G) (bisection); plain sweep q = 1.
| variant | status | c1 (a=2) | a=4 | r=4.26: 2^n / 1.6423^n / 1.30704^n | r=1 |
|---|---|---|---|---|---|
| pad, f=1 (= P10's corrected S2) | [P] (Steps A,B,D; caps settled) | 3.312 | 4.684 | 1514 / 2911 / 1.002e4 | <=50 / 82.8 / 279.6 |
| pad + (iii) | [P] | 3.284 | 4.645 | 1452 / 2849 / 9916 | <=50 / 81.8 / 271.5 |
| split, f=1 | [P] (+ Step C) | 3.220 | 4.554 | 1401 / 2753 / 9457 | <=50 / 80.7 / 263.5 |
| **split + (iii)** | [P] (+ Step E) | **3.166** ([4000,8000]: 3.174) | 4.477 | **1352 / 2683 / 9269** | <=50 / 79.7 / 256.4 |
("<=50" = bisection floor of the threshold code.) The pad row reproduces P10 exactly (DP and thresholds), which validates
the modified code path. Split+(iii) is now 0.19 below floor B (3.36), i.e. clearly beyond the uniform-width design.

## 3. Status / what is NOT proved
* Proved here (on top of [R]: P5's Lemma S2 sd bound (P10-reviewed), R3'' (P10-reviewed), Prop. 3 / Lemma 8 constructions,
  App. B pieces with the HR circle limit, P2's top Lemma 8'): Steps A-E above; in particular (ii) is settled by the
  random-placement cap argument (caps are sums over h^2 edges, never the max width) and (iii) by the covariance lemma.
* Lower-order terms dropped exactly as in P2/P5/P10's model: eps_pl (placement concentration, O(sqrt(log g/h))),
  eps'_pl N_par (third order), piece log terms c log(2 lambda) per rectangle (need all widths -> inf; tiny widths can be
  padded to lambda_0 = o(mean) at o(N) cost), rho_bar finite-lambda corrections (thresholds use P2's conservative rb^depth).
* Numerical, not closed form: base q via a 160-point geometric delta-grid (uses monotonicity of G; the interval (0,1e-5]
  assumes right-differentiability with derivative <= nearby chords); level q is closed form (w-trick).
* DP bins keep the max-mu solver per bin and ignore q (valid, maybe slightly suboptimal); unchanged P10 conventions otherwise.
* Top level: only the top grouping term uses S2 (P5's caveat; the top cluster network/Lemma 8' is P2's).
* Not attempted: q-aware pruning, per-rectangle mixed pad/split for bases, Dbar-route for plain-sweep children
  (they keep f = 1, q = 1), splitting at the top clique. All could only lower c1 further.
## Files
notes.md; dp13.py (DP, modes pad/padiii/split/splitiii); mythr13.py (thresholds, top S2q, env III=1 for f);
zlib10.py, mylib.py, mythr5.py, dp10.py (copied from P10, unchanged); sol13_R3pp_*_1.02.json; dp13_R3pp_*.out;
thr13_*.out; test_cov13.py/.out. No process left running; nothing written to /home/user/pnp-recherche.
