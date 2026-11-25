# G28-C1 notes (start 09:44 UTC). Corner of App-B-type constructions.
Tags: [P] proof sketch (unreviewed), [E] exact computation, [N] numeric, [H] heuristic, [C] counted (dimension count only).
## 0. App B as a parametrised family (paper (53)-(56))
N = p(h); pieces A_s, s in [-2h,2h]; row (R) cost R_s^row = sum_i p(s+i) ~ p(2h+s) (intermediate weights <= 2h+s),
col (C) cost ~ p(2h-s). Threshold s0 = sigma*h: R = sum_{s<s0} A_s, C = rest.  Pole-order constraint W + W' = 4h - 1.
## 1. App B corner: exact mechanism [P + N]
log p(w) = c sqrt(w) - log(4 sqrt3 w) + o(1), c = pi sqrt(2/3).  Threshold s0 = sigma h, sigma in [-2,2]:
  theta(sigma) = sqrt(2+sigma) + O(log h/sqrt h),  T(sigma) = sqrt(2-sigma) + O(log h/sqrt h)   -> theta^2+T^2 = 4.
Corner (C cheap): keep only intermediate C-weights <= w = eps*h  (sigma = 2-eps):
  eta = T = sqrt(w/h) = sqrt(eps);   theta = sqrt(4h-w)/sqrt h = sqrt(4-eps)
  => delta = 2 - sqrt(4-eta^2) = eta^2/4 + eta^4/64 + ...   EXACTLY the circle, gamma = 1/2.
WHY quadratic: it is NOT a Gaussian tail. It is the concavity of the Hardy-Ramanujan entropy log p(w) ~ c sqrt(w)
in the variable in which the pole-order constraint is ADDITIVE (W + W' = 4h - 1):
  * R's saving (expensive side, near its top W = 4h): linear in w, slope d/dW sqrt(W) at W=4h = 1/(4 sqrt h) -> delta = w/(4h);
  * C's cost (cheap side, near weight 0): c sqrt(w) -> eta = sqrt(w/h): infinite slope at 0.
  The parameter controlling the corner is the C-side weight cutoff w (equivalently s0 = 2h - w); delta ~ w, eta ~ sqrt(w).
Reparametrisation check: in "mass" units m = sqrt(level) dims are exponential (Hagedorn-like) and the constraint becomes
  m_R^2 + m_C^2 >= 4h: the circle IS the constraint.  The shape of the frontier = (growth law) o (constraint law).
Finite size (appb_curve.py, appb_curve.out; exact partition numbers, multiplicities (2h+1) included): the counted
  theta stays capped at 2 (p = N^2) for eta < ~1.9 (h=50), ~1.3 (h=500), ~0.9 (h=2000): offset
  delta_count = eta^2/4 - O(log h)/log N.  Exact ranks (appb_exact.py h=4,5,6 -> appb_exact_h*.out, over Q): true ranks
  of thresholds == counts at the corner (h=6, N=11: (121,0),(111,1),(121,2),(121,4),...), i.e. no hidden saving.
## 2. Graded two-channel (GTC) model: what controls gamma [P, counted]
Model: X = sum_k A_k, k in [0,H] (H = pole-order span); threshold k0: rank R <= e^{F_R(k0)}, osr C <= e^{F_C(H-k0)}
(F = log of cumulative dims of intermediate states, increasing).  theta = F_R(k0)/log N, T = F_C(H-k0)/log N.
 (E) Endpoints: R=X at k0=H has rank N^2 -> F_R(H) >= 2 log N; same for F_C.  The corner (2,0) is approached
     CONTINUOUSLY iff F_R(H) = 2 log N (1+o(1)) ("tight endpoint"); otherwise the counted curve crosses theta=2 at
     T* > 0 (a jump: no corner at all, every point with theta<2 has T>=T*).
 (S) Near the corner (C keeps weights <= w): delta = [F_R(H)-F_R(H-w)]/log N ~ w F_R'(H)/log N, eta = F_C(w)/log N.
     Scaling form F(w) = log N * phi(w/H):  gamma = local exponent of phi_C at 0 (phi_C(x) ~ x^gamma), PROVIDED
     phi_R'(1) > 0 (R side has an "edge", not a bulk maximum).
 Power laws phi_R = 2x^{bR}, phi_C = 2x^{bC} (tight endpoints): frontier (theta/2)^{1/bR} + (T/2)^{1/bC} = 1,
     delta ~ bR 2^{1-1/bC} eta^{1/bC}  => gamma = bC.  bR = bC = 1/2: App B circle; bR=bC=1: the line.
     Mixed bR = 1/2, bC = 1: parabola theta^2/4 + T/2 = 1, corner slope T ~ 2 delta (gamma = 1, alpha* = 1),
     balanced theta = T = sqrt5 - 1 = 1.236 (< sqrt2).  => a "Cardy R-channel + Hagedorn C-channel" identity would give step 2.
 Concavity lemma: if phi_C, phi_R concave (log-concave dims), delta <= 2w/H <= eta: consistent with Theorem C; linear
     corner <=> phi_C(x) <= C x on [eps,1] (cheap channel EXPONENTIAL in the level from level ~0, at rate ~log N/H) and
     phi_R'(1) >~ 1 (R-channel still growing at the top).
 Single-module power laws: tight endpoint F(4h) = 2F(h) with F ~ w^b  <=>  4^b = 2  <=>  b = 1/2.  So in App-B-type
     identities (pole order linear in h, both channels in the same graded module, power-law growth) the endpoint
     condition ITSELF forces b = 1/2 and hence gamma = 1/2.  Any VOA (Cardy: F ~ 2 pi sqrt(c_eff w/6), any c) is tight
     at the endpoint and has gamma = 1/2: this is the Cardy barrier seen at the corner.
## 3. Asymmetric family A1: rank-m Heisenberg (m colours), m = kappa h  [C + N]
Same App B identity (Lemma 9 is VOA locality, Lemma 10 = Shapovalov dual bases, both hold for any rank); dims
p_m(w) = [q^w] prod(1-q^n)^{-m}; log p_m(xh) ~ h G_k(x), G_k(x) = min_t [k L(t) + t x], L(t) = sum -log(1-e^{-nt}).
Small weights: G_k(x) ~ x log(e k/x)  (near-linear: the cheap side becomes "Hagedorn up to a log") -- what B2 asked for.
BUT endpoint: theta_end = G_k(4)/G_k(1) > 2 for every k > 0 (multicolor_asym.py -> multicolor_asym.out):
  k=1e-4: theta_end 2.011, jump T* = 0.21, balanced 1.419 | k=0.01: 2.068, T* 0.47, 1.441 | k=0.2: 2.204, 0.77, 1.493
  k=1: 2.349, 0.96, 1.545 | k=5: 2.563, 1.18, 1.618.   Finite h (multicolor_count.py, exact p_m): same picture.
  => counted curve has NO corner (jump), and is worse than the circle everywhere. Asymmetric index weights E != h
  (y,u in V_E): need F(2E+2h) <= F(E)+F(h) <= 2F((E+h)/2), i.e. G(4) <= 2 G(1) at scale (E+h)/2: false for all k>0.
  Restricting indices (N' < p_m(h)) only worsens the endpoint.  Counted verdict: A1 has gamma = "jump" (worse than 1/2).
### 3a. A1 in leading order vs subleading: WHY the jump [P/N]
Leading order (pure Cardy beyond the colour scale): F(w) ~ lambda w (w <~ m), F(w) ~ c sqrt(m w) (w >> m). If the large-w
part were EXACTLY sqrt, F(4h) = 2F(h) and the small-w linear part would give a GENUINE linear corner:
  eta = w/(sqrt(k) h), delta = w/(4h) for w <= k h  =>  delta = (sqrt(k)/4) eta for eta <= sqrt(k)  (slope 4/sqrt(k)).
What kills it is the subleading Meinardus term: log p_m(w) = pi sqrt(2mw/3) - ((m+3)/4) log w - log C_m + ...,
  so F(4h) - 2F(h) ~ (m/4) log(h/m) + O(m):  endpoint excess E(k) = theta_end - 2 ~ sqrt(k) ln(1/k) / (4 pi sqrt(2/3))
  (numerics: k=1e-4: .011, 1e-3: .028, 1e-2: .068 -- matches), jump T* ~ 2 sqrt(E) ~ k^{1/4} sqrt(ln(1/k)) >> sqrt(k).
  The would-be linear piece (eta <= sqrt k) lies entirely under the jump.  Counted gamma: no corner, for every k > 0.
### 3b. Exact realisability of A1 at small sizes [E]
appb_mc.py (full pieces, ranks over Q) for (h,m) = (2,2..6), (3,2), (3,3): sum_s A_s = L X verified exactly (L = lcm z_mu).
True ranks are BELOW the counts in the middle when m >~ h, e.g. h=2: (rank R, osr C) = (31,190)@20 (m=5), (43,315)@27,
(359,49)@27 (m=6), (117,25)@14 (m=4); h=3 m=3: (123,305),(335,135)@22  -- hidden cross-weight relations exist.
corner_fast.py / corner_np.py: only C-channel weights <= w0 are needed for s0 = 2h - w0:
  C(x,y,u,v) = sum_{w'<=w0} sum_{|sig|=w'} <u| at_x[h-w'] sig><sig| b_v[w'-h] y>, at = sum_{i<=w0-w'} (-1)^i a<i>;
  cross-checked against the full computation (h=2 m=6, h=3 m=3, h=6 m=1: identical).
Corner sequences s0 = 2h - w0 (corner_runs.out, corner_runs2.out; mod 1048573, exact rows; (rank R, osr C), N^2):
  1 colour h=6..10 (N=11..42): rank R = N^2 for every w0 >= 1 up to w0=6 (h=10) / w0=13 (h=7, N=15: R drops only at
     w0 >= 14 where osr C = N^2 already): App B is TRIVIAL at the corner at all computable sizes (= counts).
  m=2 h=4,5 (k=.5,.4; N=20,36): rank R = N^2 for w0=1..7.   m=3 h=4 (k=.75, N=51): N^2 for w0 = 1..4, N^2-11 at w0=5.
  m=4 h=3 (k=1.33, N=40): N^2-0,11,35,243,745 at w0=1..5 (osr 5,19,59,164,286).
  m=6,8 h=2 (k=3,4): w0=2: (574,34)@27, (1516,53)@44 -- constant-factor savings only, eta ~ 1.
  => no small-size evidence that the counted jump is filled by hidden relations; savings appear only for k >~ 1 and
     are constant factors at eta ~ 1. (Asymptotic question open: counting cannot decide it.)
  Top piece alone (w0 = 0, osr 1): N^2 - rank = C(N',2) pattern (10, 28, 36, 55, 78, 136, 231, ...) -- the F(N,1) type
  point, far from optimal N(N-1)/2.
## 4. Other asymmetric knobs (counted) [P]
 * Unequal index weights (x,v in V_h, y,u in V_E, E = lam h): span 2E+2h, endpoint theta_end = 2 sqrt(2(1+lam))/(1+sqrt lam)
   > 2 for lam != 1 (lam=2: 2.029, 4: 2.108) -> jump; and the swap is non-square: squaring by Sigma_{N1,N2} (x) Sigma_{N2,N1}
   = SWAP_{N1 N2} is a tensor product, which only worsens corners (T3).  Unequal string lengths: same.
 * Index subsets (N' < p(h)) / "typical" states only: endpoint overshoot p(4h) > N'^2 -> jump.
 * "Temperatures" (bigradings L0 + mu*J, number of parts, colours): the pole-order constraint is additive in L0 only;
   truncating in another grading is not exact unless pieces are supported there (B6 Test B: no gain). Reweighting pieces
   does not change ranks. No effect on (theta,T).
 * Lower-weight duals a (shorter pole order): need dim V_{<=h*} >= p(h) -> h* = h - o(h) (B6 idea 4). No.
 * VOA intertwiners (a in module A, b in module B): R-channel states in A(x)M, C-channel in B(x)M -- different modules,
   but all modules of one VOA share the Cardy exponent 1/2 (only constants differ: ellipse-type curves, gamma = 1/2).
 * Infinitely-generated VOAs with g(n) ~ n^a fields per weight (higher-spin free fields, MacMahon for a=1): growth
   exponent b = (a+1)/(a+2) > 1/2 in BOTH channels -> endpoint 4^b > 2 -> jump (pole order H = wt a + wt b = 2h is forced
   because the pairing a_{(2h-1)} b is the top pole).
## 5. Synthesis: span-doubling law and the edge/bulk x Cardy/Hagedorn dichotomy [P, counted]
Span-doubling. In every locality-type (App-B) identity the pairing a_{(2h-1)}b is the top pole, so H = wt a + wt b = 2h and
 the pairwise-optimal span is W + W' = 2E + H = 4h (B6 Lemma 1-D), i.e. TWICE the information scale 2h of N^2 = (dim V_h)^2.
 Endpoint tightness F(4h) = 2F(h) then holds iff growth is sqrt (Cardy): sqrt(4) = 2. Hagedorn growth (F linear) gives
 theta_end = 4 (curve theta+T ~ 4, useless); exponent b gives theta_end = 4^b. Cardy is exactly what makes VOA crossing
 tight -- and the same sqrt makes the cheap side cost eta = sqrt(w/h): gamma = 1/2 is the price of tightness.
 Any rank profile r(h) of free bosons: excess E ~ sqrt(r/h) ln(h/r)/(4 pi sqrt(2/3)) vs linear region eta <= s = sqrt(r/h):
 jump T* ~ 2 sqrt(E) ~ sqrt(s ln(1/s)) >> s. The linear part is ALWAYS hidden under the jump (numbers in 3/3a).
Dichotomy. Parametrise the corner by one control parameter w; R-saving delta ~ w^{aR}, C-cost eta ~ w^{bC}: gamma = bC/aR.
   R side "edge" (dims increase up to the end of the span, nonzero slope): aR = 1.  R side "bulk" (error = tail mass of a
   statistic whose maximum is interior, e.g. binomial): aR = 2 (Gaussian).
   C side Cardy: bC = 1/2.  C side entropic/Hagedorn (cost C(k,w) c^w or e^{lambda w}): bC = 1 (up to log).
   App B = (edge, Cardy): gamma = 1/2.   Bell/singlet family of G27-T3 = (bulk, entropic): gamma = 1/2 (with 1/log^2 loss).
   Linear corner needs (edge, Hagedorn) [or (bulk, bC = 2): windows of width sqrt k at O(1) cost -- not plausible].
   The (edge, Hagedorn) target is NOT produced by any graded/locality construction known: locality forces both channels
   to be modules of one VOA (same exponent) and span doubling forces that exponent to be 1/2.
Nested-shell reformulation of (edge, Hagedorn) in commuting models: filtration 1 = g_0 >= g_1 >= ... of commuting projectors,
 shells g_j - g_{j+1} inside Sym (j even) / Anti (j odd), tr g_j <= rho^j N^2, osr(g_j) <= c^j  =>  K_d = sum_{j<=d} (-1)^j
 (g_j - g_{j+1}) + ... has q <= sum c^j, p <= tr g_{d+1}: delta/eta = log(1/rho)/log c (LINEAR). Necessary: gauge norms Q(g_j)+Q(g_{j+1}) >= (1-rho) rho^j N (trace argument: a projector P <= Sym or <= Anti with
 osr decomposition sum G(x)G' has |tr P| = |tr(SWAP P)| = |sum tr(G G')| <= N Q), so Q >= N^{1-delta} at depth d.
 Bell model: uncertainty |supp 1^_S| >= 4^k/|S| forces c >= 1/rho; equality needs affine shells, and affine shells are
 IMPOSSIBLE already at level 1 (Q = c.d constant on S0 \ S1 with S1 affine of codim r forces Q linear on S1, so dim S1 <= k+1,
 density <= 2/N, contradicting shell 0 inside {beta=+1} of density ~1/2). Non-affine shells pay extra (open).
### 3c. A1 jump is not a finite-size artefact [N]
Exact p_m counts, min T with theta <= 2 (jump): k=0.05: h=200 .767, h=500 .685, h=inf .622; k=0.2: .825, .793, .768;
k=1: .985, .973, .965 (multicolor_count.out, multicolor_asym.out). Converges to the saddle-point limit.
Hidden-saving trend at fixed k = 1 (balanced-ish s0 = 0): N=5 true theta 1.21 vs count 1.80; N=22: 1.557 vs 1.705 --
the gap shrinks with size. Heuristic reason the jump is real for true ranks: with many colours the R-channel
"Taylor coefficients" (b_{(-1-j)}u in V_{2h+j}, dim >> N^2) determine b (x) u long before the end of the span, so
R_W reaches rank N^2 early (no "edge"): exact data k <= 1 show rank R = N^2 for small w0 (3b).
## 6. One-sentence mechanism
Edge (tight R endpoint) needs R-channel dims to reach N^2 only at the END of the span; Hagedorn C needs fast growth from
weight 0. In a locality identity both channels are the same module and the span is 2x the information scale, so the
only compatible growth is exactly sqrt (Cardy): the corner exponent of every App-B-type identity is 1/2 (counted);
faster growth (colours ~ h, higher-spin fields, large-c holomorphic VOAs with Cardy offset c/24) gives a jump instead.
Holomorphic c = 24k VOAs (h = lam k) [H]: F(w) ~ 4 pi sqrt(k(w-k)) -> theta_end = sqrt((4lam-1)/(lam-1)) > 2, after the
jump delta ~ nu eta^2 / (4(2+nu)): quadratic again.
## 7. VERDICT
gamma achieved (realisable): 1/2, App B (circle, delta = eta^2/4 + eta^4/64 + ...): PROVED (paper's costs + HR asymptotics;
 exact ranks = counts at the corner for h <= 10). Bell/singlet family (T3): 1/2 up to 1/log^2 (exact formula).
No App-B-type family with gamma > 1/2: power-law GTC => gamma = growth exponent of cheap side and tight endpoint forces 1/2
 [P, counted]; asymmetric/rank-growing A1 has a counted jump (numeric, converged), true ranks at small N show no fix [E].
Counted-only candidate with gamma = 1: "Cardy R-channel + Hagedorn C-channel" (parabola theta^2/4 + T/2 = 1); no
 realisation known; ruled out inside VOA locality identities by span doubling. Commuting analogue = nested shells with
 geometric traces and exponential osr; affine shells impossible in the Bell model [P].
## 8. Remarks
- Concavity alone does NOT force sqrt: F linear on [0,kh] then slope 1/3 is concave with F(4h)=2F(h) and would give a linear
  corner (delta <= eta/3 is the concavity bound). The obstruction is which modules exist, not counting.
- In modular-invariant settings a Hagedorn window exists only below the Casimir offset (w <~ c/24 chiral, c/12..c/6 HKS);
  the same offset makes F(4h)/F(h) = sqrt((4h - w0)/(h - w0)) > 2: Hagedorn region <=> endpoint excess <=> jump.
- Modules M for y,u (instead of V): span 2(E-Delta)+2h, same Cardy tightness iff E-Delta = h. Nothing new.
Files: appb_curve.py(.out), appb_exact.py (copied from gen26/b6), appb_exact_h{4,5,6}.out, appb_mc.py, appb_mc_h*_m*.out/json,
 corner_fast.py, corner_np.py, corner_runs.out, corner_runs2.out, multicolor_count.py(.out), mc_count2.py, multicolor_asym.py(.out).
Processes: all mine finished (10:11 UTC). (ap_sym.py processes belong to another agent.)
