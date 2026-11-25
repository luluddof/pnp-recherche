# G27-R1 review notes of G27-B2 (start 08:08 UTC)

## Reading
- Paper Prop 3 / Thm A: rectangle crossing on 4 strings of lambda bits each (all edges carry lambda bits); row split (x,y)|(u,v),
  column split (x,u)|(y,v). In a g x h grid nothing changes in the crossing: same N = 2^lambda for both splits. Pieces are
  single-N objects. Row network: g nodes, P = g h^2 lam, gP = N, c_{a,a'} <= h lam + C(h,2) 2lam <= P/g (Lemma 8(ii) ok).
  U_row = N/h + 2 sum_q c_row (straight edges between rows: h per pair). Entry cost 2^{2N/K_a}. Matches P9 Lemma U.
- Choice rows/cols is per iota (global over grid). Worst iota -> average of piece points -> conv hull. Level map
  e' = max_{conv P} min(e_a u, e_b v). OK.
- G25-P9 Pythagorean law: derived FOR THE CIRCLE (App. B costs on u^2+v^2=1). P9 itself says power change needs frontier
  strictly inside circle "for some weighting". So no contradiction with B2: B2 generalises.

## Own DP (r1_dp.py, geometric level map, rigorous downward rounding on log grid r=1.02) [08:15]
circle -> 0.500, line -> 1.000, lp1.5 -> 0.667 (exact). p2q PLATEAU 0.860, circ+lin0.3 plateau 0.633, pow0.75 0.613,
same with overheads (ovh=1: zeta=1/x, log t, entry 2^{1-xmin}). Plateau = GRID ARTIFACT: on a ratio-r grid the big child
is <= x_k/r so the small child is >= (r-1) x: kappa >= 0.02. Predicted plateau from kappa^{alpha-1} = alpha c (c=corner slope):
p2q c=2 -> 0.865, circ+lin0.3 c=6.6 -> 0.644. Matches. (B2's dp had a similar artefact, they noticed.)
## Analytic: on the corner segment v=c(1-u): Phi' = Phi_a + Phi_b/c EXACTLY (Phi=1/e). So with tiny child y:
Phi(x) >= Phi(x-y) + Phi(y)/c -> Phi linear in x (overhead-free). With overheads y ~ C log x needed -> e(x) = x^{-1} (slow factor).
=> B2's self-similar (kappa fixed) route is far from optimal; the corner gives e = x^{-1+o(1)} with small losses. To test.

## Integer-x DP allowing tiny children (r1_int.py) [08:20]
ovh=0 (pure level map): p2q gives e*x -> 1.000 exactly (Phi linear, alpha = 1 ATTAINED, child y=1). line e*x=0.5.
ovh=1 (zeta = x^-1.5, log t ~ 4 log2 x per child, straight 2^-y, entry 2^{1-y}):
  x:        1e2   1e3   1e4   1e5          (e*x)
  circle    35    128   420   1335   (slope 1/2)
  line      25    32    34.6  35.0   (slope 1)
  p2q       33    73    106   127    (local slope 1e4->1e5: 0.92, increasing)
  circ+lin0.3 35  128   319   436    (local slope 0.865, increasing)
  best small child y ~ 64 = O(log x)  -> NOT self-similar.
Heuristic: dPhi/dx ~ Phi(y)/(c (y + O(log x))) with y ~ C log x -> e(x) = L(x)/x, L(x) ~ c' L(C log x) (iterated-log growth).

## Hand checks
- Level map: e' = max_H min(e_a u, e_b v) = s* (ray exit) = min_mu h(mu)/(mu/e_a+(1-mu)/e_b) (minimax) -- OK.
- Hoelder lower bound: mu Phi_a + (1-mu)Phi_b <= ||mu||_q ||Phi||_{1/a}, so h >= ||.||_q for all mu => Phi'^{1/a} <= Phi_a^{1/a}+Phi_b^{1/a},
  x' >= x_a+x_b -> e >= (mu0 x)^{-a}. Valid for ALL trees (not only self-similar). OK.
- Achievability: self-similar e = C x^-a works with slack iff exists (mu,w): h(mu) < mu w^a+(1-mu)(1-w)^a, max_w = ||(mu,1-mu)||_{1/(1-a)}. OK.
- MISSING in B2: lower bound must hold for ANY piece family, not just chains. Proof (mine): for M = h_P(mu), threshold split
  J = {s: u_s <= t} gives realizable (t, (M - mu t)/(1-mu)) up to log S/lambda; so M >= max_u [mu u + (1-mu) T_F(u)] = h_F(mu).
  So h_P >= h_F; with chains h_P -> h_F. alpha* depends only on the concave majorant of the frontier. OK.
- alpha <= 1 is TRIVIAL (frontier contains (2,0): h(mu) >= mu); Theorem C not needed.
- B2 statement "alpha* = 1 iff inf_mu [h(mu) - max(mu,1-mu)] = 0" is VACUOUS as written (=0 at mu=1 always).
- B2 "step 2 <=> finite slope": only <=. alpha*=1 <=> gamma=1 with gamma := liminf_{s->0} log v_hull(s)/log s (exact, Legendre
  argument: h(1-t)-(1-t) = max_s[t v(s) - (1-t)s] = o(t^q) for all q iff gamma = 1). E.g. T = d log(1/d): alpha* = 1 not attained.
- Finite slope: STRONGER than B2 says: on the corner segment Phi' = Phi_a + Phi_b/c exactly -> e = L(x)/x with slowly growing L
  (tiny children y ~ log x), not just alpha -> 1 via kappa ~ c^{-1/(1-alpha)}.

## Explicit chain, brute-force small child (r1_chain.py, chain.out) [08:30]
Corner slope c (H inside v <= c(1-u)): Phi' >= Phi_a + Phi_b/c for ALL e_a,e_b (checked numerically on 1e5 random pairs, p2q, line).
Chain: big child = previous level, small child = BRUTE FORCE (e_b=1/2) with y = 2log2 x + log2 c + 4, t = x^3 (zeta ~ 1/x),
x' = x + y + 6 log2 x, all Lemma-8(ii) overheads: e*x/log2 x -> 3.9 (c=1), 7.9 (c=2), 26 (c=6.6), 200 (c=50), flat to x=1e9.
=> finite corner slope c gives e(x) <= ~4c log2(x)/x: time 2^{O(c m loglog m / log m)} with a TRIVIAL small child;
recursive small children do better (int DP: e*x grows slower than log x). B2's "slow approach kappa ~ c^{-1/(1-alpha)}" is an
artefact of restricting to self-similar splits.
## Lemma U in extreme unbalance: Hoeffding term of Lemma 8 is 2 zeta N (relative to N, not U): column side needs
zeta_b << e_a/(c e_b) ~ 1/x^2 -> t_b ~ x^4 polylog, log t = O(log x) = o(x). Entry 2/K_b << e_a^2/c -> x_b >= 2 log2 x + O(1).
Straight 1/h, 1/g fine. c_{j,j'} <= P/h holds for every iota. Single-N pieces suffice (all strings lambda bits). OK.

## p2q integer DP with overheads to 1e6 (int1e6_p2q.out) [08:39]
e*x: 1e4 106, 1e5 127, 1e6 143.6 (decade slope 0.946, still rising); circle at 1e6: 4224. Small child y=135 ~ 7 log2 x.

## Rate refinements (analytic)
- Corner modulus: v_hull(s) = s*omega(s). Tiny child: gain Phi_b/omega(Phi_b/Phi_a) -> e(x) <~ omega(1/x) log x / x.
  finite slope: e <= ~4c log2 x/x (proved by explicit chain modulo design-level overheads), heuristic c^{O(log* x)}/x with recursive
  children; omega = log(1/s) (gamma=1, infinite slope): e ~ (log x)^2/x; gamma<1: x^{-gamma} log x (consistent with alpha* >= gamma).
- Exact e = Theta(1/x) only shown for c = 1 (the line). With c>1 the half-plane bound Phi' = Phi_a + Phi_b/c is tight, overhead
  forces y >= 2log2 x: L(x) ~ c L(C log x).
## Reconciliation with G25-P9: P9's Pythagoras (Prop P(i)) uses the HR asymptotics of App-B widths = circle only; P9 itself says
 a power change needs a frontier inside the circle "for some weighting". No contradiction; B2 = correct generalisation.
 Brief's "4 rho = sup(theta+T)" / "line => step 2, only route" is the balanced special case.
## Not checked: finite_chain.py entropy claims, targets.py numbers, Bell-diagonal remark, OpenAI-side Q2, route list R4-R12.
## Processes: all mine killed/finished (08:39).
