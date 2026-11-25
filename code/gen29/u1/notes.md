# G29-U1 notes (independent adversarial review of Lemma U / alpha* / corner criterion)
Start 12:49 UTC. Not yet read: gen27/r1 (first review) -- read only in step 3.

## Step 1. Own derivation of one unbalanced level (from the paper text, Prop 3 + Lemma 8 + Def 7 + Thm A(b))
Level: (K', lam)-clique, K' = g*h nodes on a g x h grid; N = K'^2 lam = g^2 h^2 lam.
Rectangles |R| = C(g,2)C(h,2) ~ g^2h^2/4; 4|R|lam/N = (1-1/g)(1-1/h) <= 1.
Row network: g nodes, setting = h(gh-1)lam <= P_r := g h^2 lam = N/g bits; c_{a,a'} = h lam + sum_{q in rows a,a'} c^row
  <= h lam + C(h,2) 2lam <= h^2 lam = P_r/g  (Lemma 8(ii) hypothesis OK). U_row = g(g-1)h lam + 2 sum_q c^row ~ N(1/h + ubar).
Column network: h nodes, P_c = g^2 h lam = N/h, c_{j,j'} <= g^2 lam = P_c/h OK. U_col ~ N(1/g + vbar).
(ubar, vbar) = average over rectangles of (c^row_i, c^col_i)/(2 lam): any point of conv(piece set) (dense as |R| grows);
 worst iota is uniform. Prop 3 text: "this part of that proof works for every g >= 2" -> nothing uses g = h.
Lemma 8(ii) on row network with K_a groups: needs K_a | g, 2 zeta_a^2 g >= K_a ln(2 g K_a)  => g = t_a K_a,
  t_a >= ln(2gK_a)/(2 zeta_a^2). L_a = ceil((U_row + 2 zeta_a g P_r)/K_a^2), gP_r = N; entry 2^{2gP_r/K_a} = 2^{2N/K_a}.
Column: same with h = t_b K_b, entry 2^{2 h P_c/K_b} = 2^{2N/K_b}.
Per-iota exponent / N:
  min( e_a (ubar + 1/h + 2 zeta_a) + 2/K_a ,  e_b (vbar + 1/g + 2 zeta_b) + 2/K_b ) + O(K'^2 log(lam+B))/N
#iota = |I|^{|R|} = (2lam)^{c K'^2/4}: fits Thm A(b)'s Xi^{C0 K'^2} form. Choice row/col per iota needs computable
upper bounds on e_a, e_b (fine: schedule numbers). Memory: not an issue (time-only statement, Lemma 8 recomputes entries).
x' = log K' = x_a + x_b + log t_a + log t_b.
=> e' = max_{(u,v) in H} min(e_a u, e_b v) * (1 + O(zeta + 2^-x)) + O(2^-x_min). Lemma U: CONFIRMED (my derivation).
Constraints found: (1) the grouping factor t >= ln(..)/zeta^2 on BOTH sides (adds O(log x) to x per level, so the
small side needs x_b >> log x); (2) entry term 2/K_b needs x_b > log(1/e') ~ alpha log x; (3) along a long big-side
chain (depth ~ x/x_b) need sum of zeta_a <= O(1) -> zeta_a ~ x_b/x -> log t_a ~ 3 log x (still O(log x)).
(4) child lam is huge (L >= m^{7/9}) at every level, so o(lam) slack in the piece costs is harmless.
(5) Square grid NOT needed; children need NOT be comparable in size.

Frontier -> pieces: chain of realizable (R_s, C_s), R_0 = 0, R_S = X; A_s = R_s - R_{s-1} = C_{s-1} - C_s:
 row rank <= p_s + p_{s-1}, col rank <= q_{s-1} + q_s -> piece costs (theta_s/2, T_{s-1}/2) + O(1/lam): staircase
 corners above the curve; fine chain -> hull = region under concave majorant of T(theta). (Integrality/computability
 of a rank-p decomposition as an integer factorization through log p bits is assumed: design level.)
Swap symmetry: X(x,y,u,v)=[x=v][y=u] is invariant under y<->u which exchanges the two splits -> frontier symmetric.

Duality (own): max_H min(a u, b v) = min_nu h(nu)/(nu/a + (1-nu)/b) (minimax, H convex). Psi = 1/e:
 Psi' = max_nu (nu Psi_a + (1-nu) Psi_b)/h(nu).
 Upper half (power <= alpha if h >= ||.||_q all nu, q=1/(1-alpha)): Hoelder -> Psi'^{1/alpha} <= Psi_a^{1/alpha}+Psi_b^{1/alpha}
   and x' >= x_a + x_b -> Psi <= (C x)^alpha.   [needs leaves with Psi^{1/alpha}/x bounded: brute force Psi=2 OK]
 Lower half: if h(nu0) < ||(nu0,1-nu0)||_q, self-similar split s/(1-s) = (nu0/(1-nu0))^{1/(1-alpha)} gives margin
   ||.||_q/h(nu0) = 1+eps > 1, absorbs (1+o(1)) overheads and the O(log x) additive x loss -> power >= alpha.
 Corner slope c (v <= c(1-u)): at nu = c/(1+c) h(nu) = nu -> Psi' >= Psi_a + Psi_b/c (exactly, no Hoelder).
   With brute-force small child Psi_b = 2 at x_b ~ C log x: slope dPsi/dx ~ 2/(c C log x) -> e(x) = O(c log x/x).
   With recursive small children of size x_b ~ (log x)^2: F(x) := x/Psi(x) satisfies F(x) ~ c F(x_b) -> F ~ c^{log* x}
   (to be checked numerically) -> the stated rate looks correct but NOT tight.
 gamma: alpha(nu) -> gamma as nu -> 1 (h = nu + C eps^{1/(1-gamma)}, ||.||_q = nu + eps^q/q). alpha* = 1 iff some
   nu in (0,1) has h(nu) = max(nu,1-nu) (finite slope) OR alpha(nu) -> 1 at the corner. CAVEAT: gamma = 1 without a
   finite slope (e.g. T ~ (2-theta) log(1/(2-theta))) gives alpha* = 1 only as a supremum: time 2^{m/(log m)^{1-o(1)}},
   not 2^{O(m loglog m/log m)}.  So "step 2 <=> gamma = 1" needs the reading "finite slope" for the strong form.

## Step 2. Own simulation (sim.py DP on integer x grid with honest overheads; chain.py ODE big-side chain; formula2.py)
Honest overheads: zeta menu, t = ln(2gK)/(2 zeta^2) both sides, entry 2^(1-x), 1/h, 1/g; leaves = brute force e = 1/2.
NOTE (hidden constraint found while coding): at the corner point the thin side has U_col ~ v N with v ~ e_a/e_b << 1,
 but Lemma 8(ii)'s grouping error is 2 zeta N (relative to N, not U): needs zeta_b << e_a/e_b ~ x_b/x, so
 log t_b ~ 2 log(x/x_b) + log x. Still O(log x): no effect on alpha, only on constants / minimal x_b.
DP X=3e4 (out in buxxzcup0 output copied below), local alpha over [1e4,3e4]:
  circle 0.497 (sqrt(x) e -> 3.91), line 0.995 (x e -> 28.7), mono (p^2q=N^4) 0.911 rising,
  circle+corner slope 4: 0.872 rising (circle alone 0.497).  x e / log2 x at 1e4,3e4: mono 6.4, 6.3; corner4 12.4, 12.8
  -> ratio ~2 = ratio of slopes c (4 vs 2): e ~ C c log x / x in this range.
chain.py (mono, x up to 1e10, x_b <= 5% of x): recursive small children: x e/log2 x = 6.12 (1e5) -> 4.16 (1e10), x e
  101 -> 138 (slowly growing); brute-force small children only: x e/log2 x ~ 6.5-6.8 flat (= Theta(c log x/x) exactly).
  => rate O(c log x/x) CONFIRMED as upper bound; with recursive small children it is not tight (my F ~ c^{log* x}
  heuristic: x e grows much slower than log x). Numerical breakdown of chain.py beyond x ~ 1e10 (not chased).
formula2.py (precise, symmetric families, nu in [1/2,1)): circle 0.5 everywhere; line 1; mono & corner c=2,4,10: alpha(nu)=1
  at nu = c/(1+c) up to float cap (0.996-0.998); alpha(1/2) = log2(1/rho) (0.585 mono, 0.5 circle+corner);
  power corner gamma=0.75: alpha(nu) -> 0.7456 at nu = 1-1e-12 (-> gamma, slowly).

## Step 3. Read gen27/r1/notes.md (13:05 UTC) -- comparison
Agreements: level map; Lemma U incl. extreme unbalance; hidden constraint "Hoeffding error 2 zeta N relative to N" on the
 thin side (R1 found it too); alpha <= 1 trivial (endpoint (2,0)); Hoelder lower half valid for all trees; self-similar
 upper half; on a slope-c corner Psi' = Psi_a + Psi_b/c exactly; brute-force small child -> e = Theta(c log x/x)
 (R1: e x/log2 x ~ 7.9 for c=2; mine ~6.5-6.8 for mono c=2 -- overhead conventions differ); recursive small children
 beat c log x/x (R1 int DP to 1e6; my chain ODE to 1e10: 6.12 -> 4.16); B2's "step 2 <=> finite slope" only <=.
Correction to MY notes: the entry term 2/K_b must be << the per-level GAIN e_a^2/(c e_b), not << e', so x_b >= 2 log2 x
 + O(1) (R1 is right; DP picks x_b ~ 49 at x = 3e4).
R1 extra I verified by hand: h_P >= h_F for ANY piece family (threshold split J = {u_s <= t}) -> lower half applies to
 every piece system, not only chains. Correct.
Lemma U (P9) constant: entry term is 2/K_a or 2/K_b (Lemma 8: 2^{2gP/K} = 2^{2N/K}), P9 writes 1/min(K_a,K_b): factor 2,
 lower order only.
DISAGREEMENT with R1: R1 says alpha*=1 <=> gamma=1 is EXACT with gamma := liminf_{s->0} log v_hull(s)/log s.
 Counterexample (oscillating hull): vertices (s_k, s_k^g), g<1, s_k = s_{k-1}^{M_k}, M_k -> inf, plus (0,0): concave,
 v <= s^g everywhere and = at vertices, so liminf = g < 1. Excess(t) := h(1-t)-(1-t) = max_j (t s_j^g - s_j) (vertices);
 at t = s_{k-1}^{1-g}/2 all j <= k-1 give negative values, so excess(t) <= t s_k^g = t (2t)^{M_k g/(1-g)}, smaller
 than t^q/q for any fixed q once M_k is large -> alpha(nu) -> 1 along nu_k -> 1, alpha* = 1 while gamma_liminf = g.
 Correct general statement: alpha* >= gamma_liminf (so gamma_liminf = 1 => alpha* = 1); alpha* = 1 <=> some nu in (0,1)
 with h(nu) = max(nu,1-nu) (finite slope) or liminf_{t->0} log excess(t)/log t = +inf. Equivalence with gamma holds when
 log v/log s has a limit (regular corner) -- all physically relevant cases. Pathological, but R1's "exact" is too strong.

## Final runs (13:20 UTC)
chain.py circ_corner4 (circle + slope-4 corner): x e/log2 x = 12.7 (1e5) -> 10.4 (1e10); ~2x the mono (c=2) values -> linear in c.
sim.py circ_pow0.75 (gamma = 0.75) X=1e4: local alpha 0.42, 0.57, 0.62, 0.67 (rising toward <= 0.75; formula alpha(nu)->0.75).
Verdicts: Lemma U CONFIRMED (constant: entry 2/K not 1/min K); unbalanced level cost CONFIRMED (+ thin-side zeta, x_b >= 2log2 x);
alpha* lower half CONFIRMED (any tree; any piece family via R1's h_P >= h_F); upper half CONFIRMED (design level);
step2 <=> gamma=1 CONFIRMED WITH CORRECTION (regular corners; liminf-gamma fails for oscillating hulls; B2's finite-slope
version only <=); rate O(c log x/x) CONFIRMED (upper bound, tight only with brute-force small children).
No processes left running. Nothing written to /home/user/pnp-recherche.
