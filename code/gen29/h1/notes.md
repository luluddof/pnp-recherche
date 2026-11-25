# G29-H1 notes (start 12:49 UTC). Full Bell model / Walsh-Hadamard low-rank rigidity.
Tags: [L] literature (source given), [P] proof sketch, [E] exact computation, [N] numeric, [H] heuristic.
Plan: 12:55-13:25 literature; 13:25-14:30 constructions; 14:30-15:05 exact small optima; report by ~15:15.
## Reformulation [P]
phi = beta + E, |supp E| = p, |supp phi^| = q  <=>  beta = (q-Fourier-sparse) - (p-sparse). Self-dual.
Linear corner <=> error fraction eps = p/N^2 <= q^{-c} uniformly in k (sparsity poly(1/eps), independent of k).
Symmetries: z -> g z (g in O^+(2k,2), preserves Q(c,d)=c.d), and twisted translations phi'(z)=beta(t)phi(z+t)chi_{Jt}(z):
supp phi^ -> supp phi^ + Jt. So Fourier supports matter up to AGO^+(2k,2).
## 13:10 Literature [L] (WEB SEARCH SNIPPETS ONLY: arxiv.org, ar5iv, joshalman.com do not resolve in this sandbox; no paper read)
- Alman-Williams STOC'17 (arXiv 1611.05558): over any field, changing 2^{eps n} entries/row of H_n (N=2^n) drops rank below
  2^{n(1-Omega(eps^2/log(1/eps)))}. Abstract also: "non-trivial rigidity upper bounds for H_n with smaller target rank" -- exact
  small-rank statement NOT retrieved. In our coordinates their regime is eta -> 1 (rank N^{1-f}, p = N^{1+eps}); the mechanism is
  the same polynomial method (window around |x AND y| ~ n/4) as the singlet family => quadratic.
- Rank 1: R_{H_n}(1) = 2^{2n-1} - 2^{1.5n-1} (n even), 2^{2n-1} - 2^{1.5(n-1)} (n odd) (Alman et al., arXiv 2211.05217 Thm 5.3, snippet).
- Alman-Liang STOC'25 (2502.19580): R_H(c1 log N) >= ~N^2/2 over finite fields (tight); previous R(c log N) >= cN.
  Alman 2026 (2608.06592): over F_3, 1% changes cannot reduce rank to <= log^2 N/80.
- Kivva MFCS'21 (2103.05631), Alman STOC'21 (2102.11992): Kronecker products non-rigid (high-rank regime, Valiant parameters).
- Servedio-Viola, "On a special case of rigidity" (ECCC TR12-144, 2012): EXACTLY the sparse-polynomial special case of rigidity
  for IP: the sign of a polynomial with <= s monomials (+-1 / parity basis, as far as snippets say) on 2n vars disagrees with IP on
  >= Omega(1/s) fraction; they note such bounds "seem unknown in the rigidity setting"; progress on this special case => IP not in
  small AC0 o MOD2. => our "line" eps >= 1/(4q) is their Omega(1/s); beating it (delta < eta) is their open direction.
  Reverse direction [P, my remark]: a size-S depth-d AC0oMOD2 circuit for IP gives (via real probabilistic polynomials of degree
  t=(log S)^{O(d)} log(1/eps) in the parity gates; product of characters = character) a Fourier-sparse exact eps-approximant with
  q <= S^t, i.e. eps = q^{-1/polylog(S)}: corner slope -> 0 like 1/polylog(k) if S = poly(k). Not a linear corner.
- Bindu-Hatami-Karimi-Robere 2026 (arXiv 2605.01038): eps-approximate SIGN-rank of 2^m x 2^m Hadamard <= m^{O(sqrt m log(1/eps))}
  (adapting Alman-Williams) and >= Omega_eps(m). So even the sign version at constant eps is open between log N and N^{o(1)}.
- Alman-Chattopadhyay-Williams ITCS'25: probabilistic-polynomial sparsity lower bounds (AND of ORs; generic Omega(2^d) from prob.
  degree d) -- monomial basis, NOT Fourier basis; IP-relevant statement not found.
- Clean equivalence [P]: Bell (p,q) <=> exists probabilistic Fourier-q-sparse exact representation of IP_k with error eps=p/N^2
  (random twisted translation one way, averaging the other). Linear corner <=> eps <= q^{-c} uniformly in k.
## 13:20 Short proofs / remarks [P]
- Span lemma sharpened: if phi = g(Lz) with ker L = W != 0 then agreement <= 3/4 (B|_W != 0: fibre bias <= 1/2; W isotropic:
  only a 2^{-dim W} fraction of fibres has beta constant, the rest are balanced). So eps < 1/4 forces span(S) = F2^{2k}, q >= 2k.
- Line via twisted translations: F_t = phi(.+t) - beta(t) chi_{Jt} phi has spectrum <= 2q and vanishes on Z cap (Z-t)
  (density >= 1-2eps); if 2eps < 1/(2q) all F_t = 0 => phi = const*beta. So eps >= 1/(4q) (same as the line).
- Row view: row c is an approximation of chi_c by spectrum pi_2(S); per row the line is tight (1 - 2^t 1_A), but making the
  row coefficient functions g_b(c) sparse forces the trivial cost. The corner needs eps = q^{-c} with SMALL c (above the line).
- Block variant of singlet (sites = blocks of b pairs, variable y_B = IP_b(block) in {0,1}, cost 4^b, Pr = 1/2): heuristic
  optimisation over b gives delta ~ eta^2/log(1/eta) instead of eta^2/log^2 -- still quadratic (covered by S' with d = 2^b).
- [N] block_singlet.py (k=600, exact log-binomial formulas, upper bound on q): blocks b=1..6 give delta/eta^2 within ~+-30% of
  b=1 (e.g. eta=.1: .072,.079,.078,.092,.080,.089; eta=.5: .034..030). The heuristic log-factor gain is NOT visible at
  realistic k; block variant = singlet up to constants. (k=2400 run killed: too slow, not needed.)
## 13:55 Constructions, first results
- [E] singlet_front.py: singlet frontiers (q,p): k=4: (13,67),(67,13),(120,1); k=5: (16,349),(106,106),(376,16); k=6: (18,1423),
  (154,694),(694,154),(1909,19).
- greedy.py (OMP-like support growth + reweighted-L1 (HiGHS) coefficients + EXACT rational re-solve & count = rigorous points):
  k=4: (q,p) = (4,112), (10,46)  [(10,46)@16 beats singlet (13,67)]; k=5: (22,332), (40,184) ~ on singlet log-interpolation.
  (10,46)@16 has all 10 coefficients +-1/2: phi = (1/2) sum eps_i chi_{s_i}.
- Reformulation [P]: chi_s * beta = beta(Js) beta(. + Js)  =>  phi = beta on Z  <=>  (w * beta)(z) = 1 on Z, w = sum kappa_s beta(Js) delta_{Js}:
  a q-sparse measure whose convolution with the bent function is constant 1 off an eps-set. Translates of beta are pairwise
  uncorrelated, so a "generic" T gives S(z) = sum w_t beta(z+t) ~ Gaussian: needs strong algebraic structure.
## 13:15 Structure of the (q,p)=(10,46)@16 greedy optimum [E]
kappa_s = (1/2)(-1)^{Q(s)} for all 10 s (truncated Fourier series of beta, rescaled by N/2). S = A u C, A = 7 vectors pairwise
B=1 (symplectic form), C = 3 nonzero vectors of a nondegenerate plane (sum 0), A perp C.
Identity [P]: (-1)^{Q(s)} chi_s(z) = beta(z) beta(s + Jz)  (Q(s+Jz) = Q(s) + s.z + Q(z)).  Hence for UNIFORM-type
phi = (1/c) sum_{s in S} (-1)^{Q(s)} chi_s:  phi(z) = beta(z)  <=>  f_S(w) := sum_{s in S} (-1)^{Q(s+w)} = c, w = Jz.
=> "design" subproblem: a set S whose translates S+w all contain the same number (|S|+c)/2 of singular points, for all but
eps of w.  (Weighted version = kappa * beta = const off eps; the general Bell model.)
## 13:13 MAJORANA FAMILY (new) [E + P]
Fact [E, k<=7; P: standard]: for a maximal anticommuting set A (2k+1 Jordan-Wigner Majorana vectors, pairwise B=1) the bits
y_s(w) = Q(s)+B(s,w) satisfy beta(w) = F(J mod 4), J = sum_s y_s (fixed parity).  I.e. IP_k is a SYMMETRIC function of 2k+1
linear forms (the k=4 greedy optimum (10,46) is exactly this, up to O(Q)).  Family: phi = f(u), u = (J-par)/2 in {0..k},
f interpolates F on a window of d+1 consecutive u; q <= sum_{j<=d} C(2k+1,j) (product of j Majorana characters = ONE
character), p = sum_{u notin W} C(2k+1, 2u+par).  majorana_front.py (exact p, q by WHT):
 k=2: (1,6),(6,1) | k=3: (8,8) | k=4: (10,46),(46,10) | k=5: (12,232),(67,67),(232,12) | k=6: (14,1093),(92,378),(378,92)
 k=7: (16,4944),(121,1941),(576,576),(1941,121) — balanced 576@128 vs singlet 1156, HR tiling 1024.  Self-dual frontier.
Compare singlet: (10,10)@8, (13,67)@16, (106,106)@32, (154,694)@64, (1156,1156)@128.  STRICTLY better at every k tested.
Corner: still QUADRATIC [P sketch]: phi symmetric in the 2k+1 modes => characters of Majorana weight w have orbit C(2k+1,w);
q < C(2k+1,u+1) forces phi = P(J) with deg P <= u (weights identified j <-> 2k+1-j); P - (+-1 pattern with period 2 in u)
cannot vanish on u+2 consecutive u-values ((u+1)-th difference), so errors sit within distance ~u/2 of the centre of
u ~ approx N(k/2, k/8): eps >= exp(-c u^2/k)/poly => delta <= C eta^2/log^2(1/eta).
- [N] maj_asym.py k=400 (exact log-binomial masses, q upper bound): delta(eta) Majorana / singlet = 1.24 (eta=.05), 1.25 (.1),
  1.44 (.2), 1.42 (.3), 1.65 (.5), 1.79 (.8); delta/eta^2 Majorana .265,.117,.075,.064,.062,.070 (slowly decreasing at small
  eta = the log^2 factor).  Constant improvement only.
- Generalisation [P/H]: "code statistic" class: n linear forms (vectors s_i, columns of G), J = wt(Gw + c); if beta = F(J) then
  degree-d interpolation costs <= C(n,<=d) characters. Distinct nonzero s_i => the bits are pairwise independent =>
  Var J = n/4 >= (2k+1)/4 for any spanning family; Majorana (n = 2k+1, even-weight code) is the minimal-length member.
  Beating Gaussian tails would need codes of dim 2k, length n = 2^{O(eta k)} whose coset-weight distribution is concentrated
  on few weights up to mass 2^{-c eta k} AND determines beta there (few-weight codes from quadrics have n ~ N^2: trivial).
  Not found.
- Exact N=8 (exact8.py, all O^+(6,2)-orbits of supports with 0 in S, exact max hyperplane incidence, rational check):
  q=1: 28 | q=2: 28 | q=3: 28 | q=4: 24 (S = one qubit's {I,X,Z,Y}) | q=5: running in 4 parts.
## 13:25 Majorana balanced points (maj_balanced.py/.out; exact formula, = WHT count for k<=7) [E for k<=7, formula beyond]
 k : Majorana bal | singlet | HR tiling N*ceil(N/rho(N)) | N^2/8
 5 : 67 | 106 | 128 | 128      6 : 378 | 694 | 384 | 512      7 : 576 | 1156 | 1024 | 2048
 8 : 3214 | 7459 | 4096 | 8192    9 : 5036 | 12826 | 14848     11: 44552 | 143980 | 176128    13: 397594 | 1630084 | 2588672
 15: 3572224 | 18575560 | 33554432    31: 1.73e14 (exp 1.526) | 6.0e15 (1.691) | 7.2e16
 => As crossing pairs for the GENERAL model (Bell is a submodel: rank(SWAP-K)=p, osr(K)=q), the Majorana family gives the best
 known balanced points for all N = 2^k, 5 <= k <= 31 except k=4 (46 vs HR 32).  Asymptotic balanced exponent [P, entropy calc]: d = tau(2k+1): q ~ 2^{2k H(tau)}, p/N^2 ~ 2^{-2k(1-H(1/2-tau))}
 => exponents 2H(tau) and 2H(1/2-tau), balanced at tau=1/4: 2H(1/4) = 1.623 < singlet 1.792 (> sqrt2 = App B circle).
- Corner constant [P, entropy calc]: tau -> 0: eta = 2H(tau), delta = 2 - 2H(1/2-tau) ~ (4/ln2) tau^2
  => delta ~ ln2 * eta^2/ln^2(1/eta) = 0.69 eta^2/ln^2(1/eta)  vs singlet 0.96 eta^2/log2^2(1/eta) = 0.46 eta^2/ln^2(1/eta).
  Ratio 1.5 (numerics at k=400: 1.25-1.8). Still quadratic.
- Clean restatement [P]: via Majorana coordinates y in the parity class {|y| = par mod 2} of F2^{2k+1} (a hyperplane coset),
  beta = F(|y|) with F of period 4 (= +-(-1)^{C(|y|,2)} = Re(c i^{|y|})), and characters of F2^{2k} = characters chi_T of
  F2^{2k+1} modulo T ~ T^c.  So the full Bell model = "probabilistic Fourier sparsity of the symmetric MOD-4 function
  (-1)^{e_2(y)} on n = 2k+1 bits".  Symmetric approximants are provably quadratic; the question is whether NON-symmetric
  sparse approximants of this symmetric target can be exponentially better.
- Rigidity by-product [E]: Fourier-q-sparse => rank <= q, so Majorana gives R_H(16) <= 4944 = 0.30 N^2 at N=128,
  R_H(14) <= 1093 = 0.27 N^2 at N=64, R_H(12) <= 232 = 0.23 N^2 at N=32 (vs naive H_r (x) rank-1 ~ N^2/2(1-sqrt(r/N)):
  5295 at N=128, r=16; line lower bound N^2/(4r) = 256).
- [E] exact8 RESULT N=8: p_min(q) = 28, 28, 28, 24, 24 for q = 1..5 (all O^+(6,2) orbits: 1,2,6,24,110 reps; rational check
  of every improving solution).  q=6: p >= 16 by the span lemma (translate so 0 in S: span <= 5).  Known upper: Majorana (8,8),
  singlet (10,10), (28,1).  By self-duality q_min(p) = 28,28,28,24,24 for p = 1..5... i.e. points (24,4),(24,5) dual.
## 13:40 Non-symmetric perturbations of Majorana [N, heuristic optimiser only]
- maj_perturb.py k=5 d=1 (+12 greedy extra chars), maj_union.py (union of 2-3 random O(Q)-images of the Majorana frame,
  cold and warm-started IRLS): never below the single-frame p=232 at |S|=12 (cold LAD on unions: p ~ 400-510, i.e. LAD
  is a poor L0 optimiser).  No evidence either way beyond "heuristics find nothing better".
## 13:45 Weight-statistic / divisible-code framework [P reformulation + H]
General form of all "symmetric" constructions so far: n linear forms (columns of G) and an affine code c+Gw; if
beta(w) = F(wt(Gw+c)) on the bulk, interpolate F at the weights in a window: q <= C(n, <= #weights-in-window - 1).
 singlet: n=3k (site Paulis); Majorana: n=2k+1 (minimal), weights in one parity class (step 2).
 D-divisible cosets (all weights = j0 mod D, beta = (-1)^{(wt-j0)/D}): step D, so the window in u=(wt-j0)/D has
 Var(u) = n/(4D^2) (projective code => pairwise independent bits).  wt/4 mod 2 IS a quadratic form on any doubly-even code
 (wt(a+b)/4 = wt a/4 + wt b/4 - |a&b|/2), so doubly-even codes of dim 2k, length n in [4k, 8k) whose wt/4 form is IP-type
 would beat Majorana's constant (Var(u)/k = n/64k vs 1/8).  Heuristic: a LINEAR corner would follow from codes of dim 2k,
 length n ~ k*D*polylog, divisibility D ~ log^2 k/eta with (wt-j0)/D mod 2 = IP-type quadratic form [needs the divided
 weight to be QUADRATIC, which for D >= 8 requires |a&b&c| = 0 mod D/2 etc.].  RM(2,m)-subcodes (D ~ sqrt n, few weights)
 checked on paper: rank-deficiency tail 2^{-Theta(j^2)} vs cost n^{2j} => quadratic again.  Few-weight codes are long
 (exactness forces n^{#weights} >= N^2).  Lead, not tested numerically.
- 13:40 CORRECTION (divisible.py first run INVALID): (wt(c+x)-j0)/4 mod 2 is quadratic only if C is self-orthogonal AND
  triply orthogonal (|x&y&z| even for all codewords; the cubic term of the divided weight is |x&y&z| mod 2). The run used
  only self-orthogonality: u-distribution 12/52 at k=3 shows (-1)^u is not quadratic (a nondeg. form has 28 or 36 ones).
  Triple orthogonality adds ~K^2/2 linear constraints per new generator => random constructions need n >~ K^2/2 = 2k^2,
  making Var(u) = n/64 ~ k^2/32 >> Majorana's k/8.  Structured triply-orthogonal codes (RM(1,m)) are exponentially long.
  => The D=4 (and a fortiori D>=8: higher-order orthogonality) weight-statistic route gives NO improvement over Majorana
  [H, constraint counting; not a proof].  Majorana is the optimal member of the D=2 class (minimal length 2k+1).
## 13:50 Majorana identity PROOF [P]
Frame s_1..s_{2k+1} pairwise B=1, sum 0, spanning. y(w) = (B(s_i,w))_i is injective into the even-weight code of F2^{2k+1}.
T(w) := sum_i B(s_i,w) s_i satisfies T(s_j) = sum_{i!=j} s_i = s_j, so T = I and y(w).y(w') = B(w,w').  For even y,y' the polar
form of e_2 is |y||y'| - y.y' = y.y'.  Hence Q(w) + e_2(y(w)) is linear = y(w).c for some c (the s_i span), and
e_2(y+c) = e_2(y) + y.c + e_2(c).  So beta(w) = +-(-1)^{e_2(y(w)+c)} = +-(-1)^{C(|y(w)+c|,2)}: a function of the weight
(mod 4) of a coset of the even-weight code of length 2k+1.  QED.
- [E] greedy k=3: (q,p) = (6,19)@8 (rigorous exact recount): S = {0} u a Majorana frame (5 vectors) of a 4-dim nondegenerate
  subspace, kappa = (1/2)(-1)^{Q(s)}.  => N=8: p_min(6) in [16,19].  As a crossing pair (p,q) = (19,6)@8 this beats the
  project's general-model point (20,6)@8; (46,10)@16 (Majorana) dominates (62,12)@16 [G27-B1].
- [E] explicit operator checks (explicit_check.py, explicit_maj.py): K = sum kappa_s P_s(x)P_s built explicitly, rank(SWAP-K) and
  osr(K) mod 1000003 and 998244353, both Pauli conventions: (q,p) = (6,19)@8, (8,8)@8, (10,46)@16, (67,67)@32 ALL CONFIRMED.
- [L] Savicky 1994 (Europ. J. Combin. 15, 407-410): characterises symmetric bent functions on an even number of variables
  (snippet; stated in later papers as the quadratic ones, i.e. e_2 + affine). So "IP is symmetric in suitable linear
  coordinates" is classical; the 2k+1-bit parity-class (Majorana) form is what halves the effective variance: on 2k bits
  F(J) has period 4 in J (window of d+1 values covers J-range d), on the odd/even class of 2k+1 bits it has period 2 in
  u=(J-par)/2 (J-range 2d per degree d).
- [E within |kappa|<=2] milp_fixedS.py (HiGHS big-M, proven optimal, gap 0): on the FIXED Majorana degree-<=1 support the
  window solution is optimal: N=8 |S|=8 -> p=8; N=16 |S|=10 -> p=46.
- Tensor check [P]: (8,8)@8 (x) (8,8)@8 gives (q,p) = (64, 4096-56^2 = 960)@64 vs Majorana (92,378)@64: splitting frames is
  worse, as the variance argument predicts (one frame of 2k+1 modes minimises Var(u)/k).
## 13:50 Exact N=8 via orbit reps + fixed-support MILP (milp_reps.py) [E]
Rigorous coefficient bound: an optimal phi is fixed by q independent agreement equations (+-1 system), Cramer =>
|kappa| <= q^{q/2}/2^{q-1} (Hadamard bound / det divisible by 2^{q-1}); so the big-M MILP with that bound is EXACT.
q=5 reproduced: p_min = 24 (110 reps, all solves Optimal, 52 s).  q=6 running.
(Full support-free MILP is useless: q<=6 at N=8, 500 s: p=24 found, dual bound 3, gap 0.53.)
- [E, MILP-certified] q=6: 577 O^+(6,2)-orbit reps (fastreps.py, 38 s), every fixed-support MILP solved to optimality with the
  rigorous bound |kappa| <= 6^3/2^5: p_min(6) = 19 at N=8 (upper bound rationally verified + explicit operator rank check;
  lower bound relies on HiGHS optimality certificates, float tolerances).  Optimiser = 0 + Majorana frame of a 4-dim subspace.
  N=8 Bell frontier so far: q=1..3: 28 | 4,5: 24 | 6: 19 | 7: <= 19 | 8: <= 8 (Majorana) | 10: <= 8 | 28: 1 (dual: (p,q) symm.)
- q=7: 3366 reps (reps7.pkl, 161 s); cutoff MILP (p<=18?) costs ~1.4 s/rep => ~80 min: NOT run (budget). First 22 reps: p >= 19.
## FINAL SUMMARY (14:16 UTC)
LITERATURE [L, snippets only]: no upper bound on R_H(N^eta) better than the polynomial method is known to me; Servedio-Viola
 2012 is exactly the sparse-polynomial special case (line Omega(1/s) for sign agreement; beyond = open, tied to AC0oMOD2 for
 IP); approximate sign-rank of H open between Omega(log N) and log^{O(sqrt(log N) log 1/eps)} (Bindu et al. 2026);
 R_H(1) = N^2/2 - N^{1.5}/2; Alman-Liang / Alman: rigidity lower bounds only up to rank log^2 N over finite fields.
NEW: Majorana family (IP_k = period-4 symmetric function of 2k+1 linear forms; proof above): exact points
 (8,8)@8, (10,46)@16, (12,232),(67,67)@32, (92,378)@64, (576,576)@128, (834,3214)@256; explicit operator checks up to N=32.
 Best known balanced crossing pairs for N=2^k, 5<=k<=31 (beats HR tiling: 576 vs 1024 @128, 5036 vs 14848 @512).
 Corner: delta ~ ln2 eta^2/ln^2(1/eta) (1.5x singlet) -- QUADRATIC; S_{2k+1}-symmetric class provably quadratic [P sketch].
 Weight-statistic/divisible-code generalisation needs triply-orthogonal codes => long => no gain [H].
EXACT N=8: p_min(q) = 28,28,28,24,24,19 for q=1..6 (q<=5 by exhaustive incidence enumeration; q=6 by MILP per orbit rep).
VERDICT: no linear corner found; every mechanism found (singlet, block, Majorana, divisible codes) is Gaussian-tailed.
