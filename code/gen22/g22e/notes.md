# G22e — is sqrt2 optimal for vertex-algebra-type constructions?  (working notes, saved as I go)
Tags: [P] proved here (short, unreviewed), [E] exact computation, [N] numerical, [H] heuristic, [L] literature quoted from memory, [O] open.
Nothing here is a new algorithm or a new bound for 3-SAT.

## 1. Exact reduction  [P]
Data (A1)-(A4) as in G21a.  D_M(E) = dim M_{<=E}, L_M = log D_M; same for V.  Pieces A_x (x = row weight in [0,B], B = 2E0+H),
rr(A_x) <= c D_M(x), rc(A_x) <= c D_M(B-x), c <= (H+1)T'.  CERTIFIED cost K = max_{x in [0,B]} [L_M(x) + L_M(B-x)], 4 rho lambda = K.
(The theorem is about this certificate; true ranks could be smaller: escape route R5.)
 Since L_M is nondecreasing, K >= 2 L_M(B/2)  (the midpoint piece; App. B has every s, so it is present).

Three upper bounds on lambda = log N:
 (U) UPPER ANCHOR (always, no labels, no locality):  X viewed as an (u,v)|(x,y) matrix is a permutation matrix of size N^2, so
     N^2 = rc(X) <= sum_x rc(A_x) <= c(B+1) D_M(B):   2 lambda <= L_M(B) + l,  l = log(c(B+1)).
 (Lab) labels: fix x = v; the (u,y) matrix of <u|o(.)|y> is d*I_N, so the N labels are independent in M_{E0}: lambda <= L_M(E0).
 (Loc) LOWER ANCHOR in V.  CORRECTION to G21a (V'): "N <= dim V_{<=H/2}" is NOT a consequence of (A3) in general.
     What is true [P]: suppose every c_xv = a_x(H-1)b_v lies in V_{<=K} (K = 0: vacuum normalisation, App. B) and phi is a
     functional on V_{<=K} (extended by 0) with rank [phi(c_xv)] = N.  Write a = sum_m a^m, b = sum_n b^n (weight components);
     a^m(H-1)b^n has weight m+n-H, so only pairs with m+n <= S := H+K contribute.  Split them into m <= S/2 and m > S/2
     (then n < S/2): phi(c) = phi(c(P_{<=S/2} a, b)) + phi(c(P_{>S/2} a, P_{<S/2} b)), two bilinear forms of ranks
     <= dim V_{<=S/2} and <= dim V_{<S/2}.   So  N <= 2 D_V(S/2).   Locality is NOT used; only the grading.
     For module labels with non-vacuum c (K > 0) the bound degrades to S = H+K, and for K large it is void:
     e.g. H = 0, a, b in commuting tensor factors V1 (x) V2, labels in the top space W1 (x) W2 of an irreducible module, where
     o(a(-1)b) = o(a)(x)o(b) and o(V_i) = End(W_i) (Zhu algebra onto End of the top level [L]); see sect. 5 (N = 2 > D_V(0) = 1).
Reduction.  4 rho >= 2 L_M(E0 + H/2) / min( L_M(E0), log2 + L_V((H+K)/2), (L_M(B)+l)/2 ).
 For M = V and K = 0, with w := (E0 + H/2)/2 = B/4:  min(L(E0), L(H/2)) <= L(min(E0,H/2)) <= L(w), so
     4 rho >= inf_w 2 L(2w) / min( L(w) + log2, (L(4w) + l)/2 ).                                       (*)
 The same one-variable form holds for M != V under DOMINANCE  D_V(E) <= C0 D_M(E) (replace log2 by log(2 C0)).
 Without dominance the V-anchor can be at a larger argument than the M-anchor (V grows faster than M): this is the only
 place where c_eff(M) < c_eff(V)-type phenomena enter (sect. 5, R2).
 REFINEMENT (sect. 5, regrade.out): the lemma holds in EVERY admissible grading (one where a(j) maps V_m x V_n -> V_{m+n-j-1});
 H is grading-free, so one may take the anchor in one grading (e.g. conformal L0) and the cost in another.  The bound in the
 cost grading alone can be far too weak (regraded lattice, sect. 5): the honest anchor is the RANK of the pairing.

## 2. Barrier theorems  [P unless marked]
Thm A (upper anchor; asymptotic; no unitarity, no labels, no locality bound, any K).  Let (V,M) be fixed and suppose
 L_M(E) = A sqrt(E) (1+o(1)) for some A > 0 (more generally liminf L_M(E/2)/L_M(E) >= 2^{-alpha}).  Then for every
 VA-type decomposition (A1)-(A4) of the crossing of N = 2^lambda strings, certified 4 rho >= 2 sqrt2 - o(1) as lambda -> oo
 (resp. 4 rho >= 2^{2-alpha} - o(1)).
 Proof.  (U): L_M(B) >= 2 lambda - log(c(B+1)), c <= T'(H+1) <= T'(B+1).  If log B >= lambda/ log-scale then L_M(B) >> lambda and
 K >= 2L_M(B/2) >> lambda anyway; otherwise l = O(log B) = o(lambda) and B -> oo.  K >= 2 L_M(B/2) = (sqrt2 - o(1)) L_M(B)
 >= (sqrt2 - o(1)) 2 lambda.  QED.   [Elementary; holds for ANY module, twisted or not, any (half-)integral grading.]
 Remark: the paper's "growth E^alpha gives rho = 2^{alpha-1}" (alpha >= 1/2, lower anchor) and Thm A (alpha <= 1/2, upper
 anchor) together give rho >= max(2^{-alpha}, 2^{alpha-1}) >= 1/sqrt2 for power-law profiles, equality only at alpha = 1/2.
 So for a FIXED algebra the circle is forced by the single fact "log dim grows like sqrt": the c/(24u) corrections of G21a
 (which could have either sign for c < 0) are O(1/lambda) and irrelevant for the exponent alpha'.

Thm B (exact, finite size, M = V or dominance).  K = 0 (vacuum normalisation, as App. B).  Then (*) of sect. 1 holds.  Put
 g(E) = L(E)/sqrt(E).  If g is monotone (either direction) on the relevant range, 4 rho >= 2 sqrt2 - O(l/lambda):
 g nondecreasing => L(2w) >= sqrt2 L(w) (lower anchor wins);  g nonincreasing => L(4w) <= sqrt2 L(2w) (upper anchor wins).
 Failure at scale w needs  g(2w) < g(w)  AND  g(4w) > g(2w): a DIP (local minimum) of L/sqrt E.
Local-exponent lemma (G21a's (L), re-proved) [P].  For the Legendre entropy S(E) = inf_b [phi(b) + bE], phi = log Z:
 d log S / d log E = 1/2 - psi/(2S) with psi(b) := d/db[b phi(b)] = log Z - b<E>_b at the optimal b; so g_S is increasing
 where psi < 0 and decreasing where psi > 0.  (Chernoff gives L <= S exactly; L = S - O(log) needs a local CLT, so exact
 finite-size statements are checked on true dimensions, not on S.)  For tensor-power (scaling) families V^{(x)n}, n -> oo,
 L_n(eps n)/n -> s(eps) = S per copy EXACTLY (Cramér/Gärtner-Ellis for the cumulative count), and the criterion is
 r(eps) = 2 s(2eps)/min(s(eps), s(4eps)/2) >= 2 sqrt2  for all eps.  A dip needs psi > 0 at some b1 and psi < 0 at b2 < b1.
 Endpoints [P/L]:  b -> oo: psi -> log dim M_0 (= 0^- for CFT type: psi ~ g1 e^{-b D}(1 - bD)).
   b -> 0 (one-term modular regime): log Z_M = pi^2 c_eff/(6b) + sigma_M + kappa_M b, sigma_M = log S_{M,min},
   kappa_M = h_M - c/24 (+ (1/2)log(b/2pi)-type terms for eta-quotients), so psi = sigma_M + 2 kappa_M b.
   Vacuum of a unitary rational VOA: sigma = log S_00 < 0, kappa = -c/24 < 0: psi < 0 at both ends.
   Module with sigma_M < 0 < kappa_M (e.g. sl2_k L(j), j > 0): psi changes sign at b* = |sigma|/2kappa -> a dip.
 Checks of these asymptotic claims against exact dimensions [E/N, psiscan.out]: V#: psi = -2.000, -4.000, -6.000 at b = 1,2,3
 (= -cb/12, sigma = 0); E8: -0.333, -0.667 at b = .5, 1 (= -8b/12); NS fermion: -b/24 (sigma = 0); Ramond x sqrt2: +b/12;
 Ising vacuum: log(1/2) - b/24; sl2_1 L(1/2): exact sign change at b = 0.83 vs predicted |log(1/sqrt2)|/(2 * 5/24) = 0.832.
 Cardy growth L/(A sqrt E) -> 1 (slowly, log corrections) for all 12 profiles in crit.out (0.90-0.98 at E = 400; M(2,7) 0.62).

Where the hypotheses are theorems [L, references from memory]:
 * Rational, C2-cofinite VOAs of CFT type, irreducible (ordinary or finite-order twisted) modules: Zhu, JAMS 9 (1996)
   (characters span an SL2(Z)-invariant space), Dong-Li-Mason, CMP 214 (2000) (twisted case, Gamma(T)-type invariance).
   Then chi_M(e^{-b}) ~ S_{M,min} e^{pi^2 c_eff/(6b)}, c_eff = c - 24 h_min, provided S_{M,min} != 0 (true in all unitary cases,
   where min = vacuum and S_{M0} > 0; for non-unitary minimal models checked case by case, e.g. Yang-Lee).
   Tauberian step: only LOG asymptotics of the CUMULATIVE count are needed:  log Z(b) ~ C/b  <=>  log D(E) ~ 2 sqrt(CE)
   (de Bruijn's exponential Tauberian theorem / Kasahara; Bingham-Goldie-Teugels, Regular Variation, sect. 4.12).
   It needs only that D is nondecreasing, which is automatic (D is a cumulative count of nonnegative dims).  Asymptotics of
   dim M_E itself (Ingham 1941) would need monotone coefficients — NOT needed here, and false in general (Ising vacuum has
   dim V_1 = 0).  Half-integer gradings and Z2-graded superalgebras: use the character, never the supercharacter.
 * C2-cofinite non-rational (triplet W(p)): Miyamoto, Duke 122 (2004) — pseudo-traces, log-modular; log Z ~ C/b still [L].
 * Free fields (Heisenberg, fermions, symplectic fermions, lattice): Meinardus' theorem, or directly eta-quotients x theta.
 * Families with c -> oo: Thm A says nothing; Thm B needs the dip test, done numerically per family (sect. 3).
 Range: the asymptotic regime is E >> c; tensor-power families live at E ~ eps c, covered by s(eps).

## 3. Exact checks  [E dims, N Legendre]
Fixed algebras, M = V, exact cumulative dims to E = 400 (crit.out): min_w 2L(2w)/min(L(w), L(4w)/2) over all w (only exception: tiny-lambda NS n=1, below):
 Heisenberg 2.967, Heis^8 2.983, NS^8 2.870, E8 2.859, Leech 2.856, V# 2.856, Yang-Lee 2.980 (vac) / 2.947 (h=-1/5),
 M(2,7) 3.009, symplectic fermions 2.927, Ramond^8 2.917: all >= 2 sqrt2 = 2.8284, decreasing to it from above.
 (Tiny-lambda artefact: NS n=1 gives 2.0 at w = 1/2, lambda = log 2; O(log) terms dominate there.)
Scaling families (tensor powers), per-copy criterion (scaling.out, psiscan.out): no psi sign change and min r >= 2.83 for
 Heisenberg, NS/Ramond fermion, E8, Leech, V#, Yang-Lee (vac, module), M(2,7), M(2,9), symplectic fermions, twisted Heisenberg,
 Ising (vac, eps, sigma).  Minima sit at the largest eps computed: approach to 2 sqrt2 from above.
Modules, two-variable exact reduction (modscale.out): min_{u,v} 2 s_M(u+v)/min(s_M(u), s_V(v), s_M(2u+2v)/2):
 sl2 levels 1,2,4,8, all L(j): 2.88-3.09;  Ising sigma/eps 3.09/3.04;  Yang-Lee phi: 3.12;  NS with Ramond labels: 2.841.
FIRST DIP FOUND (dip.out): the profiles of sl2_k L(j), j>0, HAVE a dip; used as their own anchor they give
 min r = 2.816, 2.812, 2.811, 2.812 (k = 1, 2, 4, 8; at eps ~ 1-6 per copy) < 2 sqrt2.  In the true pair (V = L(0), M = L(j))
 the vacuum anchor rescues it (s_V < s_M: 0.44 vs 1.04 at eps = 0.1, 3.24 vs 3.47 at eps = 2): 2.90.

## 4. Abstract structure (task 3f)  [P]
 With only (A2)-(A3) on the labels (identity on labels, no growth input), the structure is equivalent to an arbitrary
 decomposition X = R + C: given rr(R) = p, rc(C) = q, put M = labels + W_row + W_col, a_x, b_v the factor maps; then
 a b + b a = id on labels is App. H (75) "directions are orders" (H = 1).  So locality + non-degeneracy alone imply only the
 line (dim M_{<=B} >= N^2/(B+1)); there is NO combinatorial barrier beyond the line.  All the force of the circle is in the
 GROWTH input (sqrt law from modularity) plus the anchor structure.  For H = 1 with all modes (Clifford relations for every
 mode) the cyclic module is a fermion Fock space: NS labels give rho ~ 1, Ramond ground-space labels give rho = 3/4 (G21a).

## 5. Escape routes (ranked by plausibility)
 R1  Anchor-grading mismatch (most concrete).  The cost is computed in one admissible grading g1, the pair count is bounded in
     another g0 (H is grading-free).  If g1 makes the profile dip while the true pair count stays large, (*) fails.
     Computation (regrade.out): A1 lattice VOA regraded by L0 - alpha(0)/2 (profile = sl2_1 L(1/2), V'_0 = 2-dim):
       loose anchor (N <= dim V'_{<=h}, the lemma in g1)  -> 2.8157 < 2 sqrt2  (ON PAPER);
       anchor in conformal grading                         -> 2.916;
       true pairing rank (sector n pairs -n; = profile of sum q^{n^2+|n|}/eta)  -> 3.030.
     So the paper-gain is an artefact of a loose label bound; but it shows exactly what is needed: a VA whose pairing rank at
     locality order H is ~ dim of a DIPPING profile.  First computation to do: search spectral-flow / charge regradings of
     lattice and affine VOAs (and orbifolds) for which pairs of sectors (q, -q) have EQUAL dimensions in the cost grading
     (pairing rank = full dim) — e.g. self-conjugate sectors, or b-c systems with symmetric charge.
 R2  Module with a dip and a vacuum that grows at least as fast near the dip scale ("anti-dominance").  For unitary rational
     sigma_M >= sigma_V (quantum dimensions >= 1), so dominance holds and this is excluded asymptotically; non-unitary models
     and twisted modules of orbifolds are the place to look (need sigma_M < 0 < kappa_M and sigma_V >= sigma_M).  The gain
     available is small: dips seen so far give 2.81 vs 2.83 (rho 0.703 vs 0.707), but in a tensor-power family it would be
     UNIFORM in n, i.e. a genuine improvement of alpha'.  First computation: scan (sigma, kappa) over minimal models M(p,q)
     and lattice orbifold twisted sectors with the one-term formula, then the exact two-variable test.
 R3  c -> oo families with non-tensor structure (W-algebras, cosets at large rank, holographic/extremal-type characters):
     only Thm B applies; needs a vacuum character with psi > 0 at intermediate b.  V# and Leech checked: psi < 0 everywhere.
 R4  Midpoint gap (no pieces near B/2): pieces at x <= B/4 and >= 3B/4 would give 4 rho -> 1 + sqrt3 = 2.73, but that is a
     cheap balanced point (rr, rc ~ N), as hard as the original problem (cf. G20b).
 R5  True ranks below window dimensions (certificate not tight); G20b found equality for App. B.
 R6  Non-sqrt growth: W_{1+oo} (plane partitions, E^{2/3}) is worse (2^{-1/3} lower anchor 0.79); infinite-dim weight
     spaces (beta-gamma): window dims infinite, the certificate is void; would need a different finite filtration [O].
 R7  Fractional locality (intertwiners, parafermions): excluded from (A2) itself (G21a).
Correction to G21a (V'): "N <= dim V_{<=H/2}" is false without vacuum normalisation: h0counter.out gives N = 2 (top dims 3x3)
 and N = 3 (4x4) solutions with H = 0, K = wt(a)+wt(b) [N exact to 1e-15, realised by Zhu-algebra surjectivity [L]];
 no solution with 2x2 tops.  Harmless: labels in a top level cost rho >= 3/4 (G21a (G)).

## 6. Verdict
 PROVED (elementary, given the sqrt law): for every FIXED vertex (super)algebra and module / finite-order twisted module with
 log dim M_{<=E} ~ A sqrt E — all rational C2-cofinite VOAs of CFT type and their modules (Zhu/DLM + de Bruijn), free-field and
 lattice algebras — every construction of type (A1)-(A4) has certified rho >= 1/sqrt2 - o(1).  No unitarity, no c_eff relation
 and no label bound are needed (only the rank of X).  So nothing at fixed c can beat sqrt2.
 For FAMILIES with c ~ E (the only place left): barrier holds iff no dip of L/sqrt E is exploitable; all 20+ tested algebras
 and modules pass; module dips exist (sl2 L(j), 2.81) but are neutralised by the vacuum anchor; the regraded lattice shows the
 loose dimension anchor is not enough and the pairing rank must be used.
 Adversarial notes: Thm B checks are finite (E <= 400, eps up to ~30) + analytic endpoints, not a proof for all eps; the
 Legendre/large-deviation identification is exact only for n -> oo tensor powers; Thm A's o(1) hides O(log B / lambda);
 the App-B-type construction (N ~ dim V_h dual pairs with H = 2h) for a general VOA is assumed, not built.
Files: chars.py crit.py scaling.py psiscan.py modscale.py dip.py regrade.py h0counter*.py and *.out (this folder).
