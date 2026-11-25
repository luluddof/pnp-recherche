# G28-C3 notes (start 10:11 UTC). Routes H (Hagedorn) and P (projector chains) to a linear corner.
Tags: [P] proof sketch, [E] exact computation, [N] numeric, [C] counted only, [H] heuristic.
## 0. Plan (10:30)
- P1: "site-pattern model" (K = function of the per-site Sym/Anti pattern; = local-U(n)-invariant K on (F^n)^{(x)k}):
  random restriction mu_b = mixture of uniform subcubes + Smolensky over R  => delta <= C_n eta^2 (generalises Thm S to
  NON-symmetric functions). Verify formula for p,q vs explicit ranks; numeric search for non-symmetric g.
- P2: G-invariant model for groups with edge-type (geometric) dimension profiles (iterated wreath products / ultrametric).
- H1: symmetric orbifold Sym^M(Heis) dims in C1's GTC counting; H2 exp. many species / matrix (large-N) Fock spaces.
## 1. STOPPED at ~10:25 UTC (coordinator pause). No computations run yet. State of the analysis:
### 1a. Theorem S' (site-pattern / local-U(n)-invariant model is quadratic) [P, unreviewed, not yet numerically checked]
Setting: N = n^k, K in (x)_i span{1, SWAP_i} (equivalently K commutes with (U_1(x)..(x)U_k)^{(x)2} for all U_i in U(n)).
Then K = g(T), T in {0,1}^k = pattern of antisymmetric sites, trace measure mu_b, b = (n-1)/(2n); SWAP = (-1)^{|T|}.
Writing g = sum_U c_U prod_{i in U} tau_i (tau_i = SWAP_i - 1/n = centred variable), osr(K) = sum_{U: c_U != 0} (n^2-1)^{|U|}
(per-site left-factor space M_n = span{1} + traceless). Hence deg g <= log_{n^2-1} q.  p = mu_b(g != (-1)^{|T|}) N^2.
Proof of the bound: mu_b = mixture: each coordinate independently fixed to 0 w.p. 1-2b, else uniform on {0,1}.
On each random subcube (m free coords, uniform measure) g restricted has degree <= d and the target is the parity of the
free coords; Smolensky over R (y_S = parity * y_{S^c} on the agreement set) gives agreement <= 2^{-m} sum_{i<=(m+d)/2} C(m,i).
So p/N^2 >= E_{m~Bin(k,2b)} P(Bin(m,1/2) > (m+d)/2) ~ exp(-d^2/(2m)).  With d = eta k ln n/ln(n^2-1):
   delta <= eta^2 * n ln n / (2 (n-1) (ln(n^2-1))^2) + o(1)     (n=2: delta <= 0.574 eta^2; n large: -> eta^2/(8 ln n)).
Covers ALL functions of the singlet pattern (no permutation symmetry needed): generalises Thm S of G27-T3 (which was S_k-
symmetric Bell) for this submodel.  Does NOT cover the full Bell model (there a character costs 1, not 3^{|U|}; restriction
argument then only gives the row-uncertainty line bound B3).  Also extends (same proof) to O(n)-local-invariant K
(site algebra span{1,SWAP,Phi}) -- to be checked.
### 1b. Reflection-product model (Route P abstraction) [P/H]
SWAP = prod_i (1-2E_i), E_i commuting projectors; K = degree-d truncation in E's: p = mass(M>d), M = sum E_i,
q <= sum_{|S|<=d} osr(E_S).  Parity constraint P(M odd) ~ 1/2 is compatible with edge distributions (a0~a1~1/2, then
geometric).  T3's spectral test (|E[beta phi]| <= ||phi^||_1/N) excludes super-geometric tails but allows geometric tails
with rho >= 1/3 (Taylor truncation).  In Bell, Fourier-sparse indicators have granularity => osr(E_i) >= 1/eps_i; uniform
small eps design forced L ~ 2k/log k and is killed by the spectral test (coset indicators have ||^||_1 = 1).
### 1c. Route H [C, from reasoning only]
Sym^M(Heis) untwisted VOA: multisets of <= M Heisenberg states; for M >~ w growth log d(w) ~ c^2 w/(4 log w) (near-Hagedorn),
Cardy (c_eff = M) beyond; inside a VOA locality identity C1's span-doubling => endpoint excess => jump (expected; NOT computed).
Exp. many species / matrix large-N Fock spaces = Heisenberg with g(n) ~ e^{lambda n} generators: theta_end = 4 (C1 sect. 4).
### 1d. G-invariant model idea (not started): K invariant under diagonal conj. by G (e.g. iterated wreath S_n wr ... wr S_n,
whose orbitals on X x X have geometric = edge-type masses); p and q decompose over irreps -> small multiplicity matrices;
search block-vanishing patterns at N = 8, 9, 16.  Partition-algebra (G = S_N) case: trivial (p ~ N^2/2 or q ~ N^2 - O(N)) [P].
## 2. NEXT STEPS (precise)
 1. Script s_prime_check.py: for n=2, k<=5 build K = g(T) explicitly (L = sum SWAP_i), compare rank(SWAP-K), osr(K) (mod 2
    primes) with pattern formulas; then for k<=12 LP/L1 search of non-symmetric degree-d g under mu_{1/4}; check bound.
 2. Compute exact bound curve delta_bound(eta) at finite k vs singlet family (T3 bell_poly.py).
 3. Sym^M(Heis) dims -> C1 GTC curve (reuse gen28/c1/appb_curve.py logic).
 4. G-invariant wreath model at N=8,9,16.
Processes: none running.
