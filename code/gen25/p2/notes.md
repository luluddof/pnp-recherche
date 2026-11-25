# G25-P2: an explicit, re-optimised Theorem A (Claude Oct-2026 preprint, App. A)

Tags: [P] proof written here (unreviewed by a second agent), [N] own numerics (scripts in this folder),
[R] taken from the paper as written (untrusted; correctness of the paper's Lemma 8 / Prop. 3 / Thm B NOT re-checked).
Nothing here changes an exponent of log m. All improvements are in constants and in the threshold m*.

Files (all in this folder):
- `sweep_check.py` / `sweep_check.out`: sweep lemma vs brute force (random tables + end-to-end #SAT of random 3-CNFs).
- `thmA_explicit.py` / `thmA_table.out`: exact recursion, three top levels, all lower-order terms, thresholds.
- `asym_c1.py` / `asym_c1.out`: asymptotic constant c1. `schedule_c1.out`: level table of the c1 schedule.
- `sensitivity.out`: thresholds vs the two machine constants left free (kappa, c_ev).
- `claude_cross2.py`: copy of b4's model (rows A/B = paper as written / chain wire count, main term only).

Notation. log = log2. A 3-CNF with m clauses, n variables all occurring, W = #equality tests <= 3m - n
(the paper uses W <= km = 3m). r = m/n, so W/n <= 3r - 1. A (K,lam)-clique is the paper's Def. 3.
"Exponent" = log2(#operations).

---------------------------------------------------------------------------------------------------------
## 1. Sweep lemma (idea I2)  [P + N]

**Lemma S.** Let K >= 2, lam, B >= 1, s(K) = floor((K+1)^2/4) - 1. The value of a (K,lam)-clique with
|W_w| < 2^B is computed with at most K * 2^{lam s(K)} table queries and as many multiply-adds on integers of
fewer than lam K^2/2 + K B + 1 bits. Hence level 0 of Theorem A has e0 = s(K0)/K0^2 <= 1/4 + 1/(2K0)
instead of the paper's (K0-1)/(2K0) ~ 1/2 (time; memory 2^{lam floor(K^2/4)}, which Theorem A does not count).

*Proof.* Nodes 1..K. For 0 <= j <= K let E_j = edges inside [j], C_j = {{u,v}: u <= j < v} (the cut), and for
sigma in ({0,1}^lam)^{C_j}
   F_j(sigma) = sum over tau in ({0,1}^lam)^{E_j} of prod_{w <= j} W_w((sigma, tau) restricted to the edges at w).
(Every edge at w <= j lies in E_j or C_j, so this is well defined.) F_0 = 1 and F_K = Z, the clique value
(C_K is empty, E_K is all edges). Let D_j = {{u, j+1}: u <= j} (subset of C_j) and N_j = {{j+1, v}: v > j+1}, so
E_{j+1} = E_j + D_j and C_{j+1} = (C_j - D_j) + N_j. Splitting the product over w <= j+1 into the product over
w <= j and the factor W_{j+1}, whose edges are D_j + N_j, and the sum over E_{j+1} into sums over E_j and D_j:
   F_{j+1}(sigma') = sum_{pi in ({0,1}^lam)^{D_j}} F_j(sigma'|_{C_j - D_j}, pi) * W_{j+1}(pi, sigma'|_{N_j}).
All values of F_{j+1} are obtained by one pass over ({0,1}^lam)^{C_j + N_j} (|C_j| + |N_j| = j(K-j) + K-j-1
strings), one query and one multiply-add per term. f(j) = j(K-j) + K - j - 1 is a concave quadratic with real
maximum (K+1)^2/4 - 1 at j = (K-1)/2, so max_{0<=j<K} f(j) = s(K). Sizes: |F_j| <= 2^{lam|E_j| + jB}. QED

*Check* (`sweep_check.out`, ALL OK): (a) 30 random integer-table cliques, K in {3,4,5}, lam in {1,2}: sweep value =
brute force value, and the largest step enumerates exactly s(K) edge strings (3, 5, 8 for K = 3, 4, 5 vs
K(K-1)/2 = 3, 6, 10). (b) 24 random 3-CNFs (n <= 9, m <= 9), clauses -> K in {3,4} groups, (K,L)-clique tables
built exactly as in the paper's Lemma 1 / eq. (26) from the network of Section 6, sweep over the clique:
equals brute-force #SAT in every case. This is the whole algorithm at l = 0, end to end.

---------------------------------------------------------------------------------------------------------
## 2. Exact recursion (idea I3)  [P, modulo the paper's Prop. 3 / Lemma 8 constructions, R]

Three corrections to the paper's bookkeeping in Thm A(b), each read off its own proof:
(i) **Exact cost terms.** With U <= N/g + |R| S(lam) (|R| = C(g,2)^2 <= g^4/4, straight edges <= N/g) and
    K^2 L <= U + (zeta1 + zeta2) N + K^2, the main term is e_l (rho_bar + 1/g + zeta1 + zeta2) N, not
    e_l (rho + 3 zeta) N with 1/g <= zeta folded in.
(ii) **Entry cost 2^{N/K}, not 2^{2N/K}.** In Lemma 8 case (ii) every group has exactly g/K rows, a row setting
    has g(g^2-1)lam < g^3 lam bits, so a table entry of the (K,L)-clique is a sum of < 2^{(g/K) g^3 lam} = 2^{N/K}
    terms (the paper's 2^{2gP/K} is the case-(i) bound, groups <= 2g/K). The paper's third 1/K pays for c*lam;
    here that term is kept separately as c_ev/g^4 (negligible).
(iii) **Split Hoeffding.** Lemma 8 (ii)'s proof uses K group events (range P) and gK node-group events (range
    P/g), each bounded by exp(-2(g/K) zeta^2). Using zeta1 for the first kind, zeta2 for the second, the
    grouping exists when K e^{-2 t zeta1^2} + t K^2 e^{-2 t zeta2^2} < 1 (t = g/K) and then
    Y < (U + (zeta1 + zeta2) gP)/K^2 (same two-line argument as (39)). We take
    zeta1 = sqrt(ln(2.0001 K)/(2t)), zeta2 = sqrt(ln(2.0001 t K^2)/(2t)).

**Recursion (Theorem A'(b)).** Parameters: x0, integers t_l >= 2; K0 = 2^{x0}, g_l = t_l K_l, K_{l+1} = g_l^2,
   e0 = s(K0)/K0^2 (sweep) [or (K0-1)/(2K0), brute force],
   e_{l+1} = e_l (rho_bar + 1/g_l + zeta1_l + zeta2_l) + 1/K_l + c_ev/g_l^4.
Claim: the (K_l, lam)-clique value is computed within 2^{e_l K_l^2 lam + O_l(lam)} operations, where O_l is the
explicit sum (coded in `overhead_per_n`) of: |R| log2|I(lam_j)| (the iota's), g log K (exhaustive search for the
grouping), e_{j-1} K_{j-1}^2 (ceiling of L), c_ev lam_j (evaluating pieces), kappa log2 Xi_j (polynomial cost per
step), over the levels j <= l, with lam_{j-1} = ceil((rho_bar + 1/g + zeta1 + zeta2) g^4 lam_j / K_{j-1}^2) + 1
(L is padded up to this value, which the Lemma-1 construction allows; so lam only grows going down).
*Proof sketch:* the paper's proof of (b), with the paper's multiplicative bound (lam+B)^{C0 K^2} replaced by
additive bookkeeping: per level, #ops = #iota * [grouping search + T_{l}(L) * (entry cost) + poly], so the logs
add. The paper's C0 then appears only as kappa (one per level, additive), not multiplied by K^2.
rho_bar(lam) = S(lam)/(4 lam), S(lam) = max_s (c_s^row + c_s^col) of App. B (53)-(56):
c_s^row = ceil log sum_{i<=2h} p(s+i) <= ceil log((2h+1) p(2h+s)), etc. [N: computed with the Hardy-Ramanujan
leading term + 1 bit; not a certified bound]. rho_bar/(1/sqrt2) = 1.017 (lam = 2^10), 1.0052 (2^12),
1.00004 (2^20). At every lam the algorithm meets (>= 2^{42} in all rows below) the excess is < 1e-9.

---------------------------------------------------------------------------------------------------------
## 3. Top level (ideas I4a, I4)  [P]

The paper (case (i)): g = 3K^8 clusters, Chebyshev, K <= m^{1/9} ("a = 9").

**Lemma T4 (refined case (i), deterministic).** Network on [g] with c_ab <= P, total U <= U0. Independent
uniform chi: [g] -> [K], Y = Y_{gamma,gamma'}. Covariances: same ordered pair <= 1/K^2; pairs sharing their
first (or second) node <= 1/K^3; reversed pair or pairs sharing a node in opposite roles <= 0; disjoint 0.
So Var Y <= (1/K^2) sum c_ab^2 + (2/K^3) sum_a c_a^2 <= P U (1 + 2/K)/K^2. Chebyshev with deviation zeta' U0/K^2
and a union over K^2/2 pairs: failure <= P(1+2/K)K^4/(2 zeta'^2 U0). Group sizes: Chernoff, union over K,
<= 1/4 with delta' = sqrt(3K ln(4K)/g). If the pair failure is <= 1/4, a grouping with all Y <= (U + zeta'U0)/K^2
and groups <= (1+delta')g/K exists and exhaustive search finds it (cost g log K, counted). For the CNF network
(P = 6 ceil(m/g), U0 = 2W) this needs g ~ 2.1 K^4/zeta'^2 (r = 4.26): **a = 4** instead of 9. Cost bits:
e (2W(1+zeta') + K^2) + 3(1+delta')(m+g)/K + g log K + lower order.

**Lemma 8' (random clustering, Las Vegas).** Give each clause an independent uniform cluster label in [g],
g = tK. Test multigraph on clauses (tests inside one clause ignored): max degree <= 6, so it has a proper
edge colouring with <= 9 matchings (Shannon). Within a matching the crossing indicators for (a,b) are independent
Bernoulli(2/g^2); with mu_bar = 2W/g^2 and Bernstein per colour class (t_i = d(mu_i + mu_bar/9)):
   Pr(c_ab > (1+2d) mu_bar) <= 9 exp(-d^2 mu_bar / (18(1+d/3))).
Cluster sizes: Pr(|cl_a| > (1+dc) m/g) <= exp(-dc^2 (m/g)/(2(1+dc/3))).
**Explicit failure probability:** p_fail <= (9/2) g^2 exp(-d^2 2W/(18 g^2 (1+d/3))) + g exp(-dc^2 m/(2g(1+dc/3))).
The algorithm checks a sample in poly(m) time and resamples; with p_fail <= 1/2 the expected number of samples is
<= 2. On success put P = max((1+2d) 2W/g, 3(1+dc)m/g): then c_ab <= P/g, Lemma 8 (ii) applies (deterministic
search over equal groupings, cost g log K), L = ceil(2W(1 + (zeta1+zeta2)(1+2d))/K^2), entries 2^{3(1+dc)m/K}.
Needs m >~ g^2 log g / d^2: **a = 2**. Result: zero-error randomized (Las Vegas) algorithm.
(Correctness never depends on the sample; only the running-time bound does.)

---------------------------------------------------------------------------------------------------------
## 4. Theorem A' (explicit)

**Theorem A'.** Fix x0 and t_l as in Section 2 (the algorithm may pick, from m, the best of a finite family of
schedules). For a 3-CNF with m clauses and W tests, #SAT (hence 3-SAT) is computed
 (det.)        in 2^{2 e_l (1+zeta') W + 3(1+delta')(m+g)/K_l + g log K_l + O_l} operations  (Lemma T4), and
 (Las Vegas)   in expected 2^{2 e_l (1 + (zeta1+zeta2)(1+2d)) W + 3(1+dc) m/K_l + g log K_l + O_l + 1} operations
               (Lemma 8'), whenever the parameter conditions above hold.
Asymptotically (`asym_c1.out`, rho_bar -> 1/sqrt2, top deviation -> 0 inside each level's range of m):
   time <= 2^{(c1 + o(1)) W / sqrt(log2 m)}, c1 = limsup_l 2 e_l sqrt(a x_{l+1}):
   paper as written (a = 9, W read as 3m)     42.4
   C'  exact recursion, brute base, a = 9     30.9
   D   + sweep base, a = 9                    15.9
   E1  + Lemma T4 top, a = 4 (deterministic)  10.6
   E2  + Lemma 8' top, a = 2 (Las Vegas)       7.47
(schedule x0 = 9, zeta2_l ~ 2^{-6-l}; c1 oscillates, its liminf is c1/sqrt2: e drops by rho only once per doubling
of log K.) Level table of that schedule (`schedule_c1.out`): log K_l = 9, 49.2, 137, 320, 691, 1439, 2943, 5955,
11986; e_l = .2510, .1855, .1335, .0953, .0677, .0480, .0340, .0240, .0170.

**Thresholds** (`thmA_table.out`; least log2 m at which the full exponent, all lower-order terms included, is
<= tau n; every variable occurring). Baselines (from memory, not re-checked): 2^n; #3-SAT O(1.6423^n)
[Kutzkov, IPL 2007] (there may be later improvements I do not recall); PPSZ 3-SAT 1.30704^n [Hertli, FOCS 2011 /
SICOMP 2014; randomized, decision only; Hansen-Kaplan-Zamir-Zwick STOC 2019 and Scheder FOCS 2021 are slightly
better, which moves tau by < 1e-3].

| row | change | r=4.26: 2^n | #3SAT | PPSZ | r=1: 2^n | #3SAT | PPSZ |
|---|---|---|---|---|---|---|---|
| A | paper as written (main term only, b4 model) | 2.76e5 | 5.52e5 | 1.10e6 | 8.05e3 | 1.66e4 | 6.84e4 |
| B | + W = 3m-n (chain wire count) | 1.37e5 | 2.76e5 | 1.10e6 | 3.78e3 | 8.05e3 | 3.39e4 |
| C | + I3 exact Hoeffding schedule, x0 free, all terms counted (I8) | 8.74e4 | 1.76e5 | 4.91e5 | 2206 | 3904 | 1.32e4 |
| C' | + I3 split zeta, entry 2^{N/K} | 8.03e4 | 1.56e5 | 4.62e5 | 1967 | 3667 | 1.23e4 |
| D | + I2 sweep base | 2.10e4 | 4.26e4 | 1.20e5 | 289 | 574 | 2837 |
| E1 | + I4a Lemma T4 top (deterministic) | 9375 | 1.90e4 | 5.34e4 | 151 | 299 | 1297 |
| E2 | + I4 Lemma 8' top (Las Vegas) | 4740 | 9530 | 2.68e4 | 108 | 202 | 696 |

Effect of each idea on the headline threshold (r = 4.26 vs 2^n): B->C (I3 + I8 accounting) x0.64; C->C' x0.92;
C'->D (I2) x0.26 (= "two levels"); D->E1 (I4a) x0.45; E1->E2 (I4) x0.51. Total vs paper: 2.76e5 -> 9.4e3
(deterministic), 4.7e3 (Las Vegas), i.e. factor 29 / 58 in log2 m*. These match b4's heuristic model
(D 2.41e4, E 5.4e3) and are now backed by Sections 1-3.

**I8 (explicit c, C0).** In the additive bookkeeping the paper's C0 survives only as kappa (poly cost per step,
once per level) and c only as c_ev (pieces computed in 2^{c_ev lam} ops) plus |I| = 4h+1 and S(lam), which are
explicit. `sensitivity.out`: every threshold of rows D-E2 is unchanged (to 4 digits) for kappa in {1e2, 1e4, 1e6}
and c_ev in {20, 2e3, 2e5}: at the thresholds n >= 2^{100} and lam_top >= 2^{42}, so these terms are < 1e-20 n.
appb_exact.py was not needed: c_s^row, c_s^col are given in closed form by (54)-(55).

---------------------------------------------------------------------------------------------------------
## 5. Caveats / what is not proved here
* Everything rests on the paper's Prop. 3 / Lemma 8 constructions and Theorem B (pieces with
  rho = 1/sqrt2) [R]; I re-derived only the cost accounting, Lemma S, Lemma T4, Lemma 8' and the split-zeta
  variant. The additive-overhead form of Thm A(b) is a proof sketch (Section 2), not a line-by-line rewrite.
* rho_bar(lam) uses the Hardy-Ramanujan leading term for log p(n) (+1 bit), not a certified inequality; its
  effect at the lam's used (>= 2^42) is < 1e-9 relative.
* Decision of an arbitrary n-variable 3-CNF needs sparsification [Cal09] with non-explicit M(eps); Theorem A'
  is stated for formulas with m clauses (so for density-r inputs directly). #SAT thresholds are compared with the
  decision bound of PPSZ only because the task asks.
* Schedules are optimised over a grid (x0 in 3..20, zeta2 = z0 q^l, z0 = 2^-2..2^-10, q in .5..9); a finer
  search can only lower the numbers. At realistic sizes (log2 m <= 40) nothing here beats brute force.
* No process left running; nothing written to /home/user/pnp-recherche.
