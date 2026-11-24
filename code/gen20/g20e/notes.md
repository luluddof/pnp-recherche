# G20e — Instance-adaptive variants of the two Oct-2026 sub-exponential 3-SAT algorithms (working notes)
Tags: [P] proved here (short argument, unreviewed), [N] numerical (own code, this folder), [H] heuristic, [L] literature recalled
(not re-checked), [O] open.  Papers are taken "as written"; their correctness is not assessed.

## 0. Plan / save 0
1 list oblivious steps; 2a cut-tensor ranks (exp_a_cutrank.py); 2b expansion of random 3-CNF graphs (exp_b_expansion.py);
2c branching to destroy expansion (exp_c_branch.py + mixing-lemma bound); 3 verdicts.

## 1. Where each algorithm is oblivious (save 1)
(A) Claude (sec. 2.1-2.6, App. A/B), (B) OpenAI App. D.
 O1 (A, B) Sparsification + clause/chain graph: every variable of degree d becomes d-1 one-bit wires (M ~ 3m - n; at
    density 4.26, M ~ 11.8 n).  The wire count M, not n or a separator, is what the exponent multiplies.
 O2 (A) Random grouping into K groups + padding (sec. 2.2 "(K,lambda)-clique", Step 3, Lemma 1): every pair of groups is
    given exactly lambda = max bits; the algorithm tries all K^g groupings but only for BALANCE, never for small cuts.
    (B) Lemma D.3: derandomised uniform labelling valid for every network; loads made near-uniform on all r^{2a} pair
    classes (needed by the Cauchy/packet machinery, G19a sec. 2).
 O3 (A, B) The crossing/packet gadget is the FULL swap/identity on 2^lambda x 2^lambda (resp. 2^L) states, independent of
    which states the adjacent node tensors can actually produce (W_a are treated as arbitrary tables / queried procedures).
 O4 (A) Fixed geometry: g x g grid, all rectangles, all rows OR all columns merged, K_l = 9 K_{l-1}^16 fixed in advance;
    (B) u = Theta(log s / 2^d) valuation levels and children of M/8 wires fixed by s only.
 O5 (A, B) Uniform group sizes (<= 2g/K nodes; beta-balanced directions) and no use of community / separator structure.
 O6 (A, B) Running time is a function of (S or M, K, lambda, Delta) only: the instance enters only through these sizes.

## 2a. Structured wires  (save 2)
[P] Prop. R (restricted crossing is a smaller full swap).  Let bundle 1 join nodes a,d and bundle 2 join b,c, N = 2^lambda.
 Fix all other tensors; let V_a in K^N be the column space of the flattening of a's tensor onto its bundle-1 legs (all other
 legs and sums absorbed), V_d likewise, r_1 = rank of the pairing V_a x V_d -> K given by the identity, same r_2.
 Then the bundle can be replaced by a wire of dimension r_1 without changing Z, and the crossing restricted to
 (V_a (x) V_b) x (V_d (x) V_c) is, after a change of bases, the full swap on K^{r_1} (x) K^{r_2}.  Proof: a rank-r bilinear
 form factors as I_r through K^r; the swap of two such factorisations is the swap of K^{r1}, K^{r2}.  The same holds for
 restriction to a SET S of reachable states (span of indicator vectors of S = K^S).  Hence Theorem C (G19e) and the
 window/flattening bounds apply verbatim with N replaced by r: structured wires only replace lambda by log2 r.
 The relevant compression is therefore the CUT RANK.  For a variable bipartition (A,B) of a CNF:
   Mcut = diag(W_A) C diag(W_B),  W_A(alpha) = #extensions inside A,  C = product of crossing clauses,  Z = 1^T Mcut 1.
 Three upper bounds on rank Mcut: 2^{min(|dA|,|dB|)} (cut-conditioning / Claude's edge bits), 2^{#crossing clauses}
 (each clause 1 - [lits_A false][lits_B false] has rank <= 2 across ANY split of its variables), and Z (#models).
 WHERE THE HARDNESS GOES [P]: rank Mcut <= Z and rank Mcut = 0 iff UNSAT; so the "reachable" restriction is the decision
 problem itself.  Locally computable parts of the structure (support of W_A, clause rank 2) need the whole table W_A,
 i.e. 2^{|dA|} subproblem calls: exactly the cut-conditioning enumeration the algebraic methods are designed to avoid.

## 2b. Instance-dependent labelling / grouping  (save 3)
Expansion numbers [N] (exp_b_expansion.py -> out_b.txt; random 3-CNF, density 4.26, one instance per n, weighted graphs):
   graph   n    M     spectral LB bisection/M   KL UB bisection/M   Cheeger LB on conductance
   chain  100  1178   0.118                     0.199               0.121
   chain  200  2356   0.111                     0.204               0.115
   chain  400  4712   0.113                     0.203               0.117
   primal 100/200/400 (M = 3m co-occurrence edges): spectral LB 0.18/0.18/0.145, KL UB 0.35/0.35/0.35, Cheeger 0.29.
 So for the chain network (what both algorithms contract): 0.11 M <= min bisection <= 0.20 M, i.e. 1.3 n .. 2.4 n wires,
 conductance phi >= 0.115 (certified per instance, stable in n).  NB: a balanced cut of the CHAIN graph costs MORE than n
 bits (a variable's chain can be cut several times); cut-conditioning on variables (primal vertex separators) is the
 classical comparison, and it is <= n.
Best exponent with G19a's inequalities, as a function of phi (one level, model DF, u beta-balanced directions):
   upper (cut-conditioning, any labelling): min_i M_i, and >= phi' M for every balanced direction (3.2);
   lower for every identity (Theorem E):    phi' M/(2u) - u   (r = 2: phi' M/u - u);
   App. D (uniform labelling):              c M/sqrt(u).
 Adaptive labelling can only matter in the window [phi' M/(2u), c M/sqrt u]; the only identity known for non-uniform
 labellings is cut-conditioning (>= phi' M), so in that window NOTHING is constructed (G19a's open staircase problem).
Does 2^{phi' M} ever beat 2^{c M/sqrt u}?  Asymptotically never on expanders (phi' fixed, u -> infinity).  In any physical
 range YES, by far [H, constants are estimates]: (B) has u <= a/2 with r^{2a} <= s^{1/10}/2 (Delta = s^{9/10}, D.6), i.e.
 u <= log2(s)/80 for r = 2, and c in the tens (G18c sec.2(ii)); c/sqrt u < phi' = 0.11 needs u > 83 c^2 ~ 8*10^3, i.e.
 log2 s > 6*10^5.  Even beating brute force 2^n needs 11.8 c/sqrt u < 1, u > 139 c^2, log2 s > 10^6.  For (A) with
 rho = 11/12 and K_l = 9K^16: e_l ~ 12 (11/12)^l < 0.11 needs l >= 54 levels, log K ~ 16^54: astronomically worse.
 (A with App. A/B, rho = 1/sqrt2, log K = Theta(2^l): l >= 14, log2 m >~ 9 * 2^14 log2 64 ~ 10^6.)
 So: for every n < 2^(10^5) the instance-ADAPTIVE cut bound (and brute force, and PPSZ 2^{0.386 n}) is smaller than
 both oblivious bounds; asymptotically the oblivious ones win on expanders.  Neither fact is a new algorithm.

Small-separator classes [H, composition argument; ALREADY KNOWN baseline]:
 Classical: #SAT in 2^{O(tw)} poly; planar / bounded genus incidence graph: tw = O(sqrt n) (genus g: O(sqrt(g n))),
 so 2^{O(sqrt n)}, versus the oblivious algebraic 2^{O(n/sqrt log n)} (worse for every n).
 Adaptive hybrid H(x): cut the formula into x regions (r-division), compute each region's boundary tensor by tree-
 decomposition DP, run App. D (Theorem D.1 allows arbitrary vertex tensors of degree <= s^{9/10} computable in
 2^{A s^{9/10}}) on the region network with S boundary wires.  Exponent ~ S/sqrt(log x) + (region cost).
  * separators n^alpha (planar alpha = 1/2, d-dim grids 1-1/d): S(x) = Theta(x^{1-alpha} n^alpha), region cost
    (n/x)^alpha:  exponent >= n^alpha [ x^{1-alpha}/sqrt(log x) + x^{-alpha} ], minimised at x = O(1):  Theta(n^alpha).
    CONSTANT ONLY (the sqrt(log x) gain is eaten by the polynomial growth of the boundary).  For planar: 2^{Theta(sqrt n)}
    either way; no 2^{o(sqrt n)} from this lever.
  * separators n/f(n), f slowly varying (hereditary): an r-division pays ~ n/f per recursion level, S(x) ~ n log x / f,
    region cost (n/x)/f: exponent ~ (n/f) [ sqrt(log x) + 1/x ] >= n/f: CONSTANT ONLY (boundary grows like log x, gain
    like sqrt(log x)).
  * Non-hierarchical "expander of communities" (x communities of size L = n/x, O(1) boundary each, community graph an
    expander): S = O(x), tw = Theta(x): classical 2^{O(n/L)}; hybrid 2^{O(n/(L sqrt(log(n/L)))) + O(L)}: a sqrt(log)
    factor over tw-DP and a factor L over the oblivious algorithm.  This IS beyond a constant, but it is nothing but
    "kernelise to the S interface bits, then run the oblivious algorithm on S bits" (S = tw here); the algebraic part is
    still oblivious.  Unproved: Theorem D.1 as stated (taken at face value) + bipartite subdivision of the region network.
  * Power-law (Chung-Lu, exponent > 2) random graphs: linear treewidth w.h.p. [L, recalled: Gao 2012] -> behaves like
    the expander case; industrial instances with small tw are already easy for 2^{O(tw)} DP.
 General statement [P modulo Thm D.1]: for any partition of the variables into regions with S wires between regions and
 region treewidth <= S^{9/10}, #3-SAT is solvable in 2^{O(S/sqrt log S)} poly + sum of region DP costs.  On a graph with
 conductance phi, any partition into regions of volume <= M has S >= phi M: CONSTANT ONLY on expanders (factor <= 1/phi).

## 2c. Branch on b variables to destroy expansion, then algebra  (save 4)
[P] Wire accounting.  Fixing a variable v deletes its d_v - 1 chain wires, and each clause it satisfies disappears,
 deleting at most 2 more wires (one per other variable: the chain is re-linked around the removed clause node).  So the
 residual chain network G_B has M_B >= M - 3 sum_{v in B} d_v wires, and for every bipartition (S,T) of the residual
 clause nodes, cut_{G_B}(S,T) >= cut_G(S',T') - 3 sum_B d_v, where S',T' extend S,T by the removed clause nodes (at most
 sum_B d_v of them, so balance degrades by that amount).  Hence
     minbis(G_B) >= minbis_{slack sum_B d}(G) - 3 sum_B d_v >= phi' M - 3 d_max b   (conductance form).
 d-regular spectral version (expander mixing lemma, second eigenvalue lambda): for |S| = s,
     e_{G-B}(S, V\B\S) >= (d - lambda) s (n-s)/n - d s b/n - lambda sqrt(s b),
 positive for balanced S as long as b/n < c((d-lambda)/d)^2-ish: a CONSTANT FRACTION of the vertices must be deleted.
 To make a balanced direction cut <= eps M (needed for Hamming-1 loads with u = 1/eps, G19a (3.4), or for a cheap
 separator) one needs  sum_B d_v >= (phi' - eps) M/3,  so b >= (phi' - eps) M/(3 d_max).
 Random 3-CNF, density 4.26: degrees ~ Poisson(12.8); with phi' >= 0.11 (certified above) and the top-degree variables
 (d <= ~30 for b >= 1e-6 n): b >= ~0.11*11.8 n/90 ~ 0.014 n even for eps -> 0.  Cost 2^{Omega(n)}.
Trade-off [P modulo "the residual costs c M_B/sqrt u"]: T(b) >= b + c (M - 3 d_max b)/sqrt u; dT/db = 1 - 3 c d_max/sqrt u > 0
 once sqrt u > 3 c d_max, so the optimum is b = 0 asymptotically: no interior trade-off point, and no point beats
 n/sqrt(log n) by more than a constant.  Branching + cut-conditioning (no algebra) costs b + phi''M = Theta(n).
[N] (out_b.txt) deleting the top-degree variables, n = 400: KL bisection 955 wires (b = 0) -> 389 wires (b = 120),
 spectral certificate falls to 0.02 M_B (weak, leaves of degree 1); b + bisection = 509 < 955: branching on high-degree
 variables beats chain-cut conditioning by a CONSTANT (chains are over-counted by cuts), consistent with the classical
 picture (DPLL / resolution need 2^{Omega(n)} on random 3-CNF because restrictions of Omega(n) variables keep boundary
 expansion [L: Chvatal-Szemeredi; Ben-Sasson-Wigderson]).

## 3. Verdicts  (save 5)
 L1 Structured wires / restricted crossing (O3):  NO GAIN in the exponent beyond a constant.
    Prop. R: P|_V is the full swap of dimension r = cut-pairing rank, so Theorem C and the window bounds still apply with
    lambda -> log2 r.  Locally computable structure (clause rank 2, support pruning of W) gives a constant factor at best
    (measured log2 r / bits ~ 0.4-0.8, small n only); the global restriction (reachable states) has rank <= #models and is
    0 iff UNSAT, so it is the problem itself; computing either costs the 2^{boundary} table enumeration.
 L2 Adaptive labelling for the algebraic level (O2, O5) on expanders: NO KNOWN GAIN.  Theorem E (G19a) leaves the window
    [phi' M/(2u), c M/sqrt u] open, but the only identity known for non-uniform labellings is cut-conditioning, which costs
    >= 2^{phi' M}; phi' >= 0.11 certified for random 3-CNF chain networks (n <= 400, density 4.26).
 L3 Adaptive coarsening / separator hybrid (O1, O4, O6):
    - expanders: CONSTANT ONLY (any coarsening keeps S >= phi M);
    - polynomial separators (planar, bounded genus, grids): CONSTANT ONLY; 2^{Theta(sqrt n)} = ALREADY KNOWN tw-DP;
    - hereditary n/f(n) separators: CONSTANT ONLY (r-division boundary grows like log x, gain like sqrt(log x));
    - "expander of small communities" (S << M, non-hierarchical): GAIN beyond a constant (factor L over the oblivious
      exponent, sqrt(log) over tw-DP), but it is kernelisation + the oblivious algorithm (ALREADY KNOWN ingredients);
      unproved steps: Theorem D.1 at face value for region tensors.
 L4 Branch-then-algebra (2c): NO GAIN (b must be Omega(n) to lower expansion; T(b) increasing in b asymptotically).
 L5 Fixed constants: in every physical range (n < 2^(10^5)) the instance-adaptive cut bound, brute force and PPSZ are
    all far below both oblivious bounds; asymptotically the oblivious bounds win on expanders.  [H: constants estimated]
 Bottom line: obliviousness is not where the slack is.  On the instances that are hard (expanders) every adaptive lever
 is capped by a constant through S >= phi M or Prop. R; on instances where adaptivity helps beyond a constant, the
 help comes from classical structure (separators / kernels), and the algebraic factor merely multiplies it.

## 4. Adversarial self-review
 * Prop. R fixes all OTHER tensors; in the recursion the other tensors are themselves sums of pieces, so the reachable
   subspace changes per term.  The proposition still shows that no single restriction escapes Theorem C.
 * "Constant factor" for local structure is from n <= 20 experiments plus the trivial bounds; no asymptotic lower bound
   on rank C for random 3-CNF across balanced cuts was proved here (literature on rank / d-SDNNF lower bounds for CNF
   was not checked).  An o(1) ratio log2 rank C / bits for large n is NOT excluded by my data.  [O]
 * Expansion numbers: one instance per n; spectral bound is rigorous for the instance (eigsh numerical); KL is only an
   upper bound; the chain order is random (Claude's chain order is arbitrary too).
 * The hybrid exponents of 2b (S/sqrt(log x) + region cost) assume Theorem D.1's bound with u ~ log x, i.e. that the
   number of regions limits u; this is my reading of (D.6) (loads need S/Delta >= poly), not a theorem.
 * 2c: assumes the residual is attacked by the same O(M/sqrt u) bound; a cleverer residual algorithm is not excluded,
   but it would then also work at b = 0.
 * Constants c of App. D are G18c's rough estimate ("tens"); the 2^(10^5) threshold scales like c^2.
 * Neither paper's correctness was assessed.  Nothing here is a lower bound for 3-SAT.

## 5. Measured structure of real cut tensors [N]  (exp_a_cutrank.py -> out_a.txt, out_a2.txt, out_a3.txt;
##    exp_a2_crossrank.py -> out_a2cross.txt).  Ranks: exact mod 2^31-1 for sizes <= 256, float SVD above, cross-checked.
 Variable bisection (KL, chosen to minimise boundary) of random 3-CNF; bits = min(|dA|,|dB|) = what cut-conditioning pays.
   n  m   cross  |dA| |dB|  supp W_A / 2^dA   rank C      rank C restricted to supports   Z (#models)
   16 68   47-57  8   8     62-160 /256      6-19 (2^2.6-2^4.3)    0-4                  0-11   (density 4.26, near UNSAT)
   20 85   56-59  10  10    94-188 /1024     20-50 (2^4.3-2^5.6)   0-2                  0-24
   20 40   21     10  9-10  286-388 /1024    153-240 (2^7.3-2^7.9) 34-65 (2^5.1-2^6.0)  2091-3962 (density 2)
 Crossing-clause product alone (exp_a2_crossrank): log2 rank C / k = 0.69-0.81 for k = 4..9 at c/k = 1, 2, 4, FLAT in k.
 Reading: (i) the clause product across a cut is NOT full rank: ~0.7 bits of bond per boundary bit (clause rank 2 +
 repeated falsification patterns) -- a constant factor, no decay with k visible in this range;
 (ii) local support pruning removes 1.5-2.5 bits of 8-10 (constant fraction);
 (iii) the restriction to globally reachable states gives rank <= Z and = 0 on UNSAT instances (rows with Z = 0):
 the dramatic drops are exactly the instances where the answer is already "UNSAT/few models"; they cannot be exploited
 without solving.  All at n <= 20: small-size data, nothing asymptotic.
 Runs at n = 24-48 with boundary > 10 bits were skipped (too slow on this machine); out_a3 was stopped after 3 rows.
