# G25-B1 — brainstorm: new constructions for the crossing trade-off (Naslund P12)

Setting (ev.py convention): S = swap on K^N (x) K^N; (p,q) feasible iff S = R + C, rank R <= p, osr C <= q.
Useful restatements used below: (a) "transpose map = elementary operator of length q on a subspace of codim p";
(b) 4-party: |Omega>_{x v}|Omega>_{y u} = R + C, rank across {x,y}|{u,v} <= p, across {x,u}|{y,v} <= q (two EPR pairs,
both cuts cut both pairs); (c) App B framework: 4rho = max_s [L(2h+s)+L(2h-s)]/L(h), L = log dim V_{<=E}.
Tags: [E] exact computation here, [A] short analytic argument here (unreviewed), [H] heuristic.

## Ideas (a idea / b why / c likely failure / d tried? / e first test / f prob.)

1. Laser method / type-zeroing on tensor powers. S_2 = I + (S_2 - I) has only the two endpoint pieces (2,0),(0,2);
   S_{2^k} = sum_T R^T (x) I^{T^c}; middle types are individually cheap (N,N) but binom(k,j)-many. Restrict to a sub-swap on
   a code X' in F_2^k hoping to kill multiplicity.  c: restriction does not select pieces (pieces are indexed by intermediate
   data, not by external indices).  d: new framing.  e: DONE [E] laser.py/laser2.py, codes even5/6, Hamming, RM(1,3),
   simplex, all assignments of |T|-classes: best 10@8 (=Pauli window), 56@16, 72@16, 332@32 — worse than direct sums
   (32@16, 128@32).  General form (base pieces block-supported on external indices) remains open.  f: 1%.
2. Code loops (Griess/Parker Moufang loops from doubly-even codes) = the natural non-associative generalisation of the
   octonion loop, outside G23a's group ansatz.  e: DONE [E] codeloop.py, cocycle identities verified; composition ansatz:
   octonions recovered (p=8@8), simplex4 N=16 p=256, RM(1,4)/simplex5 N=32 p=1024, RM(1,5) N=64 p=4096 (trivial).
   Reason: e_a^2 = (-1)^{|a|/4} = +1 for weight-8 words (split form) — the anisotropy/Chevalley–Warning obstruction of G23a
   again.  KILLED.  f: 0.1%.
3. Rank-metric / Gabidulin coordinates: M_N(F_q) = sum_i F_{q^N} x^{[i]}; transpose = slot reversal i -> -i with
   Frobenius twist; elementary operators = slot shifts.  [A] over F̄_q this IS the Weyl clock/shift basis (G20c b': p+q=N^2);
   ranks are field-independent.  MRD codes as the subspace W: irrelevant (only codim matters).  KILLED.  f: 0.1%.
4. MUBs / weighted 2-designs: S + I = sum_{b=0}^{N} P_b, P_b = sum_i |e_i^b e_i^b><.|, each piece rank N, osr N (perfect
   (1,1) pieces, better than any App B piece).  [A] R-side overlaps only inside Sym^2, C-side MASAs are orthogonal =>
   (mN, (N+1-m)(N-1)+1): the trivial line p+q ~ N^2.  Twist: weighted designs whose points sit in few low-dim subspaces
   (R side) / few MASAs (C side); the design identity forces spreading.  e: exact LP over subsets of a complete MUB/SIC at
   N=4,5,7 with flint ranks.  f: 0.5%.
5. Quantum-group R-matrix at a root of unity (R = q^{H(x)H/2} sum_{n<l} c_n E^n (x) F^n has osr <= l^2...).  [A] the object
   close to S is Ř = S R, whose osr is N^2; S - Ř is diagonal on components (rank = dims with q^{c_k} != 1).  No low-osr
   object appears.  KILLED.  f: 0.2%.
6. Casimir interpolation for other (g,V): S = f(Omega_g), osr(deg d) = dim U_{<=d}(g)|_V, rank = non-interpolated
   components.  e: DONE [E/A] casimir.py, optimistic exponents with n -> inf: fermions Lambda^h(C^n) = 1.625 (= Hrubes 2H(1/4),
   sanity check), bosons Sym^h(C^n) = 1.918.  Spinors = Clifford (capped >= 1.5, G24a).  Exceptional: fixed sizes.  f: 0.3%.
7. Holographic / modular families (symmetric orbifolds, extremal VOAs, c -> inf, sparse light spectrum).  [A] Hagedorn-then-
   Cardy profile L = 4piE (E<=E*), 4pi sqrt(E E*) (E>=E*): min_h of 2L(2h)/L(h)... gives exactly sqrt2 (attained h>=E*);
   a sparser light spectrum only lowers N.  Modular invariance (HKS) forces Cardy above c/12 => >= sqrt2 for every modular
   family: the Cardy barrier extends from fixed VAs to c -> inf families.  KILLED.  f: 0.3%.
8. Convex-jump ("two-temperature") profiles: framework needs L(h)=1, L(2h)<1.41, L(3h)<1.83, L(4h)>=2 — log-dim growth that
   ACCELERATES on [h,4h] (non-concave; G20c F2 only excludes concave).  Mechanism: the two-particle products a_(-1)b of
   weight-h states nearly null up to 2h (huge relations) but nondegenerate pairing kept, e.g. W (x) V' with W a finite graded
   Frobenius/Gorenstein piece (sup-convolution of a saturating Hagedorn piece with a slow piece gives exactly the kink), non-
   modular simple quotients at collapsing levels, Heisenberg-paired jet algebras of fat points.  c: Lemma 9 (locality) may
   need a genuine VA, and non-concave growth may be incompatible with nondegenerate pairing + locality.  d: partly G20c F2,
   G24b, G23c (open "multi-channel tensor power" gap 2.799–3.030).  e: paper check of which axioms Lemmas 9–10 use; then
   profile of a toy W (x) Heisenberg (dims from Hilbert series) and of Stanley-type non-unimodal Gorenstein W.  f: 2%.
9. Catalytic decompositions: is (p,q) for S_N (x) Z (fixed catalyst Z with rank r, osr s) better than the product scaling?
   Matrix-multiplication analogues (catalytic/support rank) do gain.  c: Theorem C with spectator, and the algorithm may not
   accept a catalyst.  e: numerics (vp.py-type) at N=2,3 with Z = small structured (Omega, swap_2, random) .  f: 1%.
10. Border / degeneration families: border values 3,5,6 at N=4,5,6 are below exact ones; an epsilon-degree e family costs only
   a factor e+1.  Look for S = lim g_eps (R + C), g_eps in a GL^4 torus.  c: G24e numerics already see a jump after 8 that
   presumably includes border solutions.  e: exact epsilon-polynomial ansatz (degree 1–2) for (7,7)@7 and (9,9)@9, and
   structure analysis of the known border (3,3)@4, (5,5)@5, (6,6)@6 to extrapolate.  f: 1.5%.
11. Association-scheme (non-group) diagonal family: "sparse + scheme-Fourier-sparse" with Bose–Mesner of non-symmetric
   schemes (cyclotomic, Paley tournaments).  d: G22c (tournaments N^2 - O(N)), G20c sect. 3.  f: 0.5%.
12. Data hiding (Sym vs Lambda indistinguishable by PPT/LOCC; Eggeling–Werner): gives LOWER-bound intuition (low-osr C cannot
   separate Sym/Lambda), not a construction.  Discarded as a construction idea.

## Top 3 (probability x value / cost)
1. Idea 8 (convex-jump profile): 2%, highest value (any such family gives rho < 1/sqrt2, rho -> 2/3 possible).  First test:
   list the axioms Lemma 9/10 actually use; compute the profile of W (x) Heisenberg for a Gorenstein W with large middle.
2. Idea 10 (border family): 1.5%, cheap-moderate.  First test: exact epsilon-ansatz at N=7 (7,7) and N=9 (9..10).
3. Idea 9 (catalyst): 1%, cheap.  First test: numerical feasibility of S_2 (x) Z, S_3 (x) Z with Z = swap_2 / Omega_2.

Files: codeloop.py (idea 2), laser.py, laser2.py (idea 1), casimir.py (idea 6).  No process left running.
