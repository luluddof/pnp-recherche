# G28-O2 notes (start 12:17 UTC)
Task: origin of entropy penalty in tensor-powered pieces; inherent? exact merged ranks k=2,3.
## Reading: BRIEF, O1 notes+scripts, B2 notes/finite_chain.py, paper Lemma 2 / Prop 3 (lines 900-1110).
## 12:45 Findings so far
Q1 origin (paper Lemma 2(b)): A_i = sum_{|Q|=i} (prod_{s in Q} C_s)(prod_{s notin Q} R_s); the row factorization index is theta = (Q, theta'),
 i.e. rank SUBADDITIVITY over the binom(k,i) sub-terms of one type class; "t bits to name Q" in Lemma 2 = the multiplicity.
 Sub-terms kept as separate iota instead: entropy paid at weight 1 (number of iota, not multiplied by e) and 2^{k|R|} iota -> worse.
Spectral theory (spec.py): M_row = R_row S^-1, M_col = C_col S^-1. All tested pieces: semisimple, integer eigenvalues in {-1,0,1,2,3,4}
 (none in (0,1)). In the eigen-word basis A_i is DIAGONAL: entry = coef of T^i in prod_s (m_s + T(1-m_s)); m outside (0,1) =>
 all roots m/(m-1) > 0 => no cancellation => EXACT merged rank = #{words: a eigen-0, b eigen-1, c other; a <= i <= a+c}.
 Verified exactly (mod 2 primes) for F(2,1),F(3,1) k=2,3 and (5,3)@4 k=2,3 (942 / 1899 at k=3,i=1 predicted and found).
Lower bound (any grouping into type-homogeneous classes, 0 semisimple in M_row, M_col): words with zero-set exactly Q
 => rank_row(G) >= sum_{Q in G} (N^2-p)^{|Q|} p^{k-|Q|}, rank_col(G) >= sum_{Q in G} (N^2-q)^{k-|Q|} q^{|Q|}; Hoelder =>
 h >= h_LB = log_{N^2}(p^mu (N^2-q)^{1-mu} + (N^2-p)^mu q^{1-mu}) (up to log(#classes)/k). Merged type classes ATTAIN h_LB.
 => merging only replaces N^2 by N^2-p, N^2-q in O1's formula; entropy stays. h_LB(0)=h_LB(1)=1, slope at mu=1 < 1 strictly
 => above circle near the ends; numerically above circle on (0,1) for N=2,3,4 pieces (asym.py).
## 12:50 more
- merged_runs.out: exact k=2,3 data for (5,3)@4,(4,4)@4,(6,1)@4 (k=3 4096x4096), oct k=2, N=5,7,8 k=2 (running).
- chain_merged.py (N=4 chain (4,4)->(6,1), k=2): middle piece A_2 = R_2-R_1 has row rank 9 (ideal p2-p1 = 2), col 5 (ideal 3);
  merged class ranks are far ABOVE the 'ideal' multinomial prod (p_s-p_{s-1})^{k_s}, slightly below subadditive.
- ideal_chain.py: even the IDEAL chain formula h = log_{N^2} sum_s (p_s-p_{s-1})^mu (q_{s-1}-q_s)^(1-mu) stays above the circle
  for every chain of known points (min (h-circle)/(mu(1-mu)) = +0.46..+0.60, h(1/2) = 0.85..0.92).
- Asymptotic theorem (single semisimple seed, large N): pieces R:(a,d), C:(c,b) (log_N of row/col sub-ranks) with Theorem C
  a+b>=2, c+d>=2 and p + rank_row(C) >= N^2, q + rank_col(R) >= N^2 (so max(a,c)=max(b,d)=2) => h(mu) >= 1-mu+mu^2 > circle
  (since (1-a)^2 >= 1-2a with a = mu-mu^2). Entropy is O(1/log N) per factor: matters only for small seeds, where it is
  exactly what Theorem C requires (mult-free (1,1)@2 powers would give pq = 1 << (2/9) 4^k: impossible).
- scan_LB.py with h_LB on (p,q) only: spurious 'violators' at p,q ~ N^2 (LB too weak there: ignores eigen-1 counts /
  Theorem C on the swapped pair (rank_row C, rank_col R)). Killed; replaced by the 4-parameter asymptotic argument above.
- NON-SEMISIMPLE loophole: if M_row has Jordan blocks at 0 (im M cap ker M != 0), sub-terms of one type class share row
  space: J(x)I + I(x)J has rank 2 not 4 (2-blocks). Model: row exponent H(x) -> H(x/2) when p = N^2/2 (x = C-fraction).
  tr M_row + tr M_col = N^2 (identity), so both sides cannot be nilpotent.
## 13:10 Exact data / Jordan / final
predict.py: spectral formula matches ALL exact merged ranks (k=2,3) for semisimple seeds: F(2,1),F(3,1),(5,3),(4,4),(6,1)@4,
 oct (8,8)@8, (6,6),(7,4)@5. Mismatch only for exact_7_11_7 and exact_8_22_5 = NON-semisimple seeds (jordan.py):
 7_11_7: rank M_row^j = 49,11,10,10 (one Jordan 2-block at 0); 8_22_5: Jordan block at eigenvalue 1 of M_row.
 Also dominated MC points mc_N4_2_8,(2,9),(3,9),(4,6) have Jordan blocks; their semisimple ranks (p_s,q_s) = a better known point.
General LB (any seed): lower semicontinuity of rank under conjugation (semisimple part M_s in orbit closure) + tops =>
 h >= h_LB(p_s,q_s), p_s = N^2 - algmult_0(M_row), q_s = N^2 - algmult_0(M_col). All known seeds: min (h_LB-circle)/(mu(1-mu)) >= +0.60
 (jordan.py output). Product of different seeds: h is a log-weighted average -> no gain from mixing.
Oct check: merged rho = h_LB(1/2) = log_64(2 sqrt(8*56)) = 0.9006 vs paper Lemma 2 11/12 = 0.9167 vs mult-free log_64(sqrt(8*64)) = 0.75.
Interface: semisimple merged A_i factors through the admissible eigen-words (a<=i<=a+c), rankable in poly(k); entries products
 of rational eigenvector entries -> fits Definition 1 (row and column factorizations both exist). Gain usable but tiny.
Jordan loophole: nilcheck.py verifies 2-block formula (inclusion matrices); jd_check.py: d>=3 blocks NOT full rank on weight spaces.
 nilmodel.py: idealized nilpotent-row (2-blocks) seeds, scan over (p,q) with pq>=2N^2/9: min normalized gap +0.019..+0.032 (q=1);
 nilmodel2.py with proven constraints (q=1 => M_row semisimple-only [tr argument: ((A(x)B)S)^2 = AB(x)BA], q=2 => p>=N(N-1)/2):
 +0.53 (N=2) .. +0.14 (N=32). Never below the circle. q=1 nilpotent impossible: tr(I-M)=N^2 vs tr(AB)=+-N.
## 12:50 (clock) wrap-up
jc.py: Jordan-Chevalley semisimple part (R_s = M_s S) is NOT always a better decomposition (mc_N4_2_9: rank 1 but osr 10 > 7). Dropped.
General tops LB WITHOUT semisimplicity: with t = dim(ker M cap im M): rank_row(class i) >= #{words over (N^2-p-t, t, p-t): a <= i <= a+xi}
 (t=0: C(k,i)(N^2-p)^i p^{k-i}); Jordan part makes 'xi' positions flexible -> weaker LB (loophole only at small N).
Extreme classes i=0,k are single sub-terms (exactly multiplicative, no entropy): h >= max(mu l(p)+(1-mu) l(N^2-q), mu l(N^2-p)+(1-mu) l(q))
 (l = log_{N^2}) for EVERY seed; with Theorem C => h >= 1-mu+mu^2 - O(1/log N) > circle for huge N (trivial large-N regime;
 bites only for N^2 >~ 1e20). Small/moderate N: decided by entropy = exact merged formula (all known seeds above circle).
asym.out: merged h(1/2) table. Files: pieces.py merged.py spec.py predict.py asym.py jordan.py chain_merged.py ideal_chain.py
 nilcheck.py jd_check.py nilmodel.py nilmodel2.py jc.py; outputs merged_runs.out asym.out jordan.out ideal_chain.out.
