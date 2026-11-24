# G22a — primary literature on the swap / Hurwitz / Naslund problem (reading notes, saved as I go)

Tags: [R] = read in the paper itself (page = page number printed in the PDF / ECCC page), [P] = argued here, [O] open.
Downloads (untrusted data, read only): dl/ (a<arxivid>/p.pdf,.txt ; eccc/<year>_<nr>.pdf,.txt).

## 0. Reading status (updated at the end)
see section 9.

## 1. Naslund, arXiv:1701.04475v3 (JCTA), "The partition rank of a tensor and k-right corners in F_q^n"
Status: READ PARTLY: sections 1, 2 (incl. Problems 12, 13) and bibliography in full; sections 3-5 (right corners) skimmed.
* Def. 8 (p.5): for a collection P of non-trivial partitions of [n], P-rank = least r with F = sum of r functions each a product
  over the blocks of SOME partition in P.
* Problem 12 (p.8, verbatim content): X finite, |X| = n, k = 4, P = { {12|34}, {13|24} }. "Does delta(x1,x4) delta(x2,x3) have
  superlinear P-rank? That is, P-rank >> n^{1+eps} for some eps > 0?"  Context (p.8): counting gives tensors of P-rank >> n^2,
  "no explicit superlinear lower bounds are known for any tensor"; "a positive answer would lead to improved lower bounds for
  non-commutative circuits, see [9, Theorem 3.6]", [9] = Shpilka-Yehudayoff survey (FnT TCS 2010).  Problem 13: P-rank of
  delta_P for P not refining any partition of the family.
* Dictionary [R+P]: our X(x,y,u,v) = [x=v][y=u] with legs (x,y,u,v) = (1,2,3,4) is exactly delta(x1,x4)delta(x2,x3);
  12|34 = (x,y)|(u,v) = rr, 13|24 = (x,u)|(y,v) = rc.  So Naslund's P-rank = min{p+q : (p,q) feasible} = our P(N), EXACTLY
  (no constant).  Naslund asks only about the SUM, over an arbitrary field F (field unspecified in Problem 12).
* Technique in the paper: Lemma 11 (partition rank of a diagonal tensor = support), proof by induction, contracting a variable
  with a vector v of large support in the annihilator of the slice-rank-1 terms. This is the same "annihilate the terms of one
  type, then count rank" mechanism as our Theorem C (we use a k-dim family instead of one vector).

## 2. Karam, arXiv:2310.11356v1, "The interplay between bounded ranks of tensors arising from partitions"
Status: READ PARTLY: sections 1 (all), 2 (all), 3 (all), 4.1; 4.2 skimmed.
* Question 1.2 (p.4) = Naslund Problem 12 verbatim (R = {12|34, 13|24}, T = 1_{x1=x4} 1_{x2=x3}, R-rank >= n^{1+eps}?).
  Karam attributes the motivation as: would improve "Nisan's noncommutative formula-size lower bound ([15], Theorem 3.6)".
* The only quantitative bound stated for Q1.2 (p.4-5): r + s >= n, by fixing x1 = x4 = 1 (each term becomes rank <= 1 in
  (x2,x3), RHS = identity of rank n).  "no lower bounds of the type Omega(n^{1+eps}) are known for the R-rank of any specific
  tensor" (p.5).
* Theorem 1.4 (p.6): if P not finer than any partition of R then R-rank(Delta_P) >= n^{c(d)}, c(d) = 3^{-2^{d+1}}.
  For d = 4: c = 3^{-32}.  This is the "power lower bound ... qualitative answer to a question of Naslund" of the abstract:
  it answers Problem 13 qualitatively (power bound, not linear), NOT Problem 12 (where n is already trivial).  Earlier notes
  (G18g "I do not know its exponent") are hereby settled: exponent 3^{-32}, irrelevant for us.
* Techniques: (a) "fragmentation" (Prop. 2.5): two decompositions of the same tensor w.r.t. two families => a decomposition
  w.r.t. a finer family with length k1((k1+k2)^2+1); uses dual functions of small support (support <= r).
  (b) Example 2.3 (p.12): a {12|34}-decomposition of length r and a {23|14}-decomposition of length s of the SAME tensor =>
  tensor rank <= r^2 s^2.  (c) Case 1 of Thm 1.4 (p.17): u of support <= r+1 killing the A_i, then partition rank <= s(r+1),
  giving min(.., n^{1/2}).  (d) Prop. 4.1: {P}-rank r and {Q}-rank s => {P ^ Q}-rank <= (rs)^{|P|}.
  None gives anything beyond n for Problem 12; all are "bounded vs unbounded" (qualitative) tools.
* NB: Karam's setting is ranks of a SINGLE tensor in two families simultaneously (T has small R1-rank AND small R2-rank);
  our problem is a SUM T = R + C with R of small R1-rank and C of small R2-rank.  Different questions; his Thm 1.9 is about
  intersections, ours about sums.  Nothing in Karam on trade-offs (p, q).

## 3. Hrubes-Wigderson-Yehudayoff, ECCC TR10-021 (JAMS 24 (2011) 871-898). Status: READ PARTLY: sects 1-7 in full, sect 8 skimmed.
* S_F(k) (p.2) = least n with (sum x_i^2)(sum y_j^2) = sum z_i^2, z_i bilinear. B_F(f) (p.4) = least n with f = sum z_i z_i'.
  Remark 1.5: B <= S; alg. closed, char != 2: S <= 3B.
* ID_k = sum_{i,j} x_i y_j x_i y_j (1.4). Lemma 3.5 (p.14): width w(ID_k) = least n of terms each g[X1,X2]h[X3,X4] (split 12|34)
  or sum_i h_i[X1] g[X2,X3] h'_i[X4] (split 23|14). Tensor = delta(1,3)delta(2,4); after swapping legs 3,4 it is Naslund's tensor.
  So w(ID_N) = P(N) = min p+q EXACTLY [R+P].
* Prop 4.1 (p.14-16): (i) B(SOS_k) <= w(ID_k) (any field: BOTH term types become products of two bilinear forms);
  (ii) char != 2: w <= 4B; the proof turns each product z z' into 2 terms of each type => feasible (2B, 2B).
* Prop 3.2: circuit size C(f) >= Omega(d^-3 w(f)). Thm 1.6: C(ID_k) >= Omega(B). Thm 1.7 / Cor 5.4: B(SOS_k) >= k^{1+eps}
  => PERM_n needs 2^{Omega(n)} nc circuits; proof uses only w(ID_k) (so P(N) >= N^{1+eps} suffices, any field), via lifted
  identity LID_r and Lemma 5.2 (central -> 2^r block-central, loss 2^r, Prop 5.3). Thm 1.8: any explicit biquadratic, char != 2.
* Thm 1.9 / sect 6: explicit real biquadratic with S_R = Theta(k^2) (vector matrices; uses "unit real vectors with inner product 1
  are equal" -- fails over C). Thm 1.10 / sect 7: S_Z(k) = Omega(k^{6/5}): intercalate matrices (Prop 7.1), colour class of size
  >= k^2/n gives a full s x s intercalate submatrix, Lemma 7.4 full => Omega(s^{3/2}) colours (Prop 7.8: each ordered colour pair
  occurs in <= 2 positions, via 9x9 forbidden configuration, Lemma 7.7).
* Background stated (p.2-3): S_R(k) > k unless k in {1,2,4,8} (Hurwitz); S_R(k) >= (2-o(1))k (James; K.Y. Lam); S_Z(k) <= O(k^2/log k)
  (Hurwitz-Radon, integer, coeffs in {0,+-1}); Pfister: rational z_i give k squares for k a power of 2.
* Nothing on (p,q) trade-offs: only the sum.

## 4. Hrubes, ECCC TR24-026 (CCC'24, 12 pp.) READ IN FULL; TR24-099 (22 pp., extended) READ IN FULL.
WHAT is bounded: sigma_F(n,m) = least s with (sum_{i<=n} x_i^2)(sum_{j<=m} y_j^2) = sum_{k<=s} f_k^2, f_k bilinear, F = C, any field
 with sqrt(-1) (TR24-099 Thm 1, e.g. Q(i)), or char p > 2.  sigma(n) <= O(n^{1.62}) (exponent 1 + 2H(alpha/2), H(alpha) = 1/2,
 alpha = 0.110: 1.6147). beta_F(n) (products f_k f_k') <= O(n^{1.62}) over ANY field (TR24-026 Thm 12 sketch; TR24-099 Thm 17:
 signed squares over any field, so beta <= 2s). sigma_C <= 2 beta_C.
HOW (not recursive; one algebraic step): (1) Lemma 3/5: sigma(n,m) <= s iff H_1..H_m in F^{n x s} with H_i H_i^t = I,
 H_i H_j^t + H_j H_i^t = 0 (Hurwitz-Radon conditions). (2) t = rho(n)-1 ~ 2 log n anticommuting orthonormal antisymmetric integer
 matrices e_i (Clifford algebra); e_A = product over A, |A| = k: e_A, e_B anticommute iff |A cap B| != k mod 2 (Lemma 6/8).
 (3) (k,t)-parity representation xi: k-subsets -> F^s, <xi(A),xi(A)> = 1, <xi(A),xi(B)> = 0 when |A cap B| = k mod 2; built from a
 degree <= k/2 multilinear polynomial interpolating the parity condition on Hamming slices (polynomial method), s = C(t, <= k/2)
 (= orthonormal representation of a Kneser-type graph, Lovasz). (4) H_A = e_A (x) xi(A): m = C(t,k) matrices, s' = n s.
 Monotonicity sigma(n) <= sigma(n,m) for m >= n. Sect 5.2: tensor construction sigma(n^l,m^l) <= sigma(n,m)^l (4m)^{l/2}; from
 [8,8,8]: 8^{11l/6}.
Limits: Remark 11/15: a (k,t)-parity rep needs dim >= C(t/2,k/2). [P] With C(t,k) >= n and t ~ 2 log n this forces dim >= sqrt(n), so
 THIS framework can never beat n^{1.5}; Open Problem 1 (t-parity rep of dim 2^{(0.5+o(1))t}) would reach exactly 1.5.  N^{sqrt2}
 (App. B) is below even this conjectural target. Problem 2 (over R), Problem 3 (char 2: subquadratic nc circuit for ID_n?) open.
TR24-099 sect 6: nc circuit size of ID_n, S_{4,n}, perm_{4,n} <= O(n^{1.96}) (char != 2), black box: Lemma 22 (size(ID) = O(gamma),
 gamma = circuit size of the forms of an SOS identity over F(sqrt-1)), Lemma 24 (m = n^r copies of y computed by rectangular
 matrix multiplication: n^{omega(r)+delta}), Thm 25: size <= n^{(r+omega(r)+delta)/(r+1)}, r = 3, delta = 0.616 -> 1.954.
COROLLARY [P, derived here, novelty unchecked]: plugging App. B (sigma_{Q(i)} <= 2P_Q(N) <= N^{sqrt2+o(1)}, delta = 0.4142):
 r = 2, omega(1,1,2) < 3.2516: size^{nc}(ID_n) <= n^{1.889} (r = 3: 1.903), over any field of char 0; same for S_{4,n}, perm_{4,n}.
 (Small, side result; may be in Siemaszko/Debnath, NOT OBTAINED.)
Dictionary: Claude preprint App. G/H (p.45-47) already states P(N) between s and 4s and P(N) = O(N^1.62) from Hrubes Thm 12;
 it does NOT itself state N^{sqrt2}.

## 5. Zhang-Zhu arXiv:2605.00590v2 (11 pp.). READ IN FULL.
 [12,12,18] over every field with sqrt(-1) (Q(i), C), not over formally real fields (Yiu: no real [12,12,20]). Construction:
 12x12 vector matrix from 4x4 quaternionic sign templates (D,P,Y,Z,S) in a 3x3 block pattern + three vectors A,B,C with
 NON-PSD Gram matrix [[1,-1,1],[-1,1,1],[1,1,1]] realised by the bilinear dot product over Q(i). Coefficient criterion Prop 2.1
 (= HR conditions in vector form). Disproves Shapiro's field-independence conjecture (1984), Adem's value 26 for r = 11,12, signed
 formula 12#12 <= 17 < 21 <= 12*12. Relevance to us: none asymptotically; for products (beta, our P) positivity is never needed, so
 the "Gram realisation" freedom is automatically available to us; P(12) <= 4*18 is far worse than our exact 28 (G20b).

## 6. Other papers read
* Rossman-Zhu arXiv:2512.01227 (Dec 2025) READ PARTLY (sects 1, 2, 3, 6 full; 4-5 skimmed). KEY FIND (not in our journal before):
  PT-rank of an n^d x n^d matrix = least number of PT-basic summands (rank 1 after some partial transpose); Lemma 3.5:
  = min sum_kappa rank(N_kappa^{T kappa}). For d = 2, PT-rank(I_{n^2}) = min{rank A + rank (I-A)^Gamma} = our P(n) EXACTLY
  (our formulation (ii') of G20c). Lemma 3.13 (alg. closed): SoS <= PT-rank. Lemma 3.15 (char != 2): PT-rank <= 2^{d-1} SoS,
  proof: 2I = D + D^{T1}, D = L^T L of rank s => an s-square formula gives the feasible pair (s, s) (better than HWY's (2B,2B)).
  Prop 6.1-6.2, Thm 6.3: PT-rank(I_{n^d}) <= n^{0.81d} (two transposes) and <= n^{O(d^0.7)} (recursive nesting of SOS formulas,
  (1.62)^{log2 d}); with App. B in place of Hrubes this gives exponent d^{1/2}: the same 1/sqrt(D) law as our G18c heuristic [P].
  Lemma 3.8: rank(M^{T kappa}) <= n^{2 min(|kappa|, d-|kappa|)} rank M. Thm 1.5: PT-rank(M)/n^{d-log d+1} <= set-multilinear
  formula size; Cor 1.6: VNC1 != VNP from PT-rank n^{d-o(log d)}. Open (p.5): over GF(2), PT-rank(I_{n^d}) could be n^{d-o(log d)}.
  p.3: SoS_C(q_n) >= 2^{ceil log2 n} (Adams-Atiyah, as cited; not checked); "Omega(n log^5 n) suffices for VPnc != VNPnc" (cited to
  HWY; not in the ECCC version I read). Only the SUM is considered; no trade-off.
* Brand-Kaski-Wang arXiv:2607.02241 READ PARTLY (intro, related work): two-slice rank R_{d-2,1,1} lower-bounds commutative
  multiplicative complexity; not about Problem 12.
* Efremenko-Garg-Oliveira-Wigderson arXiv:1710.09502 READ PARTLY (sects 1.1-1.4): "rank methods" = ANY linear map L, mu = rank L(f);
  barrier c <= 2^d n^{floor(d/2)} for tensor rank. [P] Their proof scheme applied to our two families gives only the trivial
  barrier ~N^2 (the half-tensors have N^2 coefficients), so a barrier for ALL linear maps in our problem is OPEN; our Theorem B
  covers flattenings + local multiplications only. EGOW is the genre in which Theorem B lives.
* Lenzhen-Morier-Genoud-Ovsienko arXiv:1007.2337, Zhang-Huang arXiv:1705.04913: ABSTRACT + INTRO ONLY (Z_2^n-graded twisted
  (non-associative, cubic-form) algebras; Yuzvinsky-Lam-Smith [2n+2, 2^n - phi(n), 2^n]; doubling [r + rho(2^{m-1}), 2^m s, 2^m n]).
  All near Hurwitz-Radon (r ~ log n), asymptotically n^2/log n on the diagonal.
* Claude preprint (local): passages p.45-48 (App. G end, App. H) re-read.
NOT OBTAINED: Siemaszko (Zenodo 10.5281/zenodo.24024689: zenodo.org blocked by proxy; not on arXiv by author search), Debnath
 (GitHub puterfaction/hurwitz-sqrt2: access request denied), Shapiro's book, K.Y. Lam / Yiu / Adem originals, Dugger-Isaksen,
 Carmosino-Impagliazzo-Lovett-Mihajlin (not found on arXiv/ECCC quick search), HWY Canad. Math. Bull. 2013 (content = HWY sect 7).
 arXiv searches (site search, 2024-2026) for Hurwitz / sum-of-squares / partial transpose rank / Naslund / parity representation /
 Hrubes found nothing else relevant.

## 7. A. Dictionary (N = n; F field)
 (iv) Naslund P-rank = (v) Karam Q1.2 = (iii) HWY width w(ID_N) = (vii') RZ PT-rank(I_{N^2}) = P(N) := min{p+q} : all EQUAL
      exactly (definitions; HWY Lemma 3.5 p.14, RZ Lemma 3.5 p.10, Naslund p.8).  P_bal(N) := min{p:(p,p) feasible}: P/2 <= P_bal <= P.
 (ii) beta_F <= P <= 4 beta_F (HWY Prop 4.1, char != 2 for the right side); beta gives the pair (2beta, 2beta).
      sigma_F gives the pair (sigma, sigma) (RZ L.3.15, char != 2); sigma_F <= P when F alg. closed (RZ L.3.13); with sqrt(-1) in F:
      beta <= sigma <= 2 beta, so sigma <= 2P.  Over Q: P_Q <= 4 beta_Q; Q(i) costs nothing for P (G18g 1(a)).
 (i)  sigma_R >= sigma_C; gives (sigma_R, sigma_R) over R.  sigma_Z >= sigma_R.
 (iii') nc circuit size C(ID_N) >= Omega(P(N)) (HWY Prop 3.2); P(N) >= N^{1+eps} => 2^{Omega(r)} for LID_r, PERM (HWY 5.4, 1.7).
 (vi) (p,q): only the sum is visible in (i)-(v),(vii): every commutative identity gives a BALANCED pair, every pair gives
      beta <= p+q.  G20c's "commutative shadow sees only p+q" is CONFIRMED (HWY 4.1 proof; RZ 3.15) -- precise form: shadow
      data determine P up to factor 2-4, nothing about T(theta).
 Our line p q >= (2/9)N^2 translated: p+q >= 0.943 N < N: says NOTHING new in their language (weaker than trivial).
 Field: Hrubes 1.62 any field char != 2 (beta); App. B char 0 (Q); char 2: sigma trivial (=1), P_F2 open (our G9b: rho 0.811).

## 8. B. Records (Nov 2026; [R] page, (*) = depends on App. B + our G4a, external attribution unverified)
 sigma_R(n): <= O(n^2/log n) (HR; HWY p.3, Hrubes p.2); >= (2-o(1))n (HWY Thm 1.2 p.3). Unchanged; App. B does not apply over R.
 sigma_Z(n): Omega(n^{6/5}) (HWY Thm 1.10); <= O(n^2/log n).
 sigma_C(n), sigma over fields with sqrt-1: <= 2 N^{sqrt2+o(1)} (*) ; published O(n^1.62) (Hrubes Thm 1); >= ~ n..2n.
 beta_F, P(N), w(ID), PT-rank(I): <= N^{sqrt2+o(1)} over char 0 (*); O(N^1.62) any char != 2 (Hrubes); >= N (Karam p.4-5).
 nc circuit ID_n: <= n^{1.889} char 0 (*, sect 4 corollary); n^{1.96} char != 2 (Hrubes); >= Omega(n).
 trade-off curve: circle (App. B) / line (ours); not in literature.
 Is N^{sqrt2} the record? Yes for P, beta, sigma_C, PT-rank, as far as everything I read goes (RZ Dec 2025 and Zhang-Zhu Jun 2026
 still quote 1.62; the Claude preprint itself quotes 1.62). Hrubes's own framework is capped at 1.5 [P, Remark 15].

## 9. C. Technique inventory (tried by us? per PISTES 4-5, ARTICLES C)
Constructions:
 * Hurwitz-Radon / Clifford (O(n^2/log n), [rho(n),n,n]): stops at rho(n) ~ 2 log n operators. Tried (Pauli/Clifford diagonal, G20c).
 * Intercalate / Yuzvinsky monomial pairings / doubling / LMGO cubic-twisted algebras: near-HR families only. NOT TRIED (cubic,
   non-associative twisting) -- low expectation.
 * Hrubes Clifford e_A (x) parity rep: stops at the parity-rep dimension (>= sqrt n => 1.5 cap). Only its balanced point compared
   (G20c); its natural non-symmetrised (p,q) pieces NEVER computed. NOT TRIED in (p,q) form.
 * Hrubes tensor construction (11/6): = our Kronecker family (G18g), tried.
 * Zhang-Zhu Gram realisation over Q(i): irrelevant for products. NOT TRIED, low value.
 * App. B graded VA: tried extensively (G9a, G18h, G20c, G21a). Rectangular variant [H, here]: log P(N,M) ~ sqrt(log^2 N + log^2 M)
   ("Pythagorean"); splitting/monotonicity then gives exactly sqrt2 again -- no gain (rotation invariance of the circle).
 * RZ recursive nesting (multi-block): gives exponent d^{1/2} with App. B = the 1/sqrt D law (G18c). Tried in spirit (G20a).
Lower bounds:
 * Topological (Hopf-Stiefel, Adams-Atiyah K-theory, Lam, Yiu, Dugger-Isaksen): need SQUARES (nonsingular maps); give <= ~2n;
   do not apply to products/ranks. NOT TRIED; cannot give superlinear.
 * HWY integer/intercalate counting (6/5): needs monomial (one nonzero entry) structure. NOT TRIED -- candidate in a restricted
   "monomial decomposition" model, with a trade-off version.
 * HWY real Theta(k^2) explicit (positivity): not applicable over C/Q.
 * HWY central / block-central structure and lifting (Prop 5.3): NOT USED by us (relevant to multi-scale/recursion questions).
 * Karam fragmentation (qualitative, n^{3^{-32}}): not useful.
 * Naslund Lemma 11 / Karam k=1 annihilation: = our Theorem C mechanism (generalised by us).
 * EGOW rank-method barriers (any linear map): NOT TRIED; their scheme gives only ~N^2 here; general barrier open.
 * Low-degree equations / Young flattenings: tried (G21d).
 * char 2 (RZ open question; HWY transfer works in any field; Hrubes Problem 3): partly tried (G21b char p, P(3)=4 over F2).

## 10. D. Novelty
 Theorem C (p >= k(N-kq), q >= k(N-kp), pq >= (2/9)N^2, spectator form r >= kQ(n-kb)): NOT in Naslund, Karam, HWY, Hrubes x2,
 Rossman-Zhu, BKW, Zhang-Zhu, EGOW (all READ as above). Only k = 1 (p+q >= N) is in the literature (Karam p.4-5; implicit in HWY
 "M >= N" and Naslund). No paper read considers the two ranks separately at all. Unread possible sources: Shapiro's book,
 quantum-information literature on operator-Schmidt decompositions, Siemaszko/Debnath. Verdict: new as far as checked; the proof
 is an elementary extension of the k = 1 annihilation argument, so "folklore-level lemma" remains the honest description, but the
 trade-off QUESTION itself (T(theta)) is absent from the literature.
 Theorem B (corrected statement I): not found. Same genre as EGOW 2017 (barriers for rank methods), but EGOW's results do not
 imply it (different simple sets; their method gives only ~N^2 here). Should cite EGOW and position B as a restricted-class barrier.

## 11. E. Ideas
Construction below sqrt2:
 C1. Compute the exact (rr, rc) of the natural NON-symmetrised pieces of Hrubes's construction (R(i,j,i',j') = (H_j H~_{j'}^t)_{ii'},
     C = ID - R, plus App.-H "two orders" reading of e_A (x) xi(A)) for t = 5..9: does the Clifford/parity mechanism ("sign
     locality", |A cap B| parity) give a trade-off curve of a different SHAPE from the circle? First step: script for t <= 9.
 C2. Hybrid: App. B pole-order locality combined with Clifford parity signs (super VA whose modes carry e_A factors, labels in a
     parity rep). First step: abstract-theorem accounting (G21a A1-A4) with an extra Z_2^t grading; check if parity cancellations
     reduce the window dimension below dim V_{<=E}.
 C3. (cheap, likely negative) LMGO cubic-twisted group algebras as a non-associative analogue of the Pauli-diagonal family.
Lower bound above the line:
 L1. Characteristic 2: RZ flag GF(2); HWY transfer (P >= N^{1+eps} => exp. nc lower bounds) holds over any field; no
     symmetrisation upper bound exists there (Hrubes Problem 3). First step: exact P_F2(4), P_F2(5) and the F2 frontier (p,q)
     with our exhaustive/SAT tools; compare with best F2 construction (G9b).
 L2. Monomial-model trade-off: define decompositions whose row/column factors are signed partial permutations (covers HR, Pauli,
     Clifford, Hrubes e_A parts) and adapt HWY's 6/5 counting to get p q >= N^{2+eps} there. First step: definition + check which
     known constructions are monomial.
 L3. General rank-method barrier: decide whether ANY linear map L excludes a pair with pq >= N^2/4 (extend Theorem B beyond
     flattenings, EGOW-style). First step: N = 3, (p,q) = (2,2) (infeasible over F2 by G19b Thm 4): search for a general linear L
     witnessing it (LP over random/optimised L of moderate size); if none, conjecture B for all L.
