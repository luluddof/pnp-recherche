# G30-S1 notes (start 14:56 UTC). "SWAP with a side register" = one CNOT layer, unrestricted Pauli support.
Tags: [P] proof (design level, unreviewed), [E] exact, [N] numeric, [H] heuristic, [R] read.
## 14:56-15:05 reading: BRIEF, gen29 h1/r3/p3, gen28 c2/r2.
Precise problem (checked vs P3): n = 2k qubits, D = 2^n = N^2. M12|x,y> = |x+y,y>. Hadamard on x (Clifford) gives
 M12 ~ D_beta = diag((-1)^{xi.y}) = CZ^{(x)k}.  R(k,q) := min rank(D_beta - C), C in span of <= q Paulis on 2k qubits.
 Bell = C restricted to Z-type Paulis: rank = #{beta != g}, g q-Fourier-sparse => R(k,q) <= p_Bell(k,q).
 Cost = Pauli count |S| (P3: ring cost in [|S|^{3/4},|S|^2]; polynomial cost changes keep linear corners linear).
## 14:58-15:05 Task 1: P3's reduction statements, own code [E]
- sympl_check.py / sympl_check2.py: k=2 (Sp(8,2)): random A = g F_M g^-1 with B := A^-1 F12 of order 3 found 3 times in
  2130 tries; for hit 3 explicit symplectic X with B = X F_M X^-1 found (linear solve + sampling). k=4: A4 = A(+)A, B4 = B(+)B
  on site pairs {0,1},{2,3}: A4 B4 = F12(k=4), both order 3, symplectic. CONFIRMED (group level).
- clifford_check.py (operator level, k=2, 4 qubits): Clifford unitaries U_A, U_B built from A, B (twirl formula);
  M12^dag U_A U_B = phase * Pauli (residual 2e-15). So M12 = Pauli . (conj of M) . (conj of M): CONFIRMED.
  With (V1-E1)(V2-E2) = V1V2 - E1V2 - (V1-E1)E2 this gives R_{M12}(q^2) <= 2 R_M(q) (k even). 
## 15:05 LINE for the unrestricted problem [P]
Phi_a(Y) := f(a) X^a Y X^a Z^{Ja} (f = beta, B(z,a) = z.Ja the polar form) fixes D_f, preserves rank, shifts Pauli supports
by (0,Ja). From D_f = C + E = Phi_a(C) + Phi_a(E): C - Phi_a(C) has <= 2q Paulis and rank <= 2p. Operator uncertainty
rank(Y) |supp Y| >= D (Tr Y^dag Y <= rank ||Y||^2 <= rank |supp| sum|c|^2). If 4pq < D: C = Phi_a(C) for all a => supp C
invariant under all Z-shifts => C = 0 => p = D. Hence R(k,q) * q >= D/4 (q < D): same line as Bell (p_Bell q >= N^2/4).
## 15:07 numerics k=2 (sl.py: Pauli-sparse + low-rank alternating projections, hard threshold, polish on fixed support) [N]
 D=16. Bell (Z-type only) and Full (all 256 Paulis) give the SAME frontier: q=1,2,3: 6 | q=4,5: 4 | q=6..12: 1 | 16: 0.
 (Bell values agree with H1/R3 exact 6,6,6,4; Full lower bounds by the line only: 4,2,2,1,1,1.) No sign of side-register gain.
## 15:15 STRUCTURE: E = D_f - C = D_e + O,  e = f - g_0 (g_0 = Z-part of C, q_0-sparse), O = sum_{a in A, a!=0} X^a D_{s_a}.
 (a) [P] REDUCTION Full -> Bell with loss = #X-shifts:  the support graph of E has degree <= 2|A|, the diagonal is nonzero
     on supp e (>= p_Bell(q) points), an independent set I of those gives a diagonal principal minor =>
         p_Bell(q) / (2|A|+1)  <=  R(q)  <=  p_Bell(q).
     So a Full construction beating Bell by a factor N^{x} needs >= N^{x}/2 distinct X-shifts.  (Coset version: E[c+W,c+W]
     is diagonal if W meets A only in 0: loss 2^{codim W} ~ |A|, same.)  In N-units: delta_Bell <= delta_Full <= delta_Bell + eta.
 (b) [P] CODE-SPACE IDEAS ARE DEAD.
   (b1) K = V.Pi with 0 <= Pi <= 1 (Pi a projector, e.g. stabiliser code): Z-part of D_f Pi is D_{f h}, h = diag Pi in [0,1],
        so g := f h is q-sparse with |g| <= 1, and rank(V-K) = rank(1-Pi) >= Tr(1-Pi) = D - <f,g> >= D - D sqrt(q)/N
        (|f^| = 1/N, Cauchy-Schwarz).  =>  rank(V - V Pi) >= D (1 - sqrt(q/D)).  Useless at any q << D.
   (b2) Error supported on a stabiliser code (E = E Pi_S, code dim p = D/|S|): E's Pauli spectrum is S-periodic up to signs;
        V's spectrum is f^ on the Z-Lagrangian.  If S has a non-Z element, every t in L_Z with E_t = V_t forces E != 0 on
        t+S outside L_Z => >= (D-q)/2 extra Paulis.  If S is Z-type, periodicity f^(b+g) = +-f^(b) fails for half the b
        (Q*(b+g) - Q*(b) = Q*(g) + B*(b,g) depends on b) => >= D/2 mismatches.  So q >= ~D/2.  [P sketch]
   (Code spaces have dimension D/2^m: "agree with V on a code" means p >= D/2; "error on a code" is (b2).)
## 15:17 exact small-q results [P/E]
- R(k,1) exactly (r1exact.py; P^{-1}D_f = X^a D_f signed permutation, b=0 WLOG by twisted translation; real kappa):
  k=1..4: R(1) = 1, 6, 28, 120 = (D-N)/2 = p_Bell(1), attained at a=0; a != 0 gives 3D/4.  [E]
- SPAN LEMMA FOR FULL [P]: R(k,q) >= D/4 for q <= 2k (same as Bell's lemma). Proof: let Acomm = Paulis commuting with all of S
  (dim >= 2n - dim span S, n=2k). If T in Acomm has X-part a != 0: T C T^-1 = C, so E - T E T^-1 = D_f - D_{f(.+a)} has rank D/2,
  hence rank E >= D/4. Acomm has no such T only if Acomm is inside the Z-Lagrangian L_Z, i.e. dim Acomm <= n, with equality
  only if span S = L_Z^omega = L_Z (S Bell). If S has a Z-type element, a twisted translation puts 0 in S (dim span <= q-1).
- R(k,2) = p_Bell(2) [P]: S = {I, P} (after translation) with P non-Z and kappa_1 != 0: 2x2 blocks on {z,z+a} each rank >= 1
  => >= D/2; S with no Z-type: diagonal of E is f everywhere, blocks on cosets of <a1,a2> have independent sets of size 2
  (a1 != a2) or 1 per pair (a1 = a2) => >= D/2.
## 15:25 KEY IDENTIFICATION [P + E]: "SWAP with a side register" = PAULI-SPARSE APPROXIMATION OF SWAP (no side register)
F_SWAP (two k-qubit registers) and F_12 (CNOT layer) are both involutive unipotents of Sp(4k,2) with Lagrangian image of F-I and
alternating induced form of full rank on V/L (beta(x,y) = w(x,(F-I)y), beta(x,x) = 0) => conjugate. Checked: swapconj.py finds
symplectic X with X F_SWAP X^-1 = F_12 for k=1,2 (k=3 by direct sum). Cliffords with conjugate symplectic maps are related by
Clifford conjugation and a Pauli factor, both preserving R. Hence
   R(k,q) = R_SWAP(q) := min rank(SWAP_N - C),  C = combination of <= q Pauli products P1 (x) P2 (ARBITRARY P1, P2).
 Bell = the submodel P1 = P2 (commuting terms P(x)P).  Since osr(C) <= #Paulis:  F_general(N,q) <= R_SWAP(q) <= p_Bell(q).
 So the P3 "side register" is an artefact of the 3-register picture: the question is whether NON-COMMUTING (P1 != P2) Pauli
 terms help approximate SWAP.  A linear corner here would directly be a linear corner of the SWAP crossing problem (step 2
 for Claude's algorithm), not only for P_3.  Conversely no-corner-in-Pauli would NOT exclude non-Pauli (App-B-type) corners.
- swapnum.py [N]: SWAP_4 approximated by arbitrary Pauli products gives 6,6,6,4,4,1 for q=1..6: identical to CZ^{(x)2} (check
  of the identification) and to Bell.  (General osr model at N=4: 6,6,5,4,3,1 [G26-P4]: (5,3),(3,5) are NOT Pauli-sparse.)
## 15:23 numerics k=3, k=4 (heuristic, AP with hard threshold; Full = any Pauli, Bell = Z-type only) [N]
- k=3 (D=64): Full finds the Bell/Majorana points (8,8), (19,6) [the latter is pure Bell after re-polish: anal20.py];
  Full never below Bell: q=12,16 -> 8 (r=7 not reached, 15 starts), extension of Majorana (8,8) by t free Paulis:
  Full and Bell both stay at 8 for t <= 16.  (Heuristic is weak: it also misses Bell's (4,24).)
- k=4 (D=256), extension of Majorana (10,46) (extend4.py): Full q=18 -> 46, q=34 -> 38, q=46 -> <=11 (Majorana (46,10) region).
## 15:25 LEMMAS ON MECHANISMS [P]
- (L-bent) For f bent quadratic and any A: sparsity(f A) = |supp A~| >= D / sparsity(A), A~(x) = sum_s A^(s)(-1)^{Q*(s)} chi_s(x)
  ((fA)^(t) = f^(t) A~(J* t)).  Multiplying by f always costs >= D/sigma.
- (L-triangle) If E is block-diagonal over cosets of H (|H| >= 4) with rank-1 blocks u v^T and nonzero diagonal on a region R,
  then for a, a' in H: s_a(z) s_a'(z+a) s_{a+a'}(z+a+a') = e(z) e(z+a) e(z+a+a') on R.  With g_0 = 0 (e = f): the RHS is
  +-f chi (bent) => sigma^3 >= D on R = everything, i.e. cost >= 3 D^{1/3} at rank D/h: far above the line. With a Bell base
  g_0 (sigma_0 terms): RHS = f(chi + A) + B, sigma_A <= 3 sigma_0^2 => (L-bent) cost >= ~D^{1/5} if R = everything.
  BUT at the corner the region R (error set, density eps) is small and the identity only has to hold on R: no contradiction.
  So L-triangle kills "global" rank-1 blocking, not "fix the Bell errors by blocking" -- that mechanism (rank eps D / h at
  cost ~ h/eps) is exactly what a Full linear corner would need, and it requires error sets that are unions of H-cosets
  with sparse "transport" functions s_a = u(.+a) v on them; every attempt reintroduces f's flat spectrum (cf. P3 15:20).
- Rank-r blocks (r >= 2): (r+1)-minors give f.A = B with A, B products of r+1 shift functions: sigma >= D^{1/(2r+2)} only;
  compatible with the line for r ~ 1/eta. Not an obstruction.
- Symmetric Pauli constructions: Theorem S' (G28-C2) applies to ANY site-permutation-invariant K, Pauli or not => quadratic.
- k4_34.py [N]: the only Full "gain" seen (k=4, q=34: Full 38 vs Bell-heuristic 42) is a PURE BELL support (34 Z-type terms);
  on it rank 36 is feasible. => every Full solution found by the heuristics at k=2,3,4 collapses to Bell.
- R(2,5) = 4 [P]: cases as in the span lemma; the only new case (no Bell term, span S = 5) has <= 5 shifts, and 5 points never
  block all lines of PG(3,2) (blocking sets need 7), so some 2-dim W avoids A: the 4x4 principal block on a coset of W is
  diag(f) => rank >= 4.  No-Bell-element case of q=3, k=2: rank >= 6 [P] (hyperplane coset avoiding A, or A = H\0 with
  rank-1 blocks possible on <= 2 of the 4 cosets because chi_{Ja1} is non-constant on cosets).
- cornerfit.py: best known Full points (= Bell/Majorana) in N-units; corner points (q = 2k+2) k=4..7: (eta, delta) =
  (.83,.62), (.72,.43), (.64,.32), (.57,.25): power fit delta ~ 1.0 eta^2.6 (finite-k), asymptotically ln2 eta^2/ln^2(1/eta).
## 15:39 N=4 COMPLETE (Full = Bell) [P + N]
- orbits_k2q3.py: all 2,763,520 3-subsets of the 256 Paulis on 4 qubits, union-find under the symmetry group of D_f
  (twisted translations b -> b+e_i, CZ_ij conjugations, B-symplectic permutations z -> Gz): 73 orbits.
- rank_reps.py: per orbit rep, AP on fixed support (ranks 4..8, 8 starts) + Nelder-Mead minimisation of sigma_6 (12 starts):
  NO orbit reaches rank <= 5 (22 orbits reach 6, 1 reaches 8, 50 stay > 8; sigma_6 never -> 0).  => R(2,3) = 6 [N].
- So R(2,q) = 6,6,6,4,4,1(q=6..15),0(q=16) = Bell frontier exactly (q=1,2,4,5 proved, q=3 numeric-exhaustive, q>=6 trivial).
  k=3: R(1)=R(2)=28 (proved), R(q) >= 16 for q<=6 (span lemma); upper = Bell (4,24),(6,19),(8,8),(19,6),(28,1).
  k=4: R(1)=R(2)=120, R(q) >= 64 for q <= 8; upper = Majorana (10,46),(46,10) (Full heuristics: nothing better).
## inherited bounds (Full is a submodel of the general osr model, factors of norm 1) [P, by inclusion]
- Theorem C line, BN/norm bound p >= c N^2 ln N / sum|kappa| (a Full corner needs sum|kappa| >= N^delta),
  Theorem S' (site-symmetric => quadratic), R2 span corollary: if span(pi1(S) u pi2(S)) has dim r (one common Clifford U(x)U
  puts all factors on r qubits) then p >= (N^2 - N 2^r)/2 for ANY K_r (x) 1 => a Full construction must touch ~all qubits.
- Model chain at N=4 (costs: Pauli count for Bell/Full, osr for general): Bell = Full = 6,6,6,4,4,1 (q=1..6);
  general (G26-P4 exact) = 6,6,5,4,3,1: the general model's (5,3),(3,5) are not Pauli-sparse. Full is a STRICT submodel.
## 15:45 why "error blocking" fails for Majorana/symmetric bases [H]
 The only mechanism by which Full can beat Bell (rank eps D / h at extra cost ~ h/eps, i.e. at best the LINE) needs shift
 functions s_a supported (up to rank-small junk) on Err u (Err+a).  Sparse functions with small support are essentially
 affine-subspace indicators; Majorana/singlet error sets (tails |J - k| > w) contain affine subspaces only of density
 <= 2^{-(k+w)} << e^{-w^2/k} (a subcube in the tail must fix > k coordinates), so covering Err costs exponentially many shifts.
 Symmetric s_a (functions of J) must vanish on the window: degree >= d+1, cost ~ (k/d) x base cost, for a factor <= 2 in rank;
 raising the window degree is always better.  Rank-1 blocks of size >= 4 also face L-triangle on Err.
## VERDICT (G30-S1, 15:47)
1. Precise problem: R(k,q) = min rank(SWAP_N - C), C = real combination of <= q Pauli products P1 (x) P2 on 2k qubits
   (equivalently CZ^{(x)k} or the CNOT layer M12; P3's "side register" is the 3-register disguise of P1 != P2). [P+E]
   Chain: F_general(N,q) <= R(k,q) <= p_Bell(k,q); P3: P_3-Pauli linear <=> R linear; hence P_3-Pauli linear => general SWAP
   linear (step 2 directly). P_3 brings nothing beyond SWAP in the Pauli world.
2. Not decided. Proved: line R q >= D/4; R >= p_Bell/(2m+1) (m = # Bell-shifts of the non-Bell part; = P3's remark);
   span lemma R >= D/4 for q <= 2k; R(1) = R(2) = p_Bell; code-space ideas dead (K = V Pi: R >= D(1 - sqrt(q/D));
   error on a stabiliser code: q >= ~D/2); global rank-1 blocking costs >= D^{1/3}; symmetric => quadratic (S').
   Exact/numeric: N=4 Full = Bell for all q; k=3,4 heuristics never beat Bell (all hits collapse to Bell supports).
3. Intermediate statement: a Full construction beating Bell by N^x needs >= N^x/2 distinct shifts; the only mechanism
   (blocking Bell errors with sparse transport functions) can at best reach the line, and fails for all known error sets.
Files: notes.md, core.py, sl.py, target.py, sympl_check.py, sympl_check2.py, clifford_check.py, swapconj.py, swapnum.py,
 r1exact.py, majorana.py, extend.py, extend4.py, extract20.py, anal20.py, k4_34.py, sep3.py, bellfront3.py, cornerfit.py,
 orbits_k2q3.py, rank_reps.py (+ outputs *.out, *.npy, *.pkl). No process left running. Nothing written to the repository.
