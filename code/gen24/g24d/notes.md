# G24d — better exact encodings for the open small cases (monomial ST model; GF(2) N = 4)   (working notes)
Tags: [E] exact computation here, [P] short proof here (unreviewed).  Folder: scratchpad/gen24/g24d/.
Tools: python-sat 1.9.dev24 (CaDiCaL 1.9.5 incremental), pycryptosat 5.18; built in user space (ext/): CaDiCaL (git HEAD),
kissat (git HEAD), drat-trim (git HEAD).  core.py, ev.py, verify_gf.py are read-only copies of g23b/core.py, g22c/ev.py,
g21b/verify_gf.py (independent verifiers).

## 1. ST encoding (st.py)
 Variables x[a][b][k] (k < q: S-colour, else T), palette bit P[a][b], equality bit e(c1,c2) for every pair of cells in distinct
 rows and columns (e <-> same colour; 2n clauses).  Partial permutations: AMO per row/col per colour.  Intercalates: for each
 2x2, e(D) <-> e(A) and e(D) -> P(r1,c1) != P(r1,c2).  SIGN-FREE: the sign condition is a GF(2) system given the colouring;
 it is solved after each model (gf2_signs) and an inconsistent set of intercalates is blocked (CEGAR, incremental CaDiCaL).
 Sign-free UNSAT => model UNSAT.  (In every run so far the first model was signable, no CEGAR iteration needed.)
 Symmetry breaking  [P, soundness]:
  (B1) balanced (p = q): take the line (row or column; transpose) maximising max(#S,#T), swap palettes so S is the majority,
       make it row 0; s := #S(row 0) >= N/2.  Then EVERY row and column has #S <= s and #T <= s (totalizers).  Unbalanced:
       the line maximising #S, every line #S <= s.  One SAT instance per s.
  (B2) columns: row 0 = S_0..S_{s-1} | T_0..T_{N-s-1}; columns of each block sorted by column-#S (non-increasing; col 0 = max).
  (B3) rows 1..N-1 sorted by the colour index in column 0 (strict).  (B4) colours not in row 0: first-occurrence order per
       palette along (col 0, rows 1..) then rows 1.. row-major.  Canonical form reached by: transpose/swap, row 0, column
       permutation, relabel row-0 colours, relabel the new colours of column 0 as the smallest new indices, sort rows 1..,
       relabel the remaining new colours by first occurrence; each step preserves the previous ones (row perms do not change
       column counts; the relabelling of row-0 colours does not touch column 0 because of (I)).
  (I)  implied [P]: S_b (1 <= b < s) never occurs in column 0 (rows a, 0 / cols 0, b: antidiagonal mono S forces the diagonal
       (0,0),(a,b) mono T, but (0,0) = S_0).
  (K)  key-lemma caps: |S-class| <= g*(p), |T-class| <= g*(q) (sequential counters).
 Exact verification of every SAT array: core.predicate, core.identity_tensor_ok (identity on all E_uv), core.ranks (flint, both
 flattenings, R + C = X), and the g22c evaluator ev.p_of (rank of X - sum_S K (x) K^T).

## 2. Validation of st.py  [E]  (all on file in val/, drat/)
 Balanced ST values N = 2..8: (1,1)@2 SAT; (2,2)@3 UNSAT, (3,3)@3 SAT; (3,3)@4 UNSAT, (4,4)@4 SAT; (6,6)@5 UNSAT, (7,7)@5 SAT;
 (6,6)@6 UNSAT (0.2 s; was 65 s with mysat.py), (7,7)@6 SAT; (8,8)@7 SAT; (7,7)@8 UNSAT, (8,8)@8 SAT.  Every SAT array passes
 predicate + identity + flint ranks (= colour counts) + ev.p_of.
 N = 5 unbalanced frontier (q = 1..10): min p = 10,10,10,9,8,7,6,5,4,1 reproduced exactly (SAT at p, UNSAT at p-1, all s),
 = g23b's exhaustive/SAT frontier.
 Ablation (weak breaking: NOBAL NOCOLSORT NOCAPS, all s = 0..N): (6,6)@6 and (7,7)@7 still UNSAT (< 1 s).
 DRAT: (6,6)@6 and (7,7)@7 (all s) CaDiCaL proofs VERIFIED by drat-trim (drat/).  [DRAT certifies the CNF; soundness of the
 symmetry breaking is the [P] argument of section 1, independently supported by the ablation.]
 NEW: (7,7)@7 is UNSAT in the ST model (was undecided at 400 s)  =>  ST value at N = 7 is exactly 8.
 (10,10)@9 UNSAT (62 s, all s = 5..9)  =>  ST value at N = 9 is >= 11 (counting bound gives 10).

## 3. GF(2) encoding (gf2.py)  [E]
 Unknowns U (N^2 x p), V (p x N^2), G (q matrices on (x,u)), H (q matrices on (y,v)); 256 XORs of p+q AND-products
 (CryptoMiniSat native XOR + Gauss).  Symmetry breaking:
  * full mode: V and H in reduced row echelon form (zero rows last) - removes GL_p(F2) x GL_q(F2) exactly (RREF unique).
  * GNF=r mode [P]: some G_l != 0 (else R = X has rank N^2); make it G_1 and use the GL_N x GL_N symmetry
    (x,v by A,A^{-T}; y,u by B,B^{-T}: X invariant, G_1 -> A G_1 B^{-1}) to put G_1 = diag(I_r,0); then V in RREF
    (canonical form of the row space of the transformed R).  Cases r = 1..N cover everything.
 Validation: (1,1)@2 SAT, (2,1)@3 UNSAT, (3,1)@3 SAT, (2,2)@3 UNSAT (full and GNF=1,2), (4,4)@4 SAT (0.1 s; full and GNF=2);
 every SAT verified by verify_gf (own Gaussian elimination + re-expansion R + C = X, R = UV).
 Runs: full RREF mode (3,3)@4 and (5,2)@4: no answer in 600 s wall (machine shared, ~60% CPU).
 GNF + RREF(V): (3,3)@4: r = 1,2,3,4 all UNSAT (42 s, 34 s, 45 s, 57 s).   =>  (3,3) is NOT exactly feasible over F_2 at N = 4,
 so b_{GF(2)}(4) = 4 (exact, over F_2 only; extensions F_{2^k} not covered).
 (5,2)@4: r = 1..4 all UNSAT (2-4 s).  => (5,2) NOT exactly feasible over F_2 at N = 4.  (6,2)@4 SAT (rc(C) = 1, verified),
 so on the q = 2 row at N = 4 over F_2 the least p is 6 (Thm 4 gives only p >= 4).
 CERTIFICATES: the GNF CNFs (xor chains, gf2.py mode 'dimacs', files drat/gf2_433_gnf{1..4}.cnf, gf2_452_gnf{1..4}.cnf) were
 re-solved by kissat (independent of CryptoMiniSat) with DRAT proofs; drat-trim: s VERIFIED for all 8 (drat/gfdrat*.out).
 Soundness of the GNF/RREF reduction itself is the [P] argument above (not machine-checked).
 (5,3)@4 SAT in GNF r = 2 (1.4 s, verified rr 5, rc 3) - the reduction finds the known point.
 (4,3)@4: GNF r = 1 UNSAT (kissat; DRAT 170 MB, drat-trim s VERIFIED, 281 s); r = 2,3,4: see section 5 (status).

## 4. ST runs N >= 9  [E]   (runs/, hunt/; machine shared, load 8-20 on 4 cores: each job got 15-60% of a core)
 Monotonicity [P]: a sub-array of an ST array is an ST array (partial permutations, intercalates are inherited), so the
 balanced value b_ST(N) is non-decreasing in N: (p,p) infeasible at N => infeasible at every N' > N.
 N = 9:  (10,10) UNSAT (all s = 5..9; 42.5+15.8+2.6+0.7+0.1 s); CaDiCaL DRAT proofs of all 5 CNFs VERIFIED by drat-trim.
         (11,11): s = 7,8,9 UNSAT (206.7, 13.4, 2.6 s); s = 5, 6 TIMEOUT whole (600 s CaDiCaL; s = 6 also 590 s kissat).
           Cube split (cube/, env C0 = #S in column 0, CC = #S in another column; ranges complete by (B1),(B2)):
           s = 6: c0 = 3,5,6 UNSAT (8.8, 45.8, 2.1 s); c0 = 4: col1 = 3 UNSAT (234 s), col1 = 4 TIMEOUT 480 s.
           s = 5: c0 = 5 UNSAT (81.9 s); c0 = 4: col5 = 4 UNSAT (244 s), col5 = 5 TIMEOUT 480 s.
           => (11,11)@9: 2 leaf cubes of 10 remain open (TIMEOUT); not decided.
         (12,12): s = 8, 9 UNSAT (369.5, 24.4 s); s = 5,6,7 TIMEOUT 400 s.   (13,13) s = 5: TIMEOUT 400 s (kissat).
         SAT (14,14) (8x8 corner fixed + extension, ext.py, verified).  => 11 <= b_ST(9) <= 14.
 N = 10: >= 11 (monotone), SAT (15,15) (corner extension; (13,13),(14,14) UNSAT *with that corner*).  11 <= b_ST(10) <= 15.
 N = 11: >= 13 (counting), SAT (20,20) (corner extension).  13 <= b_ST(11) <= 20.
 N = 12: SAT (22,22) (corner extension, 27 s, verified; (18..21) with that corner: no answer in 60 s each).
         (18,18) TIMEOUT: s = 6,7,8,9 CaDiCaL 300 s each; s = 6 with class sizes forced >= 4 (SAT-hunting): CaDiCaL 600 s,
         kissat 590 s; Z_m-invariant ansatz (m = 2,3,6, ext.py ZSYM): no answer in ~540 s.  15 <= b_ST(12) <= 22.
 N = 16: (28,28) TIMEOUT: s = 8,9,10,11 CaDiCaL 300 s each; Z_4 ansatz no answer.  21 <= b_ST(16) <= 32 (8x8 tiling).
 STU at N = 9..12: not run (time went to the GF(2) cases and N = 9).

## 5. Table and interpretation
 b_ST(N) (least p with (p,p) in the ST model) | counting bound (q g*(p) + p g*(q) >= N^2) | N^2/8 | N^{4/3}/2 | N^{sqrt2}
  N=2: 1 | 1 | 0.5 | 1.3 | 2.7          N=7: 8 (NEW, (7,7) UNSAT + DRAT) | 7 | 6.1 | 6.7 | 15.7
  N=3: 3 | 3 | 1.1 | 2.2 | 4.7          N=8: 8 | 8 | 8.0 | 8.0 | 18.9
  N=4: 4 | 3 | 2.0 | 3.2 | 7.1          N=9: 11..14 (NEW lower 11, DRAT) | 10 | 10.1 | 9.4 | 22.4
  N=5: 7 | 5 | 3.1 | 4.3 | 9.7          N=10: 11..15 | 10 | 12.5 | 10.8 | 26.0
  N=6: 7 | 6 | 4.5 | 5.5 | 12.6         N=11: 13..20 | 13 | 15.1 | 12.2 | 29.7
                                        N=12: 15..22 ((18,18) TIMEOUT) | 15 | 18 | 13.7 | 33.6
                                        N=16: 21..32 ((28,28) TIMEOUT) | 21 | 32 | 20.2 | 50.5
 * The exact values exceed the counting bound by 0..2 for N <= 8 and by >= 1 at N = 9; equality only at N = 2, 3, 8.
 * N = 8 is exceptional: b/N^2 = 1/8 there (the Q8 array, every class of size g*(8) = 4), but b_ST(9) >= 11 > 81/8, so
   N = 9 does NOT reach the tiling ratio; b/N^2: .25 .33 .25 .28 .19 .16 .125 (>= .136 at N = 9).  Whether any N in
   {12, 16} beats 1/8 (i.e. (18,18)@12, (<=31,<=31)@16, the asked (28,28)@16) is still undecided.
 * Exponent trend: on the exact values N = 3..8 a log-log slope is ~1.0, from 8 to 9 it jumps (>= 2.7 locally); plateaus
   (7,7,8,8 at N = 5..8) and jumps are arithmetic effects (sporadic tight structures at N = 8), so no exponent (4/3, sqrt2, 2)
   can be read off at N <= 16: N^{4/3}/2 and the counting bound are within 1 of each other and of the data up to N = 8, and
   N^{sqrt2} is 2-3x above every value here; constants dominate.
 * GF(2), N = 4: (3,3) and (5,2) are NOT exactly feasible (DRAT-certified CNFs + [P] reduction); with (4,4), (5,3), (6,2)
   feasible.  So b_{F2}(4) = 4 = b_Q(4) and the q = 2 row needs p = 6 over F_2; no characteristic-2 gain at N = 4 (exact,
   field F_2 itself; F_4 etc. not covered; over Q (3,3)@4 is border-only).

## 6. Status / files
 (4,3)@4 over GF(2): r = 1 UNSAT (kissat + DRAT verified); r = 2,3,4: kissat TIMEOUT 900 s, CMS TIMEOUT 600 s (twice). Open.
 Files: st.py (ST encoder, CEGAR signs, verify), stk.py (same encoding via kissat/cadical binaries), ext.py (SAT-hunting:
 fixed corner / Z_m ansatz; never used for UNSAT claims), gf2.py (GF(2) encoder: RREF / GNF), drat/drat.sh, drat/gfdrat*.sh,
 drat/*.cnf + *.dt.txt (drat-trim logs), val/ (validation arrays sol_*.json), runs/log_* (ST runs), hunt/ (upper-bound arrays
 solext_*.json: 9 (14,14) in val/, 10 (15,15), 11 (20,20), 12 (22,22)), gf/ (GF(2) logs, dec_gf2_*.json verified decompositions).
 All background processes killed at the end.
