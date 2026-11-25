# G26-P4 working notes (MC model exact computations)
02:05 start. Known exact optima (g19b exact_4_5_3, 5_6_6, 5_7_4, 6_7_7; g24g exact_6_12_3) all checked: every G, G' is a
signed partial permutation => all in MC.  Convention K[x,y,u,v] = sum_l E_l[x,u] C_l[y,v], Swap[x,y,u,v] = [x=v][y=u].
Symmetry group of MC: (A,B) signed perms independently: E_l -> A E_l B^-1, C_l -> B C_l A^-1 (since (A(x)B)Swap(B^-1(x)A^-1)=Swap),
 plus term order, term sign, flip E<->C (conj by Swap), transpose all.
Method: GF(2) relaxation (rank_Q >= rank_GF2 for integer matrices) => GF(2)-UNSAT => MC-UNSAT. GF(2)-SAT supports ->
 sign search (GF(3) SAT with fixed support) -> exact Q rank (flint).  CEGAR: unsignable supports blocked.
02:20 Calibration (all exact, verified by verify_g19b.py = copy of g19b/verify.py: flint rank + explicit factorisation):
 (6,1)@4 SAT 1s; (5,1)@4 UNSAT; (5,3)@4 SAT (r=3); (4,4)@4 SAT; (4,3)@4 UNSAT (all r, CMS, 3s; r=3 kissat+DRAT VERIFIED);
 (3,4)@4 UNSAT (all r, 9s); (7,4)@5 SAT r=3 (CMS 48s); (6,6)@5 SAT r=3 cube 1 (C_1 = E_11), kissat 79s.
 Tools: mc.py (encoder + CMS + CEGAR + GF(3) sign search), cube.py (orbit cubes on C_1, kissat), orbits.py, dump.py, fromsol.py, ver.py
02:23 (12,3)@6 is in MC trivially: g24g exact_6_12_3.json has all E_l, C_l signed partial perms (ver.py: rr 12, rc 3, re-expansion OK).
 Its 3 terms are all rank 4 and all in the same C_1-orbit cube (r=4, cube 57: C_1 = diag(I_4,0) after normalising E_1).
 Implied constraint RSUM: sum_l |E_l||C_l| >= rank K >= N^2 - p (makes r=1 trivial; (11,3)@6 needs r>=3).
02:35 N=4 MC frontier (front4.txt; each SAT verified exactly by ver.py on spot checks): least p for q = 1..9:
 6, 6(thm), 5, 4, 3, 1, 1, 1, 1 ; (2,5) UNSAT, (3,4) UNSAT, (4,3) UNSAT; (1,6)@4 SAT exact (rr 1, rc 6 verified!).
 p = 0 needs q = 16 (osr(Swap) = 16).
 cegar.py (CEGAR with acyclic-set / minor lemmas, no UV) works but is slower than the UV encoding (N=4: 26 s vs 1 s): dropped.
 LINK (U = M at pivot columns) gives no speedup on kissat. Load on the machine ~14 (other agents): my effective CPU ~1 core.
02:45 [P] CORE LEMMA. Theta = K Swap, P_l(a,b) = (eps_l(b), gam_l(a)); dom P_l = cols(C_l) x cols(E_l) (as (a,b)), im P_l = rows(E_l) x rows(C_l).
 Fix(Theta) injects into Fix(Theta_{D'D'}), D' = dom ∩ im (w = Theta w => w supported in im, determined by w|dom, w|dom in span D').
 => rank(I - Theta) >= N^2 - |dom ∩ im|  (encoded as implied constraint CORE=1, plus RSUM).
 COROLLARY: if every E_l, C_l has rank <= r (r >= 2) then p >= N^2 - q r^2 + 1.  Proof: equality N^2 - q r^2 forces disjoint
 domains of size r^2 with dom = im = D' and Theta_{D'D'} = I; for c in dom P_l only P_l acts, so P_l fixes all of
 cols(C_l) x cols(E_l), i.e. eps_l(b1) = a1 = a2: impossible for r >= 2.  => (5,5)@5 branch r <= 2 is UNSAT (by proof).
 Z_3 cyclic ansatz (ans.py: E_l = A^l E_0 B^-l, C_l = B^l C_0 A^-l; (A,B) types up to conj and flip: (1,3),(1,33),(3,3),(3,33),(33,33))
 for (11,3)@6: all 5 UNSAT (ansatz only, <= 31 s each). (12,3)@6 known point has a 24-element support automorphism group
 (S_3 on terms, includes Z_3 of type (33,33)); ansatz calibration run on (12,3) pending.
02:50 [P] STAIRCASE LEMMA (MC form of g24g Lemma C without the min(b-1,.) cap): Z = union support of the E_l; for any
 column order beta, S = {(x,y): beta-pos(y) < t(x)}, t(x) = # leading beta-columns on which row x of Z vanishes, is acyclic
 in the cell digraph (potential beta-pos(second coord) strictly decreases along edges inside S; no loops) => rank >= sum_x t(x).
 Same for Z^T, and for C (flip) and C^T.  Checked: 400 random MC instances, 0 violations; tight on (6,6)@5, (7,4)@5 (bound = rank),
 11 vs 12 on (12,3)@6.  Encoded exactly (max over all beta) by a DP over column subsets (STAIR=1 in mc.build, ~O(2^N N p^2) clauses).
 Support-only (implied constraints only: RSUM + CORE + STAIR, no rank encoding; supp.py):
   (11,3)@6 r=3: UNSAT (1.7 s);  (5,5)@5 r=2: UNSAT (51 s)  [=> these branches are dead by lemmas alone]
03:00 DRAT (drat4.sh, pure encoding: no lemmas, only ORD + E_1/C_1 normal form + RREF): N=4 (4,3), (3,4), (2,5), (5,2), (5,1):
 all branches r=1..4 kissat UNSAT, drat-trim s VERIFIED (20/20).  [drat-trim prints \r before 's VERIFIED'.]
 PEEL (iterated support-level peeling bound) adds nothing in support-only tests.
 Full UV runs with STAIR: (11,3)@6 r=4,5,6 and (5,5)@5 r=3,4,5 kissat started 02:50-02:53 (590 s limit).
03:05 [P]+validated COMBINED BOUND (stair2.py; 600 random MC instances, 0 violations, strictly better than plain staircase
 in ~2/3): rank >= (N^2 - |D'|) + max_beta #{(x,y) in D' : pos(y) < t_beta(x)}, 4 variants (E,C),(E^T,C^T),(C,E),(C^T,E^T).
 (proof: rank(I-Theta) >= N^2 - |D'| + rank(I - Theta_{D'D'}) and the staircase set intersected with D' is acyclic in Theta_{D'D'}.)
 Tight on (5,3)@4, (3,5)@4, (6,6)@5, (7,4)@5; 11 on (12,3)@6.  Encoded as STAIR=2 (default now).
 Support-only with STAIR=2: (11,3)@6 r=4,5,6 and (5,5)@5 r=3,4,5 still satisfiable (no new dead branch).
 kissat full-branch (old encoding, no lemmas) (12,3)@6 r=4: SAT within 580 s (cnf/n6_12_3_r4.out); (5,5)@5 r=2 kissat UNSAT.
 03:06 launched STAIR=2 full branches (11,3)@6 r=4,5,6 and (5,5)@5 r=3,4,5 (cnf/t*.cnf, 590 s each).
03:10 (12,3)@6 rediscovered FROM SCRATCH by the GF(2) SAT (kissat, full r=4 branch, old encoding, < 580 s wall) + sign local
 search (1.9 s; signls.py; now stage A of mc.sign_search): mc_N6_12_3_sat.json, verified exactly (rr 12, rc 3); same type as
 g24g's point (all ranks 4, 24 support automorphisms, 12 diagonal +1 cells).
 Machine: load ~13 from other agents; my kissat jobs get ~8% CPU each (15 s CPU per 3 min wall) -> 590 s wall = ~50 s CPU.
03:23 Involutive-automorphism ansatz (AUT="A;B;sigma" in mc.build, support equalities E_sig(l)[Ax][Bu]=E_l[x][u] etc.; all lemmas on):
 calibration: (7,4)@5 with its automorphism type: SAT 6.9 s; (12,3)@6 with A,B of type 2^2, sigma=(12): SAT 91 s (mc_N6_12_3.json,
 verified).  Z_5 cyclic ansatz (5,5)@5: 3 types UNSAT (ansatz only). Involution ansatz runs: ansinv_555.txt, ansinv_6113.txt.
03:28 Flip-type order-2 ansatz (AUTF=1; WLOG A=B=I when B=A^-1): E_sig(l) = C_l, C_sig(l) = E_l.  sigma = id (E_l = C_l
 supports) finds (6,6)@5 in 28.6 s.  (9,3)@5: all 6 involution ansatz classes UNSAT (ansatz only).
 Hunting NEW exact points (general problem has only p <= 10 at (N,q)=(5,3), <= 12 at (6,4), <= 21 at (7,3)):
 flip ansatz sweeps (9,3)@5, (11,4)@6, (20,3)@7 queued; full-branch kissat (9,3)@5 r=3,4 running.
03:31 (9,3)@5 in MC: r=1,2 dead by lemmas (support-only UNSAT); r=3 kissat UNSAT (full encoding with lemmas, ~5 min wall,
 cnf/t5_9_3_r3.cnf); r=4, r=5 running.  Ansatz TRANS (C_l = E_l^T supports; observed in g24g (12,3)@6, g19b (7,7)@6,
 our (6,6)@5) implemented (TRANS=1) - (6,6)@5 r=3 not found within 100 s CPU limit.
03:43 (9,3)@5: r=4 whole branch TIMEOUT (590 s wall); r=5 whole TIMEOUT.  NONZERO=1 (all terms nonzero; a zero term reduces to
 q-1 = 2 terms, UNSAT by g24g Theorem 1 since 9 < 10) -> reruns cnf/u5_9_3_r4/r5.  Cubes r=4 (73 C_1-orbits): cube 0 (C_1=0)
 = zero term (Theorem 1), cubes 1,2 UNSAT (213 s, 165 s); cubechain.sh runs cubes 3..72 (NONZERO) in <600 s chunks.
 Flip ansatz: (5,5)@5 sigma=id,(01) UNSAT; (9,3)@5 sigma=id,(01) UNSAT (ansatz only).
03:50 cubes (9,3)@5 r=4: 3,4 UNSAT (~190 s each, 2 workers). cubechain r=4 (3..72) and r=5 (1..35) running; whole-branch
 NONZERO r=4 TIMEOUT 590 s. kissat DRAT proof generation for r=3 running (drat/t5_9_3_r3.drat).
03:57 (9,3)@5 r=3: DRAT proof (kissat, 149 MB) VERIFIED by drat-trim (173 s).  Cubes r=4: 16/72 UNSAT so far, r=5: 8/35.
04:03 (9,3)@5 r=5: all 35 nonzero cubes (1..35) kissat UNSAT (cubechain_5_9_3_r5.txt). r=4: 33/72 so far; second chain on 56..72.
04:08 *** (9,3)@5 is INFEASIBLE in MC ***: r=1,2 support-only (lemmas) UNSAT; r=3 kissat UNSAT + DRAT VERIFIED; r=4 all 72
 nonzero C_1-cubes UNSAT; r=5 all 35 nonzero cubes UNSAT; zero-term cases = (9,2)@5, UNSAT by g24g Theorem 1.  => least p
 for q=3 at N=5 in MC is 10 (= N(N-1)/2; the general problem's border value 8 is not reachable exactly inside MC).
 Next: (6,4)@5 (NONZERO valid since (6,3)@5 is UNSAT by the above): cubes r=3,4,5 (r<=2 dead by lemmas).
04:14 (5,5)@5 r=3 cube 5 alone (new encoding, NONZERO): TIMEOUT 560 s -> (5,5)@5 out of reach. (6,4)@5 cubes running:
 chains c (r=3, 1..45), d (r=4, 1..40), f (r=4, 41..72), e (r=5, 1..35); logs cubechain_5_6_4_r*_?.txt; coverage: cover2.py 5 6 4.
04:20 sanity: full cube pipeline (lemmas + NONZERO + ORD) finds (10,3)@5 with 3 nonzero terms in cube 29 (C_1 = I) in 4.5 s (sanity/).
 (6,4)@5 cubes: first TIMEOUTs at 290 s (r=4 cube 4, r=5 cube 2) -> rerun needed with longer limit; completion unlikely.
04:25 TOP1=1 (valid symmetry breaking: term 1 = the term with |E_l| = r maximising |C_l|, i.e. |E_l| = r => |C_l| <= |C_1|):
 (6,4)@5 r=4 cube 4 UNSAT in 113 s (was TIMEOUT 290 s).  Restarted all remaining (6,4)@5 cubes (125) with TOP1 (listchain.sh,
 logs cubechain_5_6_4_r*_top_*.txt).  Coverage before restart: r=3 4/45, r=4 23/72, r=5 1/35.
04:51 (6,4)@5 coverage: r=3 29/45, r=4 56/72, r=5 22/35 UNSAT (no SAT cube). Remaining cubes restarted with TL 540 s,
 5 workers (listchain2.sh; logs cubechain_5_6_4_r*_fin_a.txt). Oversubscription caused 280 s timeouts.
05:03 *** (6,4)@5 is INFEASIBLE in MC ***: r<=2 lemmas; all nonzero cubes r=3 (45), r=4 (72), r=5 (35) kissat UNSAT
 (TOP1/NONZERO/ORD); zero-term configs are (6,3)@5 configs, excluded by (9,3)@5. => least p at (N,q)=(5,4) in MC = 7.
05:08 (6,5)@5: r<=2 dead by lemmas (support-only); 4 test cubes (r=3: 19,40; r=5: 1,30) unfinished after 4 min -> not attempted
 further (155 cubes). All my processes killed.
FINAL TABLE (MC):
 N=4 least p, q=1..9: 6,6,5,4,3,1,1,1,1 (exact; UNSAT sides DRAT-verified for (5,1),(5,2),(4,3),(3,4),(2,5)).
 N=5: q=3 -> 10 ((9,3) UNSAT); q=4 -> 7 ((6,4) UNSAT, (7,4) SAT); q=6 -> <= 6 ((6,6) SAT); (5,5),(6,5): TIMEOUT (r<=2 dead).
 N=6: (12,3) SAT; (11,3) TIMEOUT (r<=3 dead by lemmas; ansatz classes UNSAT).
