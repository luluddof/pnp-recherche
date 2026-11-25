# G27-T1 notes (exact computation: k_max(h, D) of P1's degree-1 ansatz, D -> infinity)
Start 07:38 UTC. Scripts copied from gen26/p1 (+ b6 relax.py as relax_b6.py). New: dscan.py (exact rows c=n, no
projection, row-compression of [A|rhs] by exact rref, binary search of min W' per W).
## 1. Tools
- dscan.py: P1's F-route (all columns F^mu_e, exact rows c=n), exact rref row compression, binary search min W'. Reproduces P1 h=4 D=4.
- sscan.py: NEW, faster. Direct operator form (P1 direct_cols) reduced to mode spans S_r(D)=span{t^mu[r]:V_h->V_{h+r}},
  S'_r(D)=span{t^mu[r]:V_{h-r}->V_h}, |mu|<=D, split by parity (alpha->-alpha); prefix column spans in W and W'; exact mod P
  ranks (flint), no projection/sketch. D=inf: S_r = full parity Hom (density of modes).
  CROSS-CHECK: sscan == dscan for every W at h=4 D=4..11 and h=5 D=5..11 (identical min W'). D=inf == saturated D (h=4: D>=10,
  h=5: D>=13) == B6 Test A envelope (h=4: (4,5),(5,5),(6,0); h=5: (4,7),(6,6),(7,4)).  W'=-1 means kappa=0 (R alone).
## 2. Results (min W' / k per W; k = 4h-1-W-W')
h=4: D=4,5: k=3 (2 at W=5,6). D=6: k=3 central, 4 at W=0. D=7: k=4 (5 ends). D=8: 5-6. D=9: 5-7. D>=10 = inf:
     k(W=0..6) = 9,8,7,6,6,5,9  (min W' = 6,6,6,6,5,5,0).
h=5: D=5: 3 everywhere. D=6,7: 3 (4 at W=16 endpoint). D=8: 3 central (4 at W=0,15,16). D=9: 4 central (5,6 ends).
     D=10: 6 central (W=3..9). D=11: 7 central. D=12: (W=0..8) 11,10,9,8,8,7,7,8,11. D>=13 = inf: 11,10,9,8,8,7,7,8,12(W=8,W'=-1).
     Mode-span saturation of the relevant blocks at D=13 (h=5), D=10 (h=4).
## 3. Mode-span law (satw.py, h=4; exact mod P)  [E + short proof]
dim S_r(D) = p(D) exactly (r != 0) until near saturation; S_r saturates (= full Hom(V_src,V_tgt)) at D = src+tgt for r<0
(rho: V_h->V_{h+r}) and at D = src+tgt+3 for r>=0 (h=4). r=0 slightly smaller (zero modes).
PROOF of the upper bound p(D): (L_{-1}a)[r] = -(wt a - r) a[r] (a[r] = a_{wt a-1-r}), so a -> a[r] factors through
V_{<=D}/L_{-1}V_{<=D-1}, of dim 1 + p(D) - 1 (vacuum mode vanishes for r != 0). So the weight-<=D ansatz has only p(D)
independent operators per degree r, versus dim Hom = p(h)p(h+r): saturation needs p(D) >= p(h)p(h+r), i.e.
D >~ (sqrt h + sqrt(h+r))^2 ~ 4h (central r) asymptotically.  D = c h with c < 4 never saturates for large h.
## 4. Counting proxy (counting.py): trivial antidiagonal dimU+dimU' >= N^2 with dimU_W <= P(W):
balanced k_triv(h) = 6,6,8,7,8,7,7,7,7,6,6,5,4,3,2,1,1,0,-2 for h=3..21, -> -0.7h at h=40 (since W_bal ~ 4h).
Matches the D=inf envelope at h=4 (6), 5 (7-8), 6 (7, B6), 7 (8, B6).
h=6 (sscan, Wmax=W'max=18, so W=0,1 show '-' when W'>18 needed): D=6..11: k=3 (2 at W=9), = P1. D=12: 3 central (4 at W=1,18).
     D=13: 3 central (W=9 exception gone; ends 4-7). D=14: 5 central (W=5..12). D=15: 6 central. (D=16, inf: see below)
## 5. Exponent translation (expo.py)
If k = c h persisted: theta^2+T^2 = 4 - c, balanced 4rho = sqrt(2(4-c)) (c=1.14 -> 2.39, rho .598, alpha' .74; c=2 = line).
But the finite (theta,T) of the D=inf points (log dims / log N) are ~(1.75,1.75): theta^2+T^2 ~ 6.1, = the trivial point
p=q=N^2/2 (2 - log 2/log N = 1.74 at N=15). The large k there is the sqrt(W/h) asymptotic formula being wrong at small h.
## 6. VOA / analytic remarks
- D=inf: mode spans = full Hom (density; DLM A_{n,m}(V) picture), so the VOA structure survives only in the fixed b_v[r];
  the problem IS Test A = truncated two-channel OPE (s-channel weights <= W, t-channel <= W') reproducing the swap (crossing).
- Rank bounds available: Theorem C p q >= (2/9) N^2 -> only k <~ 2h (line). So no analytic proof of boundedness is possible
  with rank methods (Theorem B barrier); the evidence is computational + counting.
- Lemma 1-D (pairwise) with pole order H = h + wt a: larger D makes the pairwise bound WORSE (W+W' >= 3h + wt a - 1);
  any D-gain is a cross-relation; at h<=7 these appear exactly when the counting (trivial antidiagonal) allows them.
h=6 D=16: central (W=6..10) k=7 except W=10:6/k7; ends larger.  D=inf (full Hom): W=0..10 -> min W' 10,10,10,10,10,10,10,9,8,7,4
     = k 13,12,11,10,9,8,7,7,7,7,9 ; == B6 (7,9),(8,8),(9,7),(10,4).  So h=6: central k: 3 (D<=13), 5 (14), 6 (15), 7 (16 = inf).
## ADJ trick (kappa spans = adjoints of rho spans, Shapovalov form) validated: identical results h=4 (D=4..10,inf), h=5 (D=5..13,inf)
h=7 (Wmax=W'max=19, W=8..19): D=7,10,12,13,14: k=3 for all W=8..19 (no even-h-type exception at odd h).  [run_h7_a.log]
## hscan.py (incremental quotient scan, rho spans + adjoint) validated vs sscan at h=5 (D=9,11,inf, 7 W's)
h=7 D=15: k=3 for all W=8..19.  (run_h7_b.log)
h=7 (hscan, W=9..15, run_h7_c.log; D=16 scan lost to a timeout): D=17: k=5 (W=9..14), 6 at W=15. D=18: 7 (W=9..13).
     D=19: 7 (W=9..11), 8 W=12. D=20: 7-8. D=21: 8 (W=9..11). D=inf: W=9:10/k8, 10:9/k8, 11:7/k9, W>=12: R alone (W'=-1)
     == B6 (9,10),(10,9),(11,7).   So h=7 central k: 3 (D<=15), 5 (17), 7 (18), 8 (21, inf).
## 7. Jump location vs mode-span law
First central jump (k>3): h=4: D=7, h=5: D=9, h=6: D=14, h=7: D=17 (D=16 untested, D=15 -> 3).
Heuristic D_jump = min{D : p(D) >= p(h)^2} (mode span can fill Hom(V_h,V_h)): 9, 11, 14, 16, 19 (h=8), 22, 25, 28, 31 (h=12);
ratio D/h -> 4 (p(D) >= p(h)^2 <=> sqrt D >~ 2 sqrt h).  D=2h: p(2h)/p(h)^2 = .36 (h=10), .03 (h=30) -> below jump for h>=6
(confirmed: h=6 D=12 k=3, h=7 D=14 k=3).  D=3h: above until h~50 (p(3h)/p(h)^2 = 0.98 at h=50), below after.
D=4h: always above (ratio ~ sqrt3 h) but then k = D=inf value = counting k_triv(h) -> negative for h>=21.
## 8. STOPPED at coordinator deadline (~08:45 UTC). h=8 D=inf (hscan 8 inf 12,13,14 20) killed unfinished (see run_h8_inf.log).
Unfinished: h=8 finite D>10 and D=inf; h=7 D=16; second-prime confirmation of the jump points (all results here: one prime
16777213, exact rows, no sketch/projection; cross-checked dscan(F-route)=sscan=hscan and D=inf == B6 Test A at h=4..7).
## 9. VERDICT: gain from D is a saturation/counting effect. k stays 3 until p(D) ~ p(h)^2 (D_jump/h = 1.75,1.8,2.33,2.43 for
h=4..7, heuristic -> 4), then rises to the D=inf value = counting k_triv(h) (6,7-8,7,8 for h=4..7), whose finite (theta,T) are
the trivial point p=q~N^2/2 (theta^2+T^2 ~ 6, outside the circle). k_triv(h) decreases (<3 for h>=17, <0 for h>=21) -> bounded.
