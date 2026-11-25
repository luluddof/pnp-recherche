# G26-P1 notes (saved as I go)
## 0. Reduction (derived from paper App B (47),(48), Lemma 9 proof; B6 sect 2a)
Per x: B_v^- rho + kappa B_v^+ = [v=x] Id on V_h, for all b in V_h.
Ansatz (degree 1 in operators): rho = sum_{r<=r_R} c_r[r] (c_r in V_{<=D}), kappa = sum_{r>=r_C} d_r[r];
r_R = W-h, r_C = h-W'.  Constraint: sum_{mu,r} c_{r,mu} <b[-r] t^mu[r]> + d_{r,mu} <t^mu[r] b[-r]> = l(b) Id.
One-variable picture (omega1=1, omega2=x): f+ = sum_r <a[r]b[-r]> x^{-r-h}, f- = sum_r <b[-r]a[r]> x^{-r-h},
(1-x)^H f+- = F, H = h+|mu|, F supported on e in [-2h, h+|mu|].
=> O^-_{mu,r} sees F_e for e >= H-h-r, O^+_{mu,r} sees F_e for e <= -r-h (triangular, unit diagonal).
=> ansatz feasible iff  l (x) Id in span{F^mu_e : e not in gap_mu}, gap_mu = (W'-2h, |mu|+h-W),
   in s=-e: dropped s in [s0+h-|mu|, s0+k-1] where W=2h+s0-1, W'=2h-s0-k, k = 4h-1-W-W'.
F computable through intermediate weights <= 2h only.
## 1. Validation [E mod 16777213]
- F-route vs direct operator columns (O^-_{mu,r}, r<=W-h; O^+_{mu,r}, r>=h-W'): identical feasibility for all (W,W') with
  2h<=W+W'<=4h at h=4 (D=4 and D=6), h=3. Lemma 9(c) (sum_e F_e = (t^mu[-h]t^nu) Id) and e=0 two-route consistency hold.
- WARNING: sketching sigma' by c<n random vectors is NOT faithful (h=7: r_all 356 at c=8 vs 369 at c=15). Use c=n (exact)
  or a dense random row projection.  Scans with c<n (scan_h6/h7.txt) are superseded by scanX_*.txt.
- D (max weight of the states a in rho=sum a[r]) matters: D=2h at h=4 saturates (column span = all n^3 tensors) -> k=7.
  D -> infinity recovers the full system (modes of all states span Hom(V_h,V_{h+r})): this is B6's Test A, saturated at small h.
## 2. Results D=h, exact in sigma' (c=n), mod p
 h=4: k=3 (k=2 at W=5,6);  h=5: k=3 all W;  h=6: k=3 all W (k=2 at W=9). r_all = 62, 115, 228.
 Observation: the paper's fields are the Heisenberg VOA: a[r] = modes of Y(a,z), Q_s = locality (z1-z2)^H.
 h=7 (c=15 exact): k=3 for every W in [0,24] (W'>=0 cap beyond). r_all=369.
 h=8 (c=22 exact): k=3 for every W in [0,28] except W=13 (k=2). r_all=628. (scanX_h8.txt)
 even-h exceptions at W=2h-3 (h=4: W=5,6; h=6: W=9; h=8: W=13): window s in [-2,0] not droppable.
 h=5 weight tests: only |mu|=h monomials: k=2 (k=3 at W=5); weights {h-2..h}: k=3; D=h+2 (weights<=7): still k=3.
 h=9 (dense projection T=1500, r_all=974 < T): k=3 for W in [0,27] except W=16 (k=2). (scanX_h9.txt, scan killed at W=27 by timeout? see file)
## 3. Deficiency profile (number of functionals l NOT achievable) by k, D=h [E mod p]
 h=5: k<=3:0, k4:2, k5:3-4, k6:3-5 ... -> 7=p(5)
 h=6: k4:3, k5:6, k6:7-8, ... -> 11
 h=7: k4:3, k5:7, k6:10, k7:11-13 -> 15
 h=8: k4:4, k5:10, k6:14, k7:17, k8:18-19 -> 22
 achievable a(k,h)=p(h)-def:  h=8: 18,12,8,5 (k=4..7); h=7: 12,8,5; h=6: 8,5.  ~ p(h-k+3) (+0..3)
 lowfun.py: dim span{a[-h]|V_h : wt a <= w} = p(w) exactly (h<=10).  So for k>=4 the achievable functionals are essentially
 those of states of weight <= h-k+3 (B6 idea 4: pole order h+wt a), plus 0..3 extra. Missing ~ p(h)-p(h-k+3).
 x-basis trick (GL_N on x): bad functionals can be paid by X-slices, cost N each -> extra f*N ~ N^2 k/sqrt(h): useless.
 h=9 full scan: k=3 for all W in [0,32] except W=16 (k=2).
 h=10 (T=2600, r_all=1528): probes W=0,10,16,18..22,30,36: k=3 feasible, k=4 infeasible (deficiency 5; 9 at W=16);
   W=17 (=2h-3, even-h exception): k=3 infeasible, k=2 feasible.  (probe_h10.txt)
 Exact over Q (exactQ.py, integer ranks, no sketch): h=5 and h=6 points agree with mod-p (incl. h=6 W=9 exception).
 D-dependence: D=h+1, h+2 at h=6,7: deficiency profile identical to D=h -> k_max=3 (log_D_*.txt).
 Parity split (alpha->-alpha automorphism) validated at h=8 (r_all 343+285=628, same deficiencies). fpar.py.
 h=9 W=17 deficiency: k4:4 k5:12 k6:18 k7:22 k8:25 k9:27 -> achievable a = 26,18,12,8,5,3 (k=4..9).
 PATTERN: a(k,h) = g(h-k) depends only on h-k: g(1..5) = 5,8,12,18,26 (from h=6..9); vs p(h-k+3) = 5,7,11,15,22.
 Prediction for h=10 W=19: deficiency k4: 42-37?, k5: 42-26=16, k6: 42-18=24, k7: 42-12=30, k8: 42-8=34.
 h=10 W=19 deficiency: k4:5 k5:16 k6:24 k7:30 k8:34  == PREDICTION CONFIRMED exactly.
 h=8 D=9: same profile as D=8 (4,10,14,17,19).
 EMPIRICAL LAW (h=5..10, central W): #achievable l with gap k = min(p(h), g(h-k)),  g(m) = p(m+4) - floor((m+4)/2)
   (g(1..6)=5,8,12,18,26,37). k=4 deficiency = floor(h/2) (partitions into exactly 2 parts?).  k=5 deficiency ~ p(h)-p(h-1).
 With the x-basis trick (pay bad functionals by X-slices, N each): k=4 costs only floor(h/2)*N extra -> effective k_max = 4;
   k=5 costs ~ N^2/sqrt(h): useless.  Still constant.
 Prediction h=11: def k3:0 k4:5 k5:19 k6:30 k7:38.
 D-dependence (unsaturated): h=6 D=7..10 (r_all 303,429,560,745 of 1331), h=7 D=8,9, h=8 D=9,10: deficiency profile
   IDENTICAL to D=h. Change only at saturation (h=4, D>=7: column span = all n^3 tensors).
 h=11 running via fpar.py (parity sectors, 8 chunks each, T=2000/sector); chunk ~250 s.
 VALIDATION of the D->infinity limit: h=4, D=10,12,14 reproduces B6's full Test-A envelope exactly:
   min W' = 6 (W<=3), 5 (W=4), 5 (W=5), 0 (W=6)  [B6: (4,5),(5,5),(6,0)].  (D=8 not yet: W=6 -> 3.)
   h=5 D=10: not yet saturated (W=4 -> 9, B6: 7). D=13 running.
 OOM note: apply_mono lru_cache is the memory hog (h=12 b-ops > 5.5 GB). h=11 chunk 3 (pi=0) OOM-killed, rerun queued.
 h=5, D=13 reproduces B6's full envelope exactly: (4,7), (6,6), (7,4); also W<=3 -> 8, W=8 -> 0.  (D=10: (4,9),(6,7),(7,6).)
 h=6 D=11 (r_all 940/1331), D=12 (r_all 1206/1331): deficiency profile at W=11 still 3,6,8,8,9,11 = D=h profile.
 second prime 16777199 at h=8 (exact c=n): identical r_all=628 and verdicts (k=3 ok, k=4 no, W=13 k=2).
 h=11 (parity sectors, T=2000 each, r_all 1088+1168=2256 < T): W=21: k=2,3 feasible; deficiency k4:5 k5:19 k6:30 k7:38
   == PREDICTION (5,19,30,38) CONFIRMED EXACTLY (out of sample).  (def_h11.txt)
 h=11 probes: k=3 feasible at W=0,6,11,16,17,18,19,20,21,22,23,27,33,40; k=4 infeasible at all probed W (def 2+3=5).
 h=12 attempt: fpar12.py (float32 S, cache clearing, b-ops split by w). validated on h=8 sector 0 (r_all 343).
 h=12 ABORTED (03:15): CPU load avg ~17 from other jobs, ~35% CPU per process, chunk 2 OOM-killed; would need >2 h. Not reported.
 h=9 D=10 (r_all 1339): profile identical (4,12,18,22) to D=9.

## 4. FINAL SUMMARY (03:24 UTC)
k_max (D=h, best W): h=4..11 all = 3; k=4 infeasible at every tested W; central exceptions with k=2 only:
 h=4 W=5,6; h=6 W=9; h=8 W=13; h=9 W=16; h=10 W=17.  (h=3: k=5, tiny/saturated.)
Status: exact over Q (integer ranks) at h=5,6 points; mod 16777213 elsewhere (16777199 also at h=8); h<=8 exact in all
indices (no sketch); h=9,10 (fbig) and h=11 (fpar, parity sectors) via dense random row projection, r_all < T.
Law: #reachable l at gap k = min(p(h), p(h-k+4) - floor((h-k+4)/2)); fitted h=5..9, confirmed exactly h=10, h=11.
D-independence until total saturation: h=6 D<=12, h=7 D<=9, h=8 D<=10, h=9 D=10. D->inf = B6 Test A (checked h=4,5).
Verdict: k bounded (3; 4 with X-slice payment of floor(h/2) functionals) -> no change of theta within this ansatz.
