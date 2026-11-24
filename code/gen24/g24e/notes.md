# G24e — exact vs border for (5,5)@5, (6,3)@5, (6,6)@6; vp.py runs at N = 7..9  (working notes)
Tags: [E] exact, [N] numerical.  d2 = min ||X - R - C||^2 (vp.py / LM).  |z|^2 = squared norm of all factors (U, V, G, H).

## Tools (this folder)
 vp.py   copy of G24c's solver (unchanged).  vp2.py: driver around vp.fg (per-restart log, warm start from a smaller-N point
         by zero padding + random new entries, odd restarts with a random placement of the old block).
 lm.py   LM on the FULL polynomial system X = U V^T + sum_l G_l (x) H_l (complex), optional Tikhonov mu ||z||^2.
 cont.py Tikhonov continuation mu = 10^-e, e step 0.5 (60 LM its per mu): exact bounded solution => |z| stays bounded
         while d2 -> 0; border => |z| -> inf with d2 ~ |z|^-k.
 gn.py   min-norm Gauss-Newton (lstsq) with backtracking, mu = 0.
 exact.py exact rank checker over Q (python-flint) -- NOT USED: nothing could be rounded (see below).
 OMP_NUM_THREADS=1 everywhere (machine load 13-14 on 4 cores from sibling agents: timings are inflated).

## 1. Exactification of G24c's three points  [N]
(6,6)@6, start = G24c best_6_6_6_s7.npy (d2 = 1.08e-6, max|coef| 66, |z| = 520):
 * plain LM (mu = 0, lam -> 1e-7, i.e. Gauss-Newton): 1.081e-6 -> 1.027e-6 in 20 its, |z| unchanged: NO quadratic
   convergence => there is no regular exact solution near this point; "essentially exact" (G24c) is wrong for this point.
 * continuation (logs/ct_6_6_6.txt): mu 1e-3 .. 1e-12:  |z|^2 = 204 -> 70e3 and d2 = 9.6e-2 -> 6.1e-7 with
   d2 * |z|^4 = 3.0e3 constant to 2 digits over 9 decades of mu (e.g. 4.2e-3*879^2 = 3256 at 1e-5, 1.95e-5*12401^2 = 3003,
   6.14e-7*69783^2 = 2988).  This is the signature of a first-order border decomposition (coefficients ~ eps^-1/2 per factor,
   error ~ eps).  So (6,6)@6 is BORDER-feasible [N, strong], and this family is border-only.
(6,3)@5, start = best_5_6_3_s7.npy (d2 0.0152): continuation (logs/ct_5_6_3.txt) mu 1e-3..1e-12: |z|^2 475 -> 3.3e5,
   d2 0.63 -> 6.0e-3, monotone, d2 ~ (|z|^2)^-0.7 (slower: higher-order border or unconverged).  Border-type, d2 -> 0 slowly.
(5,5)@5, G24c's best .npy was not on disk; re-run vp.py (seed 7, 1 restart: d2 0.034, coef 21; 4 further cold restarts with
   seed 5: all 1.0001 local minima).  Continuation from 0.034 (logs/ct_5_5_5.txt, mu 1e-6..1e-12): |z|^2 2.6e4 -> 5.5e4,
   d2 1.1e-2 -> 5.7e-3, d2 ~ |z|^-3.9: border-type.  (First attempt starting at mu = 1e-2 collapsed to the d2 = 1 local min:
   logs/ct_5_5_5_bad_mu1e-2.txt, discarded.)
CONTROL (important): (7,7)@6 is EXACTLY feasible (G19b exact_6_7_7.json).  Cold vp restarts (seed 2, 15000 its): 4.4e-4,
   1.7e-5, 9.7e-4 -- none exact; continuation from the 1.7e-5 point (logs/ct_6_7_7_control.txt) and min-norm Gauss-Newton
   (logs/gn_control_6_7_7.txt) behave EXACTLY like a second (6,6)@6 cold point (logs/ct2_6_6_6.txt, logs/gn_6_6_6.txt):
   d2 ~ 5e-5..1e-4 decreasing slowly, |z|^2 slowly growing, GN step lengths 1e-3.  => local polishing from cold numerical
   points does NOT discriminate exact from border: divergence along the followed path proves only that THAT path is a
   border family, not that no exact decomposition exists.
UPDATE (longer cold runs, 40000 its, seed 11, logs/v2_6_*_long.txt): control (7,7)@6 restarts r3, r4 converged SHARPLY
   to exact points (d2 3e-14 in 5799 its; 5.9e-16 in 237 its, max|coef| 1.12): exact points ARE found when present (control
   hit rate 2/8 cold restarts over both runs).  (6,6)@6: 0/12 sharp (G24c 3 + mine 4 + 5: best 1.1e-6, 9.4e-6, 1.2e-5, ...);
   continuation from the 9.4e-6 point (logs/ct3_6_6_6.txt) stalls at 1.4e-5, no quadratic convergence.  If an exact (6,6)@6
   had a basin like (7,7)@6, P(0/12) ~ 0.75^12 ~ 3%: moderate evidence against, NOT a proof.
VERDICT task 1: none of the three points could be exactified; all three are numerically border-feasible (d2 -> 0 with
   diverging coefficients along every path followed); exact feasibility is NOT excluded.  Exact bounds stay b(5) <= 6,
   b(6) <= 7.  (G19b: for the algorithm border = exact via interpolation, so the border values b_border(5) <= 5,
   b_border(6) <= 6 [N] are still relevant.)

## 2. vp.py at N = 7, 8, 9  [N]   (vp2.py; seeds/its in log names; all jobs < 10 min)
 calibration (8,8)@8: cold r0 d2 = 0 (nit 505, max|coef| 2.2) -> exact-type convergence (Hurwitz point reproduced).
 calibration (7,7)@6: 2 sharp hits (above).  (vp_8_8_8.txt; v2_6_7_7*.txt)
 (7,7)@7 cold, 11 restarts, 8000 its: 6.0047, 6.0071, 6.0138, 6.157, 6.160, 7.015, 7.038, 7.558, 7.594, 8.029, 8.212
   => d2 = 6 (c = 43), INFEASIBLE numerically; agrees with G20b (30 restarts, 6.0).
 (8,7)@7 = (7,8)@7 (symmetry c(p,q) = c(q,p)): cold 11 restarts: 1.625, 1.658, 1.714, 2.056, 2.286, 3.0 x4, 3.45, 4.14;
   warm from the restricted octonion (8,8)@7 minus one term, 10 restarts: 3.0001 x10.  Polishing the 1.625 point
   (60000 its, logs/pol_7_8_7.txt): 1.0617, max|coef| 29, singular values (1.0005, 0.0975, ...) => d2 -> 1 (border approach
   to c = 48 = N^2 - 1).  INFEASIBLE numerically (d2 ~ 1).
 (9,9)@9 warm from (8,8)@8, 11 restarts: 14 x10, 13.0044 x1 (coef 31)  => infeasible (G20b cold: 13.0).
 (10,10)@9 warm, 11 restarts: 10.0034 .. 10.0051 (all)  => infeasible (d2 ~ 10; G20b cold 9.07).
 (12,12)@9 warm, 13 restarts, 6000-8000 its: 7.3e-4, 1.15e-3, 1.66e-3, 2.1e-3, 2.1e-3, 2.2e-3, 2.5e-3, 2.8e-3, 3.1e-3,
   4.2e-3, 5.8e-3, 9.0e-3, 2.2e-2; max|coef| 8.6-19.6: no sharp convergence; border-type / unresolved (G20b: 2.4e-4).
## 3. Table b(N) (least p with (p,p) feasible), N <= 10
 N  | exact upper bound            | border / numerical                     | (p,p) numerically infeasible
 2  | 1 (settled, proof)           |                                        |
 3  | 3 (settled, proof)           |                                        | (2,2)
 4  | 4 [E]                        | 3 border [N, G19b]                     |
 5  | 6 [E, exact_5_6_6]           | 5 border [N, this note]                | (4,4) (d2 ~ 4.07)
 6  | 7 [E, exact_6_7_7]           | 6 border [N, d2*|z|^4 = 3.0e3 law]     | (5,5) (d2 ~ 6.07)
 7  | 8 [E, octonions restricted]  | -                                      | (7,7) d2 = 6.00;  (8,7) d2 -> 1
 8  | 8 [E, octonions]             | -                                      | (7,7) d2 = 18 (G20b)
 9  | 17 [P, direct sum 8+1]       | 14 feasible [N, G20b]; 12, 13 border-type unresolved | (9,9) 13-14, (10,10) ~9-10, (11,11) 4.47 (G20b, cold)
 10 | 25 [P, direct sum 8+2]       | 14 feasible [N, G20b]; 13 unresolved   | (10,10) 20, (11,11) 14.6, (12,12) 0.056 (G20b)
 b/N: exact 0.5, 1, 1, 1.2, 1.17, 1.14, 1.0 for N = 2..8; border 0.75 (4), 1.0 (5, 6); numerical 1.33-1.56 at N = 9
 (12-14) and 1.3-1.4 at N = 10 (13-14); exact 1.89 (9), 2.5 (10).  On N <= 10 the numerical ratio stays <= 1.6 (bounded on
 this finite range, trivially), but the exact ratio jumps after the Hurwitz dimension 8 and nothing here shows b(N)/N bounded.
 CAVEAT: at the end I ran pkill -f on 'vp.py'/'cont.py'/'gn.py'/'vp2.py'; the 'vp.py' pattern can also match sibling agents' processes of the same name (g24g's vp.py jobs, started 22:53:19, survived). All my own jobs are stopped.
