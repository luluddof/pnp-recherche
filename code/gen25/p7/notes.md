# G25-P7 — test of the SPREAD CONJECTURE (gen25/b2 notes, idea 1d)
Tags: [E] computed here (floating point unless stated; ranks of App B exact via flint), [P] short proof here (unreviewed), [H] heuristic.
L = max(||Fl1 R||_op / sqrt q, ||Fl2 C||_op / sqrt p);  ineq = N^2 / (p||R|| + q||C||)  (<= 1 proved, 1 = tight).

## 0. Exact group and a G-INVARIANT lower bound for L   [P + E]
G = GL_N x GL_N acting by T -> (a(x)b) T (b^-1 (x) a^-1) on the 4-tensor (a on x,v; b on y,u). Fixes F, preserves rank Fl1
and rank Fl2 (C = sum G (x) G' -> sum (aGb^-1)(x)(bG'a^-1)); checked numerically (act() in spread.py).
KEY FACT: with S1 = Fl1(F), S2 = Fl2(F):  Fl1(R') S1 = (a(x)b) (Fl1(R) S1) (a(x)b)^-1  and  Fl2(C') S2^T is a similarity
of Fl2(C) S2^T (by a (x) b^-T). So rho(Fl1 R S1) and rho(Fl2 C S2^T) are G-INVARIANT and
   min over the G-orbit of L  >=  Linv := max(rho(Fl1 R S1)/sqrt q, rho(Fl2 C S2^T)/sqrt p).
(The proved inequality is the trace of this: |<F,R>| = |tr(Fl1 R S1)| <= p rho(Fl1 R S1).)  [checked numerically: invariant]

## 1. Table (spread.py, load.py, table1.py -> table1.out; appb2.py/appb1.py -> appb_*.json; jh.py -> jh.out)  [E]
 family           N   (p,q)        ||R||op  ||C||op  L      Linv   ineq
 g19 (5,3)@4      4   (5,3)        3        3        1.73   1.73   0.67
 g19 (7,4)@5      5   (7,4)        3        3.16     1.50   1.50   0.74
 g19 (6,6)@5      5   (6,6)        3.16     3        1.29   1.22   0.68
 g19/mono (7,7)@6 6   (7,7)        3        3        1.13   1.13   0.86
 g24g (12,3)@6    6   (12,3)       3        4        1.73   1.73   0.75
 mono (8,8)@7     7   (8,8)        4        4        1.41   1.41   0.77
 Q8 / CD (8,8)@8  8   (8,8)        4        4        1.41   1.41   1.00
 CD               16  (100,16)     8        8        2.00   2.00   0.28
 CD               32  (494,32)     16       16       2.83   2.83   0.12     (||R|| = N/2 = q: L = sqrt(N)/2)
 Hrubes k=2       8/16/32 (56,56)/(128,128)/(288,288): ||R|| = ||C|| = 1.83/2.54/3.67, L = 0.24/0.22/0.22 = Linv, ineq .31/.39/.48
 JH (all splits)  8..2048: L in [0.06, 0.83] always; e.g. N=256 (4096,4096) ||R||=22.8 L=0.36; N=2048 (47104,516096) L=0.58,
                  endpoints (N, ~N^2): ||R|| = N/16 = 128 at N=2048, L = 0.83; ineq falls to 0.03-0.06 at N = 2048.
 JH: ||Fl1 R|| = ||F1||, ||Fl2 C|| = ||F2|| (N x N Johnson matrices), Linv = L: checked on full tensors N = 8,16,32 (jh_check.out).
 App B (h = 2..7, N = p(h) = 2,3,5,7,11,15), threshold with p ~ q (appb_summary.out):
   N=2 (4,4):    ||R||=22.3   L=11.2     Linv=10.4      ineq 2.5e-2
   N=3 (9,9):    ||R||=617    L=206      Linv=97.6      ineq 9.0e-4
   N=5 (25,25):  ||R||=2.9e4  L=5.9e3    Linv=1.6e3     ineq 2.0e-5
   N=7 (49,48):  ||R||=1.9e6  L=2.8e5    Linv=2.2e4     ineq 3.0e-7
   N=11 (103,106): ||R||=7.6e7 L=7.4e6   Linv=5.3e5     ineq 8.9e-9
   N=15 (225,225): ||R||=9.4e9 L=6.3e8   Linv=6.8e6     ineq 6.1e-11
   local log-log slope of Linv (balanced): 5.5, 5.4, 7.9, 7.0, 8.2 (rising: superpolynomial, consistent with coefficients
   ~ exp(h log h), N = p(h) => log L ~ (log N)^2 loglog N [H]).  Even the CHEAPEST threshold (single piece, p = 1) has
   Linv = 2.7, 9, 34, 125, 461, 1720 (slopes 2.9, 2.6, 3.9, 2.9, 4.2).  max over thresholds 1.8e8 at N = 15.
   (h = 7 needed integer accumulation: float residual of sum A_s exceeded 1e-9.)
 CAVEAT: every App B point at N <= 15 is in the trivial regime (max(p,q) ~ N^2), since App B is capped at N^2 up to N ~ 256+.

## 2. Minimisation (spread.minimise: Riemannian descent on (a,b) via expm, Schatten-24 + soft-max surrogate, analytic
## gradient checked by finite differences; plus nz product rank-one moves Z = al(x)be(x)ga(x)de, (p,q) -> (p+nz,q+nz))  [E]
 mintest.out / minappb.out:
   small exact + Q8 + CD8: min over G = L0 = Linv (already G-optimal; (6,6)@5: 1.29 -> 1.24, Linv 1.22).
   with 1 / 2 moves: 1.73 -> 1.50 -> 1.34 etc.: the decrease is EXACTLY sqrt(q/(q+nz)) (3/sqrt4, 3/sqrt5): the op norms do
   not move; only the counts grow.  No norm reduction from mass moves found.
   App B: N=2: 7.5 -> 6.45 (= Linv);  N=3: 75 -> 73 (Linv 65);  N=5: 1784 -> 486 (Linv 390);  N=7: 7.4e4 -> 1.07e4 (Linv 6.0e3);
   +1 move: 5.8, 68, 466, 1.05e4.  The G-descent gains a factor <= 7 and stalls above the invariant floor.

## 3. Reading / verdict
 * Where the decomposition is near-optimal and small (N <= 8, incl. the tight Q8), L in [1.13, 1.73] and is G-optimal; the
   monomial model satisfies spread by construction (||R|| = max class size <= sqrt(2q)).  Hrubes/JH satisfy spread (L < 1) up
   to N = 2048, but their points are far from the conjectured frontier (ineq 0.03-0.5), so they do not stress it.
 * CD: a whole G-orbit with L = sqrt(N)/2 (Linv = L exactly): spread FAILS for that decomposition by N^{1/2}, yet the conjecture
   survives there because (p,q) ~ (0.5 N^2, N) is dominated (up to O(1)) by the trivial point (N(N-1)/2, 1): R = F - I,
   C = I(x)I, ||R|| = 2, ||Fl2 C|| = N, L ~ 2.  => spread is NOT a property of decompositions, only (maybe) of points; the
   "one good representative" clause is essential, and it cannot be tested on one family.
 * App B: every G-equivalent decomposition has L >= Linv, growing like N^{5..8} (rising slope) on N <= 15; the proved inequality
   is met with slack 1e-11 (||R|| ~ 1e10 vs the required N^{0.586}/2 ~ 2): the mechanism is massive cancellation in tr(R S1)
   (eigenvalues of R S1 of size 1e6..1e8 summing to <F,R> <= N^2).  But at these N App B's points are trivial, so this refutes only
   "App B is spread up to G" (rigorously, modulo floating point), not the conjecture.
 * VERDICT: NOT REFUTED, NOT SUPPORTED where it matters.  No explicit family with L >> 1 for every decomposition of the same
   (p,q) exists in the data: wherever L >> 1 (CD, App B) the point is dominated by a trivial spread decomposition at the tested N.
   The only constructions reaching the circle near theta = sqrt2 (App B; JH cannot, by g24a (L2)) are maximally non-spread, so
   the conjecture there would require a new cancellation-free decomposition -- no evidence for it.  Combined with the known
   implication (balanced 4/3 => HWY non-commutative lower bounds) the prior stays low.
 * WEAK SPREAD ||R|| <= q^{1/2+d}, ||Fl2 C|| <= p^{1/2+d} [P, one line]: N^2 <= p q^{1/2+d} + q p^{1/2+d} gives at the balanced point
   theta >= 2/(3/2 + d) > 1 for every d < 1/2 -- still beyond the line, so still implies P(N) >= N^{1+eps}: exactly as hard.
   Only a version RESTRICTED to q <= N^beta (unbalanced; no HWY consequence) is "safe": it gives theta >= 2 - beta(1/2 + d)
   (beyond the line theta >= 2 - beta for d < 1/2).  Data in that regime: JH endpoints q = N: ||R|| <= 0.83 sqrt q; (12,3)@6:
   ||R|| = sqrt3 sqrt q; CD q = N: ||R|| = q (d = 1/2, no gain, but dominated point).  Nothing tests it at q = N^beta, beta < 1.
Files: spread.py (analysis, invariant bound, minimiser), load.py, table1.py/.out/.json, jh.py, jh.out/.json, jh_check.out,
 appb1.py, appb2.py, appb_h*_N*.json, appb_small.out, appb_big.out, appb_h7.out, appb_summary.out, mintest.py/.out,
 minappb.py/.out.  Copies (unchanged) of ev.py, cd.py, hr.py, cliff.py, johnson.py, appb_exact.py and three g22b shape files.
