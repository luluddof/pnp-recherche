# G26-P14 adversarial review of G26-P13 (notes as I go)
- 04:27 start. Read P13 notes, P10 notes, paper Def 1-4, Prop 3 (lines 990-1062).
- Analytic pass of Step E proof: inclusion-exclusion m(a,b) correct (count (K-2)(K-3)), Jensen on row/col sums, CS -> eps_K correct;
  martingale differences w.r.t. common filtration -> Var(Ye-Ye') = sum E(D-D')^2 correct. Equal means by symmetry. Proof OK.
- Step C algebra: [x=v][y=u] = [x'=v'][y=u][x''=v'']; x'' at (a,j), v'' at (a',j'): separated by both the row cut {x,y}|{u,v}
  and the column cut {x,u}|{y,v} of (33) -> obs (a),(b) give +(l1-l2) bits to both directions. Psi_{a,a'} in Prop 3 is a product of
  rectangle-local factors, so inserting [x''=v''] is legitimate. Cost arithmetic E_w A_q + S_w Delta_q = E_w T_q + (S_w-E_w)Delta_q OK.
- 04:3x enum14.py/.out: EXACT enumeration over all (K!)^t stratified groupings (K=4..6, t=2..3; networks unif/bern/rank1/star/
  block/sparse). Var(Y12-Y34)/(VarY12+VarY34) max = 2.000 (K=4), 1.667 (K=5), 1.500 (K=6): = 1 + 2/(K-2), attained by star
  networks; always <= 1+eps_K (8, 4, 2.83). climb14.out (hill-climb on c): K=5 1.6649, K=6 1.4072, K=8(t=1) 1.3331 -> the
  conjectured exact worst case 1+2/(K-2) ~ 1+2/K; eps_K ~ 4/K is the right order (lemma is valid, factor-2 loose).
  Means equal (|diff| ~1e-15). => Step E CONFIRMED. Note: the ratio genuinely exceeds 1 (by ~2/K), so (1+eps) is needed.
- qcheck14.py/.out: level q. Exact model chord sup_d (cost(d)-cost(0))/(d cost(0)) with cost(d) = exact max-min (closed-form
  crossing) vs P13's closed form (S_w-E_w+O_w+ent)/e: 400 random (ea,eb,al,be,ent): max(q_exact - q_P13) = 3e-10, mean gap 0
  -> P13's level q equals the right-derivative at 0 (envelope theorem) and the chord sup is attained there. CONFIRMED.
  (First attempt with a 40001-point theta grid gave spurious q_exact up to 2.9: grid error of a kinked max-min at d ~ 1e-6.)
- basecheck14.py/.out: row-sweep base q, fine delta grid 1e-10..1e3 (3000 pts): sup chord <= P13's q in all 24 cases
  (P13 conservative by ~12% from its geometric-grid chord bound); limit sqrt2-1 = 0.4142 for large g. CONFIRMED (conservative).
- sweep14.py/.out: plain clique sweep with per-edge widths; P13 asserts q = 1 ("cut_j <= s c' + beta sum|lam-c|") without proof
  (false for a fixed ordering; needs the best ordering). Brute force over all orderings + hill-climb, K=4..7: worst
  q_needed = K(K-1)/(2s) - 1 (star) = 0.20/0.25/0.36/0.40 < 1. Plausible, unproved for general K. DP sensitivity: Q0=3 run.
- Paper Def 7 / Lemma 8 (lines 1515-1560, 1626-1650): c_row, c_col <= 2 lambda assumption; split rectangle c_row = c_row_i(l2) +
  (l1-l2) <= l1 + l2 <= 2 max: compatible. Lemma 8 (ii) needs c_ab <= P/g: Step B's M(a,a') bound. OK.
- griderr14.py: base e (P10/P13 20001-pt theta grid + 1e-7) under-estimates the exact max-min by <= 3.8e-5 relative (kink).
  Inherited from dp10; effect on c1 <= ~1e-4. Minor.
- Step B analysis: g^2 E M(a,a') = T + Dbar (1+O(1/g)) verified (P(edge not in same row) = h(g-1)/(K'-1)). Concentration: M depends on
  the 2h nodes in rows a,a' (sampling w/o replacement), bounded differences 2h lam_cap per draw -> dev O(h^{3/2} sqrt(log g)) lam_cap;
  relative eps ~ sqrt(log g/h) * lam_cap/lam_bar with lam_cap/lam_bar = N_par/T ~ 1/cos(th) ~ 1.2-1.5. Sound ASYMPTOTICALLY; but
  eps_pl is dropped in the numerics although the optimal trees contain levels/bases with g,h ~ 2^7..2^10 where sqrt(ln g/h) ~ 0.2.
  Only overhead terms (O_w + ent) scale with N_cap, so the effect is second order -> sensitivity runs EPC=1,3 (levels only).
- dp14.py = copy of dp13.py + env options Q0 (q of plain sweeps), EPC (multiply al, be, ent of each level by 1+EPC sqrt(ln(max g,h)/min g,h)).
- 04:48 dp14 (fixed epl overflow bug of my first attempt) results, a=2 [1000..16000 windows]:
  base (= P13 splitiii): 3.158 / 3.164 / 3.174 / 3.166 ; a=4 ... 4.477 -> REPRODUCES P13 exactly.
  Q0=3 (plain sweeps charged q=3): identical c1 (only solvers with x in [30,100) change, by <= 3.7%) -> plain-sweep q irrelevant.
  EPC=1 / EPC=3 (level overheads x (1+C sqrt(ln max(g,h)/min(g,h)))): identical c1 to 3 decimals (rel de <= 4e-5 for x >= 1000;
  <= 1% for x in [30,100)). Level grids in the relevant trees are large, so eps_pl at levels is negligible.
- thr14_p13sol.out: mythr13 on P13's sol json (III=1, top S2q) reproduces 1352 / 2683 / 9269 (r=4.26) and <=50/79.72/256.4 (r=1).
- 04:50 launched EPB (base caps also inflated) runs and thresholds on my own sol json.
- mcstar14.out: star network MC, K=8/16/32: ratio 1.335/1.142/1.065 = 1+2/(K-2) (vs 1+eps_K = 2.0/1.34/1.14). Consistent.
## Verdicts (draft)
E CONFIRMED (proof re-derived; exact enumeration; worst case 1+2/(K-2) <= 1+eps_K).
C CONFIRMED (valid inside Prop 3 / Lemma 8 general-pieces version; cost E_w T_q + (S_w-E_w) Delta_q exact).
B CONFIRMED WITH CAVEAT (asymptotic; eps_pl dropped; sensitivity negligible).
D CONFIRMED (level q = exact derivative; base q conservative; plain-sweep q=1 unproved but irrelevant).
c1 split+(iii) = 3.166 [8000,16000] reproduced; limsup over [1000,16000] = 3.174; thresholds reproduced.
- 04:52 EPB (row-sweep BASE caps also inflated: N_cap -> (1+eps) N_cap, eps = C sqrt(ln max(g,h)/min(g,h))) + EPC same C:
  C=1: c1 a=2 = 3.169/3.184/3.185/3.177 ; a=4 4.493.   C=3: 3.179/3.195/3.192/3.192 ; a=4 4.515.
  => the dropped eps_pl is NOT negligible at the base level (bases have g,h ~ 2^8..2^12, fixed as x -> inf). Crude McDiarmid
  constant: dev ~ sqrt(8 h^3 ln g) lam_cap vs mean h^2 lam_bar -> C ~ sqrt(8)*(lam_cap/lam_bar ~1.3) ~ 3.7. So c1 = 3.166 is the
  eps_pl -> 0 value; with the crude constant c1 ~ 3.19-3.20. (A Bernstein-type bound could lower C.)
- 04:55 C=3.7 (crude McDiarmid constant), levels+bases: c1 a=2 = 3.191/3.210/3.212/3.205 ; a=4 4.532.
  EPB=1 alone (levels unperturbed) = EPC=EPB=1 exactly -> the whole effect comes from the row-sweep BASES.
  Thresholds r=4.26 (2^n / 1.6423^n / 1.30704^n): C=0 1352/2683/9269 (= P13, reproduced on my own sol json);
  C=1 1368/2685/9275 ; C=3 1393/2745/9292.
## Final verdicts
E CONFIRMED. C CONFIRMED. D CONFIRMED (plain-sweep q=1 unproved but irrelevant to c1).
B CONFIRMED WITH CORRECTION: argument sound asymptotically for levels, but eps_pl is not negligible at the bases
  (fixed g,h ~ 2^8..2^12): c1(split+(iii)) = 3.166 only in the eps_pl -> 0 idealisation; with the placement-concentration
  term kept, c1 = 3.177 (C=1) .. 3.192 (C=3) .. 3.205 (C=3.7, crude McDiarmid) at [8000,16000]; limsup windows 3.185/3.195/3.212.
No process left running.
