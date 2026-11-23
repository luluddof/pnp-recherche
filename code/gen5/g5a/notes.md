# G5a: non-product (expander) designs for the crossing expansion

## 1. Open question (Limits paper, §9 "The gap, and what an improved scheme must change")
"Narrowing the gap needs one of two new ingredients. From above, a design on an expander with sqrt(n) groups,
O(sqrt n) interaction width, Taylor order O(sqrt n) and O(1) or polylogarithmically many rows per measured group
would bring the scheme to 2^{Õ(sqrt n)}. From below, one must show that measurements cannot absorb the cut."
Object = Def. 5.1(c) design at every point t in X: partition of the variables into groups S_{t,h} (s_t ~ sqrt n),
measurement subspaces R_{t,h} (dim polylog), with P_{e,v}(t)|_{S_{t,h}} in R_{t,h} for every foreign group h;
q_t ~ sqrt n forces |X| >= 3nR/q_t ~ 3 sqrt(n) R points.

## 2. Consequences in the clause-product language
- Random 3-CNF: any partition into sqrt n-sized groups leaves Theta(n) foreign occurrences per point;
  roots of c_o give sum_t |Z_t| <= MR = sum_t q_t, so only O(sqrt n) per point can be cut;
  Theta(n) must be ABSORBED by Õ(sqrt n) measurement dimensions.
- Rank vanishing lemma (generalises the paper's 2x2 minor remark): for a fixed set B of occurrences with distinct
  variables, M_B(T) = [P_{o_i}(T)_{v_j}] has M_B(0)=I, deg det <= |B|R, so sum_t nullity(M_B(t)) <= |B|R.
  A group with r measurement rows absorbing k>r of B forces nullity >= k-r at t.
  Counting with this lemma does NOT give a contradiction (e.g. 1-dim "sum" measurements only need group size
  s <= O(q)); so it does not refute the design: OPEN in this route.

## 3. Matrix (OpenAI/G4c) route: exact formalisation
Partition system phi: V -> prod [q_i] (D directions), s = Hamming metric, crossing C with s(x,Cx)=D.
deg_t J <= d, J(0)=C, J(t_i)=K_i+E (K_i pi_i-block-diagonal, common E). Gain g = d/D.
Lemma A (exact): feasible with error E iff E = C on all far pairs s(x,y)>d (near pairs are free: interpolate).
  => r*(d) = min-rank completion of C on far pattern.
Lemma B (lower bound): r*(d) >= alpha(H_d), H_d: a~a' iff s(a,Ca')<=d or s(a',Ca)<=d
  >= largest anticode of diameter < D-d (triangle inequality)  >= K/(1+N_d), N_d = max #points within d.
Corollary C: if phi(V) has minimum distance > d (good code / expander-like / 2-design where every pair is
  separated in >= D-1 classes), then r*(d) >= K/2 (binary antipodal) — gain must exceed relative min distance.
  Product cube: anticode = Hamming ball, gain ~ sqrt(2 ln(1/delta)/log K) (G4c).
So the interpolation mechanism needs label pairs separated by FEW directions (locality), the opposite of expansion.

## 4. Experiment (crossing_design.py, search_codes.py): see table_all.txt, search_k4.txt, search_k5.txt
Alpha LB and translation-invariant UB coincide in almost all rows => r*(d) essentially pinned.
