# G14b — Limits of the tensor-contraction route (Appendix D) of OpenAI's ETH preprint

Status: working notes, level P (proofs written, not reviewed) for Theorems A, B; level E for the rest.
Source actually read: the preprint itself (dl/paper.pdf, dl/paper.txt; OpenAI, "The Exponential Time
Hypothesis is False", 25 Oct 2026, from raw.githubusercontent.com/openai/math, family 378), App. D
p.53-64 (paper.txt lines 3007-3655).  pdftotext lost superscripts: I read "L = 2a" as L = 2^a and
"Darr = r2a" as D_arr = r^{2a} (forced by the proof of Lemma D.3: gamma ranges over pairs of labels).
I did NOT check the correctness of App. D; everything below takes it at face value.

## 1. Reconstruction (quotes are from App. D)

Object. "A Boolean tensor network is a bipartite multigraph G in which each edge carries a single
bit. Each vertex v has ... a tensor T_v"; Z = sum_x prod_v T_v(x|delta(v)).  Thm D.1: <= Cs edges, max
degree <= s^{9/10}, entries in (Z/2^{P0})[V]/(V^l): Z computable in (s+L)^{O(1)} 2^{O(s/sqrt(log s))}, poly space.
#k-SAT at m <= Dn: clause vertices + equality trees, every edge subdivided, s = n+m, ring Z/2^{n+1}.

One call (M wires, depth d).  Parameters (D.4): lambda = ceil(log2 s), H = ceil(s/floor(sqrt(lambda))) (leaf
size), r = 64, eps = 1/(128*1024) = 2^-17, u = u_d = floor(eps*lambda/2^d) "retained valuation levels",
a = a_d = 2*sum_{j<=d} u_j, L = 2^a, D_arr = r^{2a}.  Coefficient ring: O = Z_(2)[U]/(f), deg f = 4*lambda,
f irreducible mod 2 (a 2-adic unramified extension; only used to have >= s^4 interpolation points with
unit differences), A_p = O[V]/(V^l, 2^p).
 (i) Packing (Lemma D.3): vertices get labels in [r]^a (conditional expectations, pair loads nearly
     uniform: sum_gamma (X_gamma - M/D_arr)^2 <= 2 M Delta).  Edges are put in <= 3M/L "packets of L slots";
     in a packet the unordered pair of symbols is fixed in every coordinate and "bit i of the slot
     number specifies their orientation in coordinate i".
 (ii) Cauchy deformation (D.7): S = {-L..L-1}, C(T) = prod_{t in S}(T-t),
     W_xy(T) = C(T) / (C'(0) (T-(y-x))) inside a packet, 0 between packets.  "W(0) = I".  deg = 2L-1.
     F(T) = sum_{X,Y} A(Y)B(X) prod_x l_{Y_x}((W(T)X)_x)  (D.8), deg F <= M(2L-1), F(0) = Z.
     [Checked exactly, C3: W(t) = c_t * (shift by t); a surviving pair x -> x+t differs in bit v2(t)
     and agrees below.]  Hence "this choice of coordinate makes W(t) local": with i = v2(t), group = symbol
     in coordinate i, W(t) joins only ports of the same group.  NO error term at the nodes.
 (iii) Two-adic scaling: alpha = 2^{a+1}; Lemma D.4: W(t+alpha z) in O[z]; every non-local entry is
     divisible by 2.
 (iv) Tapered Hermite recovery (Lemma D.5): h = ceil(3M/u), multiplicity q_t = h*max(0,u-v2(t)),
     D_t = L/2^{v2(t)} (= number of nodes of that valuation); F(0) mod 2^p from F(t+alpha z) mod
     (z^{q_t}, 2^{p_t}), p_t = p + 8hD_t + 8M(a+1).  Loss (D.10):
     v2(E_t(0)) = D * sum_{j>=1} (j/2^j)(q(i+j)-q(i)) >= -6Dh.  "The taper accounts for the factor 1/u":
     constant multiplicity Theta(M) on the retained levels would lose Theta(M D_t).
     [(D.10) checked exactly for arbitrary valuation-dependent profiles, C1/C1b.]
 (v) Modes (D.12)-(D.14): columns of a packet are cut into tiles y = y0 + 2^{i+1} j; for a row crossing
     to a tile, 1/(T-y+x) is expanded geometrically in 2^{i+1} j/(T-y0+x); mode Z = (packet, tile, b),
     weight w_Z = b+1, L_Z(X) = sum_j j^b X_y, "the coefficient of Z is divisible by 2^{w_Z}",
     "at most (3M/L) 2^{i+1} = 6M/D variables of any one weight".  Only tuples with
     sum_Z w_Z beta_Z < p_t are kept.  Their number is exp(O(sqrt(M p_t/D))) = exp(O_C(M/sqrt u))  (D.15).
 (vi) Recursion (D.16): F(t+alpha zeta) = sum_beta prod_{g=1}^{r} Z_{g,beta} (mod 2^{p_t}); the child of
     group g keeps the vertices of group g and gets one new "packet tensor" (<= 2L ports) per packet;
     <= 6M/r < M/8 edges.  M <= Cs 8^-d, u_d = Theta(lambda/2^d): exponents sum geometrically to
     O(s/sqrt(lambda)); leaves (M <= H) by brute force 2^H.
WHERE sqrt(log) COMES FROM: a partition count.  The budget is the precision p_t ~ M D_t/u (set by the
2-adic loss of the recovery), the supply is 6M/D_t modes per unit weight; #tuples =
exp(Theta(sqrt(budget * supply))) = exp(Theta(M/sqrt u)), and u <= a = O(log M) because the L = 2^a
slots of a packet and the r^{2a} pair classes must be populated ((D.6): M >= 2 D_arr^2 Delta).
This is NOT the mechanism of the main text / O'Donnell (degree d = sqrt(l log l) parity approximation,
low-rank error, derivative order h = nd/l << n).  Here the degree (2L-1) exceeds the number of
directions, the error at the nodes is zero, and the derivative order is Theta(M) (not small): the
price of Taylor data is paid 2-adically instead of combinatorially.
Not reconstructed / not checked: correctness of D.3-D.5 and (D.16); Section D.3 (arbitrary density,
wires indexing monomials of the parity lift) was only skimmed; bit-complexity claims.

## 2. Model

Wires e = 1..M with a row-end bit Y_e (on row vertex) and a column-end bit X_e (on column vertex);
identity tensor I(Y,X) = prod_e [Y_e = X_e].  Directions i = 1..u: partitions pi_i of the vertices.
A direction-i term is tau(Y,X) = prod_{g in pi_i} tau_g(end bits on vertices of g)  (a product of
disjoint child networks).
 Model DF (direction-factored, 2-adic):  I == sum_{k<=N} kappa_k tau^(k)  entrywise mod 2^p (p >= 1),
   kappa_k in Frac(O) arbitrary (2-adically large allowed), each tau^(k) a term of some direction.
   Cost >= N.  [App. D's root identity is of this form: k = (t, zeta, beta), direction v2(t),
   kappa = (Hermite weight)*(Lagrange weight in zeta).  The reduce-and-lift steps only add multiples of
   2^{p_t}, which become multiples of 2^p after multiplication by kappa.]
 Model VT (valuation-tapered, = App. D with free parameters):  grid S = {-L..L-1}, L = 2^a, 1 <= u <= a,
   any multiplicities q_t = q(v2(t)) >= 0 with q(i) = 0 for i >= u and sum_t q_t > M(2L-1); recovery by
   the Hermite/CRT formula (D.9), so p_t >= -v2(E_t(0)); at a node of valuation i there are
   n_i = nu*M/D_t modes of each weight (nu <= 6; nu >= const on instances with populated tiles) and all
   tuples with sum w beta < p_t are enumerated.  Cost >= max_t #tuples.
 Relation to G9d.  Unrolled, App. D is a signed (2-adically weighted) sum of tree-contracted networks,
   so it lies in the 2-ADIC CLOSURE of class G (congruence mod 2^p with unbounded coefficients, the
   analogue of border decompositions), not in G over Q: the mode series (D.12) is not a finite identity
   over Q.  The harmonic law's hypothesis FAILS: gadgets = packets, pieces = per-packet mode monomials
   beta_P, and a global term is a product of one piece per packet, but (a) all packets of a term share
   the node t (hence the direction), (b) the pieces are coupled by the global budget
   sum_P wt(beta_P) < p_t, (c) the piece family of a gadget is not poly-size (it is graded by weight).
   The law's CONCLUSION is not contradicted (M/sqrt u >= M/u, u = O(log M)).

## 3. Results

THEOREM A (flattening survives 2-adic approximation; model DF).  Let C1, C2 be disjoint vertex sets
such that no group of any pi_i meets both, and let w = #wires with one end in C1 and the other in C2.
Then N >= 2^w.
Proof.  Fix every end bit not on a C1-C2 wire to 0.  I restricts to a 2^w x 2^w permutation matrix P
(rows: bits of the ends in C1, columns: ends in C2).  In a term of direction i every group meets at
most one of C1, C2, so the restricted term is f(row index)*g(column index): rank <= 1 over the field
Frac(O).  So P + 2^p E = sum_k kappa_k (rank <= 1) with E over O.  det(P + 2^p E) == +-1 mod 2O, and
2O is a proper ideal, so the left side has rank 2^w.  QED
Corollary A1 (App. D's own labeling).  Take halves S_i of [r], C1 = {labels with symbol in S_i for all
i<u}, C2 = the opposite orthant.  The pair classes joining them are a fraction 2*4^-u of D_arr; by
Cauchy-Schwarz and Lemma D.3's variance bound, w >= 2*4^-u M (1 - 2^u/(sqrt2 r^a)).  So every DF
identity on that labeling has N >= 2^{(2-o(1)) M/4^u}.  With product labels [r]^u and populated pair
classes (r^{2u} <= M) this is 2^{Omega(M^{1-1/log2 r})} (r = 64: 2^{Omega(M^{5/6})}); it is the G9d
"4^{1-d}" floor, now valid mod 2^p.  Consequence: TWO-ADIC VALUATIONS BY THEMSELVES DO NOT ESCAPE
FLATTENING RANK.  What they escape is Dias-Lima's (H2)/(H3) (generic-degree exact recovery, cost
exponential in the Taylor order).

THEOREM B (the taper is optimal up to a constant; model VT).  Put
Lambda_i = -v2(E_t(0))/D_t = sum_{j=1}^{a-i} (j/2^j)(q(i)-q(i+j)).  For every admissible profile,
max_{i<u} Lambda_i >= (2-1/L) M/(4u), hence some node has p_t >= M D_t/(4u)(2-1/L) and the scheme
enumerates at least 2^{sqrt(n p_t) - n} >= 2^{(M/2) sqrt(nu/u) - nu M/D_t} tuples there.  Since u <= a <=
log2 M, the cost is 2^{Omega(M/sqrt(log M))}: within its own design App. D cannot beat its stated bound
by more than the constant.
Proof.  sigma_i := sum_{j<=a-i} j/2^j in [1/2,2).  Lambda_i <= Lambda gives sigma_i q_i <= Lambda +
sum_{j=1}^{u-1-i} (j/2^j) q_{i+j}.  Induction downward on i: q_i <= 2 Lambda (u-i).  (i = u-1:
q <= Lambda/sigma <= 2 Lambda.  Step: sigma_i q_i <= Lambda + 2 Lambda[(u-i) sigma_i - sum_{j<=u-1-i} j^2/2^j]
<= 2 Lambda (u-i) sigma_i since the last sum is >= 1/2.)  Then M(2L-1) < sum_t q_t = L sum_i q_i/2^i <=
2 Lambda L sum_i (u-i)/2^i <= 4 Lambda L u.  Tuples: with W = floor(sqrt(p/n)) the 0/1 exponent vectors
on the nW modes of weight <= W have weight <= nW(W+1)/2 < p (p > n), so there are >= 2^{nW} >=
2^{sqrt(np)-n} tuples; n p_t >= nu M^2 (2-1/L)/(4u).  QED
Extra hypothesis used in the tuple count: p_t > n, i.e. D_t^2 > 4 u nu (true in App. D, where a >= 2u so
D_t >= 2^u; in the bare model assume a >= u + log2(u) + 3).  Also L = 2^a <= M by (D.6).
Numerics (code/checks_out.txt): LP optimum of max Lambda_i is ~ 5.9 M/u (a = 2u = 60) against 18 M/u for
the paper's taper: at most a factor ~3 in the precision, sqrt 3 in the exponent, can be gained.

MAP of the three techniques against App. D.
 - Witness determinant (Dias-Lima Lemma 6.2, our G6d): needs "few easy points or small Taylor order".
   Here det W has degree L(2L-1), rank drops sum to ~L^2: no contradiction.  The technique's charge
   binom(M, q_t) is exactly what the mode expansion avoids (q_t ~ M, yet exp(M/sqrt u) terms).  ESCAPES.
 - Dual certificate / positivity (G11b): needs deg < #directions - 1 and non-zero per-point errors;
   here deg = 2L-1 >= a and the errors are 0.  VACUOUS.
 - Flattening ranks: APPLY (Thm A) but give only M/4^u per level.
 - Harmonic law: hypothesis fails (shared node + global weight budget).
The escaping feature is therefore not the valuation as such but the pair
 {exact high-degree deformation with zero error at 2L-1 nodes spread over u directions}
 + {Taylor data of order Theta(M) priced by a 2-adic weight budget, i.e. budget-coupled, non-product terms}.

## 4. Next statement
Direction-rank problem (2-adic Naslund problem in u directions): for product labels [r]^u, is
N_DF >= 2^{c M/u} (or even 2^{c M/sqrt u})?  Proved: 2^{2M/4^u}.  Achieved (if App. D is right):
2^{O(M/sqrt u)}.  A proof of 2^{Omega(M/u)} for entangled terms would extend G6b's rho_d >= 1/d beyond
product pieces and give 2^{Omega(M/log M)} for every one-level scheme of this shape.  Smallest test:
u = 2, r = 2 is exactly N <= P(N) <= N^{sqrt 2}.

## 5. Adversarial self-review
 - Thm A is one level only and counts terms; children could in principle be cheap.  It needs the terms
   to factor over ALL groups of one direction (true for (D.16)).
 - Thm B is a statement about App. D's design: recovery by formula (D.9), valuation-dependent
   multiplicities, enumeration of all tuples of (D.14), nu bounded below.  UNPROVED: that p_t > -v2(E_t(0))
   is necessary for every recovery rule (sketch: G = prod_{q != t0} ((T-q)/(-q))^{q_q} has zero data
   outside t0 and data of valuation -v2(E_{t0}(0)) at t0, but it only has degree <= M(2L-1) for tight
   profiles, and its integrality at the other nodes was not checked).  UNPROVED: that tuples with zero
   coefficient cannot be skipped cheaply.
 - "2-adic closure of G" is my formalisation; G9d itself is level E.
 - Superscripts in (D.4) were inferred.  App. D's correctness is not assessed.
