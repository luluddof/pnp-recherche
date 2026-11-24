# G17c — Setting of Article 2 checked against the OpenAI preprint itself

Agent report (AI, same model family as the authors; not an independent human review).
P = OpenAI preprint "The Exponential Time Hypothesis Is False" (not copied into this repository);
OD = O'Donnell's exposition. The agent's work folder was named `g18b` (numbering mistake at launch).
Script: `check_preprint_local.py` (exact check, Q ≤ 6, that P's local construction satisfies the Lean
hypotheses and bound).

## Q1. MATCHES after translation; DIFFERS in presentation
P has no permanent and no A(t). It deforms replacement rows: "P⃗_{f,v}(0)=δ_v, deg ⩽ R" (Def. 3.1,
p.21-22), with H(T)=Σ_x ∏Φ_f(⟨P⃗_{f,v}(T),x⟩).
- Local object (§4.1-4.3): B=([r]^Q)², D=r^{2Q}, M0=V0·f_I(W0), V0=∏Π_s, and rows (4.10)
  P⃗_u(T)=δ_u(I−M0)+Σ M0[u,y]φ_{u,y}(T)δ_y with φ=∏_{s∈E_{u,y}}(1−T/s).
- Shape (4.14): "P⃗_{f,v}(T)=P⃗_{u_{f,v}}(T)π_{c(f,v)}ι", i.e. L(I_h⊗P^B)R′.
- P^B(0)=I, not X^{⊗ℓ}; relabelling columns by V0 gives J(0)=X^{⊗}. Niceness is "a coefficient can
  survive only if a_w(t)=d_u(t)".
- Translation: ℓ↔Q=⌊λ/128ν⌋, q↔r=⌊√λ/ν²⌋, d↔R−1 with R=4⌈√λ⌉, nodes t=1..Q, one per coordinate.
- Common error: "C^B is one integer matrix, independent of u,t"; rank c^B ⩽ D·Pr[Bin(Q,(r−1)/2r)∉I]
  ⩽ D/(3r²λ).
- B includes equal pairs, so it is a mixture of cubes {0,1}^j, j⩽Q, with Q−j direction-less nodes.
- Checked exactly, for Q⩽6: P's local construction satisfies the Lean hypotheses and bound.

## Q2. DIFFERS for the error charge
P's cost (2.11)/Prop. 3.4 is Σ_{j⩽q}2^j C(n,j)·(2D+1)^{2τn}·2^{2n/r}, with q=⌊3nR/Q⌋.
- Error: enumerated measurement tuples, (2D+1)^{r·h·c^B} by (4.16), not n^{Θ(rank)}. "A rank bound alone
  would not control the number of possible measurements" (p.22). So the exponent is
  ρ·r·(hD)·ln(2D+1), with ln D=2Q ln r.
- H1: holds for the construction. Def. 3.1 itself allows arbitrary per-node partitions S_{t,g} and
  per-node C_{t,g}.
- H2: block and Taylor terms match. The error term matches only at P's parameters (log₂D≈λ/128);
  n^{O(r)} is OD's permanent DP.
- H3: holds, and more strongly ("2(b+w)D² ⩽ M", n⩽M).
- H4: holds (QJ>d_H=dR, d⩾n; parity form d=k(m_test+2(ℓ−1))).
- OD's simplifications: permanent/A(t), base n, q=log¹⁰n, error ℓ^{−100}, unequal pairs only.

## Q3. CONSISTENT
- Ours (lower bound in the framework): Ω(n(loglog n)²/√log n).
- P interpolation: "O_k(ρ(n)), ρ(n)=n(loglog n)²/√log n" (p.2, (2.12), Prop. 6.4).
- P tensor contraction (Thm 1.1/App. D): O_k(n loglog n/√log n), and O_{k,D}(n/√log n) for m⩽Dn.

App. D is also an interpolation: W(0)=I, degree 2L−1, exact niceness at node t in coordinate v₂(t),
several nodes per coordinate. Our bound is vacuous there (N⩽d). It escapes through 2-adic "modes"
("2^{O(M/√u)}") and recursion, so H2 fails. Hence (loglog)² is not intrinsic to coordinate-nice
interpolation.

## Q4. NOT FOUND in P (novelty)
P has only the heuristic "Taking R much larger than √Q puts most basis vectors inside the interval".
OD §5.1 covers only the exact case: "exact niceness forces degree at least ℓ … The low-rank allowance is
what lets us lower the degree." Dias–Lima prove 2^{Ω(√n)} and say the gap to 2^{O(ρ(n))} "remains open".
The novelty claim stands.

## Q5. Freedoms not covered by our theorem
- arbitrary per-node partitions and measurements (Def. 3.1);
- cost between 2^{rank} and (2D+1)^{rows}; by the agent's estimate (unverified), the 2^{Ω(rank)} reading
  gives only Ω(n√(loglog n/log n));
- rank collapse through ι;
- product-degree reduction over algebras A (§5);
- 2-adic precision, tapered Hermite and recursion (App. D);
- equal-pair orbits.

## Q6. Framework X (for a future statement "the √log n loss is optimal within framework X")
- (a) Field F; J∈F[t]^{2^ℓ×2^ℓ}, deg⩽d, J(0)=X^{⊗ℓ} (P: P^B·V0 on an all-unequal swap orbit).
- (b) N distinct nonzero nodes, each with one coordinate c(k); J(τ_k)=K_k+E_k, K_k[x,y]=0 if x_c≠y_c.
- (c) Rows L(I_C⊗J)R′, M=|C|2^ℓ⩾n, (r(r−1))^ℓ⩽M, groups ⩽2n/r.
- (d) Node cost 2^{Ω(n/r)}·B^{Ω(Mρ_k)}·(n/h_k)^{Ω(h_k)}, with ln B=Ω(ℓ) for P or ln n for OD.
- (e) Hermite recovery, Σh_k>deg H, using only deg H⩽(⩾n)·d.
- (f) Single level over a field: no recursion, no 2-adic truncation.

Also to be stated: rank is charged per copy; direction-less nodes are dropped (needs r→∞); Lean covers
only the rank inequality; the asymptotics and Theorem 4 are paper proofs; Theorem 4 must be re-proved
with P's base (the agent's sketch says it survives, unverified).
