# G17d — Adversarial review of G14a §3–§4 (output-level theorem, reduced degree) + Lean

Agent report (AI, same model family as the prover; not an independent human review). The Lean file was
re-compiled and checked by the coordinator and promoted to `lean/PnpLean/Article2/OutputLevel.lean`.
The agent's work folder was named `g18c` (numbering mistake at launch).

## CLAIM 1 (output-level theorem): CONFIRMED as stated, but the minor hypothesis is essential in general
Exact statement: F a field, A(t) ∈ F[t]^{n×n} arbitrary with entry degree ≤ d; labels α, β: [n]→[q]^ℓ
(need not be injective); τ_1..τ_ℓ distinct nonzero; A(τ_k) = B_k + E_k with B_k[u,v] = 0 unless
α(u)_k = β(v)_k. For fixed S_k ⊆ [q], let w(u) = #{k: α(u)_k ∉ S_k}, w′(v) = #{k: β(v)_k ∈ S_k}, and
2m+d+1 ≤ ℓ.

- **Cleanest correct form (no genericity):** Σ_k rank E_k ≥ rank A(0)[{w ≤ m}, {w′ ≤ m}].
- **One-line proof:** entrywise Lagrange interpolation of f_u·A_uv·g_v (prefactors with f(0)=g(0)=1)
  gives A(0)[R,C] = Σ_k L_k(0)·D_f(τ_k)·E_k·D_g(τ_k).
- The notes' version follows: a nonsingular A(0)[G,πG] gives ≥ |G|, and averaging over uniform S gives
  |dom π|·p² (needs α(u)_k ≠ β(πu)_k).
- **Counterexamples without the minor hypothesis:**
  - A ≡ all-ones 64×64, q=2, ℓ=6, d=0, π = complement, E_k = A: Σ rank = 6 < 7.56.
  - With every nonzero satisfying (iii): ℓ=1, five rows labelled 0, five columns labelled 1, all-ones:
    rank 1 < 1.25.
- **Is the hypothesis honest?** Yes: the bound is exactly a rank bound on submatrices of the input A(0),
  and low-rank inputs do escape. For the actual matrices it is removable at constant cost:
  - Δ-sparse A with a matching π has D′ ⊆ dom π, |D′| ≥ |dom π|/(2Δ−1), with A[D′,πD′] diagonal (greedy
    on the conflict graph; paper argument, not in Lean). So Σ rank ≥ |dom π|·p²/(2Δ−1) unconditionally.
  - In O'Donnell's App. B the W-cycle matching gives a monomial matrix, so it is nonsingular.
- **Consequence** d ≥ (1−o(1))√(½ℓ ln ℓ), exponent Ω(n(loglog n)²/√log n): re-derived; the reviewer's
  constant is ≈1/4 versus 1/√2.

## CLAIM 2: CONFIRMED
- **(i)** O'Donnell's text has deg P ≤ 4Rd (B.11) with R = 2r₀ and r₀ = Cn·(#rates), so 4R ≥ n. Nuance:
  the "4R crossing entries of one cover" count needs every factor to have a 3-variable monomial;
  otherwise it is between 2R and 4R, a constant factor.
- **(ii)** deg A_uv ≤ p_u+q_v: same interpolation with budgets w(u)+p_u ≤ M1, w′(v)+q_v ≤ M2,
  M1+M2 = ℓ−1. The Markov step and floor handling check out, under the same minor hypothesis.

## Numerics (own exact code mod 5, 7, 101; q=2,3; ℓ ≤ 4)
- Averaging identity E_S|G_S| = n·p² holds exactly in all cases.
- Adversarial test: random elements of the full solution space of paths with prescribed low-rank E_k,
  including non-uniform degrees and non-injective labels. 655 runs, 0 violations; the bound is attained
  in about 30% of runs.
- A refined graded bound, Σ rank E_k ≥ Σ_{i=0}^{ℓ−1−d} rank A0[{w≤i},{w′≤ℓ−1−d−i}], was also checked in
  the same runs with 0 violations (paper derivation only).

## Lean (`lean/PnpLean/Article2/OutputLevel.lean`, 269 lines, no sorry, standard axioms only)
- `eval_zero_eq_sum`: the Lagrange identity above.
- `rank_eval_zero_le`: (A.map (eval 0)).rank ≤ ∑ k, (E k).rank, for any prefactors and any N injective
  points.
- `rank_submatrix_le_sum_rank`: label version with S, m₁+m₂+d < ℓ and natDegree ≤ d+p u+q v; covers
  Claim 1 and Claim 2(ii).
- `card_le_sum_rank`: det A(0)[G,πG] ≠ 0 ⇒ |G| ≤ Σ rank E_k.
- `card_le_sum_rank_of_diagonal`: the diagonal-minor case.

## Remaining (paper only)
The averaging count over S, the Markov step, the greedy 1/(2Δ−1) lemma, the refined graded bound, and the
asymptotics.
