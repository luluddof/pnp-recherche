import PnpLean.Article1.WitnessDet
import PnpLean.Article1.Anticonc
import PnpLean.Article2.PerPoint
import PnpLean.Article2.QdPlus2
import PnpLean.Article2.AnyField
import PnpLean.Article2.MultiPoint
import PnpLean.Article2.OutputLevel
import PnpLean.Article3.DFRank
import PnpLean.Outils.SwapRank
import PnpLean.Outils.Intercalate

/-!
# PnpLean.Outils — reusable tools

A single entry point for the statements of this project that do not depend on SAT or on the
articles, and that other people may want to reuse.  Every name below is an `alias` of a theorem
proved elsewhere in this project: the statement is exactly the statement of the original.

Rules for this file:
* no `sorry`, `admit`, `axiom`, `native_decide`;
* every theorem depends only on `propext`, `Classical.choice`, `Quot.sound`
  (see the `#print axioms` lines at the end);
* each tool was searched for in Mathlib (November 2026) and not found.  This is a search in
  Mathlib only: the mathematics is classical, and other libraries were not searched.

See `lean/OUTILS.md` for the statements in plain words, the limits, and what is NOT included.
-/

namespace PnpOutils

/-! ## 1. Polynomial matrices: kernel dimensions at points versus the degree of the determinant -/

/-- **Kernel dimension divides the determinant.**  For a square matrix `Φ` over `ℚ[X]` and a
point `t`, `(X - t) ^ dim ker Φ(t)` divides `det Φ`.  (Generalises "geometric multiplicity ≤
algebraic multiplicity", which is the case `Φ = X • 1 - A`.) -/
alias pow_dvd_det_of_nullity := WitnessDet.pow_dvd_det_of_nullity

/-- **Sum of nullities, normalised form.**  If every entry of `Φ` has degree `≤ R` and
`Φ(0) = 1`, then `∑ t ∈ T, dim ker Φ(t) ≤ n * R` for every finite set of points `T`. -/
alias sum_nullity_le := WitnessDet.sum_nullity_le

/-- **Sum of nullities, general form.**  The same with the weaker hypothesis `det Φ(0) ≠ 0`. -/
alias sum_nullity_le_of_det_ne_zero := PerPoint.sum_nullity_le_of_det_ne_zero

/-- **Block deficiency bound.**  Rows and columns of `A` are split into groups by `h` and `g`.
If, for each group `γ`, the rows outside `γ` restricted to the columns of `γ` lie in a subspace
`Rγ γ`, then `∑ γ (columns of γ - rows of γ) ≤ dim ker A + ∑ γ dim (Rγ γ)`. -/
alias deficiency_le_nullity_add := WitnessDet.deficiency_le_nullity_add

/-! ## 2. Probability on finite product spaces, without measure theory -/

/-- The constant of the anticoncentration inequality: `(1 - e⁻²) / 4`. -/
noncomputable abbrev c0 : ℝ := UAvg.c0

/-- **Anticoncentration of a sum of independent indicators, explicit constant.**
On a finite product space with the uniform measure, let `X = ∑ i 1[β i (ω i)]` and
`V = ∑ i p i (1 - p i)` (its variance).  Then for every real `a`:
`c0 * min V (√V) ≤ E |X - a|`. -/
alias bernoulli_anticonc := UAvg.bernoulli_anticonc

/-- **Independence on a finite product.**  If `g` does not depend on coordinate `a`, then
`E[f (ω a) * g ω] = E[f] * E[g]` for uniform averages. -/
alias uavg_pi_mul := UAvg.uavg_pi_mul

/-! ## 3. Polynomial paths of matrices on the Boolean cube -/

/-- **Sum of ranks, many points per direction, any field.**  `J` is a `2^ℓ × 2^ℓ` matrix of
polynomials of degree `≤ d` with `J(0) = X^{⊗ℓ}` (the complement permutation).  Points `τ k`
are distinct and nonzero, each with a direction `c k`.  If `J(τ k) = K k + E k` with
`K k x y = 0` whenever `x (c k) ≠ y (c k)`, then
`∑ x, (#points - d - 2 * μ x) ≤ ∑ k, rank (E k)`, where `μ x` counts the points `k` with
`x (c k) = 1` (truncated subtraction in `ℕ`). -/
alias sum_rank_ge_multipoint := MultiPoint.sum_rank_ge

/-- **Rank at zero from errors at `N` points, any field, any shape.**  `A` is a matrix of
polynomials, `f u` and `g v` are row and column prefactors with value `1` at `0`, and
`f u * A u v * g v` has degree `< N`.  If at each of `N` distinct points `τ k` the prefactors kill
`A(τ k) - E k` entrywise, then `rank A(0) ≤ ∑ k, rank (E k)`.  (Entrywise Lagrange
interpolation.) -/
alias rank_eval_zero_le := OutputLevel.rank_eval_zero_le

/-- **Rank of a submatrix of `A(0)` from label-diagonal parts plus errors.**  Rows and columns
carry labels in `Q^ℓ`; `A(τ k) = B k + E k` with `B k u v = 0` when the `k`-th labels differ;
entry degrees `≤ d + p u + q v`.  For symbol sets `S k` and budgets `m₁ + m₂ + d < ℓ`, the
submatrix of `A(0)` on the rows and columns within budget has rank `≤ ∑ k, rank (E k)`. -/
alias rank_submatrix_le_sum_rank := OutputLevel.rank_submatrix_le_sum_rank

/-- **Sum of ranks, one point per direction, two-parameter form, any field.** -/
alias sum_rank_ge_N2 := AnyField.sum_rank_ge_N2

/-- **Block determinant form** (over `ℚ`): with directions in a set `S` of coordinates,
`(|T| - d) * 2 ^ (ℓ - |S|) ≤ ∑ t ∈ T, rank (E t)`. -/
alias sum_rank_ge_block := PerPoint.g6d_determinant_lemma

/-- **Four bad edges on the cube `Q_{d+2}`.**  Vertex polynomials of degree `≤ d` with value
`1` at `0`, and distinct nonzero points `τ k`, one per direction: at least 4 edges are bad. -/
alias four_bad_edges := QdPlus2.qd2_lemma

/-- Number of points of a subcube: `2 ^ (ℓ - |S|)`. -/
alias card_subcube := PerPoint.card_subcube

/-! ## 4. Rank trade-offs: the swap operator, and direction-factored identities -/

/-- **Rank versus operator-Schmidt rank for the swap with a spectator, any field.**
`swapW K m n Q` is the swap `E ⊗ G ⊗ W → G ⊗ E ⊗ W` (dimensions `m, n, Q`).  If
`A + B = swapW`, `rank A ≤ r` and the realignment of `B` has rank `≤ b`, then for every `k ≤ m`:
`k*Q*n ≤ r + k*k*Q*b` and `k*n*Q ≤ b*Q + k*k*r`; the mirror statements hold for `k ≤ n`; and
`b < m*n → Q ≤ r`. -/
alias swap_rank_tradeoff := SwapRank.swap_rank_osr_tradeoff

/-- **Product form:** for `m = n = M` and `1 ≤ b < M²`: `2*Q*M² ≤ 9*(r*b)`. -/
alias swap_rank_product := SwapRank.swap_rank_osr_product

/-- **Number of terms of a direction-factored identity.**  If a nonsingular `2^W × 2^W` matrix
is a sum of `N` terms, each a product over the groups of one direction of tables that read only
the bits of their group, then `2^|W| ≤ ∑ k, 2^(c k)` with `c = ∑ g, min (e₁ g) (e₂ g)`. -/
alias df_flattening_min := DFRank.df_flattening_min

/-- **A matrix congruent to `1` modulo a proper ideal has nonzero determinant.** -/
alias det_ne_zero_of_sub_one_mem := DFRank.det_ne_zero_of_sub_one_mem

/-- **Intercalate matrices with two palettes (key lemma).**  A colouring of `[N]²` whose colour
classes are partial permutations and in which every intercalate is bichromatic (palette `true` on
one diagonal, `false` on the other): a class of size `h` forces at least `h(h-1)/2` used colours of
the other palette. -/
alias intercalate_key_lemma := Intercalate.key_lemma

/-- **Counting:** `N² ≤ q·g*(p) + p·g*(q)` with `p`, `q` the numbers of used colours of each
palette and `g*(t)` the largest `g ≤ N` with `g(g-1)/2 ≤ t`. -/
alias intercalate_counting := Intercalate.counting

/-- **Trade-off (real form):** if `q ≤ p` then `N² ≤ 2·p·(√(2q) + 1)`. -/
alias intercalate_tradeoff := Intercalate.tradeoff_real

end PnpOutils

#print axioms PnpOutils.pow_dvd_det_of_nullity
#print axioms PnpOutils.sum_nullity_le
#print axioms PnpOutils.sum_nullity_le_of_det_ne_zero
#print axioms PnpOutils.deficiency_le_nullity_add
#print axioms PnpOutils.bernoulli_anticonc
#print axioms PnpOutils.uavg_pi_mul
#print axioms PnpOutils.sum_rank_ge_multipoint
#print axioms PnpOutils.sum_rank_ge_N2
#print axioms PnpOutils.sum_rank_ge_block
#print axioms PnpOutils.four_bad_edges
#print axioms PnpOutils.card_subcube
#print axioms PnpOutils.rank_eval_zero_le
#print axioms PnpOutils.rank_submatrix_le_sum_rank
#print axioms PnpOutils.swap_rank_tradeoff
#print axioms PnpOutils.swap_rank_product
#print axioms PnpOutils.df_flattening_min
#print axioms PnpOutils.det_ne_zero_of_sub_one_mem
#print axioms PnpOutils.intercalate_key_lemma
#print axioms PnpOutils.intercalate_counting
#print axioms PnpOutils.intercalate_tradeoff
