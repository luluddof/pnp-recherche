import PnpLean.Article1.WitnessDet
import PnpLean.Article1.Anticonc
import PnpLean.Article2.PerPoint
import PnpLean.Article2.QdPlus2
import PnpLean.Article2.AnyField
import PnpLean.Article2.MultiPoint

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
