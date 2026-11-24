import Mathlib

/-!
# Output-level rank bound (G14a §3 and §4(ii)), deterministic form — G18c

Setting. `A(t)` is an ARBITRARY polynomial matrix over a field `F` (no factorisation
`L (I ⊗ J) R'` is assumed), rows `Row`, columns `Col`, with labels `α : Row → Fin ℓ → Q`,
`β : Col → Fin ℓ → Q` (not necessarily injective). At `ℓ` distinct nonzero points,
`A(τ k) = B k + E k` with `B k u v = 0` unless `α u k = β v k` ("block-diagonal by the `k`-th
symbol"). Entry degrees: `deg A u v ≤ d + p u + q v` (`p = q = 0` is the uniform case of §3,
`d = 0` is the potential case of §4(ii)).

Fix symbol sets `S k ⊆ Q` and budgets `m₁ + m₂ + d < ℓ`. Put
`w u = #{k : α u k ∉ S k}`, `w' v = #{k : β v k ∈ S k}`. Then for all families of rows `r` with
`w (r i) + p (r i) ≤ m₁` and of columns `c` with `w' (c j) + q (c j) ≤ m₂`:

  `rank (A(0)[r, c]) ≤ ∑ k, rank (E k)`      (`rank_submatrix_le_sum_rank`).

No genericity hypothesis: the nonsingular-minor hypothesis of the notes is only used to evaluate
the left-hand side (`card_le_sum_rank`: if `det A(0)[G, π G] ≠ 0` then `|G| ≤ ∑ rank E k`;
`card_le_sum_rank_of_diagonal`: the case where `A(0)[G, π G]` is diagonal).

Proof. Entrywise Lagrange interpolation of `f u * A u v * g v`, where
`f u = ∏_{k : α u k ∉ S k} (1 - X / τ k)`, `g v = ∏_{k : β v k ∈ S k} (1 - X / τ k)`:
this polynomial has degree `< ℓ`, value `A(0) u v` at `0`, and value
`f u (τ k) * E k u v * g v (τ k)` at `τ k` (the `B k` part is killed). Hence

  `A(0)[r, c] = ∑ k, L_k(0) • (diag (f(τ k)) * E k[r, c] * diag (g(τ k)))`   (`eval_zero_eq_sum`),

and rank is subadditive (`rank_eval_zero_le`, valid for any prefactors `f`, `g` with
`f(0) = g(0) = 1`, any number `N` of points, any placement).
-/

open Polynomial Matrix Module

namespace OutputLevel

variable {F : Type*} [Field F]

/-! ### Rank bookkeeping -/

section Rank

variable {m n : Type*} [Fintype m] [Fintype n]

omit [Fintype m] in
theorem rank_add_le (A B : Matrix m n F) : (A + B).rank ≤ A.rank + B.rank := by
  rw [← finrank_range_mulVecLin, ← finrank_range_mulVecLin, ← finrank_range_mulVecLin,
    Matrix.mulVecLin_add]
  refine le_trans (Submodule.finrank_mono ?_) (Submodule.finrank_add_le_finrank_add_finrank _ _)
  rintro _ ⟨x, rfl⟩
  exact Submodule.add_mem_sup (LinearMap.mem_range_self _ x) (LinearMap.mem_range_self _ x)

omit [Fintype m] in
theorem rank_sum_le {ι : Type*} (s : Finset ι) (M : ι → Matrix m n F) :
    (∑ k ∈ s, M k).rank ≤ ∑ k ∈ s, (M k).rank := by
  classical
  induction s using Finset.induction_on with
  | empty => simp
  | insert a s ha ih =>
    rw [Finset.sum_insert ha, Finset.sum_insert ha]
    exact (rank_add_le _ _).trans (Nat.add_le_add_left ih _)

omit [Fintype m] in
theorem rank_smul_le (c : F) (A : Matrix m n F) : (c • A).rank ≤ A.rank := by
  classical
  have h : c • A = A * (c • (1 : Matrix n n F)) := by rw [Matrix.mul_smul, Matrix.mul_one]
  rw [h]
  exact Matrix.rank_mul_le_left _ _

end Rank

/-! ### The core: entrywise Lagrange interpolation -/

/-- **Core identity.** If `f u * A u v * g v` has degree `< N`, `f u (0) = g v (0) = 1`, and at
each of the `N` distinct points `τ k` the matrix `A(τ k) - E k` is killed by the prefactors, then
`A(0) u v = ∑ k, L_k(0) * (f u (τ k) * E k u v * g v (τ k))`, `L_k` the Lagrange basis. -/
theorem eval_zero_eq_sum {R C : Type*} {N : ℕ}
    (A : Matrix R C F[X]) (f : R → F[X]) (g : C → F[X])
    (τ : Fin N → F) (hτ : Function.Injective τ) (E : Fin N → Matrix R C F)
    (hdeg : ∀ u v, (f u * A u v * g v).natDegree < N)
    (hf : ∀ u, (f u).eval 0 = 1) (hg : ∀ v, (g v).eval 0 = 1)
    (hvan : ∀ k u v,
      (f u).eval (τ k) * ((A u v).eval (τ k) - E k u v) * (g v).eval (τ k) = 0)
    (u : R) (v : C) :
    (A u v).eval 0 = ∑ k, (Lagrange.basis Finset.univ τ k).eval 0 *
      ((f u).eval (τ k) * E k u v * (g v).eval (τ k)) := by
  have hlt : (f u * A u v * g v).degree < (Finset.univ : Finset (Fin N)).card := by
    rw [Finset.card_univ, Fintype.card_fin]
    exact lt_of_le_of_lt degree_le_natDegree (by exact_mod_cast hdeg u v)
  have h := Lagrange.eq_interpolate (s := Finset.univ) (v := τ) hτ.injOn hlt
  have h0 := congrArg (eval 0) h
  rw [Lagrange.interpolate_apply] at h0
  simp only [eval_mul, eval_finsetSum, eval_C, hf, hg, one_mul, mul_one] at h0
  rw [h0]
  refine Finset.sum_congr rfl fun k _ => ?_
  have := hvan k u v
  linear_combination (eval 0 (Lagrange.basis Finset.univ τ k)) * this

/-- **Core rank bound** (any prefactors, any number of points). -/
theorem rank_eval_zero_le {R C : Type*} [Fintype R] [Fintype C] {N : ℕ}
    (A : Matrix R C F[X]) (f : R → F[X]) (g : C → F[X])
    (τ : Fin N → F) (hτ : Function.Injective τ) (E : Fin N → Matrix R C F)
    (hdeg : ∀ u v, (f u * A u v * g v).natDegree < N)
    (hf : ∀ u, (f u).eval 0 = 1) (hg : ∀ v, (g v).eval 0 = 1)
    (hvan : ∀ k u v,
      (f u).eval (τ k) * ((A u v).eval (τ k) - E k u v) * (g v).eval (τ k) = 0) :
    (A.map (eval 0)).rank ≤ ∑ k, (E k).rank := by
  classical
  have hEq : A.map (eval 0) = ∑ k, (Lagrange.basis Finset.univ τ k).eval 0 •
      (Matrix.diagonal (fun u => (f u).eval (τ k)) * E k *
        Matrix.diagonal (fun v => (g v).eval (τ k))) := by
    ext u v
    rw [Matrix.map_apply, eval_zero_eq_sum A f g τ hτ E hdeg hf hg hvan u v, Matrix.sum_apply]
    refine Finset.sum_congr rfl fun k _ => ?_
    simp [Matrix.smul_apply, Matrix.diagonal_mul, Matrix.mul_diagonal]
  rw [hEq]
  refine (rank_sum_le _ _).trans (Finset.sum_le_sum fun k _ => ?_)
  refine (rank_smul_le _ _).trans ?_
  exact (Matrix.rank_mul_le_left _ _).trans (Matrix.rank_mul_le_right _ _)

/-! ### Labels and symbol sets -/

variable {ℓ : ℕ}

/-- The prefactor `∏_{k ∈ T} (1 - X / τ k)`. -/
noncomputable def pref (τ : Fin ℓ → F) (T : Finset (Fin ℓ)) : F[X] :=
  ∏ k ∈ T, (1 - C (τ k)⁻¹ * X)

theorem pref_eval_zero (τ : Fin ℓ → F) (T : Finset (Fin ℓ)) : (pref τ T).eval 0 = 1 := by
  simp [pref, eval_prod]

theorem pref_eval_of_mem (τ : Fin ℓ → F) (T : Finset (Fin ℓ)) (k : Fin ℓ) (hk : τ k ≠ 0)
    (hkT : k ∈ T) : (pref τ T).eval (τ k) = 0 := by
  rw [pref, eval_prod]
  refine Finset.prod_eq_zero hkT ?_
  simp [hk]

theorem natDegree_pref_le (τ : Fin ℓ → F) (T : Finset (Fin ℓ)) :
    (pref τ T).natDegree ≤ T.card := by
  rw [pref]
  refine (natDegree_prod_le _ _).trans ?_
  have : ∀ k ∈ T, (1 - C (τ k)⁻¹ * X : F[X]).natDegree ≤ 1 := by
    intro k _
    compute_degree
  simpa using Finset.sum_le_card_nsmul _ _ 1 this

variable {Q : Type*} [DecidableEq Q]

/-- The coordinates where a row label is NOT in the symbol set. -/
def rowBad (S : Fin ℓ → Finset Q) (a : Fin ℓ → Q) : Finset (Fin ℓ) :=
  Finset.univ.filter (fun k => a k ∉ S k)

/-- The coordinates where a column label IS in the symbol set. -/
def colBad (S : Fin ℓ → Finset Q) (b : Fin ℓ → Q) : Finset (Fin ℓ) :=
  Finset.univ.filter (fun k => b k ∈ S k)

/-- **Output-level bound, fixed symbol sets, no genericity hypothesis.**
`rank A(0)[r, c] ≤ ∑ k, rank (E k)` for rows `r i` with `w + p ≤ m₁` and columns `c j` with
`w' + q ≤ m₂`, where `m₁ + m₂ + d < ℓ` and `deg A u v ≤ d + p u + q v`. -/
theorem rank_submatrix_le_sum_rank
    {Row Col ι κ : Type*} [Fintype Row] [Fintype Col] [Fintype ι] [Fintype κ]
    (A : Matrix Row Col F[X]) (d : ℕ) (p : Row → ℕ) (q : Col → ℕ)
    (hd : ∀ u v, (A u v).natDegree ≤ d + p u + q v)
    (α : Row → Fin ℓ → Q) (β : Col → Fin ℓ → Q)
    (τ : Fin ℓ → F) (hτ : Function.Injective τ) (hτ0 : ∀ k, τ k ≠ 0)
    (B E : Fin ℓ → Matrix Row Col F)
    (hdec : ∀ k, A.map (eval (τ k)) = B k + E k)
    (hB : ∀ k u v, α u k ≠ β v k → B k u v = 0)
    (S : Fin ℓ → Finset Q) (m₁ m₂ : ℕ) (hm : m₁ + m₂ + d < ℓ)
    (r : ι → Row) (c : κ → Col)
    (hr : ∀ i, (rowBad S (α (r i))).card + p (r i) ≤ m₁)
    (hc : ∀ j, (colBad S (β (c j))).card + q (c j) ≤ m₂) :
    ((A.map (eval 0)).submatrix r c).rank ≤ ∑ k, (E k).rank := by
  have key := rank_eval_zero_le (A.submatrix r c)
    (fun i => pref τ (rowBad S (α (r i)))) (fun j => pref τ (colBad S (β (c j)))) τ hτ
    (fun k => (E k).submatrix r c) ?_ (fun i => pref_eval_zero τ _) (fun j => pref_eval_zero τ _) ?_
  · have hsub : (A.map (eval 0)).submatrix r c = (A.submatrix r c).map (eval 0) := rfl
    rw [hsub]
    exact key.trans (Finset.sum_le_sum fun k _ => Matrix.rank_submatrix_le _ _ _)
  · intro i j
    have h1 := natDegree_pref_le τ (rowBad S (α (r i)))
    have h2 := natDegree_pref_le τ (colBad S (β (c j)))
    have h3 := hd (r i) (c j)
    have h4 := hr i
    have h5 := hc j
    have h6 : (pref τ (rowBad S (α (r i))) * A.submatrix r c i j
        * pref τ (colBad S (β (c j)))).natDegree
        ≤ (rowBad S (α (r i))).card + (d + p (r i) + q (c j)) + (colBad S (β (c j))).card :=
      natDegree_mul_le_of_le (natDegree_mul_le_of_le h1 h3) h2
    omega
  · intro k i j
    by_cases h1 : α (r i) k ∈ S k
    · by_cases h2 : β (c j) k ∈ S k
      · have hz : (pref τ (colBad S (β (c j)))).eval (τ k) = 0 :=
          pref_eval_of_mem τ _ k (hτ0 k) (by simp [colBad, h2])
        simp only [hz, mul_zero]
      · have hne : α (r i) k ≠ β (c j) k := fun h => h2 (h ▸ h1)
        have hJ : (A (r i) (c j)).eval (τ k) = B k (r i) (c j) + E k (r i) (c j) := by
          have := congrFun (congrFun (hdec k) (r i)) (c j)
          simpa using this
        simp only [Matrix.submatrix_apply]
        rw [hJ, hB k _ _ hne]
        ring
    · have hz : (pref τ (rowBad S (α (r i)))).eval (τ k) = 0 :=
        pref_eval_of_mem τ _ k (hτ0 k) (by simp [rowBad, h1])
      simp only [hz, zero_mul]

/-- **G14a §3 / §4(ii), deterministic part** (fixed symbol sets `S`): if the minor
`A(0)[G, π G]` is nonsingular and every `u ∈ G` satisfies the two budget conditions, then
`|G| ≤ ∑ k, rank (E k)`. -/
theorem card_le_sum_rank
    {Row Col : Type*} [Fintype Row] [Fintype Col] [DecidableEq Row]
    (A : Matrix Row Col F[X]) (d : ℕ) (p : Row → ℕ) (q : Col → ℕ)
    (hd : ∀ u v, (A u v).natDegree ≤ d + p u + q v)
    (α : Row → Fin ℓ → Q) (β : Col → Fin ℓ → Q)
    (τ : Fin ℓ → F) (hτ : Function.Injective τ) (hτ0 : ∀ k, τ k ≠ 0)
    (B E : Fin ℓ → Matrix Row Col F)
    (hdec : ∀ k, A.map (eval (τ k)) = B k + E k)
    (hB : ∀ k u v, α u k ≠ β v k → B k u v = 0)
    (S : Fin ℓ → Finset Q) (m₁ m₂ : ℕ) (hm : m₁ + m₂ + d < ℓ)
    (G : Finset Row) (π : Row → Col)
    (hG : ∀ u ∈ G, (rowBad S (α u)).card + p u ≤ m₁ ∧
      (colBad S (β (π u))).card + q (π u) ≤ m₂)
    (hns : ((A.map (eval 0)).submatrix (fun u : G => (u : Row)) (fun u : G => π u)).det ≠ 0) :
    G.card ≤ ∑ k, (E k).rank := by
  have h := rank_submatrix_le_sum_rank A d p q hd α β τ hτ hτ0 B E hdec hB S m₁ m₂ hm
    (fun u : G => (u : Row)) (fun u : G => π u) (fun u => (hG u u.2).1) (fun u => (hG u u.2).2)
  rwa [Matrix.rank_of_isUnit _ ((Matrix.isUnit_iff_isUnit_det _).mpr (IsUnit.mk0 _ hns)),
    Fintype.card_coe] at h

/-- Special case with no genericity at all: `A(0)[G, π G]` is diagonal with nonzero diagonal
(e.g. `G` an independent set of the conflict graph of a sparse `A(0)`; a sparse matrix with a
matching `π` in its support always has such a `G ⊆ dom π` of size `≥ |dom π| / (2Δ - 1)`). -/
theorem card_le_sum_rank_of_diagonal
    {Row Col : Type*} [Fintype Row] [Fintype Col] [DecidableEq Row]
    (A : Matrix Row Col F[X]) (d : ℕ) (p : Row → ℕ) (q : Col → ℕ)
    (hd : ∀ u v, (A u v).natDegree ≤ d + p u + q v)
    (α : Row → Fin ℓ → Q) (β : Col → Fin ℓ → Q)
    (τ : Fin ℓ → F) (hτ : Function.Injective τ) (hτ0 : ∀ k, τ k ≠ 0)
    (B E : Fin ℓ → Matrix Row Col F)
    (hdec : ∀ k, A.map (eval (τ k)) = B k + E k)
    (hB : ∀ k u v, α u k ≠ β v k → B k u v = 0)
    (S : Fin ℓ → Finset Q) (m₁ m₂ : ℕ) (hm : m₁ + m₂ + d < ℓ)
    (G : Finset Row) (π : Row → Col)
    (hG : ∀ u ∈ G, (rowBad S (α u)).card + p u ≤ m₁ ∧
      (colBad S (β (π u))).card + q (π u) ≤ m₂)
    (hdiag : ∀ u ∈ G, (A u (π u)).eval 0 ≠ 0)
    (hoff : ∀ u ∈ G, ∀ u' ∈ G, u ≠ u' → (A u (π u')).eval 0 = 0) :
    G.card ≤ ∑ k, (E k).rank := by
  refine card_le_sum_rank A d p q hd α β τ hτ hτ0 B E hdec hB S m₁ m₂ hm G π hG ?_
  have hEq : (A.map (eval 0)).submatrix (fun u : G => (u : Row)) (fun u : G => π u)
      = Matrix.diagonal (fun u : G => (A u (π u)).eval 0) := by
    ext u u'
    by_cases h : u = u'
    · subst h
      simp
    · have hne : (u : Row) ≠ (u' : Row) := fun h' => h (Subtype.ext h')
      simp [Matrix.diagonal_apply_ne _ h, hoff u u.2 u' u'.2 hne]
  rw [hEq, Matrix.det_diagonal]
  exact Finset.prod_ne_zero_iff.mpr fun u _ => hdiag u u.2

end OutputLevel

#print axioms OutputLevel.eval_zero_eq_sum
#print axioms OutputLevel.rank_eval_zero_le
#print axioms OutputLevel.rank_submatrix_le_sum_rank
#print axioms OutputLevel.card_le_sum_rank
#print axioms OutputLevel.card_le_sum_rank_of_diagonal
