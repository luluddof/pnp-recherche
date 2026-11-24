import PnpLean.Article1.WitnessDet

/-!
# Per-point errors in the interpolation framework: the G6d determinant lemma

Abstract core: a square polynomial matrix `M` over `ℚ[X]` with entries of degree `≤ d` whose value
at `0` is invertible satisfies `∑_{t ∈ T} nullity (M(t)) ≤ m * d`; hence
`∑_{t ∈ T} rank (M(t)) ≥ (|T| - d) * m`.

Concrete form (`g6d_determinant_lemma`): the lemma of `code/gen6/g6d/notes.md` §2 on the cube
`Fin ℓ → Bool`, with `J(0) = X^{⊗ℓ}` and weakly nice `K_t`.
-/

open Polynomial Matrix Module

namespace PerPoint

open WitnessDet

section Abstract

variable {m : ℕ}

theorem det_ne_zero_of_eval_det_ne_zero (M : Matrix (Fin m) (Fin m) ℚ[X]) (t₀ : ℚ)
    (h0 : (evalMat M t₀).det ≠ 0) : M.det ≠ 0 := by
  intro h
  apply h0
  rw [← eval_det, h, eval_zero]

/-- **Abstract determinant lemma.** Entries of degree `≤ d`, `det M(0) ≠ 0`: the nullities of
`M(t)` over any finite set of points sum to at most `m * d`. -/
theorem sum_nullity_le_of_det_ne_zero (M : Matrix (Fin m) (Fin m) ℚ[X]) (d : ℕ)
    (hd : ∀ i j, (M i j).natDegree ≤ d) (h0 : (evalMat M 0).det ≠ 0) (T : Finset ℚ) :
    ∑ t ∈ T, nullity (evalMat M t) ≤ m * d := by
  have hdvd : ∏ t ∈ T, (X - C t) ^ nullity (evalMat M t) ∣ M.det := by
    apply Finset.prod_dvd_of_coprime
    · intro a _ b _ hab
      exact ((pairwise_coprime_X_sub_C (s := id) Function.injective_id) hab).pow
    · intro t _
      exact pow_dvd_det_of_nullity M t
  have hdeg := natDegree_le_of_dvd hdvd (det_ne_zero_of_eval_det_ne_zero M 0 h0)
  rw [natDegree_prod_of_monic _ _ (fun t _ => (monic_X_sub_C t).pow _)] at hdeg
  simp only [natDegree_pow, natDegree_X_sub_C, mul_one] at hdeg
  exact hdeg.trans (natDegree_det_le M d hd)

/-- Rank form of the abstract determinant lemma. -/
theorem sum_rank_ge (M : Matrix (Fin m) (Fin m) ℚ[X]) (d : ℕ)
    (hd : ∀ i j, (M i j).natDegree ≤ d) (h0 : (evalMat M 0).det ≠ 0) (T : Finset ℚ) :
    (T.card - d) * m ≤ ∑ t ∈ T, (evalMat M t).rank := by
  have h1 := sum_nullity_le_of_det_ne_zero M d hd h0 T
  have h2 : ∑ t ∈ T, ((evalMat M t).rank + nullity (evalMat M t)) = T.card * m := by
    simp [nullity_add_rank]
  rw [Finset.sum_add_distrib] at h2
  have h3 : (T.card - d) * m = T.card * m - m * d := by rw [Nat.sub_mul, Nat.mul_comm d m]
  omega

/-- The rank of a submatrix (arbitrary row/column selections) is at most the rank. -/
theorem rank_submatrix_le' {ι κ : Type*} [Fintype ι] [Fintype κ] [DecidableEq ι] [DecidableEq κ]
    {p q : ℕ} (A : Matrix ι κ ℚ) (r : Fin p → ι) (c : Fin q → κ) :
    (A.submatrix r c).rank ≤ A.rank := by
  have h : A.submatrix r c
      = (1 : Matrix ι ι ℚ).submatrix r id * A * (1 : Matrix κ κ ℚ).submatrix id c := by
    ext i j
    simp [Matrix.mul_apply, Matrix.one_apply]
  rw [h]
  exact (rank_mul_le_left _ _).trans (rank_mul_le_right _ _)

/-- **Block form.** `J` a polynomial matrix with entries of degree `≤ d`; `r, c` select an `m × m`
block whose value at `0` is invertible. If at every `t ∈ T`, `J(t) = K_t + E_t` with `K_t`
vanishing on the block, then `∑_{t ∈ T} rank E_t ≥ (|T| - d) * m`. -/
theorem sum_rank_error_ge {ι : Type*} [Fintype ι] [DecidableEq ι]
    (J : Matrix ι ι ℚ[X]) (d : ℕ) (hd : ∀ x y, (J x y).natDegree ≤ d)
    (r c : Fin m → ι) (h0 : ((J.submatrix r c).map (eval 0)).det ≠ 0)
    (T : Finset ℚ) (K E : ℚ → Matrix ι ι ℚ)
    (hdec : ∀ t ∈ T, J.map (eval t) = K t + E t)
    (hK : ∀ t ∈ T, ∀ i j, K t (r i) (c j) = 0) :
    (T.card - d) * m ≤ ∑ t ∈ T, (E t).rank := by
  refine (sum_rank_ge (J.submatrix r c) d (fun i j => hd _ _) h0 T).trans ?_
  refine Finset.sum_le_sum fun t ht => ?_
  have h : evalMat (J.submatrix r c) t = (E t).submatrix r c := by
    ext i j
    have := congrFun (congrFun (hdec t ht) (r i)) (c j)
    simp only [Matrix.map_apply, Matrix.add_apply, hK t ht i j, zero_add] at this
    simp [evalMat, this]
  rw [h]
  exact rank_submatrix_le' _ _ _

end Abstract

section Cube

/-- `X^{⊗ℓ}` on the cube `Fin ℓ → Bool`: entry `(x,y)` is `∏ᵢ X[xᵢ,yᵢ]`, `X = [[0,1],[1,0]]`. -/
def Xtensor (ℓ : ℕ) : Matrix (Fin ℓ → Bool) (Fin ℓ → Bool) ℚ :=
  fun x y => ∏ i, (if x i = y i then 0 else 1)

theorem Xtensor_apply (ℓ : ℕ) (x y : Fin ℓ → Bool) :
    Xtensor ℓ x y = if y = (fun i => !x i) then 1 else 0 := by
  unfold Xtensor
  split_ifs with h
  · subst h
    simp
  · have : ∃ i, x i = y i := by
      by_contra hne
      push Not at hne
      apply h
      funext i
      have := hne i
      cases hx : x i <;> cases hy : y i <;> simp_all
    obtain ⟨i, hi⟩ := this
    exact Finset.prod_eq_zero (Finset.mem_univ i) (by simp [hi])

/-- Number of points of the cube with prescribed values on `S`. -/
theorem card_subcube (ℓ : ℕ) (S : Finset (Fin ℓ)) (p : Fin ℓ → Bool) :
    Fintype.card {x : Fin ℓ → Bool // ∀ i ∈ S, x i = p i} = 2 ^ (ℓ - S.card) := by
  rw [Fintype.card_subtype]
  have h : (Finset.univ.filter fun x : Fin ℓ → Bool => ∀ i ∈ S, x i = p i)
      = Fintype.piFinset (fun i => if i ∈ S then ({p i} : Finset Bool) else Finset.univ) := by
    ext x
    simp only [Finset.mem_filter, Finset.mem_univ, true_and, Fintype.mem_piFinset]
    constructor
    · intro hx i
      by_cases hi : i ∈ S
      · simp [hi, hx i hi]
      · simp [hi]
    · intro hx i hi
      have := hx i
      simpa [hi] using this
  rw [h, Fintype.card_piFinset]
  have h2 : ∀ i : Fin ℓ, (if i ∈ S then ({p i} : Finset Bool) else Finset.univ).card
      = if i ∈ S then 1 else 2 := by
    intro i
    split_ifs <;> simp
  simp only [h2]
  rw [Finset.prod_ite, Finset.prod_const_one, one_mul, Finset.prod_const]
  congr 1
  rw [Finset.filter_not, Finset.filter_mem_eq_inter, Finset.univ_inter,
    Finset.card_sdiff_of_subset (Finset.subset_univ S), Finset.card_univ, Fintype.card_fin]

/-- **G6d determinant lemma** (general, non-commutative `J`).
`J(t)` is a `2^ℓ × 2^ℓ` polynomial matrix with entries of degree `≤ d` and `J(0) = X^{⊗ℓ}`.
`T` is a finite set of (distinct) points, each `t ∈ T` with a direction `dir t ∈ S`, and
`J(t) = K_t + E_t` with `K_t` weakly `dir t`-nice (`K_t[x,y] = 0` when `x` and `y` differ in
coordinate `dir t`). Then `∑_{t ∈ T} rank E_t ≥ (|T| - d) * 2^(ℓ - |S|)`. -/
theorem g6d_determinant_lemma (ℓ d : ℕ) (J : Matrix (Fin ℓ → Bool) (Fin ℓ → Bool) ℚ[X])
    (hd : ∀ x y, (J x y).natDegree ≤ d) (h0 : J.map (eval 0) = Xtensor ℓ)
    (S : Finset (Fin ℓ)) (T : Finset ℚ) (dir : ℚ → Fin ℓ) (hdir : ∀ t ∈ T, dir t ∈ S)
    (K E : ℚ → Matrix (Fin ℓ → Bool) (Fin ℓ → Bool) ℚ)
    (hdec : ∀ t ∈ T, J.map (eval t) = K t + E t)
    (hnice : ∀ t ∈ T, ∀ x y, x (dir t) ≠ y (dir t) → K t x y = 0) :
    (T.card - d) * 2 ^ (ℓ - S.card) ≤ ∑ t ∈ T, (E t).rank := by
  classical
  set A := {x : Fin ℓ → Bool // ∀ i ∈ S, x i = false} with hA
  have hcard : Fintype.card A = 2 ^ (ℓ - S.card) := card_subcube ℓ S (fun _ => false)
  rw [← hcard]
  let e : Fin (Fintype.card A) ≃ A := (Fintype.equivFin A).symm
  let r : Fin (Fintype.card A) → (Fin ℓ → Bool) := fun i => (e i).1
  let c : Fin (Fintype.card A) → (Fin ℓ → Bool) := fun i => fun j => !(r i j)
  have hr : Function.Injective r := fun i j h => e.injective (Subtype.ext h)
  refine sum_rank_error_ge J d hd r c ?_ T K E hdec ?_
  · have : (J.submatrix r c).map (eval 0) = 1 := by
      ext i j
      have h := congrFun (congrFun h0 (r i)) (c j)
      simp only [Matrix.map_apply] at h
      simp only [Matrix.map_apply, Matrix.submatrix_apply, h, Xtensor_apply, Matrix.one_apply]
      congr 1
      apply propext
      constructor
      · intro hc
        apply hr
        funext k
        have := congrFun hc k
        simp only [c] at this
        cases h1 : r i k <;> cases h2 : r j k <;> simp_all
      · rintro rfl
        rfl
    rw [this, det_one]
    exact one_ne_zero
  · intro t ht i j
    apply hnice t ht
    have h1 : r i (dir t) = false := (e i).2 _ (hdir t ht)
    have h2 : r j (dir t) = false := (e j).2 _ (hdir t ht)
    simp [c, h1, h2]

end Cube

end PerPoint

#print axioms PerPoint.sum_nullity_le_of_det_ne_zero
#print axioms PerPoint.sum_rank_ge
#print axioms PerPoint.sum_rank_error_ge
#print axioms PerPoint.g6d_determinant_lemma
