import Mathlib

/-!
# Witness determinant (Lemma 1 of `note_limits_s9.md`, algebraic core)
-/

open Polynomial Matrix Module

namespace WitnessDet

variable {n : ℕ}

/-- Entrywise evaluation of a polynomial matrix at `t`. -/
noncomputable def evalMat (Φ : Matrix (Fin n) (Fin n) ℚ[X]) (t : ℚ) : Matrix (Fin n) (Fin n) ℚ :=
  Φ.map (eval t)

/-- Nullity: dimension of the (right) kernel `{x | A *ᵥ x = 0}`. -/
noncomputable def nullity (A : Matrix (Fin n) (Fin n) ℚ) : ℕ :=
  Module.finrank ℚ (LinearMap.ker A.mulVecLin)

theorem eval_det (Φ : Matrix (Fin n) (Fin n) ℚ[X]) (t : ℚ) :
    eval t Φ.det = (evalMat Φ t).det := by
  have := RingHom.map_det (Polynomial.evalRingHom t) Φ
  simpa [evalMat, RingHom.mapMatrix_apply] using this

theorem det_ne_zero_of_eval_zero (Φ : Matrix (Fin n) (Fin n) ℚ[X]) (h0 : evalMat Φ 0 = 1) :
    Φ.det ≠ 0 := by
  intro h
  have := eval_det Φ 0
  rw [h, h0, det_one, eval_zero] at this
  exact zero_ne_one this

theorem natDegree_det_le (Φ : Matrix (Fin n) (Fin n) ℚ[X]) (R : ℕ)
    (hR : ∀ i j, (Φ i j).natDegree ≤ R) : Φ.det.natDegree ≤ n * R := by
  rw [det_apply']
  apply natDegree_sum_le_of_forall_le
  intro σ _
  refine (natDegree_mul_le).trans ?_
  rw [natDegree_intCast, zero_add]
  refine (natDegree_prod_le _ _).trans ?_
  refine (Finset.sum_le_card_nsmul _ _ R (fun i _ => hR _ _)).trans ?_
  simp

theorem pow_dvd_det_of_nullity (Φ : Matrix (Fin n) (Fin n) ℚ[X]) (t : ℚ) :
    (X - C t) ^ nullity (evalMat Φ t) ∣ Φ.det := by
  set A := evalMat Φ t with hA
  set K := LinearMap.ker A.mulVecLin with hK
  obtain ⟨W, hW⟩ := Submodule.exists_isCompl K
  set k := finrank ℚ K
  let bK := Module.finBasis ℚ K
  let bW := Module.finBasis ℚ W
  let b : Basis (Fin k ⊕ Fin (finrank ℚ W)) ℚ (Fin n → ℚ) :=
    (bK.prod bW).map (Submodule.prodEquivOfIsCompl K W hW)
  have hbK : ∀ i, b (Sum.inl i) = (bK i : Fin n → ℚ) := by
    intro i
    simp [b, Basis.map_apply, Submodule.coe_prodEquivOfIsCompl', Basis.prod_apply]
  have hcard : Fintype.card (Fin k ⊕ Fin (finrank ℚ W)) = n := by
    rw [← Module.finrank_eq_card_basis b]; simp
  let e := (Fintype.equivFinOfCardEq hcard).symm
  let b' : Basis (Fin n) ℚ (Fin n → ℚ) := b.reindex e.symm
  have hb' : ∀ j, b' j = b (e j) := by intro j; simp [b']
  let U : Matrix (Fin n) (Fin n) ℚ := (Pi.basisFun ℚ (Fin n)).toMatrix b'
  have hU : ∀ i j, U i j = b' j i := by
    intro i j; simp [U, Basis.toMatrix_apply]
  have hUdet : IsUnit U.det := by
    have h1 := (Pi.basisFun ℚ (Fin n)).toMatrix_mul_toMatrix_flip b'
    have h2 := congrArg det h1
    rw [det_mul, det_one] at h2
    exact IsUnit.of_mul_eq_one _ h2
  -- the product matrix
  set P := Φ * U.map C with hP
  have hroot : ∀ l j, (e j).isLeft → (P l j).IsRoot t := by
    intro l j hj
    obtain ⟨i, hi⟩ : ∃ i, e j = Sum.inl i := by
      cases h : e j with
      | inl i => exact ⟨i, rfl⟩
      | inr i => rw [h] at hj; simp at hj
    have hker : A *ᵥ (bK i : Fin n → ℚ) = 0 := by
      have : A.mulVecLin (bK i : Fin n → ℚ) = 0 := (bK i).2
      rwa [Matrix.mulVecLin_apply] at this
    have hcol : ∀ m, U m j = (bK i : Fin n → ℚ) m := by
      intro m; rw [hU, hb', hi, hbK]
    have := congrFun hker l
    simp only [Matrix.mulVec, dotProduct, Pi.zero_apply] at this
    simp only [IsRoot, hP, Matrix.mul_apply, Matrix.map_apply, eval_finsetSum, eval_mul, eval_C]
    rw [← this]
    refine Finset.sum_congr rfl (fun m _ => ?_)
    rw [hcol, hA, evalMat, Matrix.map_apply]
  let d : Fin n → ℚ[X] := fun j => if (e j).isLeft then X - C t else 1
  let M : Matrix (Fin n) (Fin n) ℚ[X] := fun l j => if (e j).isLeft then P l j /ₘ (X - C t) else P l j
  have hPM : P = M * diagonal d := by
    refine Matrix.ext fun l j => ?_
    rw [mul_diagonal]
    by_cases hj : (e j).isLeft
    · simp only [M, d, hj, ite_true]
      rw [mul_comm]
      exact ((mul_divByMonic_eq_iff_isRoot).2 (hroot l j hj)).symm
    · simp [M, d, hj]
  have hprod : ∏ j, d j = (X - C t) ^ k := by
    rw [← e.symm.prod_comp]
    simp only [d, Equiv.apply_symm_apply]
    rw [Fintype.prod_sum_type]
    simp
  have hdet : Φ.det * C U.det = M.det * (X - C t) ^ k := by
    have := congrArg det hPM
    rw [hP, det_mul, det_mul, det_diagonal, hprod] at this
    rw [← this]
    congr 1
    have := RingHom.map_det (C : ℚ →+* ℚ[X]) U
    rw [this]; rfl
  have : (X - C t) ^ k ∣ Φ.det * C U.det := ⟨M.det, by rw [hdet, mul_comm]⟩
  exact (IsUnit.dvd_mul_right (hUdet.map C)).1 this

theorem sum_nullity_le (Φ : Matrix (Fin n) (Fin n) ℚ[X]) (R : ℕ)
    (hR : ∀ i j, (Φ i j).natDegree ≤ R) (h0 : evalMat Φ 0 = 1) (T : Finset ℚ) :
    ∑ t ∈ T, nullity (evalMat Φ t) ≤ n * R := by
  have hdvd : ∏ t ∈ T, (X - C t) ^ nullity (evalMat Φ t) ∣ Φ.det := by
    apply Finset.prod_dvd_of_coprime
    · intro a _ b _ hab
      exact ((pairwise_coprime_X_sub_C (s := id) Function.injective_id) hab).pow
    · intro t _
      exact pow_dvd_det_of_nullity Φ t
  have hdeg := natDegree_le_of_dvd hdvd (det_ne_zero_of_eval_zero Φ h0)
  rw [natDegree_prod_of_monic _ _ (fun t _ => (monic_X_sub_C t).pow _)] at hdeg
  simp only [natDegree_pow, natDegree_X_sub_C, mul_one] at hdeg
  exact hdeg.trans (natDegree_det_le Φ R hR)

/-! ## Part 2: rank bound at a point -/

theorem rank_add_le' {m : ℕ} (A B : Matrix (Fin m) (Fin n) ℚ) : (A + B).rank ≤ A.rank + B.rank := by
  rw [← finrank_range_mulVecLin, ← finrank_range_mulVecLin, ← finrank_range_mulVecLin,
    mulVecLin_add]
  exact (Submodule.finrank_mono (LinearMap.range_add_le _ _)).trans
    (Submodule.finrank_add_le_finrank_add_finrank _ _)

theorem rank_sum_le {ι : Type*} (s : Finset ι) (B : ι → Matrix (Fin n) (Fin n) ℚ) :
    (∑ i ∈ s, B i).rank ≤ ∑ i ∈ s, (B i).rank := by
  classical
  induction s using Finset.induction_on with
  | empty => simp [rank_zero]
  | insert a s ha ih =>
    rw [Finset.sum_insert ha, Finset.sum_insert ha]
    exact (rank_add_le' _ _).trans (Nat.add_le_add_left ih _)

theorem nullity_add_rank (A : Matrix (Fin n) (Fin n) ℚ) : A.rank + nullity A = n := by
  rw [← finrank_range_mulVecLin, nullity, LinearMap.finrank_range_add_finrank_ker]
  simp

/-- General form (no support hypothesis on `Rγ` needed). -/
theorem deficiency_le_nullity_add_general {G : Type*} [Fintype G] [DecidableEq G]
    (A : Matrix (Fin n) (Fin n) ℚ) (g h : Fin n → G) (Rγ : G → Submodule ℚ (Fin n → ℚ))
    (hrow : ∀ γ v, h v ≠ γ → (fun j => if g j = γ then A v j else 0) ∈ Rγ γ) :
    ∑ γ, ((Finset.univ.filter (fun j => g j = γ)).card
          - (Finset.univ.filter (fun v => h v = γ)).card)
      ≤ nullity A + ∑ γ, finrank ℚ (Rγ γ) := by
  set a : G → ℕ := fun γ => (Finset.univ.filter (fun j => g j = γ)).card with ha
  set I : G → ℕ := fun γ => (Finset.univ.filter (fun v => h v = γ)).card with hI
  set r : G → ℕ := fun γ => finrank ℚ (Rγ γ) with hr
  let D : G → Matrix (Fin n) (Fin n) ℚ := fun γ => diagonal (fun j => if g j = γ then 1 else 0)
  let B : G → Matrix (Fin n) (Fin n) ℚ := fun γ => A * D γ
  have hBapp : ∀ γ v j, B γ v j = if g j = γ then A v j else 0 := by
    intro γ v j; simp [B, D, mul_diagonal]
  have hsumD : ∑ γ, D γ = 1 := by
    ext i j
    simp only [D, Matrix.sum_apply, diagonal_apply]
    by_cases hij : i = j
    · subst hij; simp
    · simp [hij]
  have hsumB : ∑ γ, B γ = A := by
    simp only [B, ← Matrix.mul_sum, hsumD, Matrix.mul_one]
  -- rank of each block ≤ a γ
  have hBa : ∀ γ, (B γ).rank ≤ a γ := by
    intro γ
    refine (rank_mul_le_right _ _).trans ?_
    rw [rank_diagonal]
    simp only [ha, ne_eq, ite_eq_right_iff, one_ne_zero, imp_false, not_not]
    rw [Fintype.card_subtype]
  -- rank of each block ≤ I γ + r γ
  have hBI : ∀ γ, (B γ).rank ≤ I γ + r γ := by
    intro γ
    rw [← finrank_span_range_row]
    let f : {v // h v = γ} → (Fin n → ℚ) := fun v => (B γ).row v
    have hle : Submodule.span ℚ (Set.range (B γ).row) ≤
        Submodule.span ℚ (Set.range f) ⊔ Rγ γ := by
      rw [Submodule.span_le]
      rintro _ ⟨v, rfl⟩
      by_cases hv : h v = γ
      · exact Submodule.mem_sup_left (Submodule.subset_span ⟨⟨v, hv⟩, rfl⟩)
      · apply Submodule.mem_sup_right
        have : (B γ).row v = fun j => if g j = γ then A v j else 0 := by
          funext j; exact hBapp γ v j
        rw [this]; exact hrow γ v hv
    refine (Submodule.finrank_mono hle).trans
      ((Submodule.finrank_add_le_finrank_add_finrank _ _).trans ?_)
    refine Nat.add_le_add_right ?_ _
    have := finrank_range_le_card (R := ℚ) f
    rw [Fintype.card_subtype] at this
    exact this
  have hrank : A.rank ≤ ∑ γ, min (a γ) (I γ + r γ) := by
    rw [← hsumB]
    exact (rank_sum_le _ _).trans (Finset.sum_le_sum fun γ _ => le_min (hBa γ) (hBI γ))
  have hn : ∑ γ, a γ = n := by
    have := Finset.card_eq_sum_card_fiberwise (f := g) (s := Finset.univ) (t := Finset.univ)
      (fun _ _ => Finset.mem_coe.2 (Finset.mem_univ _))
    simpa [ha] using this.symm
  have hnr := nullity_add_rank A
  have hpt : ∀ γ, a γ - I γ ≤ (a γ - min (a γ) (I γ + r γ)) + r γ := by
    intro γ; omega
  calc ∑ γ, (a γ - I γ) ≤ ∑ γ, ((a γ - min (a γ) (I γ + r γ)) + r γ) :=
        Finset.sum_le_sum fun γ _ => hpt γ
    _ = (∑ γ, a γ - ∑ γ, min (a γ) (I γ + r γ)) + ∑ γ, r γ := by
        rw [Finset.sum_add_distrib, Finset.sum_tsub_distrib _ (fun γ _ => min_le_left _ _)]
    _ ≤ nullity A + ∑ γ, r γ := by
        refine Nat.add_le_add_right ?_ _
        omega

/-- Faithful statement of the rank bound at a point (dim ker ≥ def − ρ).
The support hypothesis `hsupp` is part of the setting but is not needed for the proof. -/
theorem deficiency_le_nullity_add {G : Type*} [Fintype G] [DecidableEq G]
    (A : Matrix (Fin n) (Fin n) ℚ) (g h : Fin n → G) (Rγ : G → Submodule ℚ (Fin n → ℚ))
    (_hsupp : ∀ γ, ∀ x ∈ Rγ γ, ∀ j, g j ≠ γ → x j = 0)
    (hrow : ∀ γ v, h v ≠ γ → (fun j => if g j = γ then A v j else 0) ∈ Rγ γ) :
    ∑ γ, ((Finset.univ.filter (fun j => g j = γ)).card
          - (Finset.univ.filter (fun v => h v = γ)).card)
      ≤ nullity A + ∑ γ, finrank ℚ (Rγ γ) :=
  deficiency_le_nullity_add_general A g h Rγ hrow

end WitnessDet

#print axioms WitnessDet.det_ne_zero_of_eval_zero
#print axioms WitnessDet.natDegree_det_le
#print axioms WitnessDet.pow_dvd_det_of_nullity
#print axioms WitnessDet.sum_nullity_le
#print axioms WitnessDet.deficiency_le_nullity_add
