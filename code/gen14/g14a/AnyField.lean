import Mathlib

/-!
# G14a: the G11b dual-polynomial bound over an ARBITRARY field

Same statement as `SqrtLog.sum_rank_ge_N`, with `ℚ` replaced by any field `F`.
The positivity step `∑ b_x(0)^2 > 0` is replaced by a free right test vector
`c = e_{¬x₁} * pref x₁` (only the LEFT vector has to solve `a(τ k) E k = 0`).
-/

open Polynomial Matrix Module

namespace AnyField

variable {ℓ : ℕ} {F : Type*} [Field F]

def weight (x : Fin ℓ → Bool) : ℕ := (Finset.univ.filter fun k => x k = true).card

def N (ℓ m : ℕ) : ℕ := ∑ x : Fin ℓ → Bool, (m + 1 - weight x)

def compl (x : Fin ℓ → Bool) : Fin ℓ → Bool := fun i => !x i

theorem compl_inj {x y : Fin ℓ → Bool} (h : compl x = compl y) : x = y := by
  funext i
  have := congrFun h i
  simpa [compl] using this

noncomputable def pref (τ : Fin ℓ → F) (x : Fin ℓ → Bool) : F[X] :=
  ∏ k ∈ Finset.univ.filter (fun k => x k = true), (1 - C (1 / τ k) * X)

theorem pref_eval_zero (τ : Fin ℓ → F) (x : Fin ℓ → Bool) : (pref τ x).eval 0 = 1 := by
  simp [pref, eval_prod]

theorem pref_ne_zero (τ : Fin ℓ → F) (x : Fin ℓ → Bool) : pref τ x ≠ 0 := by
  intro h
  have := pref_eval_zero τ x
  rw [h, eval_zero] at this
  exact zero_ne_one this

theorem pref_eval (τ : Fin ℓ → F) (x : Fin ℓ → Bool) (k : Fin ℓ) (hk : τ k ≠ 0)
    (hx : x k = true) : (pref τ x).eval (τ k) = 0 := by
  rw [pref, eval_prod]
  refine Finset.prod_eq_zero (i := k) (by simp [hx]) ?_
  simp [hk]

theorem natDegree_pref_le (τ : Fin ℓ → F) (x : Fin ℓ → Bool) :
    (pref τ x).natDegree ≤ weight x := by
  rw [pref]
  refine (natDegree_prod_le _ _).trans ?_
  have : ∀ k ∈ Finset.univ.filter (fun k => x k = true),
      (1 - C (1 / τ k) * X : F[X]).natDegree ≤ 1 := by
    intro k _
    compute_degree
  simpa [weight] using Finset.sum_le_card_nsmul _ _ 1 this

/-- A nonzero polynomial of degree `≤ m` vanishing at `τ k` for all `k` with `x k = true`
forces `weight x ≤ m`. -/
theorem weight_le_of_ne_zero (τ : Fin ℓ → F) (hτ : Function.Injective τ) (x : Fin ℓ → Bool)
    (p : F[X]) (hp : p ≠ 0) (m : ℕ) (hdeg : p.natDegree ≤ m)
    (hvan : ∀ k, x k = true → p.eval (τ k) = 0) : weight x ≤ m := by
  classical
  have hsub : (Finset.univ.filter fun k => x k = true).image τ ⊆ p.roots.toFinset := by
    intro t ht
    obtain ⟨k, hk, rfl⟩ := Finset.mem_image.mp ht
    simp only [Finset.mem_filter, Finset.mem_univ, true_and] at hk
    rw [Multiset.mem_toFinset, mem_roots hp]
    exact hvan k hk
  have h1 := Finset.card_le_card hsub
  rw [Finset.card_image_of_injective _ hτ] at h1
  have h2 := Multiset.toFinset_card_le p.roots
  have h3 := p.card_roots'
  unfold weight
  omega

/-- **Core lemma over any field.** -/
theorem core (ℓ d m : ℕ) (J : Matrix (Fin ℓ → Bool) (Fin ℓ → Bool) F[X])
    (hd : ∀ x y, (J x y).natDegree ≤ d)
    (h0 : ∀ x y, (J x y).eval 0 = if y = compl x then 1 else 0)
    (τ : Fin ℓ → F) (hτ : Function.Injective τ) (hτ0 : ∀ k, τ k ≠ 0)
    (K E : Fin ℓ → Matrix (Fin ℓ → Bool) (Fin ℓ → Bool) F)
    (hdec : ∀ k x y, (J x y).eval (τ k) = K k x y + E k x y)
    (hnice : ∀ k x y, x k ≠ y k → K k x y = 0)
    (hm : 2 * m + d + 1 ≤ ℓ)
    (a : (Fin ℓ → Bool) → F[X]) (hdeg : ∀ x, (a x).natDegree ≤ m)
    (hvan : ∀ k x, x k = true → (a x).eval (τ k) = 0)
    (hker : ∀ k y, ∑ x, E k x y * (a x).eval (τ k) = 0) : ∀ x, a x = 0 := by
  classical
  by_contra hne
  push Not at hne
  obtain ⟨x₀, hx₀⟩ := hne
  -- factor out the largest common power of `X`.
  have hP : ∃ n, ∃ x, ¬ (X : F[X]) ^ n ∣ a x := by
    refine ⟨(a x₀).natDegree + 1, x₀, fun h => ?_⟩
    have := natDegree_le_of_dvd h hx₀
    rw [natDegree_X_pow] at this
    omega
  obtain ⟨v, hv⟩ : ∃ v, Nat.find hP = v + 1 := by
    refine Nat.exists_eq_succ_of_ne_zero fun h => ?_
    obtain ⟨x, hx⟩ := Nat.find_spec hP
    rw [h] at hx
    exact hx (by simp)
  have hall : ∀ x, (X : F[X]) ^ v ∣ a x := by
    have := Nat.find_min hP (m := v) (by omega)
    push Not at this
    exact this
  obtain ⟨x₁, hx₁⟩ : ∃ x, ¬ (X : F[X]) ^ (v + 1) ∣ a x := by
    have := Nat.find_spec hP
    rwa [hv] at this
  choose b hb using hall
  have hb1 : (b x₁).eval 0 ≠ 0 := by
    intro h
    apply hx₁
    rw [← coeff_zero_eq_eval_zero, ← X_dvd_iff] at h
    rw [hb x₁, pow_succ]
    exact mul_dvd_mul_left _ h
  have hax₁ : a x₁ ≠ 0 := by
    intro h
    apply hx₁
    rw [h]
    exact dvd_zero _
  have hw : weight x₁ ≤ m :=
    weight_le_of_ne_zero τ hτ x₁ (a x₁) hax₁ m (hdeg x₁) (fun k hk => hvan k x₁ hk)
  -- the test polynomial, with the free right vector `c = e_{¬x₁} * pref x₁`
  set Φ : F[X] := ∑ x, a x * J x (compl x₁) * pref τ x₁ with hΦ
  have hΦdeg : Φ.natDegree ≤ 2 * m + d := by
    refine natDegree_sum_le_of_forall_le _ _ fun x _ => ?_
    have := natDegree_mul_le_of_le (natDegree_mul_le_of_le (hdeg x) (hd x (compl x₁)))
      (natDegree_pref_le τ x₁)
    omega
  have hΦroot : ∀ k, Φ.eval (τ k) = 0 := by
    intro k
    rw [hΦ]
    simp only [eval_finsetSum, eval_mul, hdec]
    by_cases hx : x₁ k = true
    · simp [pref_eval τ x₁ k (hτ0 k) hx]
    · have hK : ∀ x, (a x).eval (τ k) * K k x (compl x₁) = 0 := by
        intro x
        by_cases hxk : x k = true
        · rw [hvan k x hxk]; ring
        · rw [hnice k x (compl x₁) (by simp [compl]; simp at hx hxk; simp [hx, hxk])]; ring
      have hE := hker k (compl x₁)
      have : ∑ x, (a x).eval (τ k) * (K k x (compl x₁) + E k x (compl x₁)) * (pref τ x₁).eval (τ k)
          = (∑ x, E k x (compl x₁) * (a x).eval (τ k)) * (pref τ x₁).eval (τ k) := by
        rw [Finset.sum_mul]
        refine Finset.sum_congr rfl fun x _ => ?_
        have := hK x
        linear_combination (pref τ x₁).eval (τ k) * this
      rw [this, hE, zero_mul]
  have hΦ0 : Φ = 0 :=
    eq_zero_of_natDegree_lt_card_of_eval_eq_zero Φ hτ hΦroot (by simp; omega)
  set Ψ : F[X] := ∑ x, b x * J x (compl x₁) * pref τ x₁ with hΨ
  have hfac : Φ = X ^ v * Ψ := by
    rw [hΦ, hΨ, Finset.mul_sum]
    refine Finset.sum_congr rfl fun x _ => ?_
    rw [hb x]
    ring
  have hΨ0 : Ψ = 0 := by
    rw [hΦ0] at hfac
    rcases mul_eq_zero.mp hfac.symm with h | h
    · exact absurd h (pow_ne_zero _ X_ne_zero)
    · exact h
  have hval : Ψ.eval 0 = (b x₁).eval 0 := by
    rw [hΨ]
    simp only [eval_finsetSum, eval_mul, h0, pref_eval_zero, mul_one]
    rw [Finset.sum_eq_single x₁]
    · simp
    · intro x _ hx
      have : compl x₁ ≠ compl x := fun h => hx (compl_inj h).symm
      simp [this]
    · simp
  rw [hΨ0, eval_zero] at hval
  exact hb1 hval.symm

/-! ### The dual space `A_m` -/

noncomputable def rpoly (m : ℕ) (p : ∀ x : Fin ℓ → Bool, Fin (m + 1 - weight x) → F)
    (x : Fin ℓ → Bool) : F[X] :=
  ∑ i : Fin (m + 1 - weight x), C (p x i) * X ^ (i : ℕ)

noncomputable def avec (τ : Fin ℓ → F) (m : ℕ)
    (p : ∀ x : Fin ℓ → Bool, Fin (m + 1 - weight x) → F) (x : Fin ℓ → Bool) : F[X] :=
  pref τ x * rpoly m p x

theorem natDegree_rpoly_le (m : ℕ) (p : ∀ x : Fin ℓ → Bool, Fin (m + 1 - weight x) → F)
    (x : Fin ℓ → Bool) : (rpoly m p x).natDegree ≤ m - weight x := by
  refine natDegree_sum_le_of_forall_le _ _ fun i _ => ?_
  have := natDegree_C_mul_X_pow_le (p x i) (i : ℕ)
  have := i.2
  omega

theorem natDegree_avec_le (τ : Fin ℓ → F) (m : ℕ)
    (p : ∀ x : Fin ℓ → Bool, Fin (m + 1 - weight x) → F) (x : Fin ℓ → Bool) :
    (avec τ m p x).natDegree ≤ m := by
  by_cases h : m + 1 - weight x = 0
  · have : rpoly m p x = 0 := Finset.sum_eq_zero fun i _ => absurd i.2 (by omega)
    simp [avec, this]
  · have := natDegree_mul_le_of_le (natDegree_pref_le τ x) (natDegree_rpoly_le m p x)
    rw [avec]
    omega

theorem coeff_rpoly (m : ℕ) (p : ∀ x : Fin ℓ → Bool, Fin (m + 1 - weight x) → F)
    (x : Fin ℓ → Bool) (i : Fin (m + 1 - weight x)) : (rpoly m p x).coeff i = p x i := by
  simp only [rpoly, finsetSum_coeff, coeff_C_mul_X_pow]
  rw [Finset.sum_eq_single i]
  · simp
  · intro j _ hj
    have : (i : ℕ) ≠ (j : ℕ) := fun h => hj (Fin.ext h.symm)
    simp [this]
  · simp

noncomputable def evA (τ : Fin ℓ → F) (m : ℕ) (k : Fin ℓ) :
    (∀ x : Fin ℓ → Bool, Fin (m + 1 - weight x) → F) →ₗ[F] ((Fin ℓ → Bool) → F) where
  toFun p x := (avec τ m p x).eval (τ k)
  map_add' p q := by
    funext x
    simp [avec, rpoly, eval_finsetSum, add_mul, Finset.sum_add_distrib, mul_add]
  map_smul' c p := by
    funext x
    simp [avec, rpoly, eval_finsetSum, Finset.mul_sum]
    refine Finset.sum_congr rfl fun i _ => ?_
    ring

theorem finrank_params (ℓ m : ℕ) :
    finrank F (∀ x : Fin ℓ → Bool, Fin (m + 1 - weight x) → F) = N ℓ m := by
  rw [Module.finrank_pi_fintype]
  simp [N]

/-- **G11b theorem over an arbitrary field `F`.** -/
theorem sum_rank_ge_N (ℓ d : ℕ) (J : Matrix (Fin ℓ → Bool) (Fin ℓ → Bool) F[X])
    (hd : ∀ x y, (J x y).natDegree ≤ d)
    (h0 : ∀ x y, (J x y).eval 0 = if y = compl x then 1 else 0)
    (τ : Fin ℓ → F) (hτ : Function.Injective τ) (hτ0 : ∀ k, τ k ≠ 0)
    (K E : Fin ℓ → Matrix (Fin ℓ → Bool) (Fin ℓ → Bool) F)
    (hdec : ∀ k x y, (J x y).eval (τ k) = K k x y + E k x y)
    (hnice : ∀ k x y, x k ≠ y k → K k x y = 0)
    (m : ℕ) (hm : 2 * m + d + 1 ≤ ℓ) :
    N ℓ m ≤ ∑ k, (E k).rank := by
  classical
  by_contra hlt
  push Not at hlt
  let L : (∀ x : Fin ℓ → Bool, Fin (m + 1 - weight x) → F) →ₗ[F]
      (∀ k : Fin ℓ, LinearMap.range ((E k)ᵀ.mulVecLin)) :=
    LinearMap.pi fun k => ((E k)ᵀ.mulVecLin.rangeRestrict).comp (evA τ m k)
  have hfin : finrank F (∀ k : Fin ℓ, LinearMap.range ((E k)ᵀ.mulVecLin))
      < finrank F (∀ x : Fin ℓ → Bool, Fin (m + 1 - weight x) → F) := by
    rw [finrank_params, Module.finrank_pi_fintype]
    simpa [finrank_range_mulVecLin, rank_transpose] using hlt
  obtain ⟨p, hpker, hp⟩ := (Submodule.ne_bot_iff _).mp (LinearMap.ker_ne_bot_of_finrank_lt (f := L) hfin)
  have hker0 : ∀ k, (E k)ᵀ *ᵥ (fun x => (avec τ m p x).eval (τ k)) = 0 := by
    intro k
    have := congrArg Subtype.val (congrFun (LinearMap.mem_ker.mp hpker) k)
    simpa [L, evA, Matrix.mulVec_transpose] using this
  have hker : ∀ k y, ∑ x, E k x y * (avec τ m p x).eval (τ k) = 0 := by
    intro k y
    have := congrFun (hker0 k) y
    simpa [Matrix.mulVec, dotProduct] using this
  have hzero := core ℓ d m J hd h0 τ hτ hτ0 K E hdec hnice hm (avec τ m p)
    (natDegree_avec_le τ m p)
    (fun k x hx => by rw [avec, eval_mul, pref_eval τ x k (hτ0 k) hx, zero_mul]) hker
  apply hp
  funext x i
  have hr : rpoly m p x = 0 := by
    rcases mul_eq_zero.mp (hzero x) with h | h
    · exact absurd h (pref_ne_zero τ x)
    · exact h
  have := coeff_rpoly m p x i
  rw [hr, coeff_zero] at this
  simpa using this.symm

end AnyField

#print axioms AnyField.core
#print axioms AnyField.sum_rank_ge_N
