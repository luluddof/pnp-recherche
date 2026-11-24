import PnpLean.Article2.PerPoint

/-!
# Per-point errors: the dual polynomial bound (G11b, "sqrt(log n) loss is optimal")

Source: `code/gen11/g11b/notes.md` §1, one point per direction.

`J(t)` is a `2^ℓ × 2^ℓ` polynomial matrix with entries of degree `≤ d`, `J(0) = X^{⊗ℓ}`.
`τ : Fin ℓ → ℚ` are `ℓ` distinct nonzero points; `J(τ k) = K k + E k` with `K k` weakly `k`-nice
(`K k x y = 0` whenever `x k ≠ y k`). If `2 m + d + 1 ≤ ℓ` then

  `∑ k, rank (E k) ≥ N ℓ m := ∑ x, (m + 1 - |x|)`   (`= ∑_{j ≤ m} ∑_{i ≤ j} C(ℓ, i)`).

Proof (dual polynomial method): the space of row vectors
`a x = (∏_{k : x k = 1} (1 - X / τ k)) * r x`, `deg (a x) ≤ m`, has dimension `N ℓ m`. If
`∑ rank (E k) < N ℓ m` there is `a ≠ 0` in it with `a(τ k) ᵥ* E k = 0` for all `k`. Then
`Φ := ∑ x y, a x * J x y * a (¬y)` has degree `≤ 2m + d < ℓ` and vanishes at the `ℓ` points `τ k`,
so `Φ = 0`. But writing `a = X^v * b` with `b(0) ≠ 0`, `Φ = X^(2v) * Ψ` and
`Ψ(0) = ∑ x, (b x)(0)^2 > 0`.
-/

open Polynomial Matrix Module

namespace SqrtLog

open PerPoint

variable {ℓ : ℕ}

/-- Hamming weight of a point of the cube. -/
def weight (x : Fin ℓ → Bool) : ℕ := (Finset.univ.filter fun k => x k = true).card

/-- `N ℓ m = ∑_x max(0, m + 1 - |x|)`: the dimension of the dual space `A_m`. -/
def N (ℓ m : ℕ) : ℕ := ∑ x : Fin ℓ → Bool, (m + 1 - weight x)

/-- Complement of a point of the cube. -/
def compl (x : Fin ℓ → Bool) : Fin ℓ → Bool := fun i => !x i

@[simp] theorem compl_compl (x : Fin ℓ → Bool) : compl (compl x) = x := by
  funext i; simp [compl]

/-! ### The core (polynomial) argument -/

/-- **Core lemma.** A row vector `a` of polynomials of degree `≤ m` with `a x (τ k) = 0` whenever
`x k = true` and `a(τ k) ᵥ* E k = 0` (written `(E k)ᵀ *ᵥ a(τ k) = 0`) for all `k` is zero. -/
theorem core (ℓ d m : ℕ) (J : Matrix (Fin ℓ → Bool) (Fin ℓ → Bool) ℚ[X])
    (hd : ∀ x y, (J x y).natDegree ≤ d) (h0 : J.map (eval 0) = Xtensor ℓ)
    (τ : Fin ℓ → ℚ) (hτ : Function.Injective τ)
    (K E : Fin ℓ → Matrix (Fin ℓ → Bool) (Fin ℓ → Bool) ℚ)
    (hdec : ∀ k, J.map (eval (τ k)) = K k + E k)
    (hnice : ∀ k x y, x k ≠ y k → K k x y = 0)
    (hm : 2 * m + d + 1 ≤ ℓ)
    (a : (Fin ℓ → Bool) → ℚ[X]) (hdeg : ∀ x, (a x).natDegree ≤ m)
    (hvan : ∀ k x, x k = true → (a x).eval (τ k) = 0)
    (hker : ∀ k, (E k)ᵀ *ᵥ (fun x => (a x).eval (τ k)) = 0) : ∀ x, a x = 0 := by
  classical
  by_contra hne
  push Not at hne
  obtain ⟨x₀, hx₀⟩ := hne
  -- Step 1: `Φ = 0`.
  set Φ : ℚ[X] := ∑ x, ∑ y, a x * J x y * a (compl y) with hΦ
  have hΦdeg : Φ.natDegree ≤ 2 * m + d := by
    refine natDegree_sum_le_of_forall_le _ _ fun x _ => ?_
    refine natDegree_sum_le_of_forall_le _ _ fun y _ => ?_
    have := natDegree_mul_le_of_le (natDegree_mul_le_of_le (hdeg x) (hd x y)) (hdeg (compl y))
    omega
  have hΦroot : ∀ k, Φ.eval (τ k) = 0 := by
    intro k
    have hJ : ∀ x y, (J x y).eval (τ k) = K k x y + E k x y := by
      intro x y
      have := congrFun (congrFun (hdec k) x) y
      simpa using this
    have h1 : Φ.eval (τ k) = ∑ x, ∑ y, (a x).eval (τ k) * K k x y * (a (compl y)).eval (τ k)
        + ∑ y, (∑ x, E k x y * (a x).eval (τ k)) * (a (compl y)).eval (τ k) := by
      rw [hΦ]
      simp only [eval_finsetSum, eval_mul, hJ]
      have hE : ∑ y, (∑ x, E k x y * (a x).eval (τ k)) * (a (compl y)).eval (τ k)
          = ∑ x, ∑ y, (a x).eval (τ k) * E k x y * (a (compl y)).eval (τ k) := by
        simp only [Finset.sum_mul]
        rw [Finset.sum_comm]
        exact Finset.sum_congr rfl fun x _ => Finset.sum_congr rfl fun y _ => by ring
      rw [hE, ← Finset.sum_add_distrib]
      refine Finset.sum_congr rfl fun x _ => ?_
      rw [← Finset.sum_add_distrib]
      refine Finset.sum_congr rfl fun y _ => ?_
      ring
    rw [h1]
    have h2 : ∀ x y, (a x).eval (τ k) * K k x y * (a (compl y)).eval (τ k) = 0 := by
      intro x y
      by_cases hx : x k = true
      · rw [hvan k x hx]; ring
      · by_cases hy : y k = true
        · rw [hnice k x y (by rw [hy]; exact hx)]; ring
        · rw [hvan k (compl y) (by simp [compl, hy])]; ring
    have h3 : ∀ y, ∑ x, E k x y * (a x).eval (τ k) = 0 := by
      intro y
      have := congrFun (hker k) y
      simpa [Matrix.mulVec, dotProduct] using this
    simp [h2, h3]
  have hΦ0 : Φ = 0 :=
    eq_zero_of_natDegree_lt_card_of_eval_eq_zero Φ hτ hΦroot (by simp; omega)
  -- Step 2: factor out the largest common power of `X`.
  have hP : ∃ n, ∃ x, ¬ (X : ℚ[X]) ^ n ∣ a x := by
    refine ⟨(a x₀).natDegree + 1, x₀, fun h => ?_⟩
    have := natDegree_le_of_dvd h hx₀
    rw [natDegree_X_pow] at this
    omega
  obtain ⟨v, hv⟩ : ∃ v, Nat.find hP = v + 1 := by
    refine Nat.exists_eq_succ_of_ne_zero fun h => ?_
    obtain ⟨x, hx⟩ := Nat.find_spec hP
    rw [h] at hx
    exact hx (by simp)
  have hall : ∀ x, (X : ℚ[X]) ^ v ∣ a x := by
    have := Nat.find_min hP (m := v) (by omega)
    push Not at this
    exact this
  obtain ⟨x₁, hx₁⟩ : ∃ x, ¬ (X : ℚ[X]) ^ (v + 1) ∣ a x := by
    have := Nat.find_spec hP
    rwa [hv] at this
  choose b hb using hall
  have hb1 : (b x₁).eval 0 ≠ 0 := by
    intro h
    apply hx₁
    rw [← coeff_zero_eq_eval_zero, ← X_dvd_iff] at h
    rw [hb x₁, pow_succ]
    exact mul_dvd_mul_left _ h
  -- Step 3: `Φ = X^(2v) * Ψ`, `Ψ(0) = ∑ (b x)(0)^2`.
  set Ψ : ℚ[X] := ∑ x, ∑ y, b x * J x y * b (compl y) with hΨ
  have hfac : Φ = X ^ (2 * v) * Ψ := by
    rw [hΦ, hΨ, Finset.mul_sum]
    refine Finset.sum_congr rfl fun x _ => ?_
    rw [Finset.mul_sum]
    refine Finset.sum_congr rfl fun y _ => ?_
    rw [hb x, hb (compl y)]
    ring
  have hΨ0 : Ψ = 0 := by
    rw [hΦ0] at hfac
    rcases mul_eq_zero.mp hfac.symm with h | h
    · exact absurd h (pow_ne_zero _ X_ne_zero)
    · exact h
  have hJ0 : ∀ x y, (J x y).eval 0 = if y = compl x then 1 else 0 := by
    intro x y
    have := congrFun (congrFun h0 x) y
    rw [Xtensor_apply] at this
    exact this
  have hsq : Ψ.eval 0 = ∑ x, (b x).eval 0 * (b x).eval 0 := by
    rw [hΨ]
    simp only [eval_finsetSum, eval_mul, hJ0]
    refine Finset.sum_congr rfl fun x _ => ?_
    simp
  rw [hΨ0, eval_zero] at hsq
  have := (Finset.sum_eq_zero_iff_of_nonneg (fun x _ => mul_self_nonneg _)).mp hsq.symm x₁
    (Finset.mem_univ _)
  exact hb1 (mul_self_eq_zero.mp this)

/-! ### The dual space `A_m` and its parametrization -/

/-- The prefactor `∏_{k : x k = 1} (1 - X / τ k)`. -/
noncomputable def pref (τ : Fin ℓ → ℚ) (x : Fin ℓ → Bool) : ℚ[X] :=
  ∏ k ∈ Finset.univ.filter (fun k => x k = true), (1 - C (1 / τ k) * X)

/-- The free factor `r x`, of degree `≤ m - |x|`, given by its coefficients. -/
noncomputable def rpoly (m : ℕ) (p : ∀ x : Fin ℓ → Bool, Fin (m + 1 - weight x) → ℚ)
    (x : Fin ℓ → Bool) : ℚ[X] :=
  ∑ i : Fin (m + 1 - weight x), C (p x i) * X ^ (i : ℕ)

/-- The element of `A_m` with parameters `p`. -/
noncomputable def avec (τ : Fin ℓ → ℚ) (m : ℕ)
    (p : ∀ x : Fin ℓ → Bool, Fin (m + 1 - weight x) → ℚ) (x : Fin ℓ → Bool) : ℚ[X] :=
  pref τ x * rpoly m p x

theorem pref_eval_zero (τ : Fin ℓ → ℚ) (x : Fin ℓ → Bool) : (pref τ x).eval 0 = 1 := by
  simp [pref, eval_prod]

theorem pref_ne_zero (τ : Fin ℓ → ℚ) (x : Fin ℓ → Bool) : pref τ x ≠ 0 := by
  intro h
  have := pref_eval_zero τ x
  rw [h, eval_zero] at this
  exact zero_ne_one this

theorem pref_eval (τ : Fin ℓ → ℚ) (x : Fin ℓ → Bool) (k : Fin ℓ) (hk : τ k ≠ 0)
    (hx : x k = true) : (pref τ x).eval (τ k) = 0 := by
  rw [pref, eval_prod]
  refine Finset.prod_eq_zero (i := k) (by simp [hx]) ?_
  simp [hk]

theorem natDegree_pref_le (τ : Fin ℓ → ℚ) (x : Fin ℓ → Bool) :
    (pref τ x).natDegree ≤ weight x := by
  rw [pref]
  refine (natDegree_prod_le _ _).trans ?_
  have : ∀ k ∈ Finset.univ.filter (fun k => x k = true),
      (1 - C (1 / τ k) * X : ℚ[X]).natDegree ≤ 1 := by
    intro k _
    compute_degree
  simpa [weight] using Finset.sum_le_card_nsmul _ _ 1 this

theorem natDegree_rpoly_le (m : ℕ) (p : ∀ x : Fin ℓ → Bool, Fin (m + 1 - weight x) → ℚ)
    (x : Fin ℓ → Bool) : (rpoly m p x).natDegree ≤ m - weight x := by
  refine natDegree_sum_le_of_forall_le _ _ fun i _ => ?_
  have := natDegree_C_mul_X_pow_le (p x i) (i : ℕ)
  have := i.2
  omega

theorem natDegree_avec_le (τ : Fin ℓ → ℚ) (m : ℕ)
    (p : ∀ x : Fin ℓ → Bool, Fin (m + 1 - weight x) → ℚ) (x : Fin ℓ → Bool) :
    (avec τ m p x).natDegree ≤ m := by
  by_cases h : m + 1 - weight x = 0
  · have : rpoly m p x = 0 := Finset.sum_eq_zero fun i _ => absurd i.2 (by omega)
    simp [avec, this]
  · have := natDegree_mul_le_of_le (natDegree_pref_le τ x) (natDegree_rpoly_le m p x)
    rw [avec]
    omega

theorem coeff_rpoly (m : ℕ) (p : ∀ x : Fin ℓ → Bool, Fin (m + 1 - weight x) → ℚ)
    (x : Fin ℓ → Bool) (i : Fin (m + 1 - weight x)) : (rpoly m p x).coeff i = p x i := by
  simp only [rpoly, finsetSum_coeff, coeff_C_mul_X_pow]
  rw [Finset.sum_eq_single i]
  · simp
  · intro j _ hj
    have : (i : ℕ) ≠ (j : ℕ) := fun h => hj (Fin.ext h.symm)
    simp [this]
  · simp

/-- Evaluation of `avec` at `τ k`, as a linear map in the parameters. -/
noncomputable def evA (τ : Fin ℓ → ℚ) (m : ℕ) (k : Fin ℓ) :
    (∀ x : Fin ℓ → Bool, Fin (m + 1 - weight x) → ℚ) →ₗ[ℚ] ((Fin ℓ → Bool) → ℚ) where
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
    finrank ℚ (∀ x : Fin ℓ → Bool, Fin (m + 1 - weight x) → ℚ) = N ℓ m := by
  rw [Module.finrank_pi_fintype]
  simp [N]

/-! ### Main theorem -/

/-- **G11b theorem** (dual polynomial bound, one point per direction).
`J(t)` has entries of degree `≤ d`, `J(0) = X^{⊗ℓ}`; the `τ k` are distinct and nonzero;
`J(τ k) = K k + E k` with `K k` weakly `k`-nice. If `2 m + d + 1 ≤ ℓ` then
`∑ k, rank (E k) ≥ N ℓ m = ∑ x, (m + 1 - |x|)`. -/
theorem sum_rank_ge_N (ℓ d : ℕ) (J : Matrix (Fin ℓ → Bool) (Fin ℓ → Bool) ℚ[X])
    (hd : ∀ x y, (J x y).natDegree ≤ d) (h0 : J.map (eval 0) = Xtensor ℓ)
    (τ : Fin ℓ → ℚ) (hτ : Function.Injective τ) (hτ0 : ∀ k, τ k ≠ 0)
    (K E : Fin ℓ → Matrix (Fin ℓ → Bool) (Fin ℓ → Bool) ℚ)
    (hdec : ∀ k, J.map (eval (τ k)) = K k + E k)
    (hnice : ∀ k x y, x k ≠ y k → K k x y = 0)
    (m : ℕ) (hm : 2 * m + d + 1 ≤ ℓ) :
    N ℓ m ≤ ∑ k, (E k).rank := by
  classical
  by_contra hlt
  push Not at hlt
  let L : (∀ x : Fin ℓ → Bool, Fin (m + 1 - weight x) → ℚ) →ₗ[ℚ]
      (∀ k : Fin ℓ, LinearMap.range ((E k)ᵀ.mulVecLin)) :=
    LinearMap.pi fun k => ((E k)ᵀ.mulVecLin.rangeRestrict).comp (evA τ m k)
  have hfin : finrank ℚ (∀ k : Fin ℓ, LinearMap.range ((E k)ᵀ.mulVecLin))
      < finrank ℚ (∀ x : Fin ℓ → Bool, Fin (m + 1 - weight x) → ℚ) := by
    rw [finrank_params, Module.finrank_pi_fintype]
    simpa [finrank_range_mulVecLin, rank_transpose] using hlt
  obtain ⟨p, hpker, hp⟩ := (Submodule.ne_bot_iff _).mp (LinearMap.ker_ne_bot_of_finrank_lt (f := L) hfin)
  have hker : ∀ k, (E k)ᵀ *ᵥ (fun x => (avec τ m p x).eval (τ k)) = 0 := by
    intro k
    have := congrArg Subtype.val (congrFun (LinearMap.mem_ker.mp hpker) k)
    simpa [L, evA, Matrix.mulVec_transpose] using this
  have hzero := core ℓ d m J hd h0 τ hτ K E hdec hnice hm (avec τ m p)
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

/-! ### Optional: closed form of `N ℓ m` -/

theorem card_weight_eq (ℓ i : ℕ) :
    (Finset.univ.filter fun x : Fin ℓ → Bool => weight x = i).card = ℓ.choose i := by
  classical
  have := Finset.card_powersetCard i (Finset.univ : Finset (Fin ℓ))
  rw [Finset.card_univ, Fintype.card_fin] at this
  rw [← this]
  refine Finset.card_bij' (fun x _ => Finset.univ.filter fun k => x k = true)
    (fun s _ k => decide (k ∈ s)) ?_ ?_ ?_ ?_
  · intro x hx
    simp only [Finset.mem_filter, Finset.mem_univ, true_and] at hx
    simpa [Finset.mem_powersetCard, weight] using hx
  · intro s hs
    simp only [Finset.mem_powersetCard] at hs
    rw [Finset.mem_filter]
    refine ⟨Finset.mem_univ _, ?_⟩
    rw [← hs.2]
    unfold weight
    congr 1
    ext k
    simp
  · intro x _
    funext k
    simp
  · intro s _
    ext k
    simp

/-- `N ℓ m = ∑_{j ≤ m} ∑_{i ≤ j} C(ℓ, i)`. -/
theorem N_eq (ℓ m : ℕ) :
    N ℓ m = ∑ j ∈ Finset.range (m + 1), ∑ i ∈ Finset.range (j + 1), ℓ.choose i := by
  classical
  have h1 : ∀ w : ℕ, m + 1 - w = ∑ j ∈ Finset.range (m + 1), if w ≤ j then 1 else 0 := by
    intro w
    have : (Finset.range (m + 1)).filter (fun j => w ≤ j) = Finset.Ico w (m + 1) := by
      ext j
      simp [and_comm]
    rw [Finset.sum_boole, this]
    simp
  have h2 : ∀ j : ℕ, (Finset.univ.filter fun x : Fin ℓ → Bool => weight x ≤ j).card
      = ∑ i ∈ Finset.range (j + 1), ℓ.choose i := by
    intro j
    rw [Finset.card_eq_sum_card_fiberwise (f := weight) (t := Finset.range (j + 1))
      (fun x hx => by
        have hx' : weight x ≤ j := by simpa using hx
        simpa using Nat.lt_succ_of_le hx')]
    refine Finset.sum_congr rfl fun i hi => ?_
    rw [← card_weight_eq, Finset.filter_filter]
    congr 1
    ext x
    have := Finset.mem_range.mp hi
    simp only [Finset.mem_filter, Finset.mem_univ, true_and]
    omega
  simp only [N, h1]
  rw [Finset.sum_comm]
  refine Finset.sum_congr rfl fun j _ => ?_
  rw [Finset.sum_boole, h2]
  simp

end SqrtLog

#print axioms SqrtLog.core
#print axioms SqrtLog.sum_rank_ge_N
#print axioms SqrtLog.N_eq
