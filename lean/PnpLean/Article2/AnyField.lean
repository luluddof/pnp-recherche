import PnpLean.Article2.SqrtLog

/-!
# The dual polynomial bound over an arbitrary field, two-parameter form

Generalization of `SqrtLog.sum_rank_ge_N` (file `PnpLean/Article2/SqrtLog.lean`):

* the base field `ℚ` is replaced by an arbitrary field `F` (any characteristic);
* the single parameter `m` is replaced by two parameters `mc ≤ ma` with `ma + mc + d + 1 ≤ ℓ`.

`J(t)` is a `2^ℓ × 2^ℓ` polynomial matrix over `F` with entries of degree `≤ d`,
`J(0) = X^{⊗ℓ}`; `τ : Fin ℓ → F` are `ℓ` distinct nonzero points; `J(τ k) = K k + E k` with
`K k` weakly `k`-nice (`K k x y = 0` whenever `x k ≠ y k`). Then

  `∑ k, rank (E k) ≥ N2 ℓ ma mc := ∑_{x : |x| ≤ mc} (ma + 1 - |x|)`.

For `ma = mc = m` this is `SqrtLog.N ℓ m`, so the theorem over `ℚ` is a special case.

Proof. The dual space consists of the row vectors `a x = pref x * r x`,
`pref x = ∏_{k : x k = 1} (1 - X / τ k)`, with `r x = 0` unless `|x| ≤ mc` and `deg (a x) ≤ ma`;
it has dimension `N2 ℓ ma mc`. If `∑ rank (E k) < N2 ℓ ma mc` there is `a ≠ 0` in it with
`a(τ k) ᵥ* E k = 0` for all `k`. Write `a = X^v * b` with `b x₁ (0) ≠ 0` for some `x₁`
(necessarily `|x₁| ≤ mc`). The positivity argument of the rational proof is replaced by the
*test column* `c = pref x₁ • e_{¬x₁}`:

  `Ψ := (∑ x, b x * J x (¬x₁)) * pref x₁`

has degree `≤ ma + d + mc < ℓ`, vanishes at the `ℓ` points `τ k`, and `Ψ(0) = b x₁ (0) ≠ 0`:
contradiction.
-/

open Polynomial Matrix Module

namespace AnyField

open SqrtLog (weight N)

variable {ℓ : ℕ} {F : Type*} [Field F]

/-- `X^{⊗ℓ}` over an arbitrary field: the permutation matrix of the complement `x ↦ ¬x`. -/
def XtensorF (F : Type*) [Field F] (ℓ : ℕ) : Matrix (Fin ℓ → Bool) (Fin ℓ → Bool) F :=
  fun x y => if y = (fun i => !x i) then 1 else 0

theorem XtensorF_apply (ℓ : ℕ) (x y : Fin ℓ → Bool) :
    XtensorF F ℓ x y = if y = (fun i => !x i) then 1 else 0 := rfl

/-- Over `ℚ`, `XtensorF` is the matrix `PerPoint.Xtensor` of the existing files. -/
theorem XtensorF_rat (ℓ : ℕ) : XtensorF ℚ ℓ = PerPoint.Xtensor ℓ := by
  ext x y
  rw [PerPoint.Xtensor_apply, XtensorF_apply]

/-- Number of free coefficients at the point `x`: `ma + 1 - |x|` if `|x| ≤ mc`, else `0`. -/
def dimAt (ma mc : ℕ) (x : Fin ℓ → Bool) : ℕ :=
  if weight x ≤ mc then ma + 1 - weight x else 0

/-- `N2 ℓ ma mc = ∑_{x : |x| ≤ mc} (ma + 1 - |x|)`: the dimension of the dual space. -/
def N2 (ℓ ma mc : ℕ) : ℕ :=
  ∑ x : Fin ℓ → Bool, (if weight x ≤ mc then ma + 1 - weight x else 0)

/-- For `ma = mc = m`, `N2` is the `N` of `SqrtLog.lean`. -/
theorem N2_self (ℓ m : ℕ) : N2 ℓ m m = N ℓ m := by
  unfold N2 N
  refine Finset.sum_congr rfl fun x _ => ?_
  split_ifs <;> omega

/-! ### The prefactor -/

/-- The prefactor `∏_{k : x k = 1} (1 - X / τ k)`. -/
noncomputable def pref (τ : Fin ℓ → F) (x : Fin ℓ → Bool) : F[X] :=
  ∏ k ∈ Finset.univ.filter (fun k => x k = true), (1 - C (τ k)⁻¹ * X)

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
      (1 - C (τ k)⁻¹ * X : F[X]).natDegree ≤ 1 := by
    intro k _
    compute_degree
  simpa [weight] using Finset.sum_le_card_nsmul _ _ 1 this

/-! ### The core (polynomial) argument -/

/-- **Core lemma** (any field, two parameters). A row vector `a` of polynomials of degree `≤ ma`,
supported on the points of weight `≤ mc`, with `a x (τ k) = 0` whenever `x k = true` and
`a(τ k) ᵥ* E k = 0` (written `(E k)ᵀ *ᵥ a(τ k) = 0`) for all `k`, is zero. -/
theorem core (ℓ d ma mc : ℕ) (J : Matrix (Fin ℓ → Bool) (Fin ℓ → Bool) F[X])
    (hd : ∀ x y, (J x y).natDegree ≤ d) (h0 : J.map (eval 0) = XtensorF F ℓ)
    (τ : Fin ℓ → F) (hτ : Function.Injective τ) (hτ0 : ∀ k, τ k ≠ 0)
    (K E : Fin ℓ → Matrix (Fin ℓ → Bool) (Fin ℓ → Bool) F)
    (hdec : ∀ k, J.map (eval (τ k)) = K k + E k)
    (hnice : ∀ k x y, x k ≠ y k → K k x y = 0)
    (hm : ma + mc + d + 1 ≤ ℓ)
    (a : (Fin ℓ → Bool) → F[X]) (hdeg : ∀ x, (a x).natDegree ≤ ma)
    (hsupp : ∀ x, mc < weight x → a x = 0)
    (hvan : ∀ k x, x k = true → (a x).eval (τ k) = 0)
    (hker : ∀ k, (E k)ᵀ *ᵥ (fun x => (a x).eval (τ k)) = 0) : ∀ x, a x = 0 := by
  classical
  by_contra hne
  push Not at hne
  obtain ⟨x₀, hx₀⟩ := hne
  -- Step 1: factor out the largest common power of `X`.
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
  -- `x₁` has weight `≤ mc`, and the `b x` still have degree `≤ ma`.
  have hw : weight x₁ ≤ mc := by
    by_contra h
    push Not at h
    apply hx₁
    rw [hsupp x₁ h]
    exact dvd_zero _
  have hbdeg : ∀ x, (b x).natDegree ≤ ma := by
    intro x
    by_cases hbx : b x = 0
    · simp [hbx]
    · have := hdeg x
      rw [hb x, natDegree_mul (pow_ne_zero _ X_ne_zero) hbx] at this
      omega
  -- Step 2: the test polynomial `Ψ = (∑ x, b x * J x (¬x₁)) * pref x₁`.
  obtain ⟨S, hS⟩ : ∃ S : F[X], S = ∑ x, b x * J x (SqrtLog.compl x₁) := ⟨_, rfl⟩
  obtain ⟨Ψ, hΨ⟩ : ∃ Ψ : F[X], Ψ = S * pref τ x₁ := ⟨_, rfl⟩
  have hΨdeg : Ψ.natDegree ≤ ma + d + mc := by
    rw [hΨ]
    refine natDegree_mul_le_of_le ?_ ((natDegree_pref_le τ x₁).trans hw)
    rw [hS]
    refine natDegree_sum_le_of_forall_le _ _ fun x _ => ?_
    exact natDegree_mul_le_of_le (hbdeg x) (hd x (SqrtLog.compl x₁))
  -- Step 3: `Ψ` vanishes at every `τ k`.
  have hΨroot : ∀ k, Ψ.eval (τ k) = 0 := by
    intro k
    rw [hΨ, eval_mul]
    by_cases hk : x₁ k = true
    · rw [pref_eval τ x₁ k (hτ0 k) hk, mul_zero]
    · have hy : SqrtLog.compl x₁ k = true := by simp [SqrtLog.compl, hk]
      have hJ : ∀ x y, (J x y).eval (τ k) = K k x y + E k x y := by
        intro x y
        have := congrFun (congrFun (hdec k) x) y
        simpa using this
      have h1 : τ k ^ v * S.eval (τ k)
          = ∑ x, (a x).eval (τ k) * K k x (SqrtLog.compl x₁)
            + ∑ x, E k x (SqrtLog.compl x₁) * (a x).eval (τ k) := by
        rw [hS]
        simp only [eval_finsetSum, eval_mul, hJ, Finset.mul_sum]
        rw [← Finset.sum_add_distrib]
        refine Finset.sum_congr rfl fun x _ => ?_
        rw [hb x]
        simp only [eval_mul, eval_pow, eval_X]
        ring
      have h2 : ∀ x, (a x).eval (τ k) * K k x (SqrtLog.compl x₁) = 0 := by
        intro x
        by_cases hx : x k = true
        · rw [hvan k x hx, zero_mul]
        · rw [hnice k x (SqrtLog.compl x₁) (by rw [hy]; exact hx), mul_zero]
      have h3 : ∑ x, E k x (SqrtLog.compl x₁) * (a x).eval (τ k) = 0 := by
        have := congrFun (hker k) (SqrtLog.compl x₁)
        simpa [Matrix.mulVec, dotProduct] using this
      have h4 : τ k ^ v * S.eval (τ k) = 0 := by
        rw [h1, h3]
        simp [h2]
      rcases mul_eq_zero.mp h4 with h | h
      · exact absurd h (pow_ne_zero _ (hτ0 k))
      · rw [h, zero_mul]
  -- Step 4: `Ψ(0) = b x₁ (0) ≠ 0`.
  have hJ0 : ∀ x y, (J x y).eval 0 = if y = SqrtLog.compl x then 1 else 0 := by
    intro x y
    have := congrFun (congrFun h0 x) y
    exact this
  have hΨ0 : Ψ.eval 0 = (b x₁).eval 0 := by
    rw [hΨ, hS]
    simp only [eval_mul, eval_finsetSum, hJ0, pref_eval_zero, mul_one]
    rw [Finset.sum_eq_single x₁]
    · simp
    · intro x _ hx
      have hne : SqrtLog.compl x₁ ≠ SqrtLog.compl x := fun h => hx (by
        have := congrArg SqrtLog.compl h
        simpa using this.symm)
      simp [hne]
    · simp
  -- Step 5: a nonzero polynomial of degree `< ℓ` cannot have `ℓ` distinct roots.
  have hzero : Ψ = 0 :=
    eq_zero_of_natDegree_lt_card_of_eval_eq_zero Ψ hτ hΨroot (by simp; omega)
  rw [hzero, eval_zero] at hΨ0
  exact hb1 hΨ0.symm

/-! ### The dual space and its parametrization -/

/-- The free factor `r x`: of degree `≤ ma - |x|` if `|x| ≤ mc`, zero otherwise. -/
noncomputable def rpoly (ma mc : ℕ) (p : ∀ x : Fin ℓ → Bool, Fin (dimAt ma mc x) → F)
    (x : Fin ℓ → Bool) : F[X] :=
  ∑ i : Fin (dimAt ma mc x), C (p x i) * X ^ (i : ℕ)

/-- The element of the dual space with parameters `p`. -/
noncomputable def avec (τ : Fin ℓ → F) (ma mc : ℕ)
    (p : ∀ x : Fin ℓ → Bool, Fin (dimAt ma mc x) → F) (x : Fin ℓ → Bool) : F[X] :=
  pref τ x * rpoly ma mc p x

theorem dimAt_le (ma mc : ℕ) (x : Fin ℓ → Bool) : dimAt ma mc x ≤ ma + 1 - weight x := by
  unfold dimAt
  split_ifs <;> omega

theorem rpoly_eq_zero_of_dimAt (ma mc : ℕ) (p : ∀ x : Fin ℓ → Bool, Fin (dimAt ma mc x) → F)
    (x : Fin ℓ → Bool) (h : dimAt ma mc x = 0) : rpoly ma mc p x = 0 :=
  Finset.sum_eq_zero fun i _ => absurd i.2 (by omega)

theorem natDegree_rpoly_le (ma mc : ℕ) (p : ∀ x : Fin ℓ → Bool, Fin (dimAt ma mc x) → F)
    (x : Fin ℓ → Bool) : (rpoly ma mc p x).natDegree ≤ ma - weight x := by
  refine natDegree_sum_le_of_forall_le _ _ fun i _ => ?_
  have := natDegree_C_mul_X_pow_le (p x i) (i : ℕ)
  have := i.2
  have := dimAt_le ma mc x
  omega

theorem natDegree_avec_le (τ : Fin ℓ → F) (ma mc : ℕ)
    (p : ∀ x : Fin ℓ → Bool, Fin (dimAt ma mc x) → F) (x : Fin ℓ → Bool) :
    (avec τ ma mc p x).natDegree ≤ ma := by
  by_cases h : dimAt ma mc x = 0
  · simp [avec, rpoly_eq_zero_of_dimAt ma mc p x h]
  · have := natDegree_mul_le_of_le (natDegree_pref_le τ x) (natDegree_rpoly_le ma mc p x)
    have := dimAt_le ma mc x
    rw [avec]
    omega

theorem avec_eq_zero_of_weight (τ : Fin ℓ → F) (ma mc : ℕ)
    (p : ∀ x : Fin ℓ → Bool, Fin (dimAt ma mc x) → F) (x : Fin ℓ → Bool)
    (h : mc < weight x) : avec τ ma mc p x = 0 := by
  have hlt : ¬ weight x ≤ mc := by omega
  have : dimAt ma mc x = 0 := by simp [dimAt, hlt]
  simp [avec, rpoly_eq_zero_of_dimAt ma mc p x this]

theorem coeff_rpoly (ma mc : ℕ) (p : ∀ x : Fin ℓ → Bool, Fin (dimAt ma mc x) → F)
    (x : Fin ℓ → Bool) (i : Fin (dimAt ma mc x)) : (rpoly ma mc p x).coeff i = p x i := by
  simp only [rpoly, finsetSum_coeff, coeff_C_mul_X_pow]
  rw [Finset.sum_eq_single i]
  · simp
  · intro j _ hj
    have : (i : ℕ) ≠ (j : ℕ) := fun h => hj (Fin.ext h.symm)
    simp [this]
  · simp

/-- Evaluation of `avec` at `τ k`, as a linear map in the parameters. -/
noncomputable def evA (τ : Fin ℓ → F) (ma mc : ℕ) (k : Fin ℓ) :
    (∀ x : Fin ℓ → Bool, Fin (dimAt ma mc x) → F) →ₗ[F] ((Fin ℓ → Bool) → F) where
  toFun p x := (avec τ ma mc p x).eval (τ k)
  map_add' p q := by
    funext x
    simp [avec, rpoly, eval_finsetSum, add_mul, Finset.sum_add_distrib, mul_add]
  map_smul' c p := by
    funext x
    simp [avec, rpoly, eval_finsetSum, Finset.mul_sum]
    refine Finset.sum_congr rfl fun i _ => ?_
    ring

theorem finrank_params (ℓ ma mc : ℕ) :
    finrank F (∀ x : Fin ℓ → Bool, Fin (dimAt ma mc x) → F) = N2 ℓ ma mc := by
  rw [Module.finrank_pi_fintype]
  simp only [Module.finrank_fin_fun]
  rfl

/-! ### Main theorems -/

/-- **Dual polynomial bound over any field, two-parameter form.**
`J(t)` has entries of degree `≤ d`, `J(0) = X^{⊗ℓ}`; the `τ k` are distinct and nonzero;
`J(τ k) = K k + E k` with `K k` weakly `k`-nice. If `mc ≤ ma` and `ma + mc + d + 1 ≤ ℓ` then
`∑ k, rank (E k) ≥ N2 ℓ ma mc = ∑_{x : |x| ≤ mc} (ma + 1 - |x|)`. -/
theorem sum_rank_ge_N2 (ℓ d : ℕ) (J : Matrix (Fin ℓ → Bool) (Fin ℓ → Bool) F[X])
    (hd : ∀ x y, (J x y).natDegree ≤ d) (h0 : J.map (eval 0) = XtensorF F ℓ)
    (τ : Fin ℓ → F) (hτ : Function.Injective τ) (hτ0 : ∀ k, τ k ≠ 0)
    (K E : Fin ℓ → Matrix (Fin ℓ → Bool) (Fin ℓ → Bool) F)
    (hdec : ∀ k, J.map (eval (τ k)) = K k + E k)
    (hnice : ∀ k x y, x k ≠ y k → K k x y = 0)
    (ma mc : ℕ) (_hmc : mc ≤ ma) (hm : ma + mc + d + 1 ≤ ℓ) :
    N2 ℓ ma mc ≤ ∑ k, (E k).rank := by
  classical
  by_contra hlt
  push Not at hlt
  let L : (∀ x : Fin ℓ → Bool, Fin (dimAt ma mc x) → F) →ₗ[F]
      (∀ k : Fin ℓ, LinearMap.range ((E k)ᵀ.mulVecLin)) :=
    LinearMap.pi fun k => ((E k)ᵀ.mulVecLin.rangeRestrict).comp (evA τ ma mc k)
  have hfin : finrank F (∀ k : Fin ℓ, LinearMap.range ((E k)ᵀ.mulVecLin))
      < finrank F (∀ x : Fin ℓ → Bool, Fin (dimAt ma mc x) → F) := by
    rw [finrank_params, Module.finrank_pi_fintype]
    simpa [finrank_range_mulVecLin, rank_transpose] using hlt
  obtain ⟨p, hpker, hp⟩ :=
    (Submodule.ne_bot_iff _).mp (LinearMap.ker_ne_bot_of_finrank_lt (f := L) hfin)
  have hker : ∀ k, (E k)ᵀ *ᵥ (fun x => (avec τ ma mc p x).eval (τ k)) = 0 := by
    intro k
    have := congrArg Subtype.val (congrFun (LinearMap.mem_ker.mp hpker) k)
    simpa [L, evA, Matrix.mulVec_transpose] using this
  have hzero := core ℓ d ma mc J hd h0 τ hτ hτ0 K E hdec hnice hm (avec τ ma mc p)
    (natDegree_avec_le τ ma mc p) (avec_eq_zero_of_weight τ ma mc p)
    (fun k x hx => by rw [avec, eval_mul, pref_eval τ x k (hτ0 k) hx, zero_mul]) hker
  apply hp
  funext x i
  have hr : rpoly ma mc p x = 0 := by
    rcases mul_eq_zero.mp (hzero x) with h | h
    · exact absurd h (pref_ne_zero τ x)
    · exact h
  have := coeff_rpoly ma mc p x i
  rw [hr, coeff_zero] at this
  simpa using this.symm

/-- **Corollary: the one-parameter bound over any field.** If `2 m + d + 1 ≤ ℓ` then
`∑ k, rank (E k) ≥ N ℓ m = ∑ x, (m + 1 - |x|)` (the `N` of `SqrtLog.lean`). -/
theorem sum_rank_ge_N (ℓ d : ℕ) (J : Matrix (Fin ℓ → Bool) (Fin ℓ → Bool) F[X])
    (hd : ∀ x y, (J x y).natDegree ≤ d) (h0 : J.map (eval 0) = XtensorF F ℓ)
    (τ : Fin ℓ → F) (hτ : Function.Injective τ) (hτ0 : ∀ k, τ k ≠ 0)
    (K E : Fin ℓ → Matrix (Fin ℓ → Bool) (Fin ℓ → Bool) F)
    (hdec : ∀ k, J.map (eval (τ k)) = K k + E k)
    (hnice : ∀ k x y, x k ≠ y k → K k x y = 0)
    (m : ℕ) (hm : 2 * m + d + 1 ≤ ℓ) :
    N ℓ m ≤ ∑ k, (E k).rank := by
  rw [← N2_self]
  exact sum_rank_ge_N2 ℓ d J hd h0 τ hτ hτ0 K E hdec hnice m m le_rfl (by omega)

/-- Sanity check: the theorem `SqrtLog.sum_rank_ge_N` over `ℚ` is the case `F = ℚ`. -/
theorem sum_rank_ge_N_rat (ℓ d : ℕ) (J : Matrix (Fin ℓ → Bool) (Fin ℓ → Bool) ℚ[X])
    (hd : ∀ x y, (J x y).natDegree ≤ d) (h0 : J.map (eval 0) = PerPoint.Xtensor ℓ)
    (τ : Fin ℓ → ℚ) (hτ : Function.Injective τ) (hτ0 : ∀ k, τ k ≠ 0)
    (K E : Fin ℓ → Matrix (Fin ℓ → Bool) (Fin ℓ → Bool) ℚ)
    (hdec : ∀ k, J.map (eval (τ k)) = K k + E k)
    (hnice : ∀ k x y, x k ≠ y k → K k x y = 0)
    (m : ℕ) (hm : 2 * m + d + 1 ≤ ℓ) :
    N ℓ m ≤ ∑ k, (E k).rank :=
  sum_rank_ge_N ℓ d J hd (h0.trans (XtensorF_rat ℓ).symm) τ hτ hτ0 K E hdec hnice m hm

end AnyField

#print axioms AnyField.core
#print axioms AnyField.sum_rank_ge_N2
#print axioms AnyField.sum_rank_ge_N
#print axioms AnyField.sum_rank_ge_N_rat
