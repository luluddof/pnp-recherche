import PnpLean.Article2.AnyField

/-!
# The multi-point dual bound over an arbitrary field (factor-2 form)

Setting. `J(t)` is a `2^ℓ × 2^ℓ` polynomial matrix over a field `F` with entries of degree `≤ d`
and `J(0) = X^{⊗ℓ}`. The points are indexed by a finite type `ι` (`N = Fintype.card ι`):
`τ : ι → F` injective and nonzero, and each point `k` has a direction `c k : Fin ℓ`
(several points, or none, per direction). `J(τ k) = K k + E k` with `K k` weakly `c k`-nice.
For `x : Fin ℓ → Bool`, `μ x = #{k : x (c k) = 1}`. Then

  `∑ k, rank (E k) ≥ ∑ x, max(0, N - d - 2 μ x)`     (`MultiPoint.sum_rank_ge`).

Only "half" of the niceness is used: `K k x y = 0` whenever `x (c k) = 0` and `y (c k) = 1`
(`MultiPoint.sum_rank_ge_half`); equivalently the bound holds for the `01` off-nice parts of the
`J(τ k)` themselves (`MultiPoint.sum_rank_block01_ge`).

Proof (single test vector, as in `AnyField.core`, but with a different dual space).
Let `D = N - d - 1` and `pref x = ∏_{k : x (c k) = 1} (1 - X / τ k)` (degree `≤ μ x`).
The dual space `W` consists of the row vectors

  `a x = pref x * X ^ (μ x) * r x`,  `deg (r x) ≤ D - 2 μ x`  (`r x = 0` if `2 μ x > D`),

so `deg (a x) ≤ D` and `dim W = ∑ x, max(0, D + 1 - 2 μ x)`. If `∑ rank (E k) < dim W` there is
`a ≠ 0` in `W` with `a(τ k) ᵥ* E k = 0` for all `k`. Write `a = X^v * b` with `b x₁ (0) ≠ 0`.
Since `X ^ (μ x₁) ∣ a x₁` and `X ^ (v+1) ∤ a x₁`, we get `μ x₁ ≤ v`. The test polynomial

  `Ψ := (∑ x, b x * J x (¬x₁)) * pref x₁`

has degree `≤ (D - v) + d + μ x₁ ≤ D + d < N`, vanishes at the `N` points `τ k`, and
`Ψ(0) = b x₁ (0) ≠ 0`: contradiction.

(This is the filtration argument of the notes `code/gen14/g14d/notes.md`, Theorem 1, in closed
form: `W` is a complement of the space `V_D` of annihilators, because the lowest-order
coefficient of an element of `V_D` of `X`-adic valuation `v` is supported on `v < μ x ≤ D - v`.)

The weaker "one-parameter" form `∑ x, (m + 1 - μ x) ≤ ∑ k, rank (E k)` for `2 m + d + 1 ≤ N`
is `MultiPoint.sum_rank_ge_simple`; one point per direction is `MultiPoint.sum_rank_ge_one`.
-/

open Polynomial Matrix Module

namespace MultiPoint

open AnyField (XtensorF)

variable {ℓ : ℕ} {F : Type*} [Field F] {ι : Type*} [Fintype ι]

/-- `μ x`: the number of points `k` whose direction `c k` satisfies `x (c k) = 1`. -/
def mu (c : ι → Fin ℓ) (x : Fin ℓ → Bool) : ℕ :=
  (Finset.univ.filter fun k => x (c k) = true).card

/-- One point per direction: `μ` is the Hamming weight. -/
theorem mu_id (x : Fin ℓ → Bool) : mu (fun k : Fin ℓ => k) x = SqrtLog.weight x := rfl

/-! ### The prefactor -/

/-- The prefactor `∏_{k : x (c k) = 1} (1 - X / τ k)`. -/
noncomputable def pref (τ : ι → F) (c : ι → Fin ℓ) (x : Fin ℓ → Bool) : F[X] :=
  ∏ k ∈ Finset.univ.filter (fun k => x (c k) = true), (1 - C (τ k)⁻¹ * X)

theorem pref_eval_zero (τ : ι → F) (c : ι → Fin ℓ) (x : Fin ℓ → Bool) :
    (pref τ c x).eval 0 = 1 := by
  simp [pref, eval_prod]

theorem pref_ne_zero (τ : ι → F) (c : ι → Fin ℓ) (x : Fin ℓ → Bool) : pref τ c x ≠ 0 := by
  intro h
  have := pref_eval_zero τ c x
  rw [h, eval_zero] at this
  exact zero_ne_one this

theorem pref_eval (τ : ι → F) (c : ι → Fin ℓ) (x : Fin ℓ → Bool) (k : ι) (hk : τ k ≠ 0)
    (hx : x (c k) = true) : (pref τ c x).eval (τ k) = 0 := by
  rw [pref, eval_prod]
  refine Finset.prod_eq_zero (i := k) (by simp [hx]) ?_
  simp [hk]

theorem natDegree_pref_le (τ : ι → F) (c : ι → Fin ℓ) (x : Fin ℓ → Bool) :
    (pref τ c x).natDegree ≤ mu c x := by
  rw [pref]
  refine (natDegree_prod_le _ _).trans ?_
  have : ∀ k ∈ Finset.univ.filter (fun k => x (c k) = true),
      (1 - C (τ k)⁻¹ * X : F[X]).natDegree ≤ 1 := by
    intro k _
    compute_degree
  simpa [mu] using Finset.sum_le_card_nsmul _ _ 1 this

/-! ### The core (polynomial) argument -/

/-- **Core lemma** (any field, any placement of the points). Let `D + d < N`. A row vector `a` of
polynomials of degree `≤ D` with `X ^ (μ x) ∣ a x`, `a x (τ k) = 0` whenever `x (c k) = true`,
and `a(τ k) ᵥ* E k = 0` (written `(E k)ᵀ *ᵥ a(τ k) = 0`) for all `k`, is zero.
Only the `01` half of the niceness of `K k` is used. -/
theorem core (ℓ d D : ℕ) (J : Matrix (Fin ℓ → Bool) (Fin ℓ → Bool) F[X])
    (hd : ∀ x y, (J x y).natDegree ≤ d) (h0 : J.map (eval 0) = XtensorF F ℓ)
    (τ : ι → F) (hτ : Function.Injective τ) (hτ0 : ∀ k, τ k ≠ 0) (c : ι → Fin ℓ)
    (K E : ι → Matrix (Fin ℓ → Bool) (Fin ℓ → Bool) F)
    (hdec : ∀ k, J.map (eval (τ k)) = K k + E k)
    (hnice : ∀ k x y, x (c k) = false → y (c k) = true → K k x y = 0)
    (hD : D + d < Fintype.card ι)
    (a : (Fin ℓ → Bool) → F[X]) (hdeg : ∀ x, (a x).natDegree ≤ D)
    (hdiv : ∀ x, (X : F[X]) ^ mu c x ∣ a x)
    (hvan : ∀ k x, x (c k) = true → (a x).eval (τ k) = 0)
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
  -- The key new point: `μ x₁ ≤ v`, because `X ^ (μ x₁) ∣ a x₁` but `X ^ (v + 1) ∤ a x₁`.
  have hμ : mu c x₁ ≤ v := by
    by_contra h
    push Not at h
    exact hx₁ ((pow_dvd_pow (X : F[X]) (by omega : v + 1 ≤ mu c x₁)).trans (hdiv x₁))
  -- `v ≤ D` and the `b x` have degree `≤ D - v`.
  have hbne : b x₁ ≠ 0 := by
    intro h
    rw [h, eval_zero] at hb1
    exact hb1 rfl
  have hvD : v ≤ D := by
    have := hdeg x₁
    rw [hb x₁, natDegree_mul (pow_ne_zero _ X_ne_zero) hbne, natDegree_X_pow] at this
    omega
  have hbdeg : ∀ x, (b x).natDegree ≤ D - v := by
    intro x
    by_cases hbx : b x = 0
    · simp [hbx]
    · have := hdeg x
      rw [hb x, natDegree_mul (pow_ne_zero _ X_ne_zero) hbx, natDegree_X_pow] at this
      omega
  -- Step 2: the test polynomial `Ψ = (∑ x, b x * J x (¬x₁)) * pref x₁`.
  obtain ⟨S, hS⟩ : ∃ S : F[X], S = ∑ x, b x * J x (SqrtLog.compl x₁) := ⟨_, rfl⟩
  obtain ⟨Ψ, hΨ⟩ : ∃ Ψ : F[X], Ψ = S * pref τ c x₁ := ⟨_, rfl⟩
  have hΨdeg : Ψ.natDegree ≤ (D - v) + d + mu c x₁ := by
    rw [hΨ]
    refine natDegree_mul_le_of_le ?_ (natDegree_pref_le τ c x₁)
    rw [hS]
    refine natDegree_sum_le_of_forall_le _ _ fun x _ => ?_
    exact natDegree_mul_le_of_le (hbdeg x) (hd x (SqrtLog.compl x₁))
  -- Step 3: `Ψ` vanishes at every `τ k`.
  have hΨroot : ∀ k, Ψ.eval (τ k) = 0 := by
    intro k
    rw [hΨ, eval_mul]
    by_cases hk : x₁ (c k) = true
    · rw [pref_eval τ c x₁ k (hτ0 k) hk, mul_zero]
    · have hy : SqrtLog.compl x₁ (c k) = true := by simp [SqrtLog.compl, hk]
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
        by_cases hx : x (c k) = true
        · rw [hvan k x hx, zero_mul]
        · rw [hnice k x (SqrtLog.compl x₁) (by simpa using hx) hy, mul_zero]
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
  -- Step 5: a nonzero polynomial of degree `< N` cannot have `N` distinct roots.
  have hzero : Ψ = 0 :=
    eq_zero_of_natDegree_lt_card_of_eval_eq_zero Ψ hτ hΨroot (by omega)
  rw [hzero, eval_zero] at hΨ0
  exact hb1 hΨ0.symm

/-! ### The dual space and its parametrization -/

/-- Number of free coefficients at the point `x`: `max(0, D + 1 - 2 μ x)`. -/
def dimAt (D : ℕ) (c : ι → Fin ℓ) (x : Fin ℓ → Bool) : ℕ := D + 1 - 2 * mu c x

/-- The free factor `r x`, of degree `< dimAt D c x`. -/
noncomputable def rpoly (D : ℕ) (c : ι → Fin ℓ)
    (p : ∀ x : Fin ℓ → Bool, Fin (dimAt D c x) → F) (x : Fin ℓ → Bool) : F[X] :=
  ∑ i : Fin (dimAt D c x), C (p x i) * X ^ (i : ℕ)

/-- The element `a x = pref x * (X ^ (μ x) * r x)` of the dual space with parameters `p`. -/
noncomputable def avec (τ : ι → F) (D : ℕ) (c : ι → Fin ℓ)
    (p : ∀ x : Fin ℓ → Bool, Fin (dimAt D c x) → F) (x : Fin ℓ → Bool) : F[X] :=
  pref τ c x * (X ^ mu c x * rpoly D c p x)

theorem rpoly_eq_zero_of_dimAt (D : ℕ) (c : ι → Fin ℓ)
    (p : ∀ x : Fin ℓ → Bool, Fin (dimAt D c x) → F)
    (x : Fin ℓ → Bool) (h : dimAt D c x = 0) : rpoly D c p x = 0 :=
  Finset.sum_eq_zero fun i _ => absurd i.2 (by omega)

theorem natDegree_rpoly_le (D : ℕ) (c : ι → Fin ℓ)
    (p : ∀ x : Fin ℓ → Bool, Fin (dimAt D c x) → F) (x : Fin ℓ → Bool) :
    (rpoly D c p x).natDegree ≤ D - 2 * mu c x := by
  refine natDegree_sum_le_of_forall_le _ _ fun i _ => ?_
  have := natDegree_C_mul_X_pow_le (p x i) (i : ℕ)
  have := i.2
  unfold dimAt at this
  omega

theorem natDegree_avec_le (τ : ι → F) (D : ℕ) (c : ι → Fin ℓ)
    (p : ∀ x : Fin ℓ → Bool, Fin (dimAt D c x) → F) (x : Fin ℓ → Bool) :
    (avec τ D c p x).natDegree ≤ D := by
  by_cases h : dimAt D c x = 0
  · simp [avec, rpoly_eq_zero_of_dimAt D c p x h]
  · have h1 : (X ^ mu c x * rpoly D c p x : F[X]).natDegree ≤ mu c x + (D - 2 * mu c x) :=
      natDegree_mul_le_of_le (natDegree_X_pow_le _) (natDegree_rpoly_le D c p x)
    have h2 := natDegree_mul_le_of_le (natDegree_pref_le τ c x) h1
    unfold dimAt at h
    rw [avec]
    omega

theorem coeff_rpoly (D : ℕ) (c : ι → Fin ℓ) (p : ∀ x : Fin ℓ → Bool, Fin (dimAt D c x) → F)
    (x : Fin ℓ → Bool) (i : Fin (dimAt D c x)) : (rpoly D c p x).coeff i = p x i := by
  simp only [rpoly, finsetSum_coeff, coeff_C_mul_X_pow]
  rw [Finset.sum_eq_single i]
  · simp
  · intro j _ hj
    have : (i : ℕ) ≠ (j : ℕ) := fun h => hj (Fin.ext h.symm)
    simp [this]
  · simp

/-- Evaluation of `avec` at `τ k`, as a linear map in the parameters. -/
noncomputable def evA (τ : ι → F) (D : ℕ) (c : ι → Fin ℓ) (k : ι) :
    (∀ x : Fin ℓ → Bool, Fin (dimAt D c x) → F) →ₗ[F] ((Fin ℓ → Bool) → F) where
  toFun p x := (avec τ D c p x).eval (τ k)
  map_add' p q := by
    funext x
    simp [avec, rpoly, eval_finsetSum, add_mul, Finset.sum_add_distrib, mul_add]
  map_smul' a p := by
    funext x
    simp [avec, rpoly, eval_finsetSum, Finset.mul_sum]
    refine Finset.sum_congr rfl fun i _ => ?_
    ring

theorem finrank_params (D : ℕ) (c : ι → Fin ℓ) :
    finrank F (∀ x : Fin ℓ → Bool, Fin (dimAt D c x) → F) = ∑ x, dimAt D c x := by
  rw [Module.finrank_pi_fintype]
  simp only [Module.finrank_fin_fun]

/-! ### Main theorems -/

/-- **Multi-point dual bound, factor-2 form, any field, any placement** (half-niceness version).
`J(t)` has entries of degree `≤ d`, `J(0) = X^{⊗ℓ}`; the points `τ k` (`k : ι`,
`N = Fintype.card ι`) are distinct and nonzero, the point `k` has direction `c k`;
`J(τ k) = K k + E k` with `K k x y = 0` whenever `x (c k) = 0` and `y (c k) = 1`. Then
`∑ k, rank (E k) ≥ ∑ x, max(0, N - d - 2 μ x)` (natural subtraction). -/
theorem sum_rank_ge_half (ℓ d : ℕ) (J : Matrix (Fin ℓ → Bool) (Fin ℓ → Bool) F[X])
    (hd : ∀ x y, (J x y).natDegree ≤ d) (h0 : J.map (eval 0) = XtensorF F ℓ)
    (τ : ι → F) (hτ : Function.Injective τ) (hτ0 : ∀ k, τ k ≠ 0) (c : ι → Fin ℓ)
    (K E : ι → Matrix (Fin ℓ → Bool) (Fin ℓ → Bool) F)
    (hdec : ∀ k, J.map (eval (τ k)) = K k + E k)
    (hnice : ∀ k x y, x (c k) = false → y (c k) = true → K k x y = 0) :
    ∑ x : Fin ℓ → Bool, (Fintype.card ι - d - 2 * mu c x) ≤ ∑ k, (E k).rank := by
  classical
  by_cases hNd : Fintype.card ι ≤ d
  · have : ∑ x : Fin ℓ → Bool, (Fintype.card ι - d - 2 * mu c x) = 0 :=
      Finset.sum_eq_zero fun x _ => by omega
    omega
  push Not at hNd
  obtain ⟨D, hD⟩ : ∃ D, D = Fintype.card ι - d - 1 := ⟨_, rfl⟩
  have hsum : ∑ x : Fin ℓ → Bool, (Fintype.card ι - d - 2 * mu c x) = ∑ x, dimAt D c x :=
    Finset.sum_congr rfl fun x _ => by unfold dimAt; omega
  rw [hsum]
  by_contra hlt
  push Not at hlt
  let L : (∀ x : Fin ℓ → Bool, Fin (dimAt D c x) → F) →ₗ[F]
      (∀ k : ι, LinearMap.range ((E k)ᵀ.mulVecLin)) :=
    LinearMap.pi fun k => ((E k)ᵀ.mulVecLin.rangeRestrict).comp (evA τ D c k)
  have hfin : finrank F (∀ k : ι, LinearMap.range ((E k)ᵀ.mulVecLin))
      < finrank F (∀ x : Fin ℓ → Bool, Fin (dimAt D c x) → F) := by
    rw [finrank_params, Module.finrank_pi_fintype]
    simpa [finrank_range_mulVecLin, rank_transpose] using hlt
  obtain ⟨p, hpker, hp⟩ :=
    (Submodule.ne_bot_iff _).mp (LinearMap.ker_ne_bot_of_finrank_lt (f := L) hfin)
  have hker : ∀ k, (E k)ᵀ *ᵥ (fun x => (avec τ D c p x).eval (τ k)) = 0 := by
    intro k
    have := congrArg Subtype.val (congrFun (LinearMap.mem_ker.mp hpker) k)
    simpa [L, evA, Matrix.mulVec_transpose] using this
  have hzero := core ℓ d D J hd h0 τ hτ hτ0 c K E hdec hnice (by omega) (avec τ D c p)
    (natDegree_avec_le τ D c p)
    (fun x => by
      rw [avec]
      exact Dvd.dvd.mul_left (Dvd.intro _ rfl) _)
    (fun k x hx => by rw [avec, eval_mul, pref_eval τ c x k (hτ0 k) hx, zero_mul]) hker
  apply hp
  funext x i
  have hr : rpoly D c p x = 0 := by
    rcases mul_eq_zero.mp (hzero x) with h | h
    · exact absurd h (pref_ne_zero τ c x)
    · rcases mul_eq_zero.mp h with h' | h'
      · exact absurd h' (pow_ne_zero _ X_ne_zero)
      · exact h'
  have := coeff_rpoly D c p x i
  rw [hr, coeff_zero] at this
  simpa using this.symm

/-- **Multi-point dual bound, factor-2 form** (Theorem 1 of `code/gen14/g14d/notes.md`).
Same as `sum_rank_ge_half` with the usual weak niceness `K k x y = 0` whenever
`x (c k) ≠ y (c k)`. -/
theorem sum_rank_ge (ℓ d : ℕ) (J : Matrix (Fin ℓ → Bool) (Fin ℓ → Bool) F[X])
    (hd : ∀ x y, (J x y).natDegree ≤ d) (h0 : J.map (eval 0) = XtensorF F ℓ)
    (τ : ι → F) (hτ : Function.Injective τ) (hτ0 : ∀ k, τ k ≠ 0) (c : ι → Fin ℓ)
    (K E : ι → Matrix (Fin ℓ → Bool) (Fin ℓ → Bool) F)
    (hdec : ∀ k, J.map (eval (τ k)) = K k + E k)
    (hnice : ∀ k x y, x (c k) ≠ y (c k) → K k x y = 0) :
    ∑ x : Fin ℓ → Bool, (Fintype.card ι - d - 2 * mu c x) ≤ ∑ k, (E k).rank :=
  sum_rank_ge_half ℓ d J hd h0 τ hτ hτ0 c K E hdec
    fun k x y hx hy => hnice k x y (by rw [hx, hy]; decide)

/-- The `01` off-nice part of a matrix in direction `i`: the entries `(x, y)` with `x i = 0`,
`y i = 1`, all other entries replaced by `0`. -/
def block01 (i : Fin ℓ) (M : Matrix (Fin ℓ → Bool) (Fin ℓ → Bool) F) :
    Matrix (Fin ℓ → Bool) (Fin ℓ → Bool) F :=
  fun x y => if x i = false ∧ y i = true then M x y else 0

/-- **Hypothesis-free form**: the bound holds for the `01` off-nice parts of the `J(τ k)`. -/
theorem sum_rank_block01_ge (ℓ d : ℕ) (J : Matrix (Fin ℓ → Bool) (Fin ℓ → Bool) F[X])
    (hd : ∀ x y, (J x y).natDegree ≤ d) (h0 : J.map (eval 0) = XtensorF F ℓ)
    (τ : ι → F) (hτ : Function.Injective τ) (hτ0 : ∀ k, τ k ≠ 0) (c : ι → Fin ℓ) :
    ∑ x : Fin ℓ → Bool, (Fintype.card ι - d - 2 * mu c x)
      ≤ ∑ k, (block01 (c k) (J.map (eval (τ k)))).rank :=
  sum_rank_ge_half ℓ d J hd h0 τ hτ hτ0 c
    (fun k => J.map (eval (τ k)) - block01 (c k) (J.map (eval (τ k))))
    (fun k => block01 (c k) (J.map (eval (τ k))))
    (fun k => by simp)
    (fun k x y hx hy => by simp [block01, hx, hy])

/-- **Simple (one-parameter) multi-point form**: if `2 m + d + 1 ≤ N` then
`∑ k, rank (E k) ≥ ∑ x, max(0, m + 1 - μ x)`. Termwise weaker than `sum_rank_ge`. -/
theorem sum_rank_ge_simple (ℓ d : ℕ) (J : Matrix (Fin ℓ → Bool) (Fin ℓ → Bool) F[X])
    (hd : ∀ x y, (J x y).natDegree ≤ d) (h0 : J.map (eval 0) = XtensorF F ℓ)
    (τ : ι → F) (hτ : Function.Injective τ) (hτ0 : ∀ k, τ k ≠ 0) (c : ι → Fin ℓ)
    (K E : ι → Matrix (Fin ℓ → Bool) (Fin ℓ → Bool) F)
    (hdec : ∀ k, J.map (eval (τ k)) = K k + E k)
    (hnice : ∀ k x y, x (c k) ≠ y (c k) → K k x y = 0)
    (m : ℕ) (hm : 2 * m + d + 1 ≤ Fintype.card ι) :
    ∑ x : Fin ℓ → Bool, (m + 1 - mu c x) ≤ ∑ k, (E k).rank :=
  (Finset.sum_le_sum fun x _ => by omega).trans
    (sum_rank_ge ℓ d J hd h0 τ hτ hτ0 c K E hdec hnice)

/-- **One point per direction** (`ι = Fin ℓ`, `c = id`):
`∑ k, rank (E k) ≥ ∑ x, max(0, ℓ - d - 2 |x|)`, twice the bound `AnyField.sum_rank_ge_N`
(up to parity). -/
theorem sum_rank_ge_one (ℓ d : ℕ) (J : Matrix (Fin ℓ → Bool) (Fin ℓ → Bool) F[X])
    (hd : ∀ x y, (J x y).natDegree ≤ d) (h0 : J.map (eval 0) = XtensorF F ℓ)
    (τ : Fin ℓ → F) (hτ : Function.Injective τ) (hτ0 : ∀ k, τ k ≠ 0)
    (K E : Fin ℓ → Matrix (Fin ℓ → Bool) (Fin ℓ → Bool) F)
    (hdec : ∀ k, J.map (eval (τ k)) = K k + E k)
    (hnice : ∀ k x y, x k ≠ y k → K k x y = 0) :
    ∑ x : Fin ℓ → Bool, (ℓ - d - 2 * SqrtLog.weight x) ≤ ∑ k, (E k).rank := by
  have := sum_rank_ge ℓ d J hd h0 τ hτ hτ0 (fun k => k) K E hdec hnice
  simpa [mu_id] using this

/-! ### Theorem 2: the placement of the points gives no loophole -/

/-- Truncated subtraction is superadditive: `|s| * n - ∑ a ≤ ∑ (n - a)`. -/
theorem card_mul_sub_sum_le {α : Type*} [DecidableEq α] (s : Finset α) (n : ℕ) (a : α → ℕ) :
    s.card * n - ∑ j ∈ s, a j ≤ ∑ j ∈ s, (n - a j) := by
  induction s using Finset.induction_on with
  | empty => simp
  | insert i s hi ih =>
    rw [Finset.sum_insert hi, Finset.sum_insert hi, Finset.card_insert_of_notMem hi]
    have : (s.card + 1) * n = s.card * n + n := by ring
    omega

/-- Cyclic shift of the cube `{0,1}^(L+1)`. -/
def shift {L : ℕ} (j : Fin (L + 1)) : (Fin (L + 1) → Bool) ≃ (Fin (L + 1) → Bool) where
  toFun x := fun i => x (i + j)
  invFun x := fun i => x (i - j)
  left_inv x := by funext i; simp
  right_inv x := by funext i; simp

theorem sum_mu_shift {L : ℕ} (c : ι → Fin (L + 1)) (x : Fin (L + 1) → Bool) :
    ∑ j, mu c (shift j x) = Fintype.card ι * SqrtLog.weight x := by
  classical
  have h1 : ∀ j, mu c (shift j x) = ∑ k : ι, if x (c k + j) = true then 1 else 0 := by
    intro j
    rw [mu, Finset.card_filter]
    rfl
  have h2 : ∀ k : ι, (∑ j : Fin (L + 1), if x (c k + j) = true then 1 else 0)
      = SqrtLog.weight x := by
    intro k
    rw [SqrtLog.weight, Finset.card_filter]
    exact Equiv.sum_comp (Equiv.addLeft (c k)) (fun j => if x j = true then 1 else 0)
  simp_rw [h1]
  rw [Finset.sum_comm]
  simp_rw [h2]
  simp

/-- **Theorem 2 (combinatorial part).** For every placement `c` of the points,
`ℓ * ∑ x, max(0, n - 2 μ x) ≥ ∑ x, max(0, ℓ n - 2 N |x|)`: the bound of Theorem 1 is at least
its value at the (real) balanced placement `N / ℓ` points per direction. -/
theorem placement_le (c : ι → Fin ℓ) (n : ℕ) :
    ∑ x : Fin ℓ → Bool, (ℓ * n - 2 * (Fintype.card ι * SqrtLog.weight x))
      ≤ ℓ * ∑ x : Fin ℓ → Bool, (n - 2 * mu c x) := by
  classical
  cases ℓ with
  | zero => simp
  | succ L =>
    have h1 : (L + 1) * ∑ x : Fin (L + 1) → Bool, (n - 2 * mu c x)
        = ∑ x : Fin (L + 1) → Bool, ∑ j : Fin (L + 1), (n - 2 * mu c (shift j x)) := by
      rw [Finset.sum_comm]
      have : ∀ j : Fin (L + 1), ∑ x : Fin (L + 1) → Bool, (n - 2 * mu c (shift j x))
          = ∑ x : Fin (L + 1) → Bool, (n - 2 * mu c x) :=
        fun j => Equiv.sum_comp (shift j) (fun x => n - 2 * mu c x)
      simp [this]
    rw [h1]
    refine Finset.sum_le_sum fun x _ => ?_
    have h2 := card_mul_sub_sum_le (Finset.univ : Finset (Fin (L + 1))) n
      (fun j => 2 * mu c (shift j x))
    rw [← Finset.mul_sum, sum_mu_shift] at h2
    simpa using h2

/-- **Theorem 2 (no loophole from the placement).** Whatever the placement `c` of the `N` points
on the `ℓ` directions, `ℓ * ∑ k, rank (E k) ≥ ∑ x, max(0, ℓ (N - d) - 2 N |x|)`; dividing by
`ℓ N 2^ℓ` this is `(1/N) ∑ k, ρ k ≥ (1/ℓ) E[(S_ℓ - d ℓ / N)_+]`. -/
theorem sum_rank_ge_balanced (ℓ d : ℕ) (J : Matrix (Fin ℓ → Bool) (Fin ℓ → Bool) F[X])
    (hd : ∀ x y, (J x y).natDegree ≤ d) (h0 : J.map (eval 0) = XtensorF F ℓ)
    (τ : ι → F) (hτ : Function.Injective τ) (hτ0 : ∀ k, τ k ≠ 0) (c : ι → Fin ℓ)
    (K E : ι → Matrix (Fin ℓ → Bool) (Fin ℓ → Bool) F)
    (hdec : ∀ k, J.map (eval (τ k)) = K k + E k)
    (hnice : ∀ k x y, x (c k) ≠ y (c k) → K k x y = 0) :
    ∑ x : Fin ℓ → Bool, (ℓ * (Fintype.card ι - d) - 2 * (Fintype.card ι * SqrtLog.weight x))
      ≤ ℓ * ∑ k, (E k).rank :=
  (placement_le c (Fintype.card ι - d)).trans
    (Nat.mul_le_mul_left ℓ (sum_rank_ge ℓ d J hd h0 τ hτ hτ0 c K E hdec hnice))

end MultiPoint

#print axioms MultiPoint.core
#print axioms MultiPoint.sum_rank_ge_half
#print axioms MultiPoint.sum_rank_ge
#print axioms MultiPoint.sum_rank_block01_ge
#print axioms MultiPoint.sum_rank_ge_simple
#print axioms MultiPoint.sum_rank_ge_one
#print axioms MultiPoint.placement_le
#print axioms MultiPoint.sum_rank_ge_balanced
