import PnpLean.Article2.EnCours.MultiPoint
open Polynomial Matrix Module
namespace MultiPoint
open AnyField (XtensorF)
variable {ℓ : ℕ} {F : Type*} [Field F] {ι : Type*} [Fintype ι]

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
    simp [mu, shift, Finset.card_filter]
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
