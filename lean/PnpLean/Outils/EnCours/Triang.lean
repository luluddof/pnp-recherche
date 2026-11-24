import Mathlib.LinearAlgebra.Matrix.Rank
import Mathlib.LinearAlgebra.Matrix.Kronecker
import Mathlib.LinearAlgebra.Matrix.Transvection
import Mathlib.RingTheory.Nilpotent.Basic
import Mathlib.Algebra.BigOperators.Intervals
import Mathlib.Algebra.BigOperators.Fin
import Mathlib.Order.Interval.Finset.Fin

open Matrix
open scoped Kronecker

namespace Triang

variable {K : Type*} [Field K]

/-- `Swap (a,b) (c,d) = [a = d ∧ b = c]`, i.e. `e_x ⊗ e_y ↦ e_y ⊗ e_x`. -/
def Swap (K : Type*) [Field K] (N : ℕ) : Matrix (Fin N × Fin N) (Fin N × Fin N) K :=
  Matrix.of fun a b => if a.1 = b.2 ∧ a.2 = b.1 then 1 else 0

variable {N : ℕ}

theorem mul_swap_apply (M : Matrix (Fin N × Fin N) (Fin N × Fin N) K) (p q : Fin N × Fin N) :
    (M * Swap K N) p q = M p (q.2, q.1) := by
  rw [Matrix.mul_apply, Finset.sum_eq_single (q.2, q.1)]
  · simp [Swap]
  · intro r _ hr
    have : ¬ (r.1 = q.2 ∧ r.2 = q.1) := by
      rintro ⟨h1, h2⟩; exact hr (Prod.ext h1 h2)
    simp [Swap, this]
  · simp

theorem swap_mul_apply (M : Matrix (Fin N × Fin N) (Fin N × Fin N) K) (p q : Fin N × Fin N) :
    (Swap K N * M) p q = M (p.2, p.1) q := by
  rw [Matrix.mul_apply, Finset.sum_eq_single (p.2, p.1)]
  · simp [Swap]
  · intro r _ hr
    have : ¬ (p.1 = r.2 ∧ p.2 = r.1) := by
      rintro ⟨h1, h2⟩; exact hr (Prod.ext h2.symm h1.symm)
    simp [Swap, this]
  · simp

theorem swap_mul_swap : Swap K N * Swap K N = 1 := by
  ext p q
  rw [mul_swap_apply]
  simp only [Swap, Matrix.of_apply, Matrix.one_apply]
  by_cases h : p = q
  · subst h; simp
  · have : ¬ (p.1 = q.1 ∧ p.2 = q.2) := fun ⟨h1, h2⟩ => h (Prod.ext h1 h2)
    simp [this, h]

theorem kron_mul_swap (A B : Matrix (Fin N) (Fin N) K) :
    (A ⊗ₖ B) * Swap K N = Swap K N * (B ⊗ₖ A) := by
  ext p q
  rw [mul_swap_apply, swap_mul_apply]
  simp [Matrix.kroneckerMap_apply, mul_comm]

/-- The "lower pairs" `{(b, j) : j < b}`. -/
abbrev Om (N : ℕ) := {p : Fin N × Fin N // p.2 < p.1}

theorem card_Om : Fintype.card (Om N) = N * (N - 1) / 2 := by
  rw [Fintype.card_subtype, Finset.card_eq_sum_ones, Finset.sum_filter, Fintype.sum_prod_type]
  have h : ∀ i : Fin N, (∑ j : Fin N, if j < i then 1 else 0) = (i : ℕ) := by
    intro i
    rw [← Finset.sum_filter, ← Finset.card_eq_sum_ones]
    have : (Finset.univ.filter fun j : Fin N => j < i) = Finset.Iio i := by
      ext j; simp
    rw [this, Fin.card_Iio]
  simp_rw [h]
  rw [Fin.sum_univ_eq_sum_range (fun i => i) N, Finset.sum_range_id]

/-- Core of Lemma T: upper-triangular `E_l`, arbitrary `C_l`. -/
theorem card_Om_le_rank_of_upper {ι : Type*} [Fintype ι]
    (E C : ι → Matrix (Fin N) (Fin N) K)
    (hE : ∀ l i j, j < i → E l i j = 0) :
    Fintype.card (Om N) ≤ (Swap K N - ∑ l, E l ⊗ₖ C l).rank := by
  classical
  set Kp := ∑ l, E l ⊗ₖ C l with hKp
  set Θ := Kp * Swap K N with hΘ
  have h1 : (Swap K N - Kp) * Swap K N = 1 - Θ := by
    rw [Matrix.sub_mul, swap_mul_swap]
  have hΘapp : ∀ a b c d : Fin N, Θ (a, b) (c, d) = ∑ l, E l a d * C l b c := by
    intro a b c d
    rw [hΘ, mul_swap_apply, hKp, Matrix.sum_apply]
    simp [Matrix.kroneckerMap_apply]
  let T : Matrix (Om N) (Om N) K := Θ.submatrix Subtype.val Subtype.val
  have hT : ∀ p q : Om N, T p q ≠ 0 → p.1.2 < q.1.2 := by
    intro p q hpq
    by_contra hlt
    push_neg at hlt
    apply hpq
    show Θ p.1 q.1 = 0
    obtain ⟨⟨a, b⟩, hab⟩ := p
    obtain ⟨⟨c, d⟩, hcd⟩ := q
    simp only at hab hcd hlt
    rw [hΘapp]
    refine Finset.sum_eq_zero fun l _ => ?_
    rw [hE l a d (lt_of_le_of_lt hlt hab), zero_mul]
  have hpow : ∀ k : ℕ, ∀ p q : Om N, (T ^ k) p q ≠ 0 → (p.1.2 : ℕ) + k ≤ q.1.2 := by
    intro k
    induction k with
    | zero =>
      intro p q h
      rw [pow_zero, Matrix.one_apply] at h
      by_cases hpq : p = q
      · subst hpq; simp
      · simp [hpq] at h
    | succ k ih =>
      intro p q h
      rw [pow_succ, Matrix.mul_apply] at h
      obtain ⟨r, _, hr⟩ := Finset.exists_ne_zero_of_sum_ne_zero h
      have h1 := ih p r (left_ne_zero_of_mul hr)
      have h2 := hT r q (right_ne_zero_of_mul hr)
      have : (r.1.2 : ℕ) < q.1.2 := h2
      omega
  have hnil : IsNilpotent T := by
    refine ⟨N, ?_⟩
    ext p q
    by_contra h
    have := hpow N p q h
    have := q.1.2.isLt
    simp only [Matrix.zero_apply] at *
    omega
  have hsub : (1 - Θ).submatrix Subtype.val Subtype.val = 1 - T := by
    ext p q
    simp only [Matrix.submatrix_apply, Matrix.sub_apply, Matrix.one_apply, T]
    by_cases hpq : p = q
    · subst hpq; simp
    · have : p.1 ≠ q.1 := fun h => hpq (Subtype.ext h)
      simp [hpq, this]
  calc Fintype.card (Om N) = (1 - T).rank := (Matrix.rank_of_isUnit _ hnil.isUnit_one_sub).symm
    _ = ((1 - Θ).submatrix Subtype.val Subtype.val).rank := by rw [hsub]
    _ ≤ (1 - Θ).rank := Matrix.rank_submatrix_le _ _ _
    _ = ((Swap K N - Kp) * Swap K N).rank := by rw [h1]
    _ ≤ (Swap K N - Kp).rank := Matrix.rank_mul_le_left _ _

end Triang
