import Mathlib

/-!
# Uniform averages over finite types, and independence on finite product spaces

`uavg f = (∑ x, f x) / card α` is the expectation of `f` under the uniform distribution on
the finite type `α`.  On a product `∀ i, Ω i` the coordinates are independent; this is
expressed by `uavg_pi_mul` (resampling of one coordinate).
-/

open Finset

set_option linter.unusedSectionVars false

namespace UAvg

/-- Uniform average over a finite type. -/
noncomputable def uavg {α : Type*} [Fintype α] (f : α → ℝ) : ℝ :=
  (∑ x, f x) / Fintype.card α

section basic

variable {α β : Type*} [Fintype α] [Fintype β]

theorem card_pos' [Nonempty α] : (0 : ℝ) < Fintype.card α := by
  exact_mod_cast Fintype.card_pos

theorem uavg_congr {f g : α → ℝ} (h : ∀ x, f x = g x) : uavg f = uavg g := by
  unfold uavg; rw [Finset.sum_congr rfl fun x _ => h x]

theorem uavg_mono {f g : α → ℝ} (h : ∀ x, f x ≤ g x) : uavg f ≤ uavg g :=
  div_le_div_of_nonneg_right (Finset.sum_le_sum fun x _ => h x) (Nat.cast_nonneg _)

theorem uavg_nonneg {f : α → ℝ} (h : ∀ x, 0 ≤ f x) : 0 ≤ uavg f :=
  div_nonneg (Finset.sum_nonneg fun x _ => h x) (Nat.cast_nonneg _)

theorem uavg_const [Nonempty α] (c : ℝ) : uavg (fun _ : α => c) = c := by
  unfold uavg
  rw [Finset.sum_const, Finset.card_univ, nsmul_eq_mul]
  field_simp [(card_pos' (α := α)).ne']

theorem uavg_add (f g : α → ℝ) : uavg (fun x => f x + g x) = uavg f + uavg g := by
  unfold uavg; rw [Finset.sum_add_distrib, add_div]

theorem uavg_sub (f g : α → ℝ) : uavg (fun x => f x - g x) = uavg f - uavg g := by
  unfold uavg; rw [Finset.sum_sub_distrib, sub_div]

theorem uavg_const_mul (c : ℝ) (f : α → ℝ) : uavg (fun x => c * f x) = c * uavg f := by
  unfold uavg; rw [← Finset.mul_sum, mul_div_assoc]

theorem uavg_mul_const (f : α → ℝ) (c : ℝ) : uavg (fun x => f x * c) = uavg f * c := by
  unfold uavg; rw [← Finset.sum_mul, div_mul_eq_mul_div]

theorem uavg_sum {ι : Type*} (s : Finset ι) (f : ι → α → ℝ) :
    uavg (fun x => ∑ i ∈ s, f i x) = ∑ i ∈ s, uavg (f i) := by
  unfold uavg; rw [Finset.sum_comm, Finset.sum_div]

theorem uavg_equiv (e : α ≃ β) (f : β → ℝ) : uavg (fun x => f (e x)) = uavg f := by
  unfold uavg; rw [Fintype.card_congr e, Equiv.sum_comp e f]

/-- Jensen for `|·|`. -/
theorem abs_uavg_le (f : α → ℝ) : |uavg f| ≤ uavg (fun x => |f x|) := by
  unfold uavg
  rw [abs_div, Nat.abs_cast]
  exact div_le_div_of_nonneg_right (Finset.abs_sum_le_sum_abs _ _) (Nat.cast_nonneg _)

/-- Independence of the two factors of a product. -/
theorem uavg_prod_mul (f : α → ℝ) (g : β → ℝ) :
    uavg (fun p : α × β => f p.1 * g p.2) = uavg f * uavg g := by
  unfold uavg
  rw [Fintype.sum_prod_type, Fintype.card_prod, ← Finset.sum_mul_sum]
  push_cast
  rw [div_mul_div_comm]

theorem uavg_fst [Nonempty β] (f : α → ℝ) : uavg (fun p : α × β => f p.1) = uavg f := by
  have := uavg_prod_mul f (fun _ : β => (1 : ℝ))
  rw [uavg_const, mul_one] at this
  rw [← this]; exact uavg_congr fun p => (mul_one _).symm

theorem uavg_snd [Nonempty α] (g : β → ℝ) : uavg (fun p : α × β => g p.2) = uavg g := by
  have := uavg_prod_mul (fun _ : α => (1 : ℝ)) g
  rw [uavg_const, one_mul] at this
  rw [← this]; exact uavg_congr fun p => (one_mul _).symm

end basic

section pi

variable {ι : Type*} [Fintype ι] [DecidableEq ι] {Ω : ι → Type*} [∀ i, Fintype (Ω i)]
  [∀ i, Nonempty (Ω i)]

/-- Resampling one coordinate does not change the uniform sum. -/
theorem sum_resample (a : ι) (h : (∀ i, Ω i) → ℝ) :
    ∑ ω : (∀ i, Ω i), ∑ x : Ω a, h (Function.update ω a x)
      = Fintype.card (Ω a) * ∑ ω : (∀ i, Ω i), h ω := by
  let φ : (∀ i, Ω i) × Ω a → (∀ i, Ω i) × Ω a := fun p => (Function.update p.1 a p.2, p.1 a)
  have hφ : Function.Involutive φ := by
    rintro ⟨ω, x⟩
    simp [φ]
  have h1 : ∑ p : (∀ i, Ω i) × Ω a, h (φ p).1 = ∑ p : (∀ i, Ω i) × Ω a, h p.1 :=
    Equiv.sum_comp (hφ.toPerm φ) (fun p => h p.1)
  rw [Fintype.sum_prod_type, Fintype.sum_prod_type] at h1
  rw [h1]
  simp only [Finset.sum_const, Finset.card_univ, nsmul_eq_mul]
  rw [Finset.mul_sum]

/-- **Independence**: a function of the coordinate `a` is independent of any function that
does not depend on the coordinate `a`. -/
theorem uavg_pi_mul (a : ι) (f : Ω a → ℝ) (g : (∀ i, Ω i) → ℝ)
    (hg : ∀ ω x, g (Function.update ω a x) = g ω) :
    uavg (fun ω : (∀ i, Ω i) => f (ω a) * g ω) = uavg f * uavg g := by
  have h := sum_resample a (fun ω : (∀ i, Ω i) => f (ω a) * g ω)
  simp only [Function.update_self, hg] at h
  have h2 : ∑ ω : (∀ i, Ω i), ∑ x : Ω a, f x * g ω = (∑ x, f x) * ∑ ω, g ω := by
    rw [Finset.sum_mul_sum, Finset.sum_comm]
  rw [h2] at h
  unfold uavg
  have hc : (Fintype.card (Ω a) : ℝ) ≠ 0 := (card_pos' (α := Ω a)).ne'
  have hc' : (Fintype.card (∀ i, Ω i) : ℝ) ≠ 0 := (card_pos' (α := ∀ i, Ω i)).ne'
  rw [div_mul_div_comm, div_eq_div_iff hc' (mul_ne_zero hc hc')]
  rw [h]; ring

/-- A function of one coordinate has its marginal average. -/
theorem uavg_pi_coord (a : ι) (f : Ω a → ℝ) :
    uavg (fun ω : (∀ i, Ω i) => f (ω a)) = uavg f := by
  have := uavg_pi_mul a f (fun _ => (1 : ℝ)) (fun _ _ => rfl)
  rw [uavg_const, mul_one] at this
  rw [← this]; exact uavg_congr fun ω => (mul_one _).symm

end pi

end UAvg
