import PnpLean.Article1.EnCours.UAvg

/-!
# Anticoncentration of sums of independent Bernoulli variables (finite form)

Lemma "Anticoncentration" of the paper: if `X = ∑ i, 1[β i (ω i)]` on the uniform product
space `∀ i, Ω i` (so that `X` is a sum of independent Bernoulli variables with parameters
`p i`), `V = ∑ i, p i (1 - p i)` and `a : ℝ`, then
`𝔼 |X - a| ≥ c0 · min (V, √V)` with `c0 = (1 - e^{-2}) / 4`.
-/

open Finset

set_option linter.unusedSectionVars false

namespace UAvg

/-- Cauchy–Schwarz for uniform averages. -/
theorem uavg_cs {α : Type*} [Fintype α] {r f g : α → ℝ} (hf : ∀ x, 0 ≤ f x)
    (hg : ∀ x, 0 ≤ g x) (h : ∀ x, r x ^ 2 ≤ f x * g x) : uavg r ^ 2 ≤ uavg f * uavg g := by
  unfold uavg
  rw [div_pow, div_mul_div_comm, ← sq]
  exact div_le_div_of_nonneg_right
    (Finset.sum_sq_le_sum_mul_sum_of_sq_le_mul _ (fun x _ => hf x) (fun x _ => hg x)
      (fun x _ => h x)) (sq_nonneg _)

/-! ## Sums of independent symmetric `{-1,0,1}`-valued variables -/

section sym

variable {ι : Type*} [Fintype ι] [DecidableEq ι] {Θ : ι → Type*} [∀ i, Fintype (Θ i)]
  [∀ i, Nonempty (Θ i)]

/-- Partial sum `D_s = ∑_{i ∈ s} ξ_i`. -/
def D (ξ : ∀ i, Θ i → ℤ) (s : Finset ι) (ω : ∀ i, Θ i) : ℤ := ∑ i ∈ s, ξ i (ω i)

/-- `ξ_i` takes values in `{-1,0,1}` and `P(ξ_i = 1) = P(ξ_i = -1) = τ_i / 2 ≤ 1/3`. -/
structure SymHyp (ξ : ∀ i, Θ i → ℤ) (τ : ι → ℝ) : Prop where
  hval : ∀ i x, ξ i x = 0 ∨ ξ i x = 1 ∨ ξ i x = -1
  hpos : ∀ i, uavg (fun x => if ξ i x = 1 then (1 : ℝ) else 0) = τ i / 2
  hneg : ∀ i, uavg (fun x => if ξ i x = -1 then (1 : ℝ) else 0) = τ i / 2
  hτ : ∀ i, τ i ≤ 2 / 3

variable {ξ : ∀ i, Θ i → ℤ} {τ : ι → ℝ}

theorem D_update (s : Finset ι) (a : ι) (ha : a ∉ s) (ω : ∀ i, Θ i) (x : Θ a) :
    D ξ s (Function.update ω a x) = D ξ s ω := by
  unfold D
  refine Finset.sum_congr rfl fun i hi => ?_
  have : i ≠ a := fun h => ha (h ▸ hi)
  rw [Function.update_of_ne this]

theorem D_insert (s : Finset ι) (a : ι) (ha : a ∉ s) (ω : ∀ i, Θ i) :
    D ξ (insert a s) ω = ξ a (ω a) + D ξ s ω := by
  unfold D; rw [Finset.sum_insert ha]

namespace SymHyp

variable (H : SymHyp ξ τ)
include H

theorem τ_nonneg (i : ι) : 0 ≤ τ i := by
  have := H.hpos i
  have h0 : 0 ≤ uavg (fun x => if ξ i x = 1 then (1 : ℝ) else 0) :=
    uavg_nonneg fun x => by split_ifs <;> norm_num
  linarith

theorem pzero (i : ι) : uavg (fun x => if ξ i x = 0 then (1 : ℝ) else 0) = 1 - τ i := by
  have e : ∀ x, (if ξ i x = 0 then (1 : ℝ) else 0) =
      1 - (if ξ i x = 1 then (1 : ℝ) else 0) - (if ξ i x = -1 then (1 : ℝ) else 0) := by
    intro x
    rcases H.hval i x with h | h | h <;> simp [h]
  rw [uavg_congr e, uavg_sub, uavg_sub, uavg_const, H.hpos, H.hneg]; ring

theorem mom_odd (i : ι) (k : ℕ) (hk : Odd k) : uavg (fun x => (ξ i x : ℝ) ^ k) = 0 := by
  have e : ∀ x, (ξ i x : ℝ) ^ k =
      (if ξ i x = 1 then (1 : ℝ) else 0) - (if ξ i x = -1 then (1 : ℝ) else 0) := by
    intro x
    rcases H.hval i x with h | h | h
    · simp [h, hk.pos.ne']
    · simp [h]
    · simp [h, hk.neg_one_pow]
  rw [uavg_congr e, uavg_sub, H.hpos, H.hneg]; ring

theorem mom_even (i : ι) (k : ℕ) (hk : Even k) (hk0 : k ≠ 0) :
    uavg (fun x => (ξ i x : ℝ) ^ k) = τ i := by
  have e : ∀ x, (ξ i x : ℝ) ^ k =
      (if ξ i x = 1 then (1 : ℝ) else 0) + (if ξ i x = -1 then (1 : ℝ) else 0) := by
    intro x
    rcases H.hval i x with h | h | h
    · simp [h, hk0]
    · simp [h]
    · simp [h, hk.neg_one_pow]
  rw [uavg_congr e, uavg_add, H.hpos, H.hneg]; ring

omit H in
/-- Independence of `ξ_a` and `D_s` for `a ∉ s`, for monomials. -/
theorem indep_pow (s : Finset ι) (a : ι) (ha : a ∉ s) (k j : ℕ) :
    uavg (fun ω : (∀ i, Θ i) => (ξ a (ω a) : ℝ) ^ k * (D ξ s ω : ℝ) ^ j)
      = uavg (fun x => (ξ a x : ℝ) ^ k) * uavg (fun ω : (∀ i, Θ i) => (D ξ s ω : ℝ) ^ j) :=
  uavg_pi_mul a (fun x => (ξ a x : ℝ) ^ k) (fun ω => (D ξ s ω : ℝ) ^ j)
    (fun ω x => by rw [D_update s a ha])

omit H in
/-- Independence of `ξ_a` and `D_s` for `a ∉ s`, for indicators. -/
theorem indep_ind (s : Finset ι) (a : ι) (ha : a ∉ s) (k j : ℤ) :
    uavg (fun ω : (∀ i, Θ i) => (if ξ a (ω a) = k then (1 : ℝ) else 0) *
        (if D ξ s ω = j then (1 : ℝ) else 0))
      = uavg (fun x => if ξ a x = k then (1 : ℝ) else 0) *
        uavg (fun ω : (∀ i, Θ i) => if D ξ s ω = j then (1 : ℝ) else 0) :=
  uavg_pi_mul a (fun x => if ξ a x = k then (1 : ℝ) else 0)
    (fun ω => if D ξ s ω = j then (1 : ℝ) else 0)
    (fun ω x => by rw [D_update s a ha])

/-- Second moment: `𝔼 D_s² = ∑_{i∈s} τ_i`. -/
theorem second_moment (s : Finset ι) :
    uavg (fun ω : (∀ i, Θ i) => (D ξ s ω : ℝ) ^ 2) = ∑ i ∈ s, τ i := by
  induction s using Finset.induction_on with
  | empty => simp [D, uavg_const]
  | insert a s ha ih =>
    have e : ∀ ω : (∀ i, Θ i), (D ξ (insert a s) ω : ℝ) ^ 2 =
        (ξ a (ω a) : ℝ) ^ 2 * (D ξ s ω : ℝ) ^ 0 + 2 * ((ξ a (ω a) : ℝ) ^ 1 * (D ξ s ω : ℝ) ^ 1)
          + (ξ a (ω a) : ℝ) ^ 0 * (D ξ s ω : ℝ) ^ 2 := by
      intro ω; rw [D_insert s a ha]; push_cast; ring
    rw [uavg_congr e, uavg_add, uavg_add, uavg_const_mul, indep_pow s a ha, indep_pow s a ha,
      indep_pow s a ha, H.mom_odd a 1 (by decide), H.mom_even a 2 (by decide) (by decide), ih,
      Finset.sum_insert ha]
    simp [uavg_const]

/-- Fourth moment: `𝔼 D_s⁴ ≤ T + 3T²`, `T = ∑_{i∈s} τ_i`. -/
theorem fourth_moment (s : Finset ι) :
    uavg (fun ω : (∀ i, Θ i) => (D ξ s ω : ℝ) ^ 4) ≤ ∑ i ∈ s, τ i + 3 * (∑ i ∈ s, τ i) ^ 2 := by
  induction s using Finset.induction_on with
  | empty => simp [D, uavg_const]
  | insert a s ha ih =>
    have e : ∀ ω : (∀ i, Θ i), (D ξ (insert a s) ω : ℝ) ^ 4 =
        (ξ a (ω a) : ℝ) ^ 4 * (D ξ s ω : ℝ) ^ 0 + 4 * ((ξ a (ω a) : ℝ) ^ 3 * (D ξ s ω : ℝ) ^ 1)
          + 6 * ((ξ a (ω a) : ℝ) ^ 2 * (D ξ s ω : ℝ) ^ 2)
          + 4 * ((ξ a (ω a) : ℝ) ^ 1 * (D ξ s ω : ℝ) ^ 3)
          + (ξ a (ω a) : ℝ) ^ 0 * (D ξ s ω : ℝ) ^ 4 := by
      intro ω; rw [D_insert s a ha]; push_cast; ring
    rw [uavg_congr e, uavg_add, uavg_add, uavg_add, uavg_add, uavg_const_mul, uavg_const_mul,
      uavg_const_mul, indep_pow s a ha, indep_pow s a ha, indep_pow s a ha, indep_pow s a ha,
      indep_pow s a ha, H.mom_odd a 1 (by decide), H.mom_odd a 3 (by decide),
      H.mom_even a 2 (by decide) (by decide), H.mom_even a 4 (by decide) (by decide),
      H.second_moment s, Finset.sum_insert ha]
    have hτ := H.τ_nonneg a
    have hT : 0 ≤ ∑ i ∈ s, τ i := Finset.sum_nonneg fun i _ => H.τ_nonneg i
    simp only [pow_zero, uavg_const, mul_one, zero_mul, mul_zero, add_zero, one_mul]
    nlinarith [sq_nonneg (τ a), mul_nonneg hτ hT]

/-- Recursion for `P(D_s = 0)`. -/
theorem pzero_insert (s : Finset ι) (a : ι) (ha : a ∉ s) :
    uavg (fun ω : (∀ i, Θ i) => if D ξ (insert a s) ω = 0 then (1 : ℝ) else 0)
      = (1 - τ a) * uavg (fun ω : (∀ i, Θ i) => if D ξ s ω = 0 then (1 : ℝ) else 0)
        + τ a / 2 * uavg (fun ω : (∀ i, Θ i) => if D ξ s ω = -1 then (1 : ℝ) else 0)
        + τ a / 2 * uavg (fun ω : (∀ i, Θ i) => if D ξ s ω = 1 then (1 : ℝ) else 0) := by
  have key : ∀ x d : ℤ, (x = 0 ∨ x = 1 ∨ x = -1) →
      (if x + d = 0 then (1 : ℝ) else 0) =
        (if x = 0 then (1 : ℝ) else 0) * (if d = 0 then (1 : ℝ) else 0)
        + (if x = 1 then (1 : ℝ) else 0) * (if d = -1 then (1 : ℝ) else 0)
        + (if x = -1 then (1 : ℝ) else 0) * (if d = 1 then (1 : ℝ) else 0) := by
    intro x d h
    have e1 : (1 + d = 0) ↔ d = -1 := by omega
    have e2 : (-1 + d = 0) ↔ d = 1 := by omega
    rcases h with rfl | rfl | rfl <;> simp [e1, e2]
  have e : ∀ ω : (∀ i, Θ i), (if D ξ (insert a s) ω = 0 then (1 : ℝ) else 0) =
      (if ξ a (ω a) = 0 then (1 : ℝ) else 0) * (if D ξ s ω = 0 then (1 : ℝ) else 0)
        + (if ξ a (ω a) = 1 then (1 : ℝ) else 0) * (if D ξ s ω = -1 then (1 : ℝ) else 0)
        + (if ξ a (ω a) = -1 then (1 : ℝ) else 0) * (if D ξ s ω = 1 then (1 : ℝ) else 0) := by
    intro ω; rw [D_insert s a ha]; exact key _ _ (H.hval a (ω a))
  rw [uavg_congr e, uavg_add, uavg_add, indep_ind s a ha, indep_ind s a ha, indep_ind s a ha,
    H.pzero, H.hpos, H.hneg]

/-- `P(D_s = 0) ≤ (1 + ∏_{i∈s} (1 - τ_i)) / 2`. -/
theorem pzero_le (s : Finset ι) :
    uavg (fun ω : (∀ i, Θ i) => if D ξ s ω = 0 then (1 : ℝ) else 0)
      ≤ (1 + ∏ i ∈ s, (1 - τ i)) / 2 := by
  induction s using Finset.induction_on with
  | empty => simp [D, uavg_const]
  | insert a s ha ih =>
    rw [H.pzero_insert s a ha, Finset.prod_insert ha]
    set P0 := uavg (fun ω : (∀ i, Θ i) => if D ξ s ω = 0 then (1 : ℝ) else 0)
    set Pm := uavg (fun ω : (∀ i, Θ i) => if D ξ s ω = -1 then (1 : ℝ) else 0)
    set Pp := uavg (fun ω : (∀ i, Θ i) => if D ξ s ω = 1 then (1 : ℝ) else 0)
    set Pr := ∏ i ∈ s, (1 - τ i)
    have hsum : P0 + Pm + Pp ≤ 1 := by
      have h1 : uavg (fun ω : (∀ i, Θ i) => (if D ξ s ω = 0 then (1 : ℝ) else 0)
          + (if D ξ s ω = -1 then (1 : ℝ) else 0) + (if D ξ s ω = 1 then (1 : ℝ) else 0))
          ≤ uavg (fun _ : (∀ i, Θ i) => (1 : ℝ)) := by
        refine uavg_mono fun ω => ?_
        by_cases h0 : D ξ s ω = 0
        · simp [h0]
        · by_cases h1 : D ξ s ω = -1
          · simp [h1]
          · by_cases h2 : D ξ s ω = 1
            · simp [h2]
            · simp [h0, h1, h2]
      rwa [uavg_add, uavg_add, uavg_const] at h1
    have hτ0 := H.τ_nonneg a
    have hτ1 := H.hτ a
    have hPr : 0 ≤ Pr := Finset.prod_nonneg fun i _ => by linarith [H.hτ i]
    have hmp : τ a / 2 * (Pm + Pp) ≤ τ a / 2 * (1 - P0) :=
      mul_le_mul_of_nonneg_left (by linarith) (by linarith)
    have h3 : (1 - 3 * τ a / 2) * P0 ≤ (1 - 3 * τ a / 2) * ((1 + Pr) / 2) :=
      mul_le_mul_of_nonneg_left ih (by linarith)
    nlinarith [mul_nonneg hτ0 hPr]

/-- Step (i): `𝔼 |D| ≥ P(D ≠ 0) ≥ (1 - ∏ (1 - τ_i)) / 2`. -/
theorem abs_lower :
    (1 - ∏ i, (1 - τ i)) / 2 ≤ uavg (fun ω : (∀ i, Θ i) => |(D ξ univ ω : ℝ)|) := by
  have h1 : uavg (fun ω : (∀ i, Θ i) => (1 : ℝ) - (if D ξ univ ω = 0 then (1 : ℝ) else 0))
      ≤ uavg (fun ω : (∀ i, Θ i) => |(D ξ univ ω : ℝ)|) := by
    refine uavg_mono fun ω => ?_
    by_cases h0 : D ξ univ ω = 0
    · simp [h0]
    · have : (1 : ℝ) ≤ |(D ξ univ ω : ℝ)| := by exact_mod_cast Int.one_le_abs h0
      simpa [h0] using this
  rw [uavg_sub, uavg_const] at h1
  have := H.pzero_le univ
  linarith

/-- Step (ii) (Hölder via two Cauchy–Schwarz): `T³ ≤ (𝔼|D|)² (T + 3T²)`. -/
theorem holder :
    (∑ i, τ i) ^ 3 ≤ uavg (fun ω : (∀ i, Θ i) => |(D ξ univ ω : ℝ)|) ^ 2 *
      (∑ i, τ i + 3 * (∑ i, τ i) ^ 2) := by
  set T := ∑ i, τ i with hT
  set A1 := uavg (fun ω : (∀ i, Θ i) => |(D ξ univ ω : ℝ)|) with hA1
  set A3 := uavg (fun ω : (∀ i, Θ i) => |(D ξ univ ω : ℝ)| ^ 3) with hA3
  set A4 := uavg (fun ω : (∀ i, Θ i) => (D ξ univ ω : ℝ) ^ 4) with hA4
  have h2 : uavg (fun ω : (∀ i, Θ i) => (D ξ univ ω : ℝ) ^ 2) = T := H.second_moment univ
  have h4 : A4 ≤ T + 3 * T ^ 2 := H.fourth_moment univ
  have hT0 : 0 ≤ T := Finset.sum_nonneg fun i _ => H.τ_nonneg i
  have hA1' : 0 ≤ A1 := uavg_nonneg fun ω => abs_nonneg _
  have hA3' : 0 ≤ A3 := uavg_nonneg fun ω => by positivity
  have hA4' : 0 ≤ A4 := uavg_nonneg fun ω => by positivity
  have cs1 : T ^ 2 ≤ A1 * A3 := by
    rw [← h2]
    refine uavg_cs (fun ω => abs_nonneg _) (fun ω => by positivity) (fun ω => ?_)
    rw [← sq_abs (D ξ univ ω : ℝ)]
    exact le_of_eq (by ring)
  have cs2 : A3 ^ 2 ≤ T * A4 := by
    rw [← h2]
    refine uavg_cs (fun ω => by positivity) (fun ω => by positivity) (fun ω => ?_)
    have : |(D ξ univ ω : ℝ)| ^ 2 = (D ξ univ ω : ℝ) ^ 2 := sq_abs _
    refine le_of_eq ?_
    calc (|(D ξ univ ω : ℝ)| ^ 3) ^ 2 = (|(D ξ univ ω : ℝ)| ^ 2) ^ 3 := by ring
      _ = ((D ξ univ ω : ℝ) ^ 2) ^ 3 := by rw [this]
      _ = (D ξ univ ω : ℝ) ^ 2 * (D ξ univ ω : ℝ) ^ 4 := by ring
  rcases hT0.eq_or_lt with h0 | hpos
  · rw [← h0]; simp
  · -- T⁴ ≤ A1² A3² ≤ A1² T A4
    have h5 : T ^ 4 ≤ A1 ^ 2 * (T * A4) := by
      calc T ^ 4 = (T ^ 2) ^ 2 := by ring
        _ ≤ (A1 * A3) ^ 2 := pow_le_pow_left₀ (by positivity) cs1 2
        _ = A1 ^ 2 * A3 ^ 2 := by ring
        _ ≤ A1 ^ 2 * (T * A4) := mul_le_mul_of_nonneg_left cs2 (by positivity)
    have h6 : T ^ 3 ≤ A1 ^ 2 * A4 := by
      have : T * T ^ 3 ≤ T * (A1 ^ 2 * A4) := by nlinarith [h5]
      exact le_of_mul_le_mul_left this hpos
    exact h6.trans (mul_le_mul_of_nonneg_left h4 (by positivity))

end SymHyp

end sym

/-! ## Bernoulli sums -/

section bernoulli

variable {ι : Type*} [Fintype ι] [DecidableEq ι] {Ω : ι → Type*} [∀ i, Fintype (Ω i)]
  [∀ i, Nonempty (Ω i)]

/-- The constant `c0 = (1 - e^{-2}) / 4` of the anticoncentration lemma. -/
noncomputable def c0 : ℝ := (1 - Real.exp (-2)) / 4

theorem c0_pos : 0 < c0 := by
  unfold c0
  have : Real.exp (-2) < 1 := by
    have := Real.exp_lt_exp.2 (show (-2 : ℝ) < 0 by norm_num)
    rwa [Real.exp_zero] at this
  linarith

theorem c0_le : c0 ≤ 1 / 4 := by
  unfold c0
  have := Real.exp_pos (-2)
  linarith

variable (β : ∀ i, Ω i → Prop) [∀ i, DecidablePred (β i)]

/-- Bernoulli parameter `p_i = P(β_i)`. -/
noncomputable def bp (i : ι) : ℝ := uavg (fun x : Ω i => if β i x then (1 : ℝ) else 0)

/-- Variance `V = ∑ p_i (1 - p_i)`. -/
noncomputable def bV : ℝ := ∑ i, bp β i * (1 - bp β i)

/-- The Bernoulli sum `X(ω) = ∑_i 1[β_i (ω_i)]`. -/
def bX (ω : ∀ i, Ω i) : ℝ := ∑ i, if β i (ω i) then (1 : ℝ) else 0

/-- Symmetrized increments `ξ_i (x, y) = 1[β_i x] - 1[β_i y]`. -/
def bξ (i : ι) (p : Ω i × Ω i) : ℤ := (if β i p.1 then 1 else 0) - (if β i p.2 then 1 else 0)

theorem bp_nonneg (i : ι) : 0 ≤ bp β i := uavg_nonneg fun x => by split_ifs <;> norm_num

theorem bp_le_one (i : ι) : bp β i ≤ 1 := by
  have : bp β i ≤ uavg (fun _ : Ω i => (1 : ℝ)) :=
    uavg_mono fun x => by split_ifs <;> norm_num
  rwa [uavg_const] at this

theorem bV_nonneg : 0 ≤ bV β :=
  Finset.sum_nonneg fun i _ => mul_nonneg (bp_nonneg β i) (by linarith [bp_le_one β i])

theorem bsym : SymHyp (Θ := fun i => Ω i × Ω i) (bξ β) (fun i => 2 * (bp β i * (1 - bp β i))) where
  hval i p := by unfold bξ; split_ifs <;> simp
  hpos i := by
    have e : ∀ p : Ω i × Ω i, (if bξ β i p = 1 then (1 : ℝ) else 0) =
        (if β i p.1 then (1 : ℝ) else 0) * (1 - (if β i p.2 then (1 : ℝ) else 0)) := by
      intro p; unfold bξ; split_ifs <;> simp_all
    have := uavg_prod_mul (fun x : Ω i => if β i x then (1 : ℝ) else 0)
      (fun y : Ω i => 1 - (if β i y then (1 : ℝ) else 0))
    rw [uavg_congr e, this, uavg_sub, uavg_const]; unfold bp; ring
  hneg i := by
    have e : ∀ p : Ω i × Ω i, (if bξ β i p = -1 then (1 : ℝ) else 0) =
        (1 - (if β i p.1 then (1 : ℝ) else 0)) * (if β i p.2 then (1 : ℝ) else 0) := by
      intro p; unfold bξ; split_ifs <;> simp_all
    have := uavg_prod_mul (fun x : Ω i => 1 - (if β i x then (1 : ℝ) else 0))
      (fun y : Ω i => if β i y then (1 : ℝ) else 0)
    rw [uavg_congr e, this, uavg_sub, uavg_const]; unfold bp; ring
  hτ i := by nlinarith [sq_nonneg (bp β i - 1 / 2)]

/-- Pairs of points of the product space = points of the product of the pair spaces. -/
def pairEquiv : (∀ i, Ω i × Ω i) ≃ (∀ i, Ω i) × (∀ i, Ω i) where
  toFun θ := (fun i => (θ i).1, fun i => (θ i).2)
  invFun p := fun i => (p.1 i, p.2 i)
  left_inv _ := rfl
  right_inv _ := rfl

/-- Symmetrization: `𝔼 |D| ≤ 2 𝔼 |X - a|`. -/
theorem symmetrization (a : ℝ) :
    uavg (fun θ : (∀ i, Ω i × Ω i) => |(D (bξ β) univ θ : ℝ)|)
      ≤ 2 * uavg (fun ω : (∀ i, Ω i) => |bX β ω - a|) := by
  have e : ∀ θ : (∀ i, Ω i × Ω i), |(D (bξ β) univ θ : ℝ)| =
      (fun p : (∀ i, Ω i) × (∀ i, Ω i) => |bX β p.1 - bX β p.2|) (pairEquiv θ) := by
    intro θ
    congr 1
    simp only [D, bξ, bX, pairEquiv, Equiv.coe_fn_mk]
    push_cast
    rw [Finset.sum_sub_distrib]
    rfl
  rw [uavg_congr e,
    uavg_equiv pairEquiv (fun p : (∀ i, Ω i) × (∀ i, Ω i) => |bX β p.1 - bX β p.2|)]
  have h1 : uavg (fun p : (∀ i, Ω i) × (∀ i, Ω i) => |bX β p.1 - bX β p.2|)
      ≤ uavg (fun p : (∀ i, Ω i) × (∀ i, Ω i) => |bX β p.1 - a| + |bX β p.2 - a|) :=
    uavg_mono fun p => by
      have := abs_sub_le (bX β p.1) a (bX β p.2)
      rwa [abs_sub_comm a (bX β p.2)] at this
  rw [uavg_add, uavg_fst (fun ω : (∀ i, Ω i) => |bX β ω - a|),
    uavg_snd (fun ω : (∀ i, Ω i) => |bX β ω - a|)] at h1
  linarith

/-- **Anticoncentration lemma** (paper, Lemma "Anticoncentration"): for a sum `X` of
independent Bernoulli variables with variance `V` and any real `a`,
`𝔼 |X - a| ≥ c0 · min (V, √V)`, `c0 = (1 - e^{-2}) / 4`. -/
theorem bernoulli_anticonc (a : ℝ) :
    c0 * min (bV β) (Real.sqrt (bV β)) ≤ uavg (fun ω : (∀ i, Ω i) => |bX β ω - a|) := by
  set E := uavg (fun ω : (∀ i, Ω i) => |bX β ω - a|) with hE
  set V := bV β with hV
  set m1 := uavg (fun θ : (∀ i, Ω i × Ω i) => |(D (bξ β) univ θ : ℝ)|) with hm1
  have H := bsym β
  have hsym : m1 ≤ 2 * E := symmetrization β a
  have hV0 : 0 ≤ V := bV_nonneg β
  have hT : ∑ i, 2 * (bp β i * (1 - bp β i)) = 2 * V := by
    rw [← Finset.mul_sum]; rfl
  have hm0 : 0 ≤ m1 := uavg_nonneg fun θ => abs_nonneg _
  have hc := c0_pos
  rcases le_total V 1 with hV1 | hV1
  · -- small variance
    have h1 := H.abs_lower
    have hprod : ∏ i, (1 - 2 * (bp β i * (1 - bp β i))) ≤ Real.exp (-(2 * V)) := by
      rw [← hT, ← Finset.sum_neg_distrib, Real.exp_sum]
      refine Finset.prod_le_prod₀ (fun i _ => ?_) (fun i _ => Real.one_sub_le_exp_neg _)
      nlinarith [sq_nonneg (bp β i - 1 / 2)]
    have hconv : Real.exp (-(2 * V)) ≤ 1 - (1 - Real.exp (-2)) * V := by
      have := convexOn_exp.2 (Set.mem_univ (0 : ℝ)) (Set.mem_univ (-2 : ℝ))
        (sub_nonneg.2 hV1) hV0 (by ring)
      simp only [smul_eq_mul, Real.exp_zero, mul_zero, zero_add, mul_one] at this
      have e : V * (-2 : ℝ) = -(2 * V) := by ring
      rw [e] at this
      linarith
    have h2 : c0 * V ≤ E := by
      unfold c0
      change (1 - ∏ i, (1 - 2 * (bp β i * (1 - bp β i)))) / 2 ≤ m1 at h1
      linarith
    exact (mul_le_mul_of_nonneg_left (min_le_left _ _) hc.le).trans h2
  · -- large variance
    have h1 := H.holder
    change (∑ i, 2 * (bp β i * (1 - bp β i))) ^ 3 ≤ m1 ^ 2 *
      (∑ i, 2 * (bp β i * (1 - bp β i)) + 3 * (∑ i, 2 * (bp β i * (1 - bp β i))) ^ 2) at h1
    rw [hT] at h1
    have hE0 : 0 ≤ E := uavg_nonneg fun ω => abs_nonneg _
    -- m1² ≥ 4V/7
    have h2 : 4 * V ≤ 7 * m1 ^ 2 := by
      have hV2 : 0 < V ^ 2 := by positivity
      have h3 : (2 * V) ^ 3 ≤ m1 ^ 2 * (14 * V ^ 2) := by
        refine h1.trans (mul_le_mul_of_nonneg_left ?_ (sq_nonneg _))
        nlinarith
      have h4 : V ^ 2 * (4 * V) ≤ V ^ 2 * (7 * m1 ^ 2) := by nlinarith
      exact le_of_mul_le_mul_left h4 hV2
    have h5 : V ≤ 7 * E ^ 2 := by
      have : m1 ^ 2 ≤ (2 * E) ^ 2 := pow_le_pow_left₀ hm0 hsym 2
      nlinarith
    have hs0 := Real.sqrt_nonneg V
    have hs2 := Real.sq_sqrt hV0
    have h6 : Real.sqrt V / 4 ≤ E := by
      by_contra h
      push Not at h
      have h4 : (0 : ℝ) ≤ Real.sqrt V / 4 := div_nonneg hs0 (by norm_num)
      nlinarith [mul_nonneg (sub_nonneg.2 h.le) (add_nonneg h4 hE0)]
    calc c0 * min V (Real.sqrt V) ≤ c0 * Real.sqrt V :=
          mul_le_mul_of_nonneg_left (min_le_right _ _) hc.le
      _ ≤ 1 / 4 * Real.sqrt V := mul_le_mul_of_nonneg_right c0_le hs0
      _ ≤ E := by linarith

end bernoulli

end UAvg

#print axioms UAvg.SymHyp.abs_lower
#print axioms UAvg.SymHyp.holder
#print axioms UAvg.symmetrization
#print axioms UAvg.bernoulli_anticonc
