import PnpLean.Article2.SqrtLog

/-!
# Consequences of the dual polynomial bound `SqrtLog.sum_rank_ge_N`

Elementary (non-asymptotic) corollaries of the main theorem of Article 2:

* lower bounds on `N ℓ m` (`N_ge_sum_choose`, `N_ge_choose`, `N_zero`, `N_one`, `N_mono`);
* pigeonhole / max / average forms (`exists_rank_ge`, `avg_le_max_rank`, `exists_density_ge`);
* small cases (`sum_rank_ge_one`, `sum_rank_ge_add_two`), matching the data `ℓ = 4, d = 1 : N = 6`;
* uniform rank bound (`sum_choose_le_mul_of_rank_le`) and its contrapositive, a degree lower bound
  (`degree_ge_of_rank_le`): if every `rank (E k) ≤ r` and `ℓ * r < ∑_{i ≤ m} C(ℓ, i)` then
  `ℓ ≤ 2 m + d`.
-/

open Polynomial Matrix Module

namespace SqrtLog

open PerPoint

/-! ### Elementary bounds on `N ℓ m` -/

/-- `N ℓ m ≥ ∑_{i ≤ m} C(ℓ, i)` (the term `j = m` of `N_eq`). -/
theorem N_ge_sum_choose (ℓ m : ℕ) : ∑ i ∈ Finset.range (m + 1), ℓ.choose i ≤ N ℓ m := by
  rw [N_eq]
  exact Finset.single_le_sum (f := fun j => ∑ i ∈ Finset.range (j + 1), ℓ.choose i)
    (fun _ _ => Nat.zero_le _) (Finset.self_mem_range_succ m)

/-- `N ℓ m ≥ C(ℓ, m)`. -/
theorem N_ge_choose (ℓ m : ℕ) : ℓ.choose m ≤ N ℓ m :=
  (Finset.single_le_sum (f := fun i => ℓ.choose i) (fun _ _ => Nat.zero_le _)
    (Finset.self_mem_range_succ m)).trans (N_ge_sum_choose ℓ m)

/-- `N ℓ m` is nondecreasing in `m`. -/
theorem N_mono (ℓ : ℕ) {m m' : ℕ} (h : m ≤ m') : N ℓ m ≤ N ℓ m' := by
  rw [N_eq, N_eq]
  exact Finset.sum_le_sum_of_subset_of_nonneg
    (Finset.range_subset_range.2 (by omega)) (fun _ _ _ => Nat.zero_le _)

/-- `N ℓ (m + 1) = N ℓ m + ∑_{i ≤ m + 1} C(ℓ, i)`. -/
theorem N_succ (ℓ m : ℕ) :
    N ℓ (m + 1) = N ℓ m + ∑ i ∈ Finset.range (m + 2), ℓ.choose i := by
  rw [N_eq, N_eq, Finset.sum_range_succ]

@[simp] theorem N_zero (ℓ : ℕ) : N ℓ 0 = 1 := by
  simp [N_eq]

/-- `N ℓ 1 = 1 + (1 + ℓ) = ℓ + 2`. -/
@[simp] theorem N_one (ℓ : ℕ) : N ℓ 1 = ℓ + 2 := by
  rw [N_succ, N_zero]
  simp [Finset.sum_range_succ]
  omega

/-- `N ℓ 2 = (ℓ + 2) + (1 + ℓ + C(ℓ, 2))`. -/
theorem N_two (ℓ : ℕ) : N ℓ 2 = 2 * ℓ + 3 + ℓ.choose 2 := by
  rw [N_succ, N_one]
  simp [Finset.sum_range_succ]
  omega

theorem N_pos (ℓ m : ℕ) : 0 < N ℓ m :=
  lt_of_lt_of_le (by simp) (N_mono ℓ (Nat.zero_le m))

/-- `N ℓ m ≥ m + 1` (each of the `m + 1` inner sums contains `C(ℓ, 0) = 1`). -/
theorem N_ge_succ (ℓ m : ℕ) : m + 1 ≤ N ℓ m := by
  induction m with
  | zero => simp
  | succ m ih =>
    rw [N_succ]
    have : 1 ≤ ∑ i ∈ Finset.range (m + 2), ℓ.choose i :=
      (by simp : 1 ≤ ℓ.choose 0).trans (Finset.single_le_sum
        (f := fun i => ℓ.choose i) (fun _ _ => Nat.zero_le _) (by simp))
    omega

/-- Small-case check: `N 4 1 = 6` (the value found by computer for `ℓ = 4`, `d = 1`). -/
example : N 4 1 = 6 := by simp

/-! ### Consequences of the main theorem -/

section Main

variable (ℓ d : ℕ) (J : Matrix (Fin ℓ → Bool) (Fin ℓ → Bool) ℚ[X])
  (τ : Fin ℓ → ℚ) (K E : Fin ℓ → Matrix (Fin ℓ → Bool) (Fin ℓ → Bool) ℚ)

/-- **Pigeonhole form.** Under the hypotheses of `sum_rank_ge_N` and `0 < ℓ`, some direction `k`
has `N ℓ m ≤ ℓ * rank (E k)`. -/
theorem exists_rank_ge (hℓ : 0 < ℓ)
    (hd : ∀ x y, (J x y).natDegree ≤ d) (h0 : J.map (eval 0) = Xtensor ℓ)
    (hτ : Function.Injective τ) (hτ0 : ∀ k, τ k ≠ 0)
    (hdec : ∀ k, J.map (eval (τ k)) = K k + E k)
    (hnice : ∀ k x y, x k ≠ y k → K k x y = 0)
    (m : ℕ) (hm : 2 * m + d + 1 ≤ ℓ) :
    ∃ k, N ℓ m ≤ ℓ * (E k).rank := by
  have hmain := sum_rank_ge_N ℓ d J hd h0 τ hτ hτ0 K E hdec hnice m hm
  have : Nonempty (Fin ℓ) := ⟨⟨0, hℓ⟩⟩
  obtain ⟨k, -, hk⟩ := Finset.exists_max_image Finset.univ (fun k => (E k).rank)
    Finset.univ_nonempty
  refine ⟨k, hmain.trans ?_⟩
  calc ∑ k', (E k').rank ≤ ∑ _k' : Fin ℓ, (E k).rank :=
        Finset.sum_le_sum fun k' _ => hk k' (Finset.mem_univ _)
    _ = ℓ * (E k).rank := by simp

/-- **Uniform rank bound.** If every `rank (E k) ≤ r` then `N ℓ m ≤ ℓ * r`. -/
theorem N_le_mul_of_rank_le
    (hd : ∀ x y, (J x y).natDegree ≤ d) (h0 : J.map (eval 0) = Xtensor ℓ)
    (hτ : Function.Injective τ) (hτ0 : ∀ k, τ k ≠ 0)
    (hdec : ∀ k, J.map (eval (τ k)) = K k + E k)
    (hnice : ∀ k x y, x k ≠ y k → K k x y = 0)
    (m : ℕ) (hm : 2 * m + d + 1 ≤ ℓ) (r : ℕ) (hr : ∀ k, (E k).rank ≤ r) :
    N ℓ m ≤ ℓ * r := by
  have hmain := sum_rank_ge_N ℓ d J hd h0 τ hτ hτ0 K E hdec hnice m hm
  calc N ℓ m ≤ ∑ k, (E k).rank := hmain
    _ ≤ ∑ _k : Fin ℓ, r := Finset.sum_le_sum fun k _ => hr k
    _ = ℓ * r := by simp

/-- **Quantitative statement toward the `√` threshold.** If every `rank (E k) ≤ r` and
`2 m + d + 1 ≤ ℓ`, then `∑_{i ≤ m} C(ℓ, i) ≤ ℓ * r`. -/
theorem sum_choose_le_mul_of_rank_le
    (hd : ∀ x y, (J x y).natDegree ≤ d) (h0 : J.map (eval 0) = Xtensor ℓ)
    (hτ : Function.Injective τ) (hτ0 : ∀ k, τ k ≠ 0)
    (hdec : ∀ k, J.map (eval (τ k)) = K k + E k)
    (hnice : ∀ k x y, x k ≠ y k → K k x y = 0)
    (m : ℕ) (hm : 2 * m + d + 1 ≤ ℓ) (r : ℕ) (hr : ∀ k, (E k).rank ≤ r) :
    ∑ i ∈ Finset.range (m + 1), ℓ.choose i ≤ ℓ * r :=
  (N_ge_sum_choose ℓ m).trans
    (N_le_mul_of_rank_le ℓ d J τ K E hd h0 hτ hτ0 hdec hnice m hm r hr)

/-- Same with a single binomial coefficient: `C(ℓ, m) ≤ ℓ * r`. -/
theorem choose_le_mul_of_rank_le
    (hd : ∀ x y, (J x y).natDegree ≤ d) (h0 : J.map (eval 0) = Xtensor ℓ)
    (hτ : Function.Injective τ) (hτ0 : ∀ k, τ k ≠ 0)
    (hdec : ∀ k, J.map (eval (τ k)) = K k + E k)
    (hnice : ∀ k x y, x k ≠ y k → K k x y = 0)
    (m : ℕ) (hm : 2 * m + d + 1 ≤ ℓ) (r : ℕ) (hr : ∀ k, (E k).rank ≤ r) :
    ℓ.choose m ≤ ℓ * r :=
  (N_ge_choose ℓ m).trans
    (N_le_mul_of_rank_le ℓ d J τ K E hd h0 hτ hτ0 hdec hnice m hm r hr)

/-- **Degree lower bound** (contrapositive form). If every `rank (E k) ≤ r` and, for some `m`,
`ℓ * r < ∑_{i ≤ m} C(ℓ, i)`, then the degree satisfies `ℓ ≤ 2 m + d`, i.e. `d ≥ ℓ - 2 m`. -/
theorem degree_ge_of_rank_le
    (hd : ∀ x y, (J x y).natDegree ≤ d) (h0 : J.map (eval 0) = Xtensor ℓ)
    (hτ : Function.Injective τ) (hτ0 : ∀ k, τ k ≠ 0)
    (hdec : ∀ k, J.map (eval (τ k)) = K k + E k)
    (hnice : ∀ k x y, x k ≠ y k → K k x y = 0)
    (m r : ℕ) (hr : ∀ k, (E k).rank ≤ r)
    (hsmall : ℓ * r < ∑ i ∈ Finset.range (m + 1), ℓ.choose i) :
    ℓ ≤ 2 * m + d := by
  by_contra hcon
  have := sum_choose_le_mul_of_rank_le ℓ d J τ K E hd h0 hτ hτ0 hdec hnice m (by omega) r hr
  omega

/-- **Exact recovery is impossible below degree `ℓ`**: if `d + 1 ≤ ℓ`, the errors cannot all
vanish: `∑ k, rank (E k) ≥ 1`. -/
theorem sum_rank_ge_one
    (hd : ∀ x y, (J x y).natDegree ≤ d) (h0 : J.map (eval 0) = Xtensor ℓ)
    (hτ : Function.Injective τ) (hτ0 : ∀ k, τ k ≠ 0)
    (hdec : ∀ k, J.map (eval (τ k)) = K k + E k)
    (hnice : ∀ k x y, x k ≠ y k → K k x y = 0)
    (hℓ : d + 1 ≤ ℓ) : 1 ≤ ∑ k, (E k).rank := by
  simpa using sum_rank_ge_N ℓ d J hd h0 τ hτ hτ0 K E hdec hnice 0 (by omega)

/-- If `d + 1 ≤ ℓ`, some error matrix `E k` is nonzero. -/
theorem exists_error_ne_zero
    (hd : ∀ x y, (J x y).natDegree ≤ d) (h0 : J.map (eval 0) = Xtensor ℓ)
    (hτ : Function.Injective τ) (hτ0 : ∀ k, τ k ≠ 0)
    (hdec : ∀ k, J.map (eval (τ k)) = K k + E k)
    (hnice : ∀ k x y, x k ≠ y k → K k x y = 0)
    (hℓ : d + 1 ≤ ℓ) : ∃ k, E k ≠ 0 := by
  by_contra hcon
  push Not at hcon
  have := sum_rank_ge_one ℓ d J τ K E hd h0 hτ hτ0 hdec hnice hℓ
  simp [hcon] at this

/-- If `d + 3 ≤ ℓ` then `∑ k, rank (E k) ≥ ℓ + 2` (`N ℓ 1 = 1 + (1 + ℓ)`). For `ℓ = 4`, `d = 1`
this is the bound `6` of the small-case data. -/
theorem sum_rank_ge_add_two
    (hd : ∀ x y, (J x y).natDegree ≤ d) (h0 : J.map (eval 0) = Xtensor ℓ)
    (hτ : Function.Injective τ) (hτ0 : ∀ k, τ k ≠ 0)
    (hdec : ∀ k, J.map (eval (τ k)) = K k + E k)
    (hnice : ∀ k x y, x k ≠ y k → K k x y = 0)
    (hℓ : d + 3 ≤ ℓ) : ℓ + 2 ≤ ∑ k, (E k).rank := by
  simpa using sum_rank_ge_N ℓ d J hd h0 τ hτ hτ0 K E hdec hnice 1 (by omega)

/-- If `d + 3 ≤ ℓ`, some direction has `rank (E k) ≥ 2`. -/
theorem exists_rank_ge_two
    (hd : ∀ x y, (J x y).natDegree ≤ d) (h0 : J.map (eval 0) = Xtensor ℓ)
    (hτ : Function.Injective τ) (hτ0 : ∀ k, τ k ≠ 0)
    (hdec : ∀ k, J.map (eval (τ k)) = K k + E k)
    (hnice : ∀ k x y, x k ≠ y k → K k x y = 0)
    (hℓ : d + 3 ≤ ℓ) : ∃ k, 2 ≤ (E k).rank := by
  by_contra hcon
  push Not at hcon
  have h1 := N_le_mul_of_rank_le ℓ d J τ K E hd h0 hτ hτ0 hdec hnice 1 (by omega) 1
    (fun k => by have := hcon k; omega)
  rw [N_one] at h1
  omega

/-- If `d + 5 ≤ ℓ` then `∑ k, rank (E k) ≥ 2 ℓ + 3 + C(ℓ, 2)`. -/
theorem sum_rank_ge_quadratic
    (hd : ∀ x y, (J x y).natDegree ≤ d) (h0 : J.map (eval 0) = Xtensor ℓ)
    (hτ : Function.Injective τ) (hτ0 : ∀ k, τ k ≠ 0)
    (hdec : ∀ k, J.map (eval (τ k)) = K k + E k)
    (hnice : ∀ k x y, x k ≠ y k → K k x y = 0)
    (hℓ : d + 5 ≤ ℓ) : 2 * ℓ + 3 + ℓ.choose 2 ≤ ∑ k, (E k).rank := by
  have := sum_rank_ge_N ℓ d J hd h0 τ hτ hτ0 K E hdec hnice 2 (by omega)
  rwa [N_two] at this

/-- **Best `m` for a given degree**: with `m = (ℓ - d - 1) / 2`,
`∑ k, rank (E k) ≥ N ℓ ((ℓ - d - 1) / 2)` as soon as `d + 1 ≤ ℓ`. -/
theorem sum_rank_ge_N_opt
    (hd : ∀ x y, (J x y).natDegree ≤ d) (h0 : J.map (eval 0) = Xtensor ℓ)
    (hτ : Function.Injective τ) (hτ0 : ∀ k, τ k ≠ 0)
    (hdec : ∀ k, J.map (eval (τ k)) = K k + E k)
    (hnice : ∀ k x y, x k ≠ y k → K k x y = 0)
    (hℓ : d + 1 ≤ ℓ) : N ℓ ((ℓ - d - 1) / 2) ≤ ∑ k, (E k).rank :=
  sum_rank_ge_N ℓ d J hd h0 τ hτ hτ0 K E hdec hnice _ (by omega)

/-! ### Rational (density) forms -/

/-- **Density form (existential).** Some direction has error density
`rank (E k) / 2^ℓ ≥ N ℓ m / (ℓ * 2^ℓ)`. -/
theorem exists_density_ge (hℓ : 0 < ℓ)
    (hd : ∀ x y, (J x y).natDegree ≤ d) (h0 : J.map (eval 0) = Xtensor ℓ)
    (hτ : Function.Injective τ) (hτ0 : ∀ k, τ k ≠ 0)
    (hdec : ∀ k, J.map (eval (τ k)) = K k + E k)
    (hnice : ∀ k x y, x k ≠ y k → K k x y = 0)
    (m : ℕ) (hm : 2 * m + d + 1 ≤ ℓ) :
    ∃ k, (N ℓ m : ℚ) / (ℓ * 2 ^ ℓ) ≤ ((E k).rank : ℚ) / 2 ^ ℓ := by
  obtain ⟨k, hk⟩ := exists_rank_ge ℓ d J τ K E hℓ hd h0 hτ hτ0 hdec hnice m hm
  refine ⟨k, ?_⟩
  have hℓ' : (0 : ℚ) < ℓ := by exact_mod_cast hℓ
  have h2 : (0 : ℚ) < 2 ^ ℓ := by positivity
  have hk' : (N ℓ m : ℚ) ≤ ℓ * (E k).rank := by exact_mod_cast hk
  rw [div_le_div_iff₀ (by positivity) h2]
  nlinarith

/-- **Max form.** `N ℓ m / (ℓ * 2^ℓ) ≤ max_k rank (E k) / 2^ℓ`. -/
theorem avg_le_max_rank (hℓ : 0 < ℓ)
    (hd : ∀ x y, (J x y).natDegree ≤ d) (h0 : J.map (eval 0) = Xtensor ℓ)
    (hτ : Function.Injective τ) (hτ0 : ∀ k, τ k ≠ 0)
    (hdec : ∀ k, J.map (eval (τ k)) = K k + E k)
    (hnice : ∀ k x y, x k ≠ y k → K k x y = 0)
    (m : ℕ) (hm : 2 * m + d + 1 ≤ ℓ) :
    (N ℓ m : ℚ) / (ℓ * 2 ^ ℓ) ≤
      Finset.univ.sup' (Finset.univ_nonempty_iff.2 ⟨(⟨0, hℓ⟩ : Fin ℓ)⟩)
        (fun k => ((E k).rank : ℚ) / 2 ^ ℓ) := by
  obtain ⟨k, hk⟩ := exists_density_ge ℓ d J τ K E hℓ hd h0 hτ hτ0 hdec hnice m hm
  exact hk.trans (Finset.le_sup' (fun k => ((E k).rank : ℚ) / 2 ^ ℓ) (Finset.mem_univ k))

/-- **Average form.** The average error density over the `ℓ` directions is at least
`N ℓ m / (ℓ * 2^ℓ)`. -/
theorem avg_density_ge (hℓ : 0 < ℓ)
    (hd : ∀ x y, (J x y).natDegree ≤ d) (h0 : J.map (eval 0) = Xtensor ℓ)
    (hτ : Function.Injective τ) (hτ0 : ∀ k, τ k ≠ 0)
    (hdec : ∀ k, J.map (eval (τ k)) = K k + E k)
    (hnice : ∀ k x y, x k ≠ y k → K k x y = 0)
    (m : ℕ) (hm : 2 * m + d + 1 ≤ ℓ) :
    (N ℓ m : ℚ) / (ℓ * 2 ^ ℓ) ≤ (∑ k, ((E k).rank : ℚ) / 2 ^ ℓ) / ℓ := by
  have hmain := sum_rank_ge_N ℓ d J hd h0 τ hτ hτ0 K E hdec hnice m hm
  have hℓ' : (0 : ℚ) < ℓ := by exact_mod_cast hℓ
  have h2 : (0 : ℚ) < 2 ^ ℓ := by positivity
  have hk' : (N ℓ m : ℚ) ≤ ∑ k, ((E k).rank : ℚ) := by exact_mod_cast hmain
  rw [← Finset.sum_div, div_div, mul_comm (2 ^ ℓ : ℚ) ℓ]
  exact div_le_div_of_nonneg_right hk' (by positivity)

/-- Density form of the `C(ℓ, m)` bound: some direction has
`rank (E k) / 2^ℓ ≥ C(ℓ, m) / (ℓ * 2^ℓ)`. -/
theorem exists_density_ge_choose (hℓ : 0 < ℓ)
    (hd : ∀ x y, (J x y).natDegree ≤ d) (h0 : J.map (eval 0) = Xtensor ℓ)
    (hτ : Function.Injective τ) (hτ0 : ∀ k, τ k ≠ 0)
    (hdec : ∀ k, J.map (eval (τ k)) = K k + E k)
    (hnice : ∀ k x y, x k ≠ y k → K k x y = 0)
    (m : ℕ) (hm : 2 * m + d + 1 ≤ ℓ) :
    ∃ k, (ℓ.choose m : ℚ) / (ℓ * 2 ^ ℓ) ≤ ((E k).rank : ℚ) / 2 ^ ℓ := by
  obtain ⟨k, hk⟩ := exists_density_ge ℓ d J τ K E hℓ hd h0 hτ hτ0 hdec hnice m hm
  refine ⟨k, le_trans ?_ hk⟩
  have hℓ' : (0 : ℚ) < ℓ := by exact_mod_cast hℓ
  exact div_le_div_of_nonneg_right (by exact_mod_cast N_ge_choose ℓ m) (by positivity)

end Main

end SqrtLog

#print axioms SqrtLog.N_ge_sum_choose
#print axioms SqrtLog.N_ge_choose
#print axioms SqrtLog.N_one
#print axioms SqrtLog.exists_rank_ge
#print axioms SqrtLog.avg_le_max_rank
#print axioms SqrtLog.avg_density_ge
#print axioms SqrtLog.sum_rank_ge_one
#print axioms SqrtLog.sum_rank_ge_add_two
#print axioms SqrtLog.sum_choose_le_mul_of_rank_le
#print axioms SqrtLog.degree_ge_of_rank_le
