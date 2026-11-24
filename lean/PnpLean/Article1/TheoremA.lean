import PnpLean.Article1.TheoremB

/-!
# Theorem A (variable-uniform deformations: no absorption), unconditional
-/

open Polynomial Matrix Module

set_option linter.unusedSectionVars false
set_option linter.style.numericProjections false

namespace ClauseScheme

variable {n : ℕ} {κ : Type*} [Fintype κ] [DecidableEq κ] {I : κ → Finset (Fin n)}
  {G : Type*} [Fintype G] [DecidableEq G]

namespace Scheme

variable (S : Scheme I G)

/-- Variable-uniform deformation: `P_{e,v} = P_v`. -/
def Uniform : Prop := ∀ o o' : Occ I, o.var = o'.var → S.P o = S.P o'

/-- Home in `Option G`: `none` for split variables (constrained in every block). -/
noncomputable def homeA (t : ℚ) (v : Fin n) : Option G := by
  classical
  exact if v ∈ S.split t then none else some (S.homeNS t v)

theorem homeA_eq_some (t : ℚ) (v : Fin n) (γ : G) :
    S.homeA t v = some γ ↔ v ∈ S.J t γ := by
  classical
  rw [S.mem_J]
  unfold homeA
  split_ifs with h <;> simp [h]

/-- Measurement spaces extended by `⊥` at `none`. -/
def RspA (t : ℚ) : Option G → Submodule ℚ (Fin n → ℚ)
  | none => ⊥
  | some γ => S.Rsp t γ

theorem rowA_mem (hU : S.Uniform) (ω : WitnessMap I) (t : ℚ) (ht : t ∈ S.pts)
    (γ' : Option G) (v : Fin n) (hv : S.homeA t v ≠ γ') :
    (fun j => if (some (S.grp t j) : Option G) = γ' then
      WitnessDet.evalMat (S.Φ ω) t v j else 0) ∈ S.RspA t γ' := by
  classical
  cases γ' with
  | none =>
    convert (S.RspA t none).zero_mem using 1
    funext j; simp
  | some γ =>
    simp only [Option.some.injEq]
    show _ ∈ S.Rsp t γ
    simp only [WitnessDet.evalMat, Matrix.map_apply, Φ, row]
    rcases hω : ω.1 v with _ | o
    · have hno : ∀ o : Occ I, o.var ≠ v := ω.2.2 v hω
      have hns : v ∉ S.split t := by
        intro hs
        obtain ⟨-, o, -, ho, -⟩ := Finset.mem_filter.1 hs
        exact hno o ho
      have hh : S.homeNS t v = S.grp t v :=
        S.homeNS_noocc t v (fun ⟨o, ho⟩ => hno o ho)
      have hv' : S.grp t v ≠ γ := by
        intro h; apply hv; unfold homeA; rw [if_neg hns, hh, h]
      simp only
      convert (S.Rsp t γ).zero_mem using 1
      funext j
      by_cases hj : j = v
      · subst hj; simp [hv']
      · simp [hj]
    · have hov : o.var = v := ω.2.1 v o hω
      simp only
      -- find an occurrence `o'` of `v` whose clause is not assigned to `γ`
      obtain ⟨o', ho', hγ⟩ : ∃ o' : Occ I, o'.var = v ∧ S.α t o'.clause ≠ γ := by
        by_cases hs : v ∈ S.split t
        · obtain ⟨-, o1, o2, h1, h2, hne⟩ := Finset.mem_filter.1 hs
          by_cases h1γ : S.α t o1.clause = γ
          · exact ⟨o2, h2, fun h => hne (h1γ.trans h.symm)⟩
          · exact ⟨o1, h1, h1γ⟩
        · refine ⟨o, hov, ?_⟩
          rw [S.homeNS_occ t v hs o hov]
          intro h; apply hv; unfold homeA; rw [if_neg hs, h]
      rw [hU o o' (hov.trans ho'.symm)]
      exact S.cond3 t ht o' γ hγ

/-- Rank bound at a point for a variable-uniform scheme:
`Σ_γ (a_γ - |J_{t,γ}|)_+ - ρ_t ≤ dim ker Φ(t)`. -/
theorem defA_sub_le_nullity (hU : S.Uniform) (ω : WitnessMap I) (t : ℚ) (ht : t ∈ S.pts) :
    (∑ γ, (S.a t γ - (S.J t γ).card)) - S.ρ t
      ≤ WitnessDet.nullity (WitnessDet.evalMat (S.Φ ω) t) := by
  classical
  have key := WitnessDet.deficiency_le_nullity_add_general (WitnessDet.evalMat (S.Φ ω) t)
    (fun j => (some (S.grp t j) : Option G)) (S.homeA t) (S.RspA t)
    (fun γ' v hv => S.rowA_mem hU ω t ht γ' v hv)
  rw [Fintype.sum_option, Fintype.sum_option] at key
  have h1 : (Finset.univ.filter (fun j => (some (S.grp t j) : Option G) = none)).card = 0 := by
    simp
  have h2 : ∀ γ, (Finset.univ.filter (fun j => (some (S.grp t j) : Option G) = some γ)).card
      = S.a t γ := by
    intro γ; simp [a]
  have h3 : ∀ γ, (Finset.univ.filter (fun v => S.homeA t v = some γ)).card = (S.J t γ).card := by
    intro γ
    congr 1
    ext v; simp [S.homeA_eq_some]
  have h4 : finrank ℚ (S.RspA t none) = 0 := by rw [RspA]; exact finrank_bot _ _
  simp only [h1, h2, h3, h4, zero_tsub, zero_add] at key
  have h5 : ∑ γ, finrank ℚ (S.RspA t (some γ)) = S.ρ t := rfl
  rw [h5] at key
  omega

/-- Determinant + Fact R for variable-uniform schemes. -/
theorem exists_pointA (hU : S.Uniform) (hM : n ≤ M I) :
    ∃ t ∈ S.pts, ∑ γ, (S.a t γ - (S.J t γ).card) < S.ρ t + S.q t + 1 := by
  obtain ⟨ω⟩ := (inferInstance : Nonempty (WitnessMap I))
  by_contra hcon
  push Not at hcon
  have h1 : ∑ t ∈ S.pts, (S.q t + 1) ≤
      ∑ t ∈ S.pts, ((∑ γ, (S.a t γ - (S.J t γ).card)) - S.ρ t) :=
    Finset.sum_le_sum fun t ht => by have := hcon t ht; omega
  have h2 : ∑ t ∈ S.pts, ((∑ γ, (S.a t γ - (S.J t γ).card)) - S.ρ t)
      ≤ ∑ t ∈ S.pts, WitnessDet.nullity (WitnessDet.evalMat (S.Φ ω) t) :=
    Finset.sum_le_sum fun t ht => S.defA_sub_le_nullity hU ω t ht
  have h3 := WitnessDet.sum_nullity_le _ _ (S.Φ_natDegree ω) (S.Φ_eval_zero ω) S.pts
  have h4 := S.factR
  have h5 : n * S.R ≤ M I * S.R := Nat.mul_le_mul_right _ hM
  omega

theorem sum_J (t : ℚ) : ∑ γ, (S.J t γ).card + (S.split t).card = n := by
  classical
  have := Finset.card_eq_sum_card_fiberwise (f := S.homeNS t)
    (s := Finset.univ.filter (fun v => v ∉ S.split t)) (t := Finset.univ)
    (fun _ _ => Finset.mem_coe.2 (Finset.mem_univ _))
  have hJ : ∀ γ, (Finset.univ.filter (fun v => v ∉ S.split t)).filter
      (fun v => S.homeNS t v = γ) = S.J t γ := by
    intro γ; ext v; simp [S.mem_J]
  simp only [hJ] at this
  rw [← this]
  have hc : Finset.univ.filter (fun v => v ∉ S.split t) = (S.split t)ᶜ := by ext; simp
  rw [hc, Finset.card_compl, Fintype.card_fin]
  have := Finset.card_le_univ (S.split t)
  rw [Fintype.card_fin] at this
  omega

/-- **Theorem A** (variable-uniform deformations: no absorption). If `M ≥ n`, the deformation
is variable-uniform and `sepLB` is a lower bound for `sep(F)`, then some `t ∈ X` satisfies
`ρ_t + q_t + 1 > min (sep(F), 2n/3 - s_t)`. -/
theorem theoremA (hU : S.Uniform) (hM : n ≤ M I) (sepLB : ℕ)
    (hsep : ∀ Y : Finset (Fin n), IsSep I Y → sepLB ≤ Y.card) :
    ∃ t ∈ S.pts, min (sepLB : ℝ) (2 * n / 3 - S.s t) < S.ρ t + S.q t + 1 := by
  classical
  obtain ⟨t, ht, hlt⟩ := S.exists_pointA hU hM
  refine ⟨t, ht, ?_⟩
  set D := ∑ γ, (S.a t γ - (S.J t γ).card) with hD
  have hDlt : (D : ℝ) < S.ρ t + S.q t + 1 := by exact_mod_cast hlt
  have hpt : ∀ γ, S.a t γ ≤ (S.a t γ - (S.J t γ).card) + (S.J t γ).card := fun γ => by omega
  have hsa := S.sum_a t
  have hsJ := S.sum_J t
  by_cases hcase : ∀ γ, 3 * (S.J t γ).card ≤ 2 * n
  · -- Case 1: the split variables form a 2/3-separator
    have hsepY : IsSep I (S.split t) := by
      intro v hv
      refine le_trans ?_ (hcase (S.homeNS t v))
      refine Nat.mul_le_mul_left _ ?_
      rw [← Set.ncard_coe_finset]
      refine Set.ncard_le_ncard (fun u hu => ?_) (Finset.finite_toSet _)
      obtain ⟨h1, h2⟩ := S.reach_homeNS t v u hv hu
      exact (S.mem_J t _ u).2 ⟨h1, h2⟩
    have h1 := hsep _ hsepY
    have h2 : ∑ γ, S.a t γ ≤ D + ∑ γ, (S.J t γ).card := by
      rw [hD, ← Finset.sum_add_distrib]; exact Finset.sum_le_sum fun γ _ => hpt γ
    have : sepLB ≤ D := by omega
    refine lt_of_le_of_lt (min_le_left _ _) (lt_of_le_of_lt ?_ hDlt)
    exact_mod_cast this
  · -- Case 2: one non-split class is larger than 2n/3
    push Not at hcase
    obtain ⟨γ₀, hγ₀⟩ := hcase
    have hs : S.a t γ₀ ≤ S.s t := Finset.le_sup (f := S.a t) (Finset.mem_univ γ₀)
    have e1 := Finset.add_sum_erase Finset.univ (S.a t) (Finset.mem_univ γ₀)
    have e2 := Finset.add_sum_erase Finset.univ (fun γ => (S.J t γ).card) (Finset.mem_univ γ₀)
    have e3 := Finset.add_sum_erase Finset.univ (fun γ => S.a t γ - (S.J t γ).card)
      (Finset.mem_univ γ₀)
    have e4 : ∑ γ ∈ Finset.univ.erase γ₀, S.a t γ ≤
        ∑ γ ∈ Finset.univ.erase γ₀, (S.a t γ - (S.J t γ).card) +
        ∑ γ ∈ Finset.univ.erase γ₀, (S.J t γ).card := by
      rw [← Finset.sum_add_distrib]; exact Finset.sum_le_sum fun γ _ => hpt γ
    have key : (S.J t γ₀).card ≤ S.s t + D := by omega
    have key' : (2 * n : ℝ) < 3 * (S.s t + D) := by
      have : 2 * n < 3 * (S.s t + D) := by omega
      exact_mod_cast this
    refine lt_of_le_of_lt (min_le_right _ _) (lt_of_le_of_lt ?_ hDlt)
    linarith

end Scheme

end ClauseScheme

#print axioms ClauseScheme.Scheme.defA_sub_le_nullity
#print axioms ClauseScheme.Scheme.exists_pointA
#print axioms ClauseScheme.Scheme.theoremA
