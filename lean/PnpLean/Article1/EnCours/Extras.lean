import PnpLean.Article1.Lemma2
import PnpLean.Article1.TheoremA

/-!
# Article 1, deterministic leftovers

## (a) The hypothesis `n ≤ M` is not needed

Paper (`resultats/article/main.tex`, proof of the witness-determinant lemma, Step 1, and the
remark after it): the rows of `Φ^ω` indexed by variables that never occur are constant, so
`deg det Φ^ω ≤ n_occ · R ≤ min(n, M) · R`, where `n_occ` is the number of variables that occur.

* `WitnessDet.natDegree_det_le_card`, `WitnessDet.sum_nullity_le_card` : the sharper degree bound.
* `ClauseScheme.nOcc`, `nOcc_le_M`, `nOcc_le_n`.
* `Scheme.lemma1_occ`, `Scheme.avg_lemma1_occ` : Lemma 1 with `n_occ · R`.
* `Scheme.exists_point'`, `Scheme.theoremB'`, `Scheme.theoremB_bounded'`,
  `Scheme.theoremB_unconditional'`, `Scheme.theoremB_bounded_unconditional'`,
  `Scheme.exists_pointA'`, `Scheme.theoremA'` : the results of `TheoremB.lean`, `Lemma2.lean`
  and `TheoremA.lean` **without** the hypothesis `hM : n ≤ M I`.
-/

open Polynomial Matrix Module

set_option linter.unusedSectionVars false
set_option linter.style.numericProjections false

namespace WitnessDet

variable {n : ℕ}

/-- **Sharper degree bound.** If the rows in `s` have entries of degree `≤ R` and the other rows
are constant, then `deg det Φ ≤ |s| · R`. -/
theorem natDegree_det_le_card (Φ : Matrix (Fin n) (Fin n) ℚ[X]) (R : ℕ) (s : Finset (Fin n))
    (hs : ∀ i ∈ s, ∀ j, (Φ i j).natDegree ≤ R) (hns : ∀ i ∉ s, ∀ j, (Φ i j).natDegree = 0) :
    Φ.det.natDegree ≤ s.card * R := by
  classical
  rw [det_apply']
  apply natDegree_sum_le_of_forall_le
  intro σ _
  refine (natDegree_mul_le).trans ?_
  rw [natDegree_intCast, zero_add]
  refine (natDegree_prod_le _ _).trans ?_
  calc ∑ i, (Φ (σ i) i).natDegree ≤ ∑ i, (if σ i ∈ s then R else 0) :=
        Finset.sum_le_sum fun i _ => by
          split_ifs with h
          · exact hs _ h _
          · exact (hns _ h _).le
    _ = ∑ j, (if j ∈ s then R else 0) := Equiv.sum_comp σ (fun j => if j ∈ s then R else 0)
    _ = s.card * R := by simp

/-- Sum of nullities with the sharper degree bound. -/
theorem sum_nullity_le_card (Φ : Matrix (Fin n) (Fin n) ℚ[X]) (R : ℕ) (s : Finset (Fin n))
    (hs : ∀ i ∈ s, ∀ j, (Φ i j).natDegree ≤ R) (hns : ∀ i ∉ s, ∀ j, (Φ i j).natDegree = 0)
    (h0 : evalMat Φ 0 = 1) (T : Finset ℚ) :
    ∑ t ∈ T, nullity (evalMat Φ t) ≤ s.card * R := by
  have hdvd : ∏ t ∈ T, (X - C t) ^ nullity (evalMat Φ t) ∣ Φ.det := by
    apply Finset.prod_dvd_of_coprime
    · intro a _ b _ hab
      exact ((pairwise_coprime_X_sub_C (s := id) Function.injective_id) hab).pow
    · intro t _
      exact pow_dvd_det_of_nullity Φ t
  have hdeg := natDegree_le_of_dvd hdvd (det_ne_zero_of_eval_zero Φ h0)
  rw [natDegree_prod_of_monic _ _ (fun t _ => (monic_X_sub_C t).pow _)] at hdeg
  simp only [natDegree_pow, natDegree_X_sub_C, mul_one] at hdeg
  exact hdeg.trans (natDegree_det_le_card Φ R s hs hns)

end WitnessDet

namespace ClauseScheme

variable {n : ℕ} {κ : Type*} [Fintype κ] [DecidableEq κ] {I : κ → Finset (Fin n)}
  {G : Type*} [Fintype G] [DecidableEq G]

/-- The variables that occur in some clause (`d_v ≥ 1`). -/
def occVars (I : κ → Finset (Fin n)) : Finset (Fin n) :=
  Finset.univ.filter (fun v => 0 < dOcc I v)

/-- `n_occ`: the number of variables that occur. -/
def nOcc (I : κ → Finset (Fin n)) : ℕ := (occVars I).card

theorem mem_occVars (v : Fin n) : v ∈ occVars I ↔ ∃ o : Occ I, o.var = v := by
  simp [occVars, dOcc, Finset.card_pos, Finset.Nonempty]

theorem nOcc_le_n : nOcc I ≤ n := by
  have := Finset.card_le_univ (occVars I)
  rwa [Fintype.card_fin] at this

theorem nOcc_le_M : nOcc I ≤ M I := by
  classical
  have h : occVars I = Finset.univ.image (fun o : Occ I => o.var) := by
    ext v
    simp [mem_occVars]
  unfold nOcc M
  rw [h]
  exact Finset.card_image_le.trans (by simp)

namespace Scheme

variable (S : Scheme I G)

theorem witness_none (ω : WitnessMap I) (v : Fin n) (hv : v ∉ occVars I) : ω.1 v = none := by
  rcases hω : ω.1 v with _ | o
  · rfl
  · exact absurd ((mem_occVars v).2 ⟨o, ω.2.1 v o hω⟩) hv

/-- Rows of never-occurring variables are constant. -/
theorem Φ_natDegree_eq_zero (ω : WitnessMap I) (v : Fin n) (hv : v ∉ occVars I) (j : Fin n) :
    (S.Φ ω v j).natDegree = 0 := by
  unfold Φ row
  rw [Scheme.witness_none ω v hv]
  simp only
  split_ifs <;> simp

/-- **Lemma 1, sharp form** (paper, Step 1): `Σ_t dim ker Φ^ω(t) ≤ n_occ · R`. -/
theorem sum_nullity_le_occ (ω : WitnessMap I) :
    ∑ t ∈ S.pts, WitnessDet.nullity (WitnessDet.evalMat (S.Φ ω) t) ≤ nOcc I * S.R :=
  WitnessDet.sum_nullity_le_card _ _ (occVars I) (fun i _ j => S.Φ_natDegree ω i j)
    (fun i hi j => S.Φ_natDegree_eq_zero ω i hi j) (S.Φ_eval_zero ω) _

/-- **Lemma 1** with `n_occ · R`: `Σ_{t∈X} (def^ω_t - ρ_t)_+ ≤ n_occ · R`. -/
theorem lemma1_occ (ω : WitnessMap I) :
    ∑ t ∈ S.pts, (S.deficiency ω t - S.ρ t) ≤ nOcc I * S.R :=
  (S.lemma1 ω).1.trans (S.sum_nullity_le_occ ω)

/-- **Lemma 1** with `min(n, M) · R`. -/
theorem lemma1_min (ω : WitnessMap I) :
    ∑ t ∈ S.pts, (S.deficiency ω t - S.ρ t) ≤ min n (M I) * S.R :=
  (S.lemma1_occ ω).trans (Nat.mul_le_mul_right _ (le_min nOcc_le_n nOcc_le_M))

/-- **Lemma 1** with `M · R`, no hypothesis on the formula. -/
theorem lemma1_M (ω : WitnessMap I) :
    ∑ t ∈ S.pts, (S.deficiency ω t - S.ρ t) ≤ M I * S.R :=
  (S.lemma1_occ ω).trans (Nat.mul_le_mul_right _ nOcc_le_M)

/-- Averaged Lemma 1, sharp form: `Σ_t (𝔼 def_t - ρ_t)_+ ≤ n_occ · R`. -/
theorem avg_lemma1_occ : ∑ t ∈ S.pts, max (S.Edef t - S.ρ t) 0 ≤ nOcc I * S.R := by
  have h1 : ∀ t, max (S.Edef t - S.ρ t) 0 ≤
      avg (fun ω => ((S.deficiency ω t - S.ρ t : ℕ) : ℝ)) := by
    intro t
    apply max_le
    · have hc : avg (fun _ : WitnessMap I => (S.ρ t : ℝ)) = S.ρ t := avg_const _
      rw [Edef, ← hc, avg, avg, avg, ← sub_div, ← Finset.sum_sub_distrib]
      refine div_le_div_of_nonneg_right (Finset.sum_le_sum fun ω _ => ?_) card_pos.le
      have : (S.deficiency ω t : ℝ) - S.ρ t ≤ ((S.deficiency ω t - S.ρ t : ℕ) : ℝ) := by
        rcases le_total (S.ρ t) (S.deficiency ω t) with h | h
        · rw [Nat.cast_sub h]
        · rw [Nat.sub_eq_zero_of_le h]; simp; exact_mod_cast h
      exact this
    · have h0 : avg (fun _ : WitnessMap I => (0 : ℝ)) = 0 := avg_const 0
      rw [← h0]
      exact avg_mono fun ω => Nat.cast_nonneg _
  calc ∑ t ∈ S.pts, max (S.Edef t - S.ρ t) 0
      ≤ ∑ t ∈ S.pts, avg (fun ω => ((S.deficiency ω t - S.ρ t : ℕ) : ℝ)) :=
        Finset.sum_le_sum fun t _ => h1 t
    _ = avg (fun ω => ∑ t ∈ S.pts, ((S.deficiency ω t - S.ρ t : ℕ) : ℝ)) := (avg_sum _ _).symm
    _ ≤ avg (fun _ => ((nOcc I * S.R : ℕ) : ℝ)) := avg_mono fun ω => by
        rw [← Nat.cast_sum]; exact_mod_cast S.lemma1_occ ω
    _ = nOcc I * S.R := by rw [avg_const]; push_cast; ring

/-- Averaged Lemma 1 with `M · R`, no hypothesis on the formula. -/
theorem avg_lemma1_M : ∑ t ∈ S.pts, max (S.Edef t - S.ρ t) 0 ≤ M I * S.R := by
  refine S.avg_lemma1_occ.trans ?_
  exact_mod_cast Nat.mul_le_mul_right S.R (nOcc_le_M (I := I))

/-- Averaging + Fact R, **without** `n ≤ M`: some point has `𝔼 def_t < ρ_t + q_t + 1`. -/
theorem exists_point' : ∃ t ∈ S.pts, S.Edef t < S.ρ t + S.q t + 1 := by
  by_contra hcon
  push Not at hcon
  have h1 : ∑ t ∈ S.pts, ((S.q t + 1 : ℕ) : ℝ) ≤ ∑ t ∈ S.pts, max (S.Edef t - S.ρ t) 0 :=
    Finset.sum_le_sum fun t ht => by
      push_cast
      exact le_max_of_le_left (by linarith [hcon t ht])
  have h2 := S.factR
  have h3 := S.avg_lemma1_M
  rw [← Nat.cast_sum] at h1
  have h5 : ((M I * S.R + 1 : ℕ) : ℝ) ≤ ((∑ t ∈ S.pts, (S.q t + 1) : ℕ) : ℝ) := by
    exact_mod_cast h2
  push_cast at h5 h1
  linarith

/-- **Theorem B** (conditional on Lemma 2), **without** `n ≤ M`. -/
theorem theoremB' (c0 W : ℝ) (hc0 : 0 ≤ c0) (hL2 : S.Lemma2Hyp c0)
    (hW : ∀ Y : Finset (Fin n), IsSep I Y → (∀ v ∈ Y, 2 ≤ dOcc I v) →
      W ≤ ∑ v ∈ Y, wt (dOcc I v)) :
    ∃ t ∈ S.pts, min (c0 / (2 * Real.sqrt 2) * W / Real.sqrt (S.s t)) (2 * n / 3 - S.s t)
      < S.ρ t + S.q t + 1 := by
  obtain ⟨t, ht, hlt⟩ := S.exists_point'
  refine ⟨t, ht, lt_of_le_of_lt ?_ hlt⟩
  have hL := hL2 t ht
  rcases S.lemma3 t with h1 | h2
  · refine min_le_of_left_le ?_
    have hw := hW _ h1 (S.split_two_le t)
    refine le_trans ?_ ((le_max_right _ _).trans hL)
    rw [Real.sqrt_mul (by norm_num : (0 : ℝ) ≤ 2)]
    have hk : 0 ≤ c0 / (2 * (Real.sqrt 2 * Real.sqrt (S.s t))) := by positivity
    calc c0 / (2 * Real.sqrt 2) * W / Real.sqrt (S.s t)
        = c0 / (2 * (Real.sqrt 2 * Real.sqrt (S.s t))) * W := by
          rw [div_mul_eq_mul_div, div_div, div_mul_eq_mul_div, mul_assoc]
      _ ≤ _ := mul_le_mul_of_nonneg_left hw hk
  · refine min_le_of_right_le ?_
    exact le_trans h2.le ((le_max_left _ _).trans hL)

/-- **Theorem B, bounded occurrences** (conditional on Lemma 2), **without** `n ≤ M`. -/
theorem theoremB_bounded' (c0 : ℝ) (hc0 : 0 ≤ c0) (hL2 : S.Lemma2Hyp c0)
    (b : ℕ) (hb : 2 ≤ b) (hdb : ∀ v, dOcc I v ≤ b) (sepLB : ℕ)
    (hsep : ∀ Y : Finset (Fin n), IsSep I Y → sepLB ≤ Y.card) :
    ∃ t ∈ S.pts, min (c0 * ((b : ℝ) - 1) / (Real.sqrt 2 * (b : ℝ) ^ 2) * sepLB /
        Real.sqrt (S.s t)) ((2 * n / 3 : ℝ) - S.s t) < (S.ρ t : ℝ) + S.q t + 1 := by
  have hW : ∀ Y : Finset (Fin n), IsSep I Y → (∀ v ∈ Y, 2 ≤ dOcc I v) →
      wt b * sepLB ≤ ∑ v ∈ Y, wt (dOcc I v) := by
    intro Y hY h2
    have hwb : 0 ≤ wt b := by
      unfold wt
      have : (2 : ℝ) ≤ b := by exact_mod_cast hb
      apply div_nonneg <;> nlinarith
    calc wt b * sepLB ≤ wt b * Y.card :=
          mul_le_mul_of_nonneg_left (by exact_mod_cast hsep Y hY) hwb
      _ = ∑ v ∈ Y, wt b := by rw [Finset.sum_const, nsmul_eq_mul, mul_comm]
      _ ≤ ∑ v ∈ Y, wt (dOcc I v) :=
          Finset.sum_le_sum fun v hv => wt_anti _ _ (h2 v hv) (hdb v)
  obtain ⟨t, ht, h⟩ := S.theoremB' c0 (wt b * sepLB) hc0 hL2 hW
  refine ⟨t, ht, ?_⟩
  convert h using 3
  have hb2 : (2 : ℝ) ≤ b := by exact_mod_cast hb
  have hb0 : (b : ℝ) ≠ 0 := by linarith
  have hs2 : Real.sqrt 2 ≠ 0 := by positivity
  unfold wt
  field_simp

/-- **Theorem B, unconditional, without `n ≤ M`**: the only remaining hypothesis is that `W`
is a lower bound for the weighted separator number. -/
theorem theoremB_unconditional' (W : ℝ)
    (hW : ∀ Y : Finset (Fin n), IsSep I Y → (∀ v ∈ Y, 2 ≤ dOcc I v) →
      W ≤ ∑ v ∈ Y, wt (dOcc I v)) :
    ∃ t ∈ S.pts, min (UAvg.c0 / (2 * Real.sqrt 2) * W / Real.sqrt (S.s t)) (2 * n / 3 - S.s t)
      < S.ρ t + S.q t + 1 :=
  S.theoremB' UAvg.c0 W UAvg.c0_pos.le S.lemma2 hW

/-- **Theorem B, bounded occurrences, unconditional, without `n ≤ M`.** -/
theorem theoremB_bounded_unconditional'
    (b : ℕ) (hb : 2 ≤ b) (hdb : ∀ v, dOcc I v ≤ b) (sepLB : ℕ)
    (hsep : ∀ Y : Finset (Fin n), IsSep I Y → sepLB ≤ Y.card) :
    ∃ t ∈ S.pts, min (UAvg.c0 * ((b : ℝ) - 1) / (Real.sqrt 2 * (b : ℝ) ^ 2) * sepLB /
        Real.sqrt (S.s t)) ((2 * n / 3 : ℝ) - S.s t) < (S.ρ t : ℝ) + S.q t + 1 :=
  S.theoremB_bounded' UAvg.c0 UAvg.c0_pos.le S.lemma2 b hb hdb sepLB hsep

/-- Determinant + Fact R for variable-uniform schemes, **without** `n ≤ M`. -/
theorem exists_pointA' (hU : S.Uniform) :
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
  have h3 := S.sum_nullity_le_occ ω
  have h4 := S.factR
  have h5 : nOcc I * S.R ≤ M I * S.R := Nat.mul_le_mul_right _ nOcc_le_M
  omega

/-- **Theorem A** (variable-uniform deformations: no absorption), **without** `n ≤ M`. -/
theorem theoremA' (hU : S.Uniform) (sepLB : ℕ)
    (hsep : ∀ Y : Finset (Fin n), IsSep I Y → sepLB ≤ Y.card) :
    ∃ t ∈ S.pts, min (sepLB : ℝ) (2 * n / 3 - S.s t) < S.ρ t + S.q t + 1 := by
  classical
  obtain ⟨t, ht, hlt⟩ := S.exists_pointA' hU
  refine ⟨t, ht, ?_⟩
  set D := ∑ γ, (S.a t γ - (S.J t γ).card) with hD
  have hDlt : (D : ℝ) < S.ρ t + S.q t + 1 := by exact_mod_cast hlt
  have hpt : ∀ γ, S.a t γ ≤ (S.a t γ - (S.J t γ).card) + (S.J t γ).card := fun γ => by omega
  have hsa := S.sum_a t
  have hsJ := S.sum_J t
  by_cases hcase : ∀ γ, 3 * (S.J t γ).card ≤ 2 * n
  · have hsepY : IsSep I (S.split t) := by
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
  · push Not at hcase
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

#print axioms WitnessDet.natDegree_det_le_card
#print axioms WitnessDet.sum_nullity_le_card
#print axioms ClauseScheme.nOcc_le_M
#print axioms ClauseScheme.Scheme.lemma1_occ
#print axioms ClauseScheme.Scheme.avg_lemma1_occ
#print axioms ClauseScheme.Scheme.exists_point'
#print axioms ClauseScheme.Scheme.theoremB_unconditional'
#print axioms ClauseScheme.Scheme.theoremB_bounded_unconditional'
#print axioms ClauseScheme.Scheme.exists_pointA'
#print axioms ClauseScheme.Scheme.theoremA'
