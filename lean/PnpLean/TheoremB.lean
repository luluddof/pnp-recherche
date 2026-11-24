import PnpLean.Scheme

/-!
# Theorem B (conditional on Lemma 2), with Lemma 3 proved

Averaging Lemma 1 over all witness maps (uniform average = expectation for the independent
uniform choice of an occurrence per variable), Fact R, Lemma 3 (proved here) and Lemma 2
(taken as an explicit hypothesis `Lemma2Hyp`) give Theorem B of `note_limits_s9.md`.
-/

open Polynomial Matrix Module

set_option linter.unusedSectionVars false
set_option linter.style.numericProjections false

namespace ClauseScheme

variable {n : ℕ} {κ : Type*} [Fintype κ] [DecidableEq κ] {I : κ → Finset (Fin n)}
  {G : Type*} [Fintype G] [DecidableEq G]

noncomputable instance : Fintype (WitnessMap I) := by
  classical
  unfold WitnessMap
  infer_instance

instance : Nonempty (WitnessMap I) := by
  classical
  refine ⟨⟨fun v => if h : ∃ o : Occ I, o.var = v then some h.choose else none, ?_, ?_⟩⟩
  · intro v o ho
    have ho' : (if h : ∃ o : Occ I, o.var = v then some h.choose else none) = some o := ho
    by_cases h : ∃ o : Occ I, o.var = v
    · rw [dif_pos h] at ho'
      cases ho'
      exact h.choose_spec
    · rw [dif_neg h] at ho'; cases ho'
  · intro v hv o ho
    have hv' : (if h : ∃ o : Occ I, o.var = v then some h.choose else none) = none := hv
    by_cases h : ∃ o : Occ I, o.var = v
    · rw [dif_pos h] at hv'; cases hv'
    · exact h ⟨o, ho⟩

/-- Uniform average over all witness maps. Since the witness maps are exactly the
choice functions `v ↦ (occurrence of v)` (or `⊥` when `d_v = 0`), this is the expectation
under the independent uniform choice of Lemma 2. -/
noncomputable def avg (f : WitnessMap I → ℝ) : ℝ :=
  (∑ ω, f ω) / Fintype.card (WitnessMap I)

theorem card_pos : (0 : ℝ) < Fintype.card (WitnessMap I) := by
  exact_mod_cast Fintype.card_pos

theorem avg_mono {f g : WitnessMap I → ℝ} (h : ∀ ω, f ω ≤ g ω) : avg f ≤ avg g :=
  div_le_div_of_nonneg_right (Finset.sum_le_sum fun ω _ => h ω) card_pos.le

theorem avg_const (c : ℝ) : avg (fun _ : WitnessMap I => c) = c := by
  unfold avg
  rw [Finset.sum_const, Finset.card_univ, nsmul_eq_mul]
  field_simp [card_pos.ne']

theorem avg_sum {ι : Type*} (s : Finset ι) (f : ι → WitnessMap I → ℝ) :
    avg (fun ω => ∑ i ∈ s, f i ω) = ∑ i ∈ s, avg (f i) := by
  unfold avg
  rw [Finset.sum_comm, Finset.sum_div]

/-- Number of occurrences `d_v`. -/
def dOcc (I : κ → Finset (Fin n)) (v : Fin n) : ℕ :=
  (Finset.univ.filter (fun o : Occ I => o.var = v)).card

/-- `w(d) = 2(d-1)/d²`. -/
noncomputable def wt (d : ℕ) : ℝ := 2 * ((d : ℝ) - 1) / (d : ℝ) ^ 2

/-- The primal graph with the vertices of `Y` isolated (i.e. `Pr(F) - Y`, plus isolated
copies of the vertices of `Y`). -/
def primalMinus (I : κ → Finset (Fin n)) (Y : Finset (Fin n)) : SimpleGraph (Fin n) where
  Adj u v := u ≠ v ∧ (∃ e, u ∈ I e ∧ v ∈ I e) ∧ u ∉ Y ∧ v ∉ Y
  symm := ⟨fun _ _ ⟨h1, ⟨e, h2, h3⟩, h4, h5⟩ => ⟨Ne.symm h1, ⟨e, h3, h2⟩, h5, h4⟩⟩
  loopless := ⟨fun _ h => h.1 rfl⟩

/-- `Y` is a `2/3`-separator: every component of `Pr(F) - Y` has at most `2n/3` vertices. -/
def IsSep (I : κ → Finset (Fin n)) (Y : Finset (Fin n)) : Prop :=
  ∀ v ∉ Y, 3 * Set.ncard {u | (primalMinus I Y).Reachable v u} ≤ 2 * n

namespace Scheme

variable (S : Scheme I G)

/-- `𝔼 def_t` (average over witness maps). -/
noncomputable def Edef (t : ℚ) : ℝ := avg (fun ω => (S.deficiency ω t : ℝ))

/-- Averaged Lemma 1: `Σ_t (𝔼 def_t - ρ_t)_+ ≤ n R`. -/
theorem avg_lemma1 : ∑ t ∈ S.pts, max (S.Edef t - S.ρ t) 0 ≤ n * S.R := by
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
    _ ≤ avg (fun _ => ((n * S.R : ℕ) : ℝ)) := avg_mono fun ω => by
        rw [← Nat.cast_sum]; exact_mod_cast S.lemma1' ω
    _ = n * S.R := by rw [avg_const]; push_cast; ring

/-- Averaging + Fact R: if `M ≥ n`, some point has `𝔼 def_t < ρ_t + q_t + 1`. -/
theorem exists_point (hM : n ≤ M I) : ∃ t ∈ S.pts, S.Edef t < S.ρ t + S.q t + 1 := by
  by_contra hcon
  push Not at hcon
  have h1 : ∑ t ∈ S.pts, ((S.q t + 1 : ℕ) : ℝ) ≤ ∑ t ∈ S.pts, max (S.Edef t - S.ρ t) 0 :=
    Finset.sum_le_sum fun t ht => by
      push_cast
      exact le_max_of_le_left (by linarith [hcon t ht])
  have h2 := S.factR
  have h3 := S.avg_lemma1
  have h4 : (n * S.R : ℝ) ≤ M I * S.R := by
    exact_mod_cast Nat.mul_le_mul_right _ hM
  rw [← Nat.cast_sum] at h1
  have h5 : ((M I * S.R + 1 : ℕ) : ℝ) ≤ ((∑ t ∈ S.pts, (S.q t + 1) : ℕ) : ℝ) := by
    exact_mod_cast h2
  push_cast at h5 h1
  linarith

/-! ## Lemma 2 data and Lemma 3 -/

/-- `π_{v,γ}`: fraction of the occurrences of `v` in clauses assigned to `γ`
(for `d_v = 0`: indicator of the own group). -/
noncomputable def π (t : ℚ) (v : Fin n) (γ : G) : ℝ :=
  if dOcc I v = 0 then (if S.grp t v = γ then 1 else 0)
  else ((Finset.univ.filter (fun o : Occ I => o.var = v ∧ S.α t o.clause = γ)).card : ℝ) /
    dOcc I v

/-- `μ_γ = Σ_v π_{v,γ} = 𝔼 |I^ω_{t,γ}|`. -/
noncomputable def μ (t : ℚ) (γ : G) : ℝ := ∑ v, S.π t v γ

/-- Split variables `Σ_t`. -/
def split (t : ℚ) : Finset (Fin n) :=
  Finset.univ.filter (fun v => ∃ o o' : Occ I, o.var = v ∧ o'.var = v ∧
    S.α t o.clause ≠ S.α t o'.clause)

/-- `s_t = max_γ |S_{t,γ}|`. -/
def s (t : ℚ) : ℕ := Finset.univ.sup (S.a t)

/-- **Lemma 2** (expected deficiency), as a hypothesis, with constant `c0`. -/
def Lemma2Hyp (c0 : ℝ) : Prop :=
  ∀ t ∈ S.pts, max ((1 / 2) * ∑ γ, |S.μ t γ - S.a t γ|)
    (c0 / (2 * Real.sqrt (2 * S.s t)) * ∑ v ∈ S.split t, wt (dOcc I v)) ≤ S.Edef t

/-- Home group of a non-split variable. -/
noncomputable def homeNS (t : ℚ) (v : Fin n) : G := by
  classical
  exact if h : ∃ o : Occ I, o.var = v then S.α t h.choose.clause else S.grp t v

theorem homeNS_occ (t : ℚ) (v : Fin n) (hv : v ∉ S.split t) (o : Occ I) (ho : o.var = v) :
    S.α t o.clause = S.homeNS t v := by
  classical
  have h : ∃ o : Occ I, o.var = v := ⟨o, ho⟩
  simp only [homeNS, dif_pos h]
  by_contra hne
  exact hv (Finset.mem_filter.2 ⟨Finset.mem_univ _, o, h.choose, ho, h.choose_spec, hne⟩)

theorem homeNS_noocc (t : ℚ) (v : Fin n) (h : ¬ ∃ o : Occ I, o.var = v) :
    S.homeNS t v = S.grp t v := by
  classical
  simp only [homeNS, dif_neg h]

theorem homeNS_adj (t : ℚ) (u w : Fin n)
    (h : (primalMinus I (S.split t)).Adj u w) : S.homeNS t u = S.homeNS t w := by
  obtain ⟨_, ⟨e, hu, hw⟩, hu', hw'⟩ := h
  rw [← S.homeNS_occ t u hu' ⟨(e, u), hu⟩ rfl, ← S.homeNS_occ t w hw' ⟨(e, w), hw⟩ rfl]
  rfl

theorem reach_homeNS (t : ℚ) (v u : Fin n) (hv : v ∉ S.split t)
    (h : (primalMinus I (S.split t)).Reachable v u) :
    u ∉ S.split t ∧ S.homeNS t u = S.homeNS t v := by
  obtain ⟨p⟩ := h
  induction p with
  | nil => exact ⟨hv, rfl⟩
  | cons hadj _ ih =>
    obtain ⟨h1, h2⟩ := ih hadj.2.2.2
    exact ⟨h1, h2.trans (S.homeNS_adj t _ _ hadj).symm⟩

/-- `J_{t,γ}`: non-split variables with home group `γ`. -/
noncomputable def J (t : ℚ) (γ : G) : Finset (Fin n) := by
  classical
  exact Finset.univ.filter (fun v => v ∉ S.split t ∧ S.homeNS t v = γ)

theorem mem_J (t : ℚ) (γ : G) (v : Fin n) :
    v ∈ S.J t γ ↔ v ∉ S.split t ∧ S.homeNS t v = γ := by
  classical
  simp [J]

theorem sum_π (t : ℚ) (v : Fin n) : ∑ γ, S.π t v γ = 1 := by
  unfold π
  by_cases h0 : dOcc I v = 0
  · simp [h0]
  · simp only [h0, if_false]
    rw [← Finset.sum_div, div_eq_one_iff_eq (by exact_mod_cast h0)]
    have := Finset.card_eq_sum_card_fiberwise (f := fun o : Occ I => S.α t o.clause)
      (s := Finset.univ.filter (fun o : Occ I => o.var = v)) (t := Finset.univ)
      (fun _ _ => Finset.mem_coe.2 (Finset.mem_univ _))
    unfold dOcc
    rw [this]
    push_cast
    refine Finset.sum_congr rfl fun γ _ => ?_
    rw [Finset.filter_filter]

theorem sum_μ (t : ℚ) : ∑ γ, S.μ t γ = n := by
  unfold μ
  rw [Finset.sum_comm]
  simp [S.sum_π]

theorem sum_a (t : ℚ) : ∑ γ, S.a t γ = n := by
  have := Finset.card_eq_sum_card_fiberwise (f := S.grp t) (s := Finset.univ) (t := Finset.univ)
    (fun _ _ => Finset.mem_coe.2 (Finset.mem_univ _))
  simpa [a] using this.symm

theorem π_nonneg (t : ℚ) (v : Fin n) (γ : G) : 0 ≤ S.π t v γ := by
  unfold π; split_ifs <;> positivity

theorem π_J (t : ℚ) (γ : G) (v : Fin n) (hv : v ∈ S.J t γ) : S.π t v γ = 1 := by
  rw [S.mem_J] at hv
  obtain ⟨hns, hh⟩ := hv
  unfold π
  by_cases h0 : dOcc I v = 0
  · have hno : ¬ ∃ o : Occ I, o.var = v := by
      rintro ⟨o, ho⟩
      have : o ∈ Finset.univ.filter (fun o : Occ I => o.var = v) := by simp [ho]
      unfold dOcc at h0
      rw [Finset.card_eq_zero] at h0
      rw [h0] at this; simp at this
    rw [S.homeNS_noocc t v hno] at hh
    simp [h0, hh]
  · simp only [h0, if_false]
    rw [div_eq_one_iff_eq (by exact_mod_cast h0)]
    unfold dOcc
    congr 2
    refine Finset.filter_congr fun o _ => ⟨fun h => h.1, fun h => ⟨h, ?_⟩⟩
    rw [S.homeNS_occ t v hns o h, hh]

theorem le_half_sum_abs {x : G → ℝ} (hx : ∑ γ, x γ = 0) (γ₀ : G) :
    x γ₀ ≤ (1 / 2) * ∑ γ, |x γ| := by
  have h1 := Finset.add_sum_erase Finset.univ x (Finset.mem_univ γ₀)
  have h2 := Finset.add_sum_erase Finset.univ (fun γ => |x γ|) (Finset.mem_univ γ₀)
  have h3 := Finset.abs_sum_le_sum_abs x (Finset.univ.erase γ₀)
  have h4 : ∑ γ ∈ Finset.univ.erase γ₀, x γ = - x γ₀ := by linarith
  rw [h4, abs_neg] at h3
  have := le_abs_self (x γ₀)
  linarith

/-- **Lemma 3** (split variables form a separator, or the loads are unbalanced). -/
theorem lemma3 (t : ℚ) :
    IsSep I (S.split t) ∨
      (2 * n / 3 : ℝ) - S.s t < (1 / 2) * ∑ γ, |S.μ t γ - S.a t γ| := by
  by_cases hcase : ∀ γ, 3 * (S.J t γ).card ≤ 2 * n
  · left
    intro v hv
    refine le_trans ?_ (hcase (S.homeNS t v))
    refine Nat.mul_le_mul_left _ ?_
    rw [← Set.ncard_coe_finset]
    refine Set.ncard_le_ncard (fun u hu => ?_) (Finset.finite_toSet _)
    obtain ⟨h1, h2⟩ := S.reach_homeNS t v u hv hu
    exact (S.mem_J t _ u).2 ⟨h1, h2⟩
  · right
    push Not at hcase
    obtain ⟨γ₀, hγ₀⟩ := hcase
    have hμ : ((S.J t γ₀).card : ℝ) ≤ S.μ t γ₀ := by
      unfold μ
      calc ((S.J t γ₀).card : ℝ) = ∑ v ∈ S.J t γ₀, S.π t v γ₀ := by
            rw [Finset.sum_congr rfl fun v hv => S.π_J t γ₀ v hv]; simp
        _ ≤ ∑ v, S.π t v γ₀ :=
            Finset.sum_le_sum_of_subset_of_nonneg (Finset.subset_univ _)
              fun v _ _ => S.π_nonneg t v γ₀
    have ha : (S.a t γ₀ : ℝ) ≤ S.s t := by
      exact_mod_cast Finset.le_sup (f := S.a t) (Finset.mem_univ γ₀)
    have hsum : ∑ γ, (S.μ t γ - S.a t γ) = 0 := by
      rw [Finset.sum_sub_distrib, S.sum_μ, ← Nat.cast_sum, S.sum_a]; ring
    have := le_half_sum_abs hsum γ₀
    have hJ : (2 * n : ℝ) < 3 * (S.J t γ₀).card := by exact_mod_cast hγ₀
    linarith

theorem split_two_le (t : ℚ) (v : Fin n) (hv : v ∈ S.split t) : 2 ≤ dOcc I v := by
  obtain ⟨-, o, o', ho, ho', hne⟩ := Finset.mem_filter.1 hv
  unfold dOcc
  refine Finset.one_lt_card.2 ⟨o, by simp [ho], o', by simp [ho'], ?_⟩
  rintro rfl; exact hne rfl

/-- **Theorem B** (conditional on Lemma 2). If `M ≥ n`, `W` is a lower bound for the
weighted separator number `W(F)`, and Lemma 2 holds with a constant `c0 ≥ 0`, then some
`t ∈ X` has `ρ_t + q_t + 1 > min (c0/(2√2) · W/√s_t, 2n/3 - s_t)`. -/
theorem theoremB (hM : n ≤ M I) (c0 W : ℝ) (hc0 : 0 ≤ c0) (hL2 : S.Lemma2Hyp c0)
    (hW : ∀ Y : Finset (Fin n), IsSep I Y → (∀ v ∈ Y, 2 ≤ dOcc I v) →
      W ≤ ∑ v ∈ Y, wt (dOcc I v)) :
    ∃ t ∈ S.pts, min (c0 / (2 * Real.sqrt 2) * W / Real.sqrt (S.s t)) (2 * n / 3 - S.s t)
      < S.ρ t + S.q t + 1 := by
  obtain ⟨t, ht, hlt⟩ := S.exists_point hM
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

end Scheme

end ClauseScheme

#print axioms ClauseScheme.Scheme.avg_lemma1
#print axioms ClauseScheme.Scheme.exists_point
#print axioms ClauseScheme.Scheme.lemma3
#print axioms ClauseScheme.Scheme.theoremB
