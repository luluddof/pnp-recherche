import PnpLean.Article1.TheoremB
import PnpLean.Article1.Anticonc

/-!
# Lemma 2 (expected deficiency) and the unconditional Theorem B

* `witnessEquiv` : witness maps = product of per-variable choices (`Choice I v`), so the
  uniform average over witness maps is the expectation under independent uniform choices.
* `avg_home_ind`, `avg_Icard` : the averaging identity `𝔼 1[h^ω_t(v) = γ] = π_{v,γ}`,
  `𝔼 |I^ω_{t,γ}| = μ_γ`.
* `half_sum_abs_le_Edef` : first half of Lemma 2.
* `anticonc_le_Edef` : second half of Lemma 2 (uses `UAvg.bernoulli_anticonc`).
* `lemma2` : `Lemma2Hyp c0` holds with `c0 = (1 - e^{-2})/4`.
* `theoremB_unconditional`, `theoremB_bounded_unconditional`.
-/

open UAvg

set_option linter.unusedSectionVars false
set_option linter.style.numericProjections false

namespace ClauseScheme

variable {n : ℕ} {κ : Type*} [Fintype κ] [DecidableEq κ] {I : κ → Finset (Fin n)}
  {G : Type*} [Fintype G] [DecidableEq G]

/-- The possible values of a witness map at the variable `v`: an occurrence of `v`, or
`none` if `v` has no occurrence. -/
def Choice (I : κ → Finset (Fin n)) (v : Fin n) :=
  {x : Option (Occ I) //
    (∀ o, x = some o → o.var = v) ∧ (x = none → ∀ o : Occ I, o.var ≠ v)}

noncomputable instance (v : Fin n) : Fintype (Choice I v) := by
  classical
  unfold Choice
  infer_instance

instance (v : Fin n) : Nonempty (Choice I v) := by
  classical
  by_cases h : ∃ o : Occ I, o.var = v
  · obtain ⟨o, ho⟩ := h
    exact ⟨⟨some o, fun o' h' => (by cases h'; exact ho), fun h' => (by cases h')⟩⟩
  · exact ⟨⟨none, fun o h' => (by cases h'), fun _ o ho => h ⟨o, ho⟩⟩⟩

/-- Witness maps are exactly the families of per-variable choices. -/
def witnessEquiv : WitnessMap I ≃ (∀ v, Choice I v) where
  toFun ω := fun v => ⟨ω.1 v, ω.2.1 v, ω.2.2 v⟩
  invFun c := ⟨fun v => (c v).1, fun v => (c v).2.1, fun v => (c v).2.2⟩
  left_inv _ := rfl
  right_inv _ := rfl

theorem avg_eq_uavg (f : WitnessMap I → ℝ) : avg f = uavg f := rfl

/-- The uniform average over witness maps is the uniform average over the product of the
per-variable choice sets (independent uniform choices). -/
theorem avg_eq_uavg_pi (F : (∀ v, Choice I v) → ℝ) :
    avg (fun ω : WitnessMap I => F (witnessEquiv ω)) = uavg F :=
  uavg_equiv witnessEquiv F

theorem natsub_cast (a b : ℕ) : ((a - b : ℕ) : ℝ) = max ((a : ℝ) - b) 0 := by
  rcases le_total b a with h | h
  · rw [Nat.cast_sub h, max_eq_left]
    exact sub_nonneg.2 (by exact_mod_cast h)
  · rw [Nat.sub_eq_zero_of_le h, max_eq_right]
    · simp
    · exact sub_nonpos.2 (by exact_mod_cast h)

theorem max_zero_eq (x : ℝ) : max x 0 = (x + |x|) / 2 := by
  rcases le_total 0 x with h | h
  · rw [max_eq_left h, abs_of_nonneg h]; ring
  · rw [max_eq_right h, abs_of_nonpos h]; ring

namespace Scheme

variable (S : Scheme I G)

/-- Home group as a function of the choice at `v` only. -/
def homeC (t : ℚ) (v : Fin n) (c : Choice I v) : G :=
  match c.1 with
  | some o => S.α t o.clause
  | none => S.grp t v

theorem home_eq (ω : WitnessMap I) (t : ℚ) (v : Fin n) :
    S.home ω t v = S.homeC t v (witnessEquiv ω v) := rfl

/-- The event "the home group of `v` is `γ`", as a predicate on the choice at `v`. -/
def homeP (t : ℚ) (γ : G) (v : Fin n) (c : Choice I v) : Prop := S.homeC t v c = γ

instance (t : ℚ) (γ : G) (v : Fin n) : DecidablePred (S.homeP t γ v) := fun c => by
  unfold homeP; infer_instance

/-- `P(home of v = γ) = π_{v,γ}` under the uniform choice at `v`. -/
theorem choice_prob (t : ℚ) (v : Fin n) (γ : G) :
    uavg (fun c : Choice I v => if S.homeP t γ v c then (1 : ℝ) else 0) = S.π t v γ := by
  classical
  unfold π
  by_cases h0 : dOcc I v = 0
  · have hno : ∀ o : Occ I, o.var ≠ v := by
      intro o ho
      have : o ∈ Finset.univ.filter (fun o : Occ I => o.var = v) := by simp [ho]
      unfold dOcc at h0
      rw [Finset.card_eq_zero] at h0
      rw [h0] at this; simp at this
    have hc : ∀ c : Choice I v, S.homeC t v c = S.grp t v := by
      rintro ⟨x, h1, h2⟩
      cases x with
      | none => rfl
      | some o => exact absurd (h1 o rfl) (hno o)
    simp only [h0, ite_true, homeP, hc]
    exact uavg_const _
  · simp only [h0, ite_false]
    obtain ⟨o₀, ho₀⟩ : ∃ o : Occ I, o.var = v := by
      by_contra hne
      push Not at hne
      apply h0
      unfold dOcc
      rw [Finset.card_eq_zero, Finset.filter_eq_empty_iff]
      intro o _; exact hne o
    let φ : {o : Occ I // o.var = v} → Choice I v := fun o =>
      ⟨some o.1, fun o' h' => (by cases h'; exact o.2), fun h' => (by cases h')⟩
    have hφ : Function.Bijective φ := by
      constructor
      · intro o o' h
        have := congrArg Subtype.val h
        simp only [φ, Option.some.injEq] at this
        exact Subtype.ext this
      · rintro ⟨x, h1, h2⟩
        cases x with
        | none => exact absurd ho₀ (h2 rfl o₀)
        | some o => exact ⟨⟨o, h1 o rfl⟩, rfl⟩
    have h1 := uavg_equiv (Equiv.ofBijective φ hφ)
      (fun c : Choice I v => if S.homeP t γ v c then (1 : ℝ) else 0)
    rw [← h1]
    have h2 : ∀ o : {o : Occ I // o.var = v},
        (if S.homeP t γ v (Equiv.ofBijective φ hφ o) then (1 : ℝ) else 0) =
          if S.α t o.1.clause = γ then 1 else 0 := fun o => rfl
    unfold uavg
    congr 1
    · rw [Finset.sum_congr rfl fun o _ => h2 o,
        ← Finset.sum_subtype (Finset.univ.filter (fun o : Occ I => o.var = v)) (by simp)
          (fun o => if S.α t o.clause = γ then (1 : ℝ) else 0),
        Finset.sum_boole, Finset.filter_filter]
    · rw [Fintype.card_subtype]; rfl

/-- **Averaging identity**: under the uniform average over all witness maps, the indicator
of `h^ω_t(v) = γ` has average `π_{v,γ}`. -/
theorem avg_home_ind (t : ℚ) (v : Fin n) (γ : G) :
    avg (fun ω : WitnessMap I => if S.home ω t v = γ then (1 : ℝ) else 0) = S.π t v γ := by
  rw [← S.choice_prob t v γ,
    ← uavg_pi_coord (Ω := Choice I) v (fun c => if S.homeP t γ v c then (1 : ℝ) else 0)]
  exact avg_eq_uavg_pi
    (fun c : (∀ v, Choice I v) => if S.homeP t γ v (c v) then (1 : ℝ) else 0)

theorem Icard_cast (ω : WitnessMap I) (t : ℚ) (γ : G) :
    (S.Icard ω t γ : ℝ) = ∑ v, if S.home ω t v = γ then (1 : ℝ) else 0 := by
  unfold Icard
  rw [Finset.card_filter]
  push_cast
  rfl

/-- **Averaging identity**: `𝔼 |I^ω_{t,γ}| = μ_γ`. -/
theorem avg_Icard (t : ℚ) (γ : G) :
    avg (fun ω : WitnessMap I => (S.Icard ω t γ : ℝ)) = S.μ t γ := by
  have e : (fun ω : WitnessMap I => (S.Icard ω t γ : ℝ)) =
      fun ω => ∑ v ∈ Finset.univ, (fun v ω => if S.home ω t v = γ then (1 : ℝ) else 0) v ω :=
    funext fun ω => S.Icard_cast ω t γ
  rw [e, avg_sum]
  unfold μ
  exact Finset.sum_congr rfl fun v _ => S.avg_home_ind t v γ

theorem sum_Icard (ω : WitnessMap I) (t : ℚ) : ∑ γ, S.Icard ω t γ = n := by
  have := Finset.card_eq_sum_card_fiberwise (f := S.home ω t) (s := Finset.univ)
    (t := Finset.univ) (fun _ _ => Finset.mem_coe.2 (Finset.mem_univ _))
  simpa [Icard] using this.symm

theorem Edef_eq_sum (t : ℚ) :
    S.Edef t = ∑ γ, avg (fun ω : WitnessMap I => ((S.a t γ - S.Icard ω t γ : ℕ) : ℝ)) := by
  unfold Edef
  rw [← avg_sum]
  unfold deficiency
  congr 1
  funext ω
  push_cast
  rfl

/-- **Lemma 2, first half**: `𝔼 def_t ≥ ½ Σ_γ |μ_γ - a_γ|`. -/
theorem half_sum_abs_le_Edef (t : ℚ) :
    (1 / 2) * ∑ γ, |S.μ t γ - S.a t γ| ≤ S.Edef t := by
  have hpt : ∀ γ, max ((S.a t γ : ℝ) - S.μ t γ) 0 ≤
      avg (fun ω : WitnessMap I => ((S.a t γ - S.Icard ω t γ : ℕ) : ℝ)) := by
    intro γ
    apply max_le
    · have h1 : (S.a t γ : ℝ) - S.μ t γ =
          avg (fun ω : WitnessMap I => (S.a t γ : ℝ) - S.Icard ω t γ) := by
        rw [avg_eq_uavg, uavg_sub, uavg_const, ← avg_eq_uavg, S.avg_Icard]
      rw [h1]
      exact avg_mono fun ω => by rw [natsub_cast]; exact le_max_left _ _
    · rw [avg_eq_uavg]; exact uavg_nonneg fun ω => Nat.cast_nonneg _
  have hsum : ∑ γ, ((S.a t γ : ℝ) - S.μ t γ) = 0 := by
    rw [Finset.sum_sub_distrib, S.sum_μ, ← Nat.cast_sum, S.sum_a]; ring
  calc (1 / 2) * ∑ γ, |S.μ t γ - S.a t γ|
      = ∑ γ, max ((S.a t γ : ℝ) - S.μ t γ) 0 := by
        have : ∀ γ, max ((S.a t γ : ℝ) - S.μ t γ) 0 =
            ((S.a t γ : ℝ) - S.μ t γ) / 2 + |S.μ t γ - S.a t γ| / 2 := by
          intro γ; rw [max_zero_eq, abs_sub_comm]; ring
        rw [Finset.sum_congr rfl fun γ _ => this γ, Finset.sum_add_distrib, ← Finset.sum_div,
          ← Finset.sum_div, hsum]
        ring
    _ ≤ ∑ γ, avg (fun ω : WitnessMap I => ((S.a t γ - S.Icard ω t γ : ℕ) : ℝ)) :=
        Finset.sum_le_sum fun γ _ => hpt γ
    _ = S.Edef t := (S.Edef_eq_sum t).symm

/-! ## Second half of Lemma 2 -/

/-- Pointwise: `def^ω_t = ½ Σ_γ | |I^ω_γ| - a_γ |`. -/
theorem deficiency_eq (ω : WitnessMap I) (t : ℚ) :
    (S.deficiency ω t : ℝ) = (1 / 2) * ∑ γ, |(S.Icard ω t γ : ℝ) - S.a t γ| := by
  have hsum : ∑ γ, ((S.a t γ : ℝ) - S.Icard ω t γ) = 0 := by
    rw [Finset.sum_sub_distrib, ← Nat.cast_sum, ← Nat.cast_sum, S.sum_a, S.sum_Icard]; ring
  have : ∀ γ, ((S.a t γ - S.Icard ω t γ : ℕ) : ℝ) =
      ((S.a t γ : ℝ) - S.Icard ω t γ) / 2 + |(S.Icard ω t γ : ℝ) - S.a t γ| / 2 := by
    intro γ; rw [natsub_cast, max_zero_eq, abs_sub_comm]; ring
  unfold deficiency
  push_cast
  rw [Finset.sum_congr rfl fun γ _ => this γ, Finset.sum_add_distrib, ← Finset.sum_div,
    ← Finset.sum_div, hsum]
  ring

/-- `V_γ = Σ_v π_{v,γ} (1 - π_{v,γ})`, the variance of `|I^ω_{t,γ}|`. -/
noncomputable def Vg (t : ℚ) (γ : G) : ℝ := ∑ v, S.π t v γ * (1 - S.π t v γ)

theorem π_le_one (t : ℚ) (v : Fin n) (γ : G) : S.π t v γ ≤ 1 := by
  rw [← S.sum_π t v]
  exact Finset.single_le_sum (fun γ _ => S.π_nonneg t v γ) (Finset.mem_univ γ)

theorem Vg_nonneg (t : ℚ) (γ : G) : 0 ≤ S.Vg t γ :=
  Finset.sum_nonneg fun v _ =>
    mul_nonneg (S.π_nonneg t v γ) (by linarith [S.π_le_one t v γ])

theorem Vg_le_μ (t : ℚ) (γ : G) : S.Vg t γ ≤ S.μ t γ :=
  Finset.sum_le_sum fun v _ => by nlinarith [S.π_nonneg t v γ, S.π_le_one t v γ]

/-- Anticoncentration for the load of one group: `|I^ω_{t,γ}|` is a sum of independent
Bernoulli variables with parameters `π_{v,γ}`. -/
theorem avg_abs_Icard_ge (t : ℚ) (γ : G) (a : ℝ) :
    c0 * min (S.Vg t γ) (Real.sqrt (S.Vg t γ)) ≤
      avg (fun ω : WitnessMap I => |(S.Icard ω t γ : ℝ) - a|) := by
  have h := bernoulli_anticonc (S.homeP t γ) a
  have hV : bV (S.homeP t γ) = S.Vg t γ :=
    Finset.sum_congr rfl fun v _ => by
      have : bp (S.homeP t γ) v = S.π t v γ := S.choice_prob t v γ
      rw [this]
  rw [hV, ← avg_eq_uavg_pi] at h
  have e : ∀ ω : WitnessMap I, bX (S.homeP t γ) (witnessEquiv ω) = (S.Icard ω t γ : ℝ) := by
    intro ω; rw [S.Icard_cast]; rfl
  simpa only [e] using h

/-- One group: `𝔼 | |I_γ| - a_γ | ≥ c0/√(2s) · V_γ`. -/
theorem one_group (t : ℚ) (hs : 1 ≤ S.s t) (γ : G) :
    c0 / Real.sqrt (2 * S.s t) * S.Vg t γ ≤
      avg (fun ω : WitnessMap I => |(S.Icard ω t γ : ℝ) - S.a t γ|) := by
  set r := Real.sqrt (2 * S.s t) with hr
  set V := S.Vg t γ with hVdef
  have hs' : (1 : ℝ) ≤ S.s t := by exact_mod_cast hs
  have hr0 : 0 ≤ r := Real.sqrt_nonneg _
  have hr2 : r ^ 2 = 2 * S.s t := Real.sq_sqrt (by linarith)
  have hr1 : 1 ≤ r := by
    by_contra h
    push Not at h
    nlinarith
  have hV0 : 0 ≤ V := S.Vg_nonneg t γ
  have hVμ : V ≤ S.μ t γ := S.Vg_le_μ t γ
  have ha : (S.a t γ : ℝ) ≤ S.s t := by
    exact_mod_cast Finset.le_sup (f := S.a t) (Finset.mem_univ γ)
  have hc := c0_pos
  rcases le_or_gt (S.μ t γ) (2 * S.s t) with hcase | hcase
  · have hanti := S.avg_abs_Icard_ge t γ (S.a t γ)
    have hsV0 := Real.sqrt_nonneg V
    have hsV2 : Real.sqrt V * Real.sqrt V = V := Real.mul_self_sqrt hV0
    have hsVr : Real.sqrt V ≤ r := Real.sqrt_le_sqrt (hVμ.trans hcase)
    have h1 : V / r ≤ V := div_le_self hV0 hr1
    have h2 : V / r ≤ Real.sqrt V := by
      rw [div_le_iff₀ (by linarith)]
      nlinarith [mul_le_mul_of_nonneg_left hsVr hsV0]
    calc c0 / r * V = c0 * (V / r) := by ring
      _ ≤ c0 * min V (Real.sqrt V) := mul_le_mul_of_nonneg_left (le_min h1 h2) hc.le
      _ ≤ _ := hanti
  · have h1 : |S.μ t γ - S.a t γ| ≤
        avg (fun ω : WitnessMap I => |(S.Icard ω t γ : ℝ) - S.a t γ|) := by
      have := abs_uavg_le (fun ω : WitnessMap I => (S.Icard ω t γ : ℝ) - S.a t γ)
      rw [uavg_sub, uavg_const, ← avg_eq_uavg, S.avg_Icard] at this
      exact this
    have h2 := le_abs_self (S.μ t γ - S.a t γ)
    have h3 : c0 / r * V ≤ c0 * V :=
      mul_le_mul_of_nonneg_right (div_le_self hc.le hr1) hV0
    have h4 : c0 * V ≤ 1 / 4 * V := mul_le_mul_of_nonneg_right c0_le hV0
    linarith

/-- A split variable contributes at least `w(d_v)` to the total variance. -/
theorem wt_le_of_split (t : ℚ) (v : Fin n) (hv : v ∈ S.split t) :
    wt (dOcc I v) ≤ ∑ γ, S.π t v γ * (1 - S.π t v γ) := by
  classical
  have hd := S.split_two_le t v hv
  obtain ⟨-, o, o', ho, ho', hne⟩ := Finset.mem_filter.1 hv
  have hd0 : dOcc I v ≠ 0 := by omega
  have hdR : (2 : ℝ) ≤ dOcc I v := by exact_mod_cast hd
  have hdpos : (0 : ℝ) < dOcc I v := by linarith
  have hge : ∀ o : Occ I, o.var = v → 1 / (dOcc I v : ℝ) ≤ S.π t v (S.α t o.clause) := by
    intro o ho
    unfold π
    simp only [hd0, ite_false]
    refine div_le_div_of_nonneg_right ?_ hdpos.le
    have : 0 < (Finset.univ.filter
        (fun o' : Occ I => o'.var = v ∧ S.α t o'.clause = S.α t o.clause)).card :=
      Finset.card_pos.2 ⟨o, by simp [ho]⟩
    exact_mod_cast this
  set γ1 := S.α t o.clause
  set γ2 := S.α t o'.clause
  set u := 1 / (dOcc I v : ℝ) with hu
  set x := S.π t v γ1
  set y := S.π t v γ2
  have hx : u ≤ x := hge o ho
  have hy : u ≤ y := hge o' ho'
  have hxy : x + y ≤ 1 := by
    rw [← S.sum_π t v, ← Finset.sum_pair hne]
    exact Finset.sum_le_sum_of_subset_of_nonneg (Finset.subset_univ _)
      fun γ _ _ => S.π_nonneg t v γ
  have hsum : x * (1 - x) + y * (1 - y) ≤ ∑ γ, S.π t v γ * (1 - S.π t v γ) := by
    rw [← Finset.sum_pair (f := fun γ => S.π t v γ * (1 - S.π t v γ)) hne]
    exact Finset.sum_le_sum_of_subset_of_nonneg (Finset.subset_univ _)
      fun γ _ _ => mul_nonneg (S.π_nonneg t v γ) (by linarith [S.π_le_one t v γ])
  have hw : wt (dOcc I v) = 2 * (u * (1 - u)) := by
    unfold wt
    rw [hu]
    field_simp
  rw [hw]
  have h1 : 0 ≤ (x - u) * (1 - u - x) := mul_nonneg (by linarith) (by linarith)
  have h2 : 0 ≤ (y - u) * (1 - u - y) := mul_nonneg (by linarith) (by linarith)
  nlinarith

theorem sum_Vg_ge (t : ℚ) : ∑ v ∈ S.split t, wt (dOcc I v) ≤ ∑ γ, S.Vg t γ := by
  unfold Vg
  rw [Finset.sum_comm]
  calc ∑ v ∈ S.split t, wt (dOcc I v)
      ≤ ∑ v ∈ S.split t, ∑ γ, S.π t v γ * (1 - S.π t v γ) :=
        Finset.sum_le_sum fun v hv => S.wt_le_of_split t v hv
    _ ≤ ∑ v, ∑ γ, S.π t v γ * (1 - S.π t v γ) :=
        Finset.sum_le_sum_of_subset_of_nonneg (Finset.subset_univ _) fun v _ _ =>
          Finset.sum_nonneg fun γ _ =>
            mul_nonneg (S.π_nonneg t v γ) (by linarith [S.π_le_one t v γ])

theorem Edef_nonneg (t : ℚ) : 0 ≤ S.Edef t := by
  unfold Edef; rw [avg_eq_uavg]; exact uavg_nonneg fun ω => Nat.cast_nonneg _

/-- **Lemma 2, second half**: `𝔼 def_t ≥ c0/(2√(2 s_t)) · Σ_{v split} w(d_v)`. -/
theorem anticonc_le_Edef (t : ℚ) :
    c0 / (2 * Real.sqrt (2 * S.s t)) * ∑ v ∈ S.split t, wt (dOcc I v) ≤ S.Edef t := by
  by_cases hs0 : S.s t = 0
  · rw [hs0]
    simpa using S.Edef_nonneg t
  · have hs : 1 ≤ S.s t := Nat.one_le_iff_ne_zero.2 hs0
    have hE : S.Edef t =
        (1 / 2) * ∑ γ, avg (fun ω : WitnessMap I => |(S.Icard ω t γ : ℝ) - S.a t γ|) := by
      unfold Edef
      rw [← avg_sum, avg_eq_uavg, avg_eq_uavg, ← uavg_const_mul]
      exact uavg_congr fun ω => S.deficiency_eq ω t
    have hk : 0 ≤ c0 / Real.sqrt (2 * S.s t) :=
      div_nonneg c0_pos.le (Real.sqrt_nonneg _)
    rw [hE]
    calc c0 / (2 * Real.sqrt (2 * S.s t)) * ∑ v ∈ S.split t, wt (dOcc I v)
        = (1 / 2) * (c0 / Real.sqrt (2 * S.s t) * ∑ v ∈ S.split t, wt (dOcc I v)) := by
          rw [mul_comm 2 (Real.sqrt _), ← div_div]; ring
      _ ≤ (1 / 2) * (c0 / Real.sqrt (2 * S.s t) * ∑ γ, S.Vg t γ) := by
          have := mul_le_mul_of_nonneg_left (S.sum_Vg_ge t) hk
          linarith
      _ = (1 / 2) * ∑ γ, c0 / Real.sqrt (2 * S.s t) * S.Vg t γ := by rw [Finset.mul_sum]
      _ ≤ _ := by
          have := Finset.sum_le_sum fun γ (_ : γ ∈ Finset.univ) => S.one_group t hs γ
          linarith

/-- **Lemma 2** (expected deficiency) of the paper, fully proved:
`Lemma2Hyp` holds with `c0 = (1 - e^{-2})/4`. -/
theorem lemma2 : S.Lemma2Hyp c0 := fun t _ =>
  max_le (S.half_sum_abs_le_Edef t) (S.anticonc_le_Edef t)

/-- **Theorem B, unconditional** (no hypothesis on Lemma 2). -/
theorem theoremB_unconditional (hM : n ≤ M I) (W : ℝ)
    (hW : ∀ Y : Finset (Fin n), IsSep I Y → (∀ v ∈ Y, 2 ≤ dOcc I v) →
      W ≤ ∑ v ∈ Y, wt (dOcc I v)) :
    ∃ t ∈ S.pts, min (c0 / (2 * Real.sqrt 2) * W / Real.sqrt (S.s t)) (2 * n / 3 - S.s t)
      < S.ρ t + S.q t + 1 :=
  S.theoremB hM c0 W c0_pos.le S.lemma2 hW

/-- **Theorem B, bounded occurrences, unconditional**. -/
theorem theoremB_bounded_unconditional (hM : n ≤ M I)
    (b : ℕ) (hb : 2 ≤ b) (hdb : ∀ v, dOcc I v ≤ b) (sepLB : ℕ)
    (hsep : ∀ Y : Finset (Fin n), IsSep I Y → sepLB ≤ Y.card) :
    ∃ t ∈ S.pts, min (c0 * ((b : ℝ) - 1) / (Real.sqrt 2 * (b : ℝ) ^ 2) * sepLB /
        Real.sqrt (S.s t)) ((2 * n / 3 : ℝ) - S.s t) < (S.ρ t : ℝ) + S.q t + 1 :=
  S.theoremB_bounded hM c0 c0_pos.le S.lemma2 b hb hdb sepLB hsep

end Scheme

end ClauseScheme

#print axioms ClauseScheme.Scheme.avg_home_ind
#print axioms ClauseScheme.Scheme.avg_Icard
#print axioms ClauseScheme.Scheme.half_sum_abs_le_Edef
#print axioms ClauseScheme.Scheme.anticonc_le_Edef
#print axioms ClauseScheme.Scheme.lemma2
#print axioms ClauseScheme.Scheme.theoremB_unconditional
#print axioms ClauseScheme.Scheme.theoremB_bounded_unconditional
