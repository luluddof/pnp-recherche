import PnpLean.WitnessDet

/-!
# Clause-product interpolation schemes (Definition 5.1 of `note_limits_s9.md`)

Deterministic part: the scheme structure, witness maps, Lemma 1, Fact R.
-/

open Polynomial Matrix Module

namespace ClauseScheme

set_option linter.unusedSectionVars false

variable {n : ℕ} {κ : Type*} [Fintype κ] [DecidableEq κ]

/-- Occurrences `O = {(e,v) : v ∈ I_e}` of a CNF whose clause `e` has variable set `I e`. -/
def Occ (I : κ → Finset (Fin n)) := {p : κ × Fin n // p.2 ∈ I p.1}

instance (I : κ → Finset (Fin n)) : Fintype (Occ I) := by
  unfold Occ; infer_instance

instance (I : κ → Finset (Fin n)) : DecidableEq (Occ I) := by
  unfold Occ; infer_instance

/-- Clause of an occurrence. -/
def Occ.clause {I : κ → Finset (Fin n)} (o : Occ I) : κ := o.1.1
/-- Variable of an occurrence. -/
def Occ.var {I : κ → Finset (Fin n)} (o : Occ I) : Fin n := o.1.2

/-- Index type of the Hermite data `(f^{(j)}(t))_{t ∈ X, j ≤ q_t}`. -/
abbrev HIdx (pts : Finset ℚ) (q : ℚ → ℕ) := Σ t : pts, Fin (q t + 1)

/-- Hermite data of a polynomial. -/
noncomputable def hermiteData (pts : Finset ℚ) (q : ℚ → ℕ) (f : ℚ[X]) : HIdx pts q → ℚ :=
  fun p => (derivative^[p.2] f).eval (p.1 : ℚ)

/-- A clause-product interpolation scheme for the CNF `I` (Definition 5.1), with groups
labelled by a fixed finite type `G` at every point. -/
structure Scheme (I : κ → Finset (Fin n)) (G : Type*) [Fintype G] [DecidableEq G] where
  /-- degree bound of the deformation -/
  R : ℕ
  /-- (a) the deformation vector `P_o(T) ∈ ℚ[T]^n` -/
  P : Occ I → Fin n → ℚ[X]
  P_deg : ∀ o j, (P o j).natDegree ≤ R
  P_zero : ∀ o, (fun j => (P o j).eval 0) = Pi.single o.var 1
  /-- (b) recovery points `X ⊂ ℚ \ {0}` -/
  pts : Finset ℚ
  zero_not_mem : (0 : ℚ) ∉ pts
  /-- Taylor orders -/
  q : ℚ → ℕ
  /-- linear functional on the Hermite data -/
  Λ : (HIdx pts q → ℚ) →ₗ[ℚ] ℚ
  recovery : ∀ f : ℚ[X], f.natDegree ≤ Fintype.card (Occ I) * R →
    Λ (hermiteData pts q f) = f.eval 0
  /-- (c) designs: the partition at `t` (variable ↦ its group) -/
  grp : ℚ → Fin n → G
  /-- clause ↦ its group -/
  α : ℚ → κ → G
  /-- measurement subspaces `𝓡_{t,γ} ⊆ ℚ^{S_{t,γ}}` -/
  Rsp : ℚ → G → Submodule ℚ (Fin n → ℚ)
  Rsp_supp : ∀ t γ, ∀ x ∈ Rsp t γ, ∀ j, grp t j ≠ γ → x j = 0
  /-- condition (3) -/
  cond3 : ∀ t ∈ pts, ∀ o : Occ I, ∀ γ, α t o.clause ≠ γ →
    (fun j => if grp t j = γ then (P o j).eval t else 0) ∈ Rsp t γ

variable {I : κ → Finset (Fin n)} {G : Type*} [Fintype G] [DecidableEq G]

/-- `M = |O|`. -/
def M (I : κ → Finset (Fin n)) : ℕ := Fintype.card (Occ I)

/-- A witness map: each variable gets one of its occurrences, or `⊥` (= `none`) if it has none. -/
def WitnessMap (I : κ → Finset (Fin n)) :=
  {ω : Fin n → Option (Occ I) //
    (∀ v o, ω v = some o → o.var = v) ∧ (∀ v, ω v = none → ∀ o : Occ I, o.var ≠ v)}

namespace Scheme

variable (S : Scheme I G)

/-- Row `v` of `Φ^ω`: `P_{ω(v)}`, or `e_v` if `ω(v) = ⊥`. -/
noncomputable def row (ω : WitnessMap I) (v : Fin n) : Fin n → ℚ[X] :=
  match ω.1 v with
  | some o => S.P o
  | none => fun j => if j = v then 1 else 0

/-- The witness matrix `Φ^ω(T)`. -/
noncomputable def Φ (ω : WitnessMap I) : Matrix (Fin n) (Fin n) ℚ[X] := fun v j => S.row ω v j

/-- Home group `h^ω_t(v)`. -/
def home (ω : WitnessMap I) (t : ℚ) (v : Fin n) : G :=
  match ω.1 v with
  | some o => S.α t o.clause
  | none => S.grp t v

/-- `a_γ = |S_{t,γ}|`. -/
def a (t : ℚ) (γ : G) : ℕ := (Finset.univ.filter (fun j => S.grp t j = γ)).card

/-- `|I^ω_{t,γ}|`. -/
def Icard (ω : WitnessMap I) (t : ℚ) (γ : G) : ℕ :=
  (Finset.univ.filter (fun v => S.home ω t v = γ)).card

/-- Deficiency `def^ω_t = Σ_γ (a_γ - |I^ω_{t,γ}|)_+` (truncated subtraction in `ℕ`). -/
def deficiency (ω : WitnessMap I) (t : ℚ) : ℕ := ∑ γ, (S.a t γ - S.Icard ω t γ)

/-- `ρ_t = Σ_γ dim 𝓡_{t,γ}`. -/
noncomputable def ρ (t : ℚ) : ℕ := ∑ γ, finrank ℚ (S.Rsp t γ)

theorem Φ_natDegree (ω : WitnessMap I) (i j : Fin n) : (S.Φ ω i j).natDegree ≤ S.R := by
  unfold Φ row
  rcases hω : ω.1 i with _ | o
  · simp only; split_ifs <;> simp
  · exact S.P_deg _ _

theorem Φ_eval_zero (ω : WitnessMap I) : WitnessDet.evalMat (S.Φ ω) 0 = 1 := by
  ext i j
  simp only [WitnessDet.evalMat, Matrix.map_apply, Φ, row]
  rcases hω : ω.1 i with _ | o
  · simp only
    by_cases hij : j = i
    · subst hij; simp
    · simp [hij, Ne.symm hij]
  · have h1 := congrFun (S.P_zero o) j
    rw [h1, ω.2.1 i o hω, Matrix.one_apply, Pi.single_apply]
    by_cases hij : j = i
    · subst hij; simp
    · simp [hij, Ne.symm hij]

theorem row_mem (ω : WitnessMap I) (t : ℚ) (ht : t ∈ S.pts) (γ : G) (v : Fin n)
    (hv : S.home ω t v ≠ γ) :
    (fun j => if S.grp t j = γ then WitnessDet.evalMat (S.Φ ω) t v j else 0) ∈ S.Rsp t γ := by
  simp only [WitnessDet.evalMat, Matrix.map_apply, Φ, row]
  unfold home at hv
  split at hv
  · rename_i o ho
    simp only [ho]
    exact S.cond3 t ht o γ hv
  · rename_i ho
    simp only [ho]
    convert (S.Rsp t γ).zero_mem using 1
    funext j
    by_cases hj : j = v
    · subst hj; simp [hv]
    · simp [hj]

/-- Rank bound at a point: `def^ω_t - ρ_t ≤ dim ker Φ^ω(t)`. -/
theorem deficiency_sub_le_nullity (ω : WitnessMap I) (t : ℚ) (ht : t ∈ S.pts) :
    S.deficiency ω t - S.ρ t ≤ WitnessDet.nullity (WitnessDet.evalMat (S.Φ ω) t) := by
  have := WitnessDet.deficiency_le_nullity_add (WitnessDet.evalMat (S.Φ ω) t) (S.grp t)
    (S.home ω t) (S.Rsp t) (S.Rsp_supp t) (fun γ v hv => S.row_mem ω t ht γ v hv)
  unfold deficiency ρ a Icard
  omega

/-- **Lemma 1** (witness determinant):
`Σ_{t∈X} (def^ω_t - ρ_t)_+ ≤ Σ_{t∈X} dim ker Φ^ω(t) ≤ n R`. -/
theorem lemma1 (ω : WitnessMap I) :
    ∑ t ∈ S.pts, (S.deficiency ω t - S.ρ t)
      ≤ ∑ t ∈ S.pts, WitnessDet.nullity (WitnessDet.evalMat (S.Φ ω) t) ∧
    ∑ t ∈ S.pts, WitnessDet.nullity (WitnessDet.evalMat (S.Φ ω) t) ≤ n * S.R :=
  ⟨Finset.sum_le_sum fun t ht => S.deficiency_sub_le_nullity ω t ht,
   WitnessDet.sum_nullity_le _ _ (S.Φ_natDegree ω) (S.Φ_eval_zero ω) _⟩

theorem lemma1' (ω : WitnessMap I) :
    ∑ t ∈ S.pts, (S.deficiency ω t - S.ρ t) ≤ n * S.R :=
  (S.lemma1 ω).1.trans (S.lemma1 ω).2

/-! ## Fact R -/

/-- **Fact R** (Hermite counting): `Σ_{t∈X} (q_t + 1) ≥ M R + 1`. -/
theorem factR : M I * S.R + 1 ≤ ∑ t ∈ S.pts, (S.q t + 1) := by
  by_contra hlt
  push Not at hlt
  set f : ℚ[X] := ∏ t ∈ S.pts, (X - C t) ^ (S.q t + 1) with hf
  have hmon : f.Monic := monic_prod_of_monic _ _ (fun t _ => (monic_X_sub_C t).pow _)
  have hf0 : f ≠ 0 := hmon.ne_zero
  have hdeg : f.natDegree = ∑ t ∈ S.pts, (S.q t + 1) := by
    rw [hf, natDegree_prod_of_monic _ _ (fun t _ => (monic_X_sub_C t).pow _)]
    simp
  have hdata : hermiteData S.pts S.q f = 0 := by
    funext ⟨⟨t, ht⟩, j⟩
    simp only [hermiteData, Pi.zero_apply]
    apply isRoot_iterate_derivative_of_lt_rootMultiplicity
    have hdvd : (X - C t) ^ (S.q t + 1) ∣ f := Finset.dvd_prod_of_mem _ ht
    have := (le_rootMultiplicity_iff hf0).2 hdvd
    have hj : (j : ℕ) < S.q t + 1 := j.2
    omega
  have hrec := S.recovery f (by unfold M at hlt; omega)
  rw [hdata, map_zero] at hrec
  have : f.eval 0 ≠ 0 := by
    rw [hf, eval_prod]
    refine Finset.prod_ne_zero_iff.2 fun t ht => ?_
    simp only [eval_pow, eval_sub, eval_X, eval_C, zero_sub]
    exact pow_ne_zero _ (neg_ne_zero.2 fun h => S.zero_not_mem (h ▸ ht))
  exact this hrec.symm

end Scheme

end ClauseScheme

#print axioms ClauseScheme.Scheme.lemma1
#print axioms ClauseScheme.Scheme.lemma1'
#print axioms ClauseScheme.Scheme.factR
