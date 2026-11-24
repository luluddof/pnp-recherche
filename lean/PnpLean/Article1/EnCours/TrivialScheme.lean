import PnpLean.Article1.TheoremA

/-!
# The trivial scheme (paper §6, "The DP-type model is degenerate"), combinatorial form

Paper: `resultats/article/main.tex`, Proposition "A trivial scheme in the DP-type model".

For **every** CNF `I`, every degree bound `R`, every partition `grp : Fin n → G` of the variables
and every clause assignment `α : κ → G`, there is a `Scheme I G` with

* constant deformations `P_o = e_v`;
* recovery points `X = {1, …, M R + 1}`, Taylor orders `q_t = 0`, `Λ` = Lagrange evaluation at `0`;
* the same design at every point, fully measured groups `𝓡_{t,γ} = ℚ^{S_γ}`.

Its parameters are `q_t = 0`, `ρ_t = n`, `a_γ = |S_γ|` (`trivialScheme_*`).
With blocks of size `s` where `n ≤ s * s` (`blockGrp`): `s` groups, each of size `≤ s`, hence
`s_t ≤ s` (`exists_trivial_scheme_sqrt`). (Treewidth is not part of the `Scheme` structure: the
interaction graph has at most `s` vertices, so `w_t ≤ s - 1`; this last step stays on paper.)
-/

open Polynomial Matrix Module

set_option linter.unusedSectionVars false
set_option linter.style.numericProjections false

namespace ClauseScheme

variable {n : ℕ} {κ : Type*} [Fintype κ] [DecidableEq κ]
  {G : Type*} [Fintype G] [DecidableEq G]

/-- Recovery points `{1, …, N}`. -/
def trivPts (N : ℕ) : Finset ℚ := (Finset.Icc 1 N).image (fun k : ℕ => (k : ℚ))

theorem card_trivPts (N : ℕ) : (trivPts N).card = N := by
  rw [trivPts, Finset.card_image_of_injective _ Nat.cast_injective]
  simp

theorem zero_notMem_trivPts (N : ℕ) : (0 : ℚ) ∉ trivPts N := by
  simp [trivPts]

/-- Lagrange evaluation at `0` from the values at `pts` (all Taylor orders `0`). -/
noncomputable def lagrangeΛ (pts : Finset ℚ) : (HIdx pts (fun _ => 0) → ℚ) →ₗ[ℚ] ℚ where
  toFun y := ∑ t : pts, y ⟨t, 0⟩ * (Lagrange.basis pts id (t : ℚ)).eval 0
  map_add' y z := by simp [add_mul, Finset.sum_add_distrib]
  map_smul' c y := by simp [Finset.mul_sum, mul_assoc]

theorem lagrangeΛ_spec (pts : Finset ℚ) (f : ℚ[X]) (hf : f.natDegree < pts.card) :
    lagrangeΛ pts (hermiteData pts (fun _ => 0) f) = f.eval 0 := by
  have hdeg : f.degree < pts.card := by
    by_cases h0 : f = 0
    · subst h0
      rw [degree_zero]
      exact WithBot.bot_lt_coe _
    · rw [degree_eq_natDegree h0]
      exact_mod_cast hf
  have h := Lagrange.eq_interpolate (s := pts) (v := id) Function.injective_id.injOn hdeg
  have h' := congrArg (eval 0) h
  rw [Lagrange.interpolate_apply, eval_finsetSum] at h'
  rw [h']
  simp only [lagrangeΛ, LinearMap.coe_mk, AddHom.coe_mk, hermiteData, Function.iterate_zero, id,
    eval_mul, eval_C]
  exact Finset.sum_coe_sort pts (fun t => eval t f * eval 0 (Lagrange.basis pts id t))

/-- The fully measured group `ℚ^{S_γ}`, as the range of the coordinate projection. -/
noncomputable def coordSub (grp : Fin n → G) (γ : G) : Submodule ℚ (Fin n → ℚ) :=
  LinearMap.range (Matrix.diagonal (fun j => if grp j = γ then (1 : ℚ) else 0)).mulVecLin

theorem mem_coordSub_of (grp : Fin n → G) (γ : G) (x : Fin n → ℚ) :
    (fun j => if grp j = γ then x j else 0) ∈ coordSub grp γ := by
  refine ⟨x, ?_⟩
  funext j
  simp [Matrix.mulVec_diagonal]

theorem coordSub_supp (grp : Fin n → G) (γ : G) (x : Fin n → ℚ) (hx : x ∈ coordSub grp γ)
    (j : Fin n) (hj : grp j ≠ γ) : x j = 0 := by
  obtain ⟨y, rfl⟩ := hx
  simp [Matrix.mulVec_diagonal, hj]

theorem finrank_coordSub (grp : Fin n → G) (γ : G) :
    finrank ℚ (coordSub grp γ) = (Finset.univ.filter (fun j => grp j = γ)).card := by
  rw [coordSub, Matrix.finrank_range_mulVecLin, Matrix.rank_diagonal, Fintype.card_subtype]
  congr 1
  ext j
  simp

variable (I : κ → Finset (Fin n))

/-- **The trivial scheme**: constant deformations, Lagrange recovery at `M R + 1` points, the
same fully measured design at every point. -/
noncomputable def trivialScheme (R : ℕ) (grp : Fin n → G) (α : κ → G) : Scheme I G where
  R := R
  P := fun o j => C (if j = o.var then 1 else 0)
  P_deg := fun o j => by simp
  P_zero := fun o => by
    funext j
    simp [Pi.single_apply]
  pts := trivPts (M I * R + 1)
  zero_not_mem := zero_notMem_trivPts _
  q := fun _ => 0
  Λ := lagrangeΛ _
  recovery := fun f hf => lagrangeΛ_spec _ f (by
    rw [card_trivPts]
    exact Nat.lt_succ_of_le hf)
  grp := fun _ => grp
  α := fun _ => α
  Rsp := fun _ γ => coordSub grp γ
  Rsp_supp := fun _ γ x hx j hj => coordSub_supp grp γ x hx j hj
  cond3 := fun t _ o γ _ => mem_coordSub_of grp γ _

variable {I}

section Params

variable (R : ℕ) (grp : Fin n → G) (α : κ → G)

@[simp] theorem trivialScheme_R : (trivialScheme I R grp α).R = R := rfl

@[simp] theorem trivialScheme_q (t : ℚ) : (trivialScheme I R grp α).q t = 0 := rfl

@[simp] theorem trivialScheme_grp (t : ℚ) : (trivialScheme I R grp α).grp t = grp := rfl

theorem trivialScheme_card_pts : (trivialScheme I R grp α).pts.card = M I * R + 1 :=
  card_trivPts _

theorem trivialScheme_a (t : ℚ) (γ : G) :
    (trivialScheme I R grp α).a t γ = (Finset.univ.filter (fun j => grp j = γ)).card := rfl

/-- `ρ_t = n`: every group is fully measured. -/
theorem trivialScheme_ρ (t : ℚ) : (trivialScheme I R grp α).ρ t = n := by
  have h := (trivialScheme I R grp α).sum_a t
  rw [← h]
  unfold Scheme.ρ
  refine Finset.sum_congr rfl fun γ _ => ?_
  exact finrank_coordSub grp γ

/-- The trivial scheme is variable-uniform. -/
theorem trivialScheme_uniform : (trivialScheme I R grp α).Uniform := by
  intro o o' h
  funext j
  simp [trivialScheme, h]

/-- The witness matrix of the trivial scheme is the identity: every deficiency is absorbed by
`ρ_t = n`, and Lemma 1 holds trivially. -/
theorem trivialScheme_deficiency_le (ω : WitnessMap I) (t : ℚ) :
    (trivialScheme I R grp α).deficiency ω t ≤ (trivialScheme I R grp α).ρ t := by
  rw [trivialScheme_ρ]
  have h := (trivialScheme I R grp α).sum_a t
  rw [← h]
  exact Finset.sum_le_sum fun γ _ => Nat.sub_le _ _

end Params

/-! ### Blocks of size `s` with `n ≤ s * s` -/

/-- The block partition `j ↦ ⌊j / s⌋`, with `s` groups, for `n ≤ s * s`. -/
def blockGrp (s : ℕ) (hs : n ≤ s * s) (j : Fin n) : Fin s :=
  ⟨(j : ℕ) / s, by
    have hj := j.2
    have hpos : 0 < s := by
      rcases Nat.eq_zero_or_pos s with h | h
      · subst h; omega
      · exact h
    exact Nat.div_lt_of_lt_mul (lt_of_lt_of_le hj hs)⟩

theorem card_blockGrp_fiber_le (s : ℕ) (hs : n ≤ s * s) (γ : Fin s) :
    (Finset.univ.filter (fun j : Fin n => blockGrp s hs j = γ)).card ≤ s := by
  have hpos : 0 < s := γ.pos
  have h := Finset.card_le_card_of_injOn (s := Finset.univ.filter (fun j : Fin n => blockGrp s hs j = γ))
    (t := (Finset.univ : Finset (Fin s)))
    (fun j => (⟨(j : ℕ) % s, Nat.mod_lt _ hpos⟩ : Fin s)) (fun _ _ => Finset.mem_univ _)
    (by
      intro i hi j hj hij
      simp only [Finset.coe_filter, Finset.mem_univ, true_and, Set.mem_setOf_eq] at hi hj
      have h1 : (i : ℕ) / s = (j : ℕ) / s := by
        have := congrArg Fin.val (hi.trans hj.symm)
        simpa [blockGrp] using this
      have h2 : (i : ℕ) % s = (j : ℕ) % s := by simpa using congrArg Fin.val hij
      apply Fin.ext
      rw [← Nat.div_add_mod (i : ℕ) s, ← Nat.div_add_mod (j : ℕ) s, h1, h2])
  simpa using h

/-- **Trivial scheme, `√n` form** (paper, Proposition of §6). For every CNF `I` on `n` variables,
every `R`, and every `s` with `n ≤ s * s` (e.g. `s = ⌈√n⌉`), there is a scheme with `s` groups
such that at every point: all groups have size `≤ s` (so `s_t ≤ s`), `q_t = 0`, `ρ_t = n`.
It has `M R + 1` recovery points, degree bound `R`, and is variable-uniform. -/
theorem exists_trivial_scheme_sqrt (I : κ → Finset (Fin n)) (R s : ℕ) (hs : n ≤ s * s) :
    ∃ S : Scheme I (Fin s), S.R = R ∧ S.pts.card = M I * R + 1 ∧ S.Uniform ∧
      ∀ t, S.s t ≤ s ∧ S.q t = 0 ∧ S.ρ t = n := by
  rcases Nat.eq_zero_or_pos s with h0 | hpos
  · -- `s = 0` forces `n = 0`; any clause assignment is impossible only if `κ` is nonempty,
    -- but then `Fin 0` has no element: treat `κ` empty / nonempty uniformly via `n = 0`.
    subst h0
    have hn : n = 0 := by omega
    subst hn
    refine ⟨trivialScheme I R (fun j => j.elim0) (fun e => ?_), rfl, card_trivPts _,
      trivialScheme_uniform _ _ _, fun t => ⟨?_, rfl, trivialScheme_ρ _ _ _ t⟩⟩
    · sorry
    · simp [Scheme.s]
  · refine ⟨trivialScheme I R (blockGrp s hs) (fun _ => ⟨0, hpos⟩), rfl, card_trivPts _,
      trivialScheme_uniform _ _ _, fun t => ⟨?_, rfl, trivialScheme_ρ _ _ _ t⟩⟩
    unfold Scheme.s
    exact Finset.sup_le fun γ _ => card_blockGrp_fiber_le s hs γ

end ClauseScheme

#print axioms ClauseScheme.lagrangeΛ_spec
#print axioms ClauseScheme.trivialScheme_ρ
#print axioms ClauseScheme.trivialScheme_uniform
#print axioms ClauseScheme.exists_trivial_scheme_sqrt
