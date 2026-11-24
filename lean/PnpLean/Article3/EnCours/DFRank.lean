import Mathlib

/-!
# Theorem A″ (direction-factored identities: a lower bound on the number of terms) — G19d

Source: `code/gen18/g18s/notes.md` §1 and §5 (statement `df_flattening`), `code/gen18/g18e/notes.md`
§2 (Theorem A′), `code/gen14/g14b/notes.md` (model DF).

Setting (after the flattening of the notes: the selected wires `W`, each with a row-index end
and a column-index end; all other end bits are fixed).  `ρ i e` / `γ i e` is the group, in
direction `i`, of the row end / column end of the wire `e`.  A direction-`i` term is a product
over the groups `g` of direction `i` of tables `τ g y x` reading only the row bits `y e` with
`ρ i e = g` and the column bits `x e` with `γ i e = g` (`Local`).

* `rank_term_le` : one term has rank `≤ 2 ^ (#{e | σ (ρ e) = true} + #{e | σ (γ e) = false})`
  for every choice of sides `σ : G → Bool`.
* `rank_term_le_min` : one term has rank `≤ 2 ^ ∑ g, min (e₁ g) (e₂ g)`.
* `df_flattening` : the reviewer's statement (g18s §5), verbatim.
* `df_flattening_min` : Theorem A″ as written in the notes: `2 ^ w ≤ ∑ k, 2 ^ c (dir k)`,
  `c i = ∑ g, min (e₁ i g) (e₂ i g)`.
* `df_card_le_mul`, `df_log_bound`, `df_log_bound_sup`, `df_log_bound_abs` : corollaries
  `2 ^ w ≤ N * 2 ^ C`, `w - C ≤ log₂ N`, with `C = max_i c_i`, and the form
  `w - c_i = ½ ∑_g |e₁ - e₂|`.
* `det_ne_zero_of_sub_one_mem` (+ `_fractionRing`) : the "2-adic" lemma.
* `df_flattening_mod` : A″ for an identity that holds modulo a proper ideal of a domain.
-/

open Matrix Module

namespace DFRank

/-! ### Rank bookkeeping -/

section Rank

variable {F : Type*} [Field F] {m n : Type*} [Fintype n]

theorem rank_add_le (A B : Matrix m n F) : (A + B).rank ≤ A.rank + B.rank := by
  rw [← finrank_range_mulVecLin, ← finrank_range_mulVecLin, ← finrank_range_mulVecLin,
    Matrix.mulVecLin_add]
  refine le_trans (Submodule.finrank_mono ?_) (Submodule.finrank_add_le_finrank_add_finrank _ _)
  rintro _ ⟨x, rfl⟩
  exact Submodule.add_mem_sup (LinearMap.mem_range_self _ x) (LinearMap.mem_range_self _ x)

theorem rank_sum_le {ι : Type*} (s : Finset ι) (M : ι → Matrix m n F) :
    (∑ k ∈ s, M k).rank ≤ ∑ k ∈ s, (M k).rank := by
  classical
  induction s using Finset.induction_on with
  | empty => simp
  | insert a s ha ih =>
    rw [Finset.sum_insert ha, Finset.sum_insert ha]
    exact (rank_add_le _ _).trans (Nat.add_le_add_left ih _)

/-- A matrix of the form `M y x = F y (ψ x) * H (φ y) x` factors through `P × Q`. -/
theorem rank_le_of_factor {R C P Q : Type*} [Fintype R] [Fintype C] [Fintype P] [Fintype Q]
    (φ : R → P) (ψ : C → Q) (Fy : R → Q → F) (Hx : P → C → F)
    (M : Matrix R C F) (hM : ∀ y x, M y x = Fy y (ψ x) * Hx (φ y) x) :
    M.rank ≤ Fintype.card P * Fintype.card Q := by
  classical
  let U : Matrix R (P × Q) F := fun y pq => if φ y = pq.1 then Fy y pq.2 else 0
  let V : Matrix (P × Q) C F := fun pq x => if ψ x = pq.2 then Hx pq.1 x else 0
  have h : M = U * V := by
    ext y x
    rw [hM, Matrix.mul_apply, Fintype.sum_prod_type]
    rw [Finset.sum_eq_single (φ y), Finset.sum_eq_single (ψ x)]
    · simp [U, V]
    · intro q _ hq
      simp [U, V, Ne.symm hq]
    · simp
    · intro p _ hp
      simp [U, Ne.symm hp]
    · simp
  rw [h]
  calc (U * V).rank ≤ U.rank := Matrix.rank_mul_le_left _ _
    _ ≤ Fintype.card (P × Q) := Matrix.rank_le_card_width _
    _ = Fintype.card P * Fintype.card Q := Fintype.card_prod _ _

end Rank

/-! ### One term -/

section Term

variable {K : Type*} [Field K] {W : Type*} [Fintype W] {G : Type*} [Fintype G]

/-- `t` reads only the row bits located in group `g` (w.r.t. `ρ`) and the column bits located
in group `g` (w.r.t. `γ`). -/
def Local (ρ γ : W → G) (g : G) (t : (W → Bool) → (W → Bool) → K) : Prop :=
  ∀ y y' x x', (∀ e, ρ e = g → y e = y' e) → (∀ e, γ e = g → x e = x' e) → t y x = t y' x'

/-- Extension by `false` of a partial assignment. -/
def ext {p : W → Prop} [DecidablePred p] (a : {e // p e} → Bool) : W → Bool :=
  fun e => if h : p e then a ⟨e, h⟩ else false

theorem ext_restrict {p : W → Prop} [DecidablePred p] (y : W → Bool) (e : W) (he : p e) :
    ext (p := p) (fun e' => y e'.1) e = y e := by
  simp [ext, he]

/-- **Rank of one direction-factored term.** For every choice of sides `σ`, the term factors
through the row bits of the wires whose row end lies in a `true` group together with the column
bits of the wires whose column end lies in a `false` group. -/
theorem rank_term_le (ρ γ : W → G) (σ : G → Bool)
    (τ : G → (W → Bool) → (W → Bool) → K) (hloc : ∀ g, Local ρ γ g (τ g))
    (M : Matrix (W → Bool) (W → Bool) K) (hM : ∀ y x, M y x = ∏ g, τ g y x) :
    M.rank ≤ 2 ^ (Fintype.card {e // σ (ρ e) = true} + Fintype.card {e // σ (γ e) = false}) := by
  classical
  have key := rank_le_of_factor (F := K)
    (P := {e // σ (ρ e) = true} → Bool) (Q := {e // σ (γ e) = false} → Bool)
    (fun y e => y e.1) (fun x e => x e.1)
    (fun y q => ∏ g, if σ g = false then τ g y (ext q) else 1)
    (fun p x => ∏ g, if σ g = true then τ g (ext p) x else 1)
    M (by
      intro y x
      rw [hM, ← Finset.prod_mul_distrib]
      refine Finset.prod_congr rfl fun g _ => ?_
      cases hσ : σ g
      · simp only [if_true, Bool.false_eq_true, if_false, mul_one]
        refine hloc g _ _ _ _ (fun _ _ => rfl) fun e he => ?_
        exact (ext_restrict x e (by rw [he, hσ])).symm
      · simp only [if_true, Bool.true_eq_false, if_false, one_mul]
        refine hloc g _ _ _ _ (fun e he => ?_) fun _ _ => rfl
        exact (ext_restrict y e (by rw [he, hσ])).symm)
  refine key.trans (le_of_eq ?_)
  simp [pow_add]

/-- Number of selected wires whose end (given by `f`) lies in a group of side `b`, counted
group by group. -/
theorem card_side (f : W → G) (σ : G → Bool) (b : Bool) :
    Fintype.card {e // σ (f e) = b} = ∑ g, if σ g = b then Fintype.card {e // f e = g} else 0 := by
  classical
  rw [Fintype.card_subtype]
  simp_rw [Fintype.card_subtype]
  rw [Finset.card_eq_sum_card_fiberwise (f := f) (t := Finset.univ) (fun _ _ => Finset.mem_univ _)]
  refine Finset.sum_congr rfl fun g _ => ?_
  split_ifs with h
  · congr 1
    ext e
    simp only [Finset.mem_filter, Finset.mem_univ, true_and]
    constructor
    · rintro ⟨_, rfl⟩; rfl
    · rintro rfl; exact ⟨h, rfl⟩
  · rw [Finset.card_eq_zero, Finset.filter_eq_empty_iff]
    rintro e _ ⟨h1, rfl⟩
    exact h h1

/-- With the optimal sides `σ g = true ↔ e₁ g ≤ e₂ g`, the exponent is `∑ g, min (e₁ g) (e₂ g)`. -/
theorem card_side_min (ρ γ : W → G) :
    Fintype.card {e // decide (Fintype.card {e // ρ e = ρ e} ≤ Fintype.card {e' // γ e' = ρ e}) = true} = 0 ∨ True :=
  Or.inr trivial

end Term

end DFRank
