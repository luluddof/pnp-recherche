import Mathlib.Tactic
import Mathlib.Data.Nat.Log
import Mathlib.LinearAlgebra.Matrix.Rank
import Mathlib.LinearAlgebra.Matrix.NonsingularInverse
import Mathlib.RingTheory.Ideal.Quotient.Basic
import Mathlib.RingTheory.Localization.FractionRing

/-!
# Theorem A″ (direction-factored identities: a lower bound on the number of terms) — G19d

Source: `code/gen18/g18s/notes.md` §1 and §5 (statement `df_flattening`), `code/gen18/g18e/notes.md`
§2 (Theorem A′), `code/gen14/g14b/notes.md` (model DF).

Setting (after the flattening of the notes: `W` = the selected wires, each with a row-index end
and a column-index end; all other end bits are fixed).  `ρ i e` / `γ i e` is the group, in
direction `i`, of the row end / column end of the wire `e`.  A direction-`i` term is a product
over the groups `g` of direction `i` of tables `τ g y x` reading only the row bits `y e` with
`ρ i e = g` and the column bits `x e` with `γ i e = g` (`Local`).  Notation of the notes:
`e₁ g = cnt ρ g`, `e₂ g = cnt γ g`, `c = cExp ρ γ = ∑ g, min (e₁ g) (e₂ g)`, `w = card W`.

* `rank_term_le` : one term (with a scalar `κ`) has rank
  `≤ 2 ^ (#{e | σ (ρ e) = true} + #{e | σ (γ e) = false})` for every choice of sides `σ`.
* `rank_term_le_min` : one term has rank `≤ 2 ^ cExp ρ γ`.
* `df_flattening` : the reviewer's statement (g18s §5), verbatim; `df_flattening_smul` : the same
  with explicit scalars `κ k`.
* `df_flattening_min` : Theorem A″ as written in the notes, `2 ^ w ≤ ∑ k, 2 ^ c (dir k)`.
* `df_card_le_mul`, `df_log_bound`, `df_log_bound_sup`, `df_log_bound_abs` : the corollaries
  `2 ^ w ≤ N * 2 ^ C`, `w - C ≤ log₂ N` (`C ≥ c (dir k)` for all `k`, e.g. `C = max_i c_i`), and
  the form `log₂ N ≥ min_i ½ ∑_g |e₁ - e₂|` (`two_mul_cExp` : `2 c + ∑_g |e₁ - e₂| = 2 w`).
* `det_ne_zero_of_sub_one_mem`, `det_map_ne_zero_of_sub_one_mem`,
  `det_map_ne_zero_of_sub_one_mem_fractionRing` : the "2-adic" lemma.
* `df_flattening_mod` : A″ for an identity that holds modulo a proper ideal (`I = (2 ^ p)`).
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

/-- A matrix of the form `M y x = Fy y (ψ x) * Hx (φ y) x` factors through `P × Q`. -/
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

/-- A square matrix with nonzero determinant that is a sum of `N` matrices. -/
theorem card_le_sum_rank {n : Type*} [Fintype n] [DecidableEq n] {N : ℕ}
    (T : Fin N → Matrix n n F) (A : Matrix n n F) (hA : A = ∑ k, T k) (hdet : A.det ≠ 0) :
    Fintype.card n ≤ ∑ k, (T k).rank := by
  have hunit : IsUnit A := (Matrix.isUnit_iff_isUnit_det A).2 (isUnit_iff_ne_zero.2 hdet)
  rw [← Matrix.rank_of_isUnit A hunit, hA]
  exact rank_sum_le _ _

end Rank

/-! ### Counting wire ends -/

section Count

variable {W : Type*} [Fintype W] {G : Type*} [Fintype G] [DecidableEq G]

/-- Number of selected wires whose end (given by `f`) lies in the group `g`
(`e₁ g = cnt ρ g`, `e₂ g = cnt γ g` in the notes). -/
def cnt (f : W → G) (g : G) : ℕ := Fintype.card {e // f e = g}

/-- The exponent `c = ∑_g min (e₁ g) (e₂ g)` of the notes. -/
def cExp (ρ γ : W → G) : ℕ := ∑ g, min (cnt ρ g) (cnt γ g)

/-- The optimal sides: `true` (factor through the row bits) iff `e₁ g ≤ e₂ g`. -/
def optSide (ρ γ : W → G) (g : G) : Bool := decide (cnt ρ g ≤ cnt γ g)

/-- Number of selected wires whose end (given by `f`) lies in a group of side `b`, counted
group by group. -/
theorem card_side (f : W → G) (σ : G → Bool) (b : Bool) :
    Fintype.card {e // σ (f e) = b} = ∑ g, if σ g = b then cnt f g else 0 := by
  classical
  have h := Finset.card_eq_sum_card_fiberwise (f := f)
    (s := Finset.univ.filter fun e => σ (f e) = b) (t := Finset.univ) (by intro x _; simp)
  rw [Fintype.card_subtype, h]
  refine Finset.sum_congr rfl fun g _ => ?_
  unfold cnt
  rw [Fintype.card_subtype]
  split_ifs with hg
  · congr 1
    ext e
    simp only [Finset.mem_filter, Finset.mem_univ, true_and]
    constructor
    · rintro ⟨_, rfl⟩; rfl
    · rintro rfl; exact ⟨hg, rfl⟩
  · rw [Finset.card_eq_zero, Finset.filter_eq_empty_iff]
    intro e he hfe
    subst hfe
    exact hg (by simpa using he)

theorem sum_cnt (f : W → G) : ∑ g, cnt f g = Fintype.card W := by
  have h := card_side f (fun _ => true) true
  simp only [ite_true] at h
  rw [← h, Fintype.card_subtype]
  simp

/-- With the optimal sides, the exponent of `rank_term_le` is `∑ g, min (e₁ g) (e₂ g)`. -/
theorem card_optSide (ρ γ : W → G) :
    Fintype.card {e // optSide ρ γ (ρ e) = true} + Fintype.card {e // optSide ρ γ (γ e) = false}
      = cExp ρ γ := by
  rw [card_side, card_side, ← Finset.sum_add_distrib, cExp]
  refine Finset.sum_congr rfl fun g _ => ?_
  by_cases h : cnt ρ g ≤ cnt γ g
  · rw [optSide, decide_eq_true h]; simp [min_eq_left h]
  · rw [optSide, decide_eq_false h]; simp [min_eq_right (not_le.1 h).le]

/-- `w - c = ½ ∑_g |e₁ g - e₂ g|`, written without division or subtraction. -/
theorem two_mul_cExp (ρ γ : W → G) :
    (2 * cExp ρ γ : ℤ) + ∑ g, |(cnt ρ g : ℤ) - cnt γ g| = 2 * Fintype.card W := by
  have h : (2 * Fintype.card W : ℤ) = ∑ g, ((cnt ρ g : ℤ) + cnt γ g) := by
    rw [Finset.sum_add_distrib, ← Nat.cast_sum, ← Nat.cast_sum, sum_cnt, sum_cnt]; ring
  rw [h, cExp, Nat.cast_sum, Finset.mul_sum, ← Finset.sum_add_distrib]
  refine Finset.sum_congr rfl fun g _ => ?_
  rcases le_total (cnt ρ g) (cnt γ g) with h | h
  · rw [min_eq_left h, abs_of_nonpos (sub_nonpos.2 (by exact_mod_cast h))]; ring
  · rw [min_eq_right h, abs_of_nonneg (sub_nonneg.2 (by exact_mod_cast h))]; ring

end Count

/-! ### One term -/

section Term

variable {K : Type*} [Field K] {W : Type*} [Fintype W] [DecidableEq W] {G : Type*} [Fintype G]

/-- `t` reads only the row bits located in group `g` (w.r.t. `ρ`) and the column bits located
in group `g` (w.r.t. `γ`). -/
def Local (ρ γ : W → G) (g : G) (t : (W → Bool) → (W → Bool) → K) : Prop :=
  ∀ y y' x x', (∀ e, ρ e = g → y e = y' e) → (∀ e, γ e = g → x e = x' e) → t y x = t y' x'

/-- Extension by `false` of a partial assignment. -/
def ext {p : W → Prop} [DecidablePred p] (a : {e // p e} → Bool) : W → Bool :=
  fun e => if h : p e then a ⟨e, h⟩ else false

omit [Fintype W] [DecidableEq W] in
theorem ext_restrict {p : W → Prop} [DecidablePred p] (y : W → Bool) (e : W) (he : p e) :
    ext (p := p) (fun e' => y e'.1) e = y e := by
  simp [ext, he]

/-- **Rank of one direction-factored term.** For every choice of sides `σ`, the term factors
through the row bits of the wires whose row end lies in a `true` group together with the column
bits of the wires whose column end lies in a `false` group. -/
theorem rank_term_le (ρ γ : W → G) (σ : G → Bool) (κ : K)
    (τ : G → (W → Bool) → (W → Bool) → K) (hloc : ∀ g, Local ρ γ g (τ g))
    (M : Matrix (W → Bool) (W → Bool) K) (hM : ∀ y x, M y x = κ * ∏ g, τ g y x) :
    M.rank ≤ 2 ^ (Fintype.card {e // σ (ρ e) = true} + Fintype.card {e // σ (γ e) = false}) := by
  classical
  have key := rank_le_of_factor (F := K)
    (P := {e // σ (ρ e) = true} → Bool) (Q := {e // σ (γ e) = false} → Bool)
    (fun y e => y e.1) (fun x e => x e.1)
    (fun y q => κ * ∏ g, if σ g = false then τ g y (ext q) else 1)
    (fun p x => ∏ g, if σ g = true then τ g (ext p) x else 1)
    M (by
      intro y x
      rw [hM, mul_assoc, ← Finset.prod_mul_distrib]
      congr 1
      refine Finset.prod_congr rfl fun g _ => ?_
      cases hσ : σ g
      · simp only [ite_true, Bool.false_eq_true, ite_false, mul_one]
        refine hloc g _ _ _ _ (fun _ _ => rfl) fun e he => ?_
        exact (ext_restrict x e (by rw [he, hσ])).symm
      · simp only [ite_true, Bool.true_eq_false, ite_false, one_mul]
        refine hloc g _ _ _ _ (fun e he => ?_) fun _ _ => rfl
        exact (ext_restrict y e (by rw [he, hσ])).symm)
  refine key.trans (le_of_eq ?_)
  simp [pow_add]

/-- One direction-factored term has rank `≤ 2 ^ ∑ g, min (e₁ g) (e₂ g)`. -/
theorem rank_term_le_min [DecidableEq G] (ρ γ : W → G) (κ : K)
    (τ : G → (W → Bool) → (W → Bool) → K) (hloc : ∀ g, Local ρ γ g (τ g))
    (M : Matrix (W → Bool) (W → Bool) K) (hM : ∀ y x, M y x = κ * ∏ g, τ g y x) :
    M.rank ≤ 2 ^ cExp ρ γ := by
  have h := rank_term_le ρ γ (optSide ρ γ) κ τ hloc M hM
  rwa [card_optSide] at h

end Term

/-! ### Theorem A″ -/

section Main

variable {K : Type*} [Field K] {W : Type*} [Fintype W] [DecidableEq W]
  {ι : Type*} {G : ι → Type*} [∀ i, Fintype (G i)]

/-- The flattening bound with explicit scalars `κ k` (the `κ_k ∈ Frac(O)` of model DF). -/
theorem df_flattening_smul
    (ρ γ : ∀ i, W → G i) (σ : ∀ i, G i → Bool) {N : ℕ} (dir : Fin N → ι) (κ : Fin N → K)
    (τ : ∀ k : Fin N, G (dir k) → (W → Bool) → (W → Bool) → K)
    (hloc : ∀ k g y y' x x', (∀ e, ρ (dir k) e = g → y e = y' e) →
      (∀ e, γ (dir k) e = g → x e = x' e) → τ k g y x = τ k g y' x')
    (A : Matrix (W → Bool) (W → Bool) K)
    (hA : ∀ y x, A y x = ∑ k, κ k * ∏ g, τ k g y x)
    (hdet : A.det ≠ 0) :
    2 ^ Fintype.card W ≤
      ∑ k, 2 ^ (Fintype.card {e // σ (dir k) (ρ (dir k) e) = true} +
        Fintype.card {e // σ (dir k) (γ (dir k) e) = false}) := by
  have h := card_le_sum_rank (fun k => Matrix.of fun y x => κ k * ∏ g, τ k g y x) A
    (by ext y x; rw [hA, Matrix.sum_apply]; rfl) hdet
  have hc : Fintype.card (W → Bool) = 2 ^ Fintype.card W := by simp
  rw [← hc]
  refine h.trans (Finset.sum_le_sum fun k _ => ?_)
  exact rank_term_le _ _ _ (κ k) (τ k) (hloc k) _ (fun _ _ => rfl)

/-- **Reviewer's statement** (`df_flattening`, g18s notes §5), verbatim up to the unused
instance `[∀ i, DecidableEq (G i)]`. -/
theorem df_flattening
    (ρ γ : ∀ i, W → G i) (σ : ∀ i, G i → Bool) {N : ℕ} (dir : Fin N → ι)
    (τ : ∀ k : Fin N, G (dir k) → (W → Bool) → (W → Bool) → K)
    (hloc : ∀ k g y y' x x', (∀ e, ρ (dir k) e = g → y e = y' e) →
      (∀ e, γ (dir k) e = g → x e = x' e) → τ k g y x = τ k g y' x')
    (A : Matrix (W → Bool) (W → Bool) K)
    (hA : ∀ y x, A y x = ∑ k, ∏ g, τ k g y x)
    (hdet : A.det ≠ 0) :
    2 ^ Fintype.card W ≤
      ∑ k, 2 ^ (Fintype.card {e // σ (dir k) (ρ (dir k) e) = true} +
        Fintype.card {e // σ (dir k) (γ (dir k) e) = false}) :=
  df_flattening_smul ρ γ σ dir (fun _ => 1) τ hloc A (by simpa using hA) hdet

variable [∀ i, DecidableEq (G i)]

/-- **Theorem A″** (g18s notes §1): `2 ^ w ≤ ∑ k, 2 ^ c (dir k)` with
`c i = ∑ g, min (e₁ i g) (e₂ i g)`. -/
theorem df_flattening_min
    (ρ γ : ∀ i, W → G i) {N : ℕ} (dir : Fin N → ι) (κ : Fin N → K)
    (τ : ∀ k : Fin N, G (dir k) → (W → Bool) → (W → Bool) → K)
    (hloc : ∀ k g y y' x x', (∀ e, ρ (dir k) e = g → y e = y' e) →
      (∀ e, γ (dir k) e = g → x e = x' e) → τ k g y x = τ k g y' x')
    (A : Matrix (W → Bool) (W → Bool) K)
    (hA : ∀ y x, A y x = ∑ k, κ k * ∏ g, τ k g y x)
    (hdet : A.det ≠ 0) :
    2 ^ Fintype.card W ≤ ∑ k, 2 ^ cExp (ρ (dir k)) (γ (dir k)) := by
  have h := df_flattening_smul ρ γ (fun i => optSide (ρ i) (γ i)) dir κ τ hloc A hA hdet
  refine h.trans (le_of_eq (Finset.sum_congr rfl fun k _ => ?_))
  exact congrArg (fun n => 2 ^ n) (card_optSide _ _)

/-! ### Corollaries -/

theorem le_mul_of_le_sum {N w C : ℕ} (c : Fin N → ℕ) (h : 2 ^ w ≤ ∑ k, 2 ^ c k)
    (hC : ∀ k, c k ≤ C) : 2 ^ w ≤ N * 2 ^ C := by
  refine h.trans ?_
  calc ∑ k, 2 ^ c k ≤ ∑ _k : Fin N, 2 ^ C :=
        Finset.sum_le_sum fun k _ => Nat.pow_le_pow_right (by norm_num) (hC k)
    _ = N * 2 ^ C := by simp

theorem sub_le_log {N w C : ℕ} (h : 2 ^ w ≤ N * 2 ^ C) : w - C ≤ Nat.log 2 N := by
  have hN : N ≠ 0 := by
    rintro rfl
    simp at h
  rw [Nat.le_log_iff_pow_le (by norm_num) hN]
  by_cases hw : C ≤ w
  · have h' : 2 ^ (w - C) * 2 ^ C ≤ N * 2 ^ C := by
      rw [← pow_add, Nat.sub_add_cancel hw]; exact h
    exact Nat.le_of_mul_le_mul_right h' (by positivity)
  · rw [Nat.sub_eq_zero_of_le (not_le.1 hw).le]
    rw [pow_zero]
    exact Nat.pos_of_ne_zero hN

/-- `2 ^ w ≤ N * 2 ^ C` whenever `c (dir k) ≤ C` for every term (e.g. `C = max_i c_i`). -/
theorem df_card_le_mul
    (ρ γ : ∀ i, W → G i) {N : ℕ} (dir : Fin N → ι) (κ : Fin N → K)
    (τ : ∀ k : Fin N, G (dir k) → (W → Bool) → (W → Bool) → K)
    (hloc : ∀ k g y y' x x', (∀ e, ρ (dir k) e = g → y e = y' e) →
      (∀ e, γ (dir k) e = g → x e = x' e) → τ k g y x = τ k g y' x')
    (A : Matrix (W → Bool) (W → Bool) K)
    (hA : ∀ y x, A y x = ∑ k, κ k * ∏ g, τ k g y x)
    (hdet : A.det ≠ 0) (C : ℕ) (hC : ∀ k, cExp (ρ (dir k)) (γ (dir k)) ≤ C) :
    2 ^ Fintype.card W ≤ N * 2 ^ C :=
  le_mul_of_le_sum _ (df_flattening_min ρ γ dir κ τ hloc A hA hdet) hC

/-- `log₂ N ≥ w - C` whenever `c (dir k) ≤ C` for every term. -/
theorem df_log_bound
    (ρ γ : ∀ i, W → G i) {N : ℕ} (dir : Fin N → ι) (κ : Fin N → K)
    (τ : ∀ k : Fin N, G (dir k) → (W → Bool) → (W → Bool) → K)
    (hloc : ∀ k g y y' x x', (∀ e, ρ (dir k) e = g → y e = y' e) →
      (∀ e, γ (dir k) e = g → x e = x' e) → τ k g y x = τ k g y' x')
    (A : Matrix (W → Bool) (W → Bool) K)
    (hA : ∀ y x, A y x = ∑ k, κ k * ∏ g, τ k g y x)
    (hdet : A.det ≠ 0) (C : ℕ) (hC : ∀ k, cExp (ρ (dir k)) (γ (dir k)) ≤ C) :
    Fintype.card W - C ≤ Nat.log 2 N :=
  sub_le_log (df_card_le_mul ρ γ dir κ τ hloc A hA hdet C hC)

/-- `log₂ N ≥ w - max_i c_i` (finitely many directions). -/
theorem df_log_bound_sup [Fintype ι]
    (ρ γ : ∀ i, W → G i) {N : ℕ} (dir : Fin N → ι) (κ : Fin N → K)
    (τ : ∀ k : Fin N, G (dir k) → (W → Bool) → (W → Bool) → K)
    (hloc : ∀ k g y y' x x', (∀ e, ρ (dir k) e = g → y e = y' e) →
      (∀ e, γ (dir k) e = g → x e = x' e) → τ k g y x = τ k g y' x')
    (A : Matrix (W → Bool) (W → Bool) K)
    (hA : ∀ y x, A y x = ∑ k, κ k * ∏ g, τ k g y x)
    (hdet : A.det ≠ 0) :
    Fintype.card W - Finset.univ.sup (fun i => cExp (ρ i) (γ i)) ≤ Nat.log 2 N :=
  df_log_bound ρ γ dir κ τ hloc A hA hdet _
    (fun k => Finset.le_sup (f := fun i => cExp (ρ i) (γ i)) (Finset.mem_univ (dir k)))

/-- `log₂ N ≥ min_i ½ ∑_g |e₁ i g - e₂ i g|`: if `2 D ≤ ∑_g |e₁ - e₂|` in the direction of
every term, then `D ≤ log₂ N`. -/
theorem df_log_bound_abs
    (ρ γ : ∀ i, W → G i) {N : ℕ} (dir : Fin N → ι) (κ : Fin N → K)
    (τ : ∀ k : Fin N, G (dir k) → (W → Bool) → (W → Bool) → K)
    (hloc : ∀ k g y y' x x', (∀ e, ρ (dir k) e = g → y e = y' e) →
      (∀ e, γ (dir k) e = g → x e = x' e) → τ k g y x = τ k g y' x')
    (A : Matrix (W → Bool) (W → Bool) K)
    (hA : ∀ y x, A y x = ∑ k, κ k * ∏ g, τ k g y x)
    (hdet : A.det ≠ 0) (D : ℕ)
    (hD : ∀ k, (2 * D : ℤ) ≤ ∑ g, |(cnt (ρ (dir k)) g : ℤ) - cnt (γ (dir k)) g|) :
    D ≤ Nat.log 2 N := by
  have hC : ∀ k, cExp (ρ (dir k)) (γ (dir k)) ≤ Fintype.card W - D ∧ D ≤ Fintype.card W := by
    intro k
    have h1 := two_mul_cExp (ρ (dir k)) (γ (dir k))
    have h2 := hD k
    omega
  have h := df_card_le_mul ρ γ dir κ τ hloc A hA hdet _ (fun k => (hC k).1)
  have hN : 0 < N := by
    rcases Nat.eq_zero_or_pos N with rfl | hN
    · simp at h
    · exact hN
  have h3 := sub_le_log h
  have h4 := (hC ⟨0, hN⟩).2
  omega

end Main

/-! ### The "2-adic" lemma -/

section Modular

variable {R : Type*} [CommRing R] {n : Type*} [Fintype n] [DecidableEq n]

/-- A square matrix congruent to `1` entrywise modulo a proper ideal has nonzero determinant. -/
theorem det_ne_zero_of_sub_one_mem (I : Ideal R) (hI : I ≠ ⊤) (B : Matrix n n R)
    (h : ∀ i j, B i j - (1 : Matrix n n R) i j ∈ I) : B.det ≠ 0 := by
  intro h0
  have hmap : (Ideal.Quotient.mk I).mapMatrix B = 1 := by
    ext i j
    have hij := (Ideal.Quotient.eq (I := I)).2 (h i j)
    simp only [RingHom.mapMatrix_apply, Matrix.map_apply]
    rw [hij]
    by_cases e : i = j <;> simp [Matrix.one_apply, e]
  have h1 : Ideal.Quotient.mk I B.det = 1 := by
    rw [RingHom.map_det, hmap, Matrix.det_one]
  rw [h0, map_zero] at h1
  apply hI
  rw [Ideal.eq_top_iff_one]
  have h2 : Ideal.Quotient.mk I (1 : R) = 0 := by rw [map_one]; exact h1.symm
  exact Ideal.Quotient.eq_zero_iff_mem.1 h2

/-- The same after an injective base change (e.g. to the fraction field). -/
theorem det_map_ne_zero_of_sub_one_mem {K : Type*} [CommRing K] [Algebra R K]
    (hinj : Function.Injective (algebraMap R K)) (I : Ideal R) (hI : I ≠ ⊤) (B : Matrix n n R)
    (h : ∀ i j, B i j - (1 : Matrix n n R) i j ∈ I) : (B.map (algebraMap R K)).det ≠ 0 := by
  have h0 := det_ne_zero_of_sub_one_mem I hI B h
  have e : (B.map (algebraMap R K)).det = algebraMap R K B.det := by
    rw [RingHom.map_det]; rfl
  rw [e]
  exact fun hz => h0 (hinj (hz.trans (map_zero _).symm))

/-- **2-adic lemma** (g18s notes §5): over an integral domain `R` with a proper ideal `I`, a
square matrix `B ≡ 1 mod I` entrywise has `det ≠ 0` in `FractionRing R`. -/
theorem det_map_ne_zero_of_sub_one_mem_fractionRing [IsDomain R] (I : Ideal R) (hI : I ≠ ⊤)
    (B : Matrix n n R) (h : ∀ i j, B i j - (1 : Matrix n n R) i j ∈ I) :
    (B.map (algebraMap R (FractionRing R))).det ≠ 0 :=
  det_map_ne_zero_of_sub_one_mem (IsFractionRing.injective R (FractionRing R)) I hI B h

end Modular

/-! ### A″ for an identity modulo a proper ideal -/

section Mod

variable {R : Type*} [CommRing R] {K : Type*} [Field K] [Algebra R K]
  {W : Type*} [Fintype W] [DecidableEq W]
  {ι : Type*} {G : ι → Type*} [∀ i, Fintype (G i)] [∀ i, DecidableEq (G i)]

/-- **Theorem A″, DF form.** `R = O`, `I = (2 ^ p)` (any proper ideal), `K = Frac(O)` (any field
into which `R` embeds).  If `∑ k, κ k * ∏ g, τ k g y x = [y = x] + E y x` with `E y x ∈ I`, then
`2 ^ w ≤ ∑ k, 2 ^ c (dir k)`. -/
theorem df_flattening_mod (hinj : Function.Injective (algebraMap R K)) (I : Ideal R) (hI : I ≠ ⊤)
    (ρ γ : ∀ i, W → G i) {N : ℕ} (dir : Fin N → ι) (κ : Fin N → K)
    (τ : ∀ k : Fin N, G (dir k) → (W → Bool) → (W → Bool) → K)
    (hloc : ∀ k g y y' x x', (∀ e, ρ (dir k) e = g → y e = y' e) →
      (∀ e, γ (dir k) e = g → x e = x' e) → τ k g y x = τ k g y' x')
    (E : Matrix (W → Bool) (W → Bool) R) (hE : ∀ y x, E y x ∈ I)
    (hid : ∀ y x, ∑ k, κ k * ∏ g, τ k g y x = (if y = x then 1 else 0) + algebraMap R K (E y x)) :
    2 ^ Fintype.card W ≤ ∑ k, 2 ^ cExp (ρ (dir k)) (γ (dir k)) := by
  have hdet := det_map_ne_zero_of_sub_one_mem hinj I hI (1 + E) (fun y x => by simpa using hE y x)
  refine df_flattening_min ρ γ dir κ τ hloc _ (fun y x => ?_) hdet
  rw [hid, Matrix.map_apply, Matrix.add_apply, map_add, Matrix.one_apply]
  split_ifs <;> simp

end Mod

end DFRank

#print axioms DFRank.rank_term_le
#print axioms DFRank.rank_term_le_min
#print axioms DFRank.df_flattening
#print axioms DFRank.df_flattening_smul
#print axioms DFRank.df_flattening_min
#print axioms DFRank.df_card_le_mul
#print axioms DFRank.df_log_bound
#print axioms DFRank.df_log_bound_sup
#print axioms DFRank.df_log_bound_abs
#print axioms DFRank.two_mul_cExp
#print axioms DFRank.det_ne_zero_of_sub_one_mem
#print axioms DFRank.det_map_ne_zero_of_sub_one_mem_fractionRing
#print axioms DFRank.df_flattening_mod
