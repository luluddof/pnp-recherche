import Mathlib.Tactic
import Mathlib.Data.Fintype.Card
import Mathlib.Data.Finset.Card
import Mathlib.Data.Nat.Choose.Basic
import Mathlib.Analysis.Real.Sqrt

/-!
# Monomial (ST intercalate) model: key lemma and counting inequality — G23e

Sources: `code/gen22/g22b/notes.md` §3 (Theorem 3 (a), (c), ST corollary) and
`code/gen23/g23b/notes.md` ("Lean candidates").

Setup: `col : Fin N → Fin N → κ` (colouring of the `N × N` array), `pal : κ → Bool` (palette).
Hypotheses: `hrow`/`hcol` (every colour class is a partial permutation) and `hint` (the
sign-free bichromatic intercalate condition: if the diagonal `(x,y),(v,u)` of a genuine `2 × 2`
submatrix is monochromatic, so is the antidiagonal `(x,u),(v,y)`, with the other palette).

* `key_lemma`  : `#(cells of colour j).choose 2 ≤ #{used c | pal c ≠ pal j}`
  (two-to-one map `((a,b),(a',b')) ↦ col a b'` on ordered pairs of distinct `j`-cells).
* `counting`   : `N² ≤ q * gstar N p + p * gstar N q`, `p = #used T (pal = true)`,
  `q = #used S (pal = false)`.
* `gstar_le_sqrt` : `gstar N t ≤ Nat.sqrt (2t) + 1`; `gstar_le_real_sqrt` : `≤ √(2t) + 1`.
* `tradeoff_real` : `q ≤ p → N² ≤ 2 p (√(2q) + 1)` (reals).
-/

set_option linter.unusedSectionVars false

namespace Intercalate

open Finset

variable {N : ℕ} {κ : Type*} [Fintype κ] [DecidableEq κ]

/-- Cells of colour `j`. -/
def cells (col : Fin N → Fin N → κ) (j : κ) : Finset (Fin N × Fin N) :=
  univ.filter fun p => col p.1 p.2 = j

/-- Used colours `c` whose palette differs from that of `j`. -/
def otherUsed (col : Fin N → Fin N → κ) (pal : κ → Bool) (j : κ) : Finset κ :=
  univ.filter fun c => (∃ a b, col a b = c) ∧ pal c ≠ pal j

/-- Used colours of palette `β`. -/
def usedPal (col : Fin N → Fin N → κ) (pal : κ → Bool) (β : Bool) : Finset κ :=
  univ.filter fun c => (∃ a b, col a b = c) ∧ pal c = β

/-- `gstar N t` = the largest `g ≤ N` with `g * (g - 1) / 2 ≤ t`. -/
def gstar (N t : ℕ) : ℕ := ((range (N + 1)).filter fun g => g * (g - 1) / 2 ≤ t).sup id

lemma le_gstar {N t g : ℕ} (hg : g ≤ N) (ht : g * (g - 1) / 2 ≤ t) : g ≤ gstar N t :=
  Finset.le_sup (f := id) (by simp [Finset.mem_filter, Finset.mem_range]; omega)

lemma gstar_spec (N t : ℕ) : gstar N t ≤ N ∧ gstar N t * (gstar N t - 1) / 2 ≤ t := by
  have hne : ((range (N + 1)).filter fun g => g * (g - 1) / 2 ≤ t).Nonempty :=
    ⟨0, by simp⟩
  obtain ⟨g, hg, hgs⟩ := Finset.exists_mem_eq_sup _ hne id
  rw [Finset.mem_filter, Finset.mem_range] at hg
  unfold gstar
  rw [hgs]
  exact ⟨by simp; omega, hg.2⟩

lemma otherUsed_eq (col : Fin N → Fin N → κ) (pal : κ → Bool) (j : κ) :
    otherUsed col pal j = usedPal col pal (!pal j) := by
  ext c
  simp only [otherUsed, usedPal, Finset.mem_filter, Finset.mem_univ, true_and]
  cases pal c <;> cases pal j <;> simp

/-- `gstar N t ≤ ⌊√(2t)⌋ + 1`. -/
theorem gstar_le_sqrt (N t : ℕ) : gstar N t ≤ Nat.sqrt (2 * t) + 1 := by
  obtain ⟨-, h⟩ := gstar_spec N t
  generalize gstar N t = g at h ⊢
  rcases Nat.lt_or_ge g 2 with hg | hg
  · omega
  · have h1 : g * (g - 1) < 2 * t + 2 := by
      have := Nat.lt_succ_of_le h
      rw [Nat.div_lt_iff_lt_mul (by norm_num)] at this
      omega
    have h2 : (g - 1) * (g - 1) ≤ 2 * t := by
      obtain ⟨k, rfl⟩ : ∃ k, g = k + 2 := ⟨g - 2, by omega⟩
      have e : k + 2 - 1 = k + 1 := by omega
      rw [e] at h1 ⊢
      nlinarith
    have := Nat.le_sqrt.2 h2
    omega

/-- Real form: `gstar N t ≤ √(2t) + 1`. -/
theorem gstar_le_real_sqrt (N t : ℕ) : (gstar N t : ℝ) ≤ Real.sqrt (2 * t) + 1 := by
  have h := gstar_le_sqrt N t
  have h' : ((Nat.sqrt (2 * t) : ℕ) : ℝ) ≤ Real.sqrt (2 * t) := by
    have := Real.nat_sqrt_le_real_sqrt (a := 2 * t)
    simpa using this
  calc (gstar N t : ℝ) ≤ ((Nat.sqrt (2 * t) + 1 : ℕ) : ℝ) := by exact_mod_cast h
    _ = ((Nat.sqrt (2 * t) : ℕ) : ℝ) + 1 := by push_cast; ring
    _ ≤ _ := by linarith

/-- `q ≤ p → q √(2p) ≤ p √(2q)` (reals, `p q : ℕ`). -/
lemma mul_sqrt_le {p q : ℕ} (hqp : q ≤ p) :
    (q : ℝ) * Real.sqrt (2 * p) ≤ (p : ℝ) * Real.sqrt (2 * q) := by
  have hq0 : (0 : ℝ) ≤ q := by positivity
  have hp0 : (0 : ℝ) ≤ p := by positivity
  have hqp' : (q : ℝ) ≤ p := by exact_mod_cast hqp
  have e1 : (q : ℝ) * Real.sqrt (2 * p) = Real.sqrt ((q : ℝ) ^ 2 * (2 * p)) := by
    rw [Real.sqrt_mul (sq_nonneg _), Real.sqrt_sq hq0]
  have e2 : (p : ℝ) * Real.sqrt (2 * q) = Real.sqrt ((p : ℝ) ^ 2 * (2 * q)) := by
    rw [Real.sqrt_mul (sq_nonneg _), Real.sqrt_sq hp0]
  rw [e1, e2]
  apply Real.sqrt_le_sqrt
  have : (q : ℝ) ^ 2 * (2 * p) = 2 * (q * p) * q := by ring
  have h2 : (p : ℝ) ^ 2 * (2 * q) = 2 * (q * p) * p := by ring
  rw [this, h2]
  exact mul_le_mul_of_nonneg_left hqp' (by positivity)

section KeyLemma

variable (col : Fin N → Fin N → κ) (pal : κ → Bool)
  (hrow : ∀ a b b', col a b = col a b' → b = b')
  (hcol : ∀ a a' b, col a b = col a' b → a = a')
  (hint : ∀ x y u v, x ≠ v → y ≠ u → col x y = col v u →
    (col x u = col v y ∧ pal (col x u) ≠ pal (col x y)))
include hrow hcol hint

/-- Two cells of colour `j` in different rows: the cell `(x,u)` has the other palette. -/
lemma pal_ne_of_cells {j : κ} {x y v u : Fin N} (hx : col x y = j) (hv : col v u = j)
    (hxv : x ≠ v) : pal (col x u) ≠ pal j := by
  have hyu : y ≠ u := by
    rintro rfl
    exact hxv (hcol _ _ _ (hx.trans hv.symm))
  have := (hint x y u v hxv hyu (hx.trans hv.symm)).2
  rwa [hx] at this

/-- The fibre lemma: the map `((a,b),(a',b')) ↦ col a b'` is two-to-one on ordered pairs of
distinct cells of colour `j`. -/
lemma fibre {j : κ} {a b a' b' c d c' d' : Fin N}
    (h1 : col a b = j) (h2 : col a' b' = j) (h12 : (a, b) ≠ (a', b'))
    (h3 : col c d = j) (h4 : col c' d' = j) (h34 : (c, d) ≠ (c', d'))
    (he : col a b' = col c d') :
    ((c, d), (c', d')) = ((a, b), (a', b')) ∨ ((c, d), (c', d')) = ((a', b'), (a, b)) := by
  have haa' : a ≠ a' := by
    rintro rfl
    exact h12 (by rw [hrow _ _ _ (h1.trans h2.symm)])
  have hcc' : c ≠ c' := by
    rintro rfl
    exact h34 (by rw [hrow _ _ _ (h3.trans h4.symm)])
  by_cases hac : a = c
  · subst hac
    have hbd' : b' = d' := hrow _ _ _ he
    subst hbd'
    have : a' = c' := hcol _ _ _ (h2.trans h4.symm)
    subst this
    have : b = d := hrow _ _ _ (h1.trans h3.symm)
    subst this
    left; rfl
  · have hbd' : b' ≠ d' := by
      rintro rfl
      exact hac (hcol _ _ _ he)
    obtain ⟨hsym, hpal⟩ := hint a b' d' c hac hbd' he
    have hp1 : pal (col a b') ≠ pal j := pal_ne_of_cells col pal hrow hcol hint h1 h2 haa'
    have hpj : pal (col a d') = pal j := by
      revert hpal hp1
      cases pal (col a d') <;> cases pal (col a b') <;> cases pal j <;> simp
    have hpj' : pal (col c b') = pal j := hsym ▸ hpj
    have hac' : a = c' := by
      by_contra hne
      exact pal_ne_of_cells col pal hrow hcol hint h1 h4 hne hpj
    have hca' : c = a' := by
      by_contra hne
      exact pal_ne_of_cells col pal hrow hcol hint h3 h2 hne hpj'
    subst hac' hca'
    have : d' = b := hrow _ _ _ (h4.trans h1.symm)
    subst this
    have : d = b' := hrow _ _ _ (h3.trans h2.symm)
    subst this
    right; rfl

/-- **Key lemma.** `#(cells of colour j) choose 2 ≤ #(used colours of the other palette)`. -/
theorem key_lemma (j : κ) :
    (cells col j).card.choose 2 ≤ (otherUsed col pal j).card := by
  classical
  set s := (cells col j).offDiag with hs
  have mem_s : ∀ e : (Fin N × Fin N) × (Fin N × Fin N), e ∈ s ↔
      col e.1.1 e.1.2 = j ∧ col e.2.1 e.2.2 = j ∧ e.1 ≠ e.2 := by
    intro e; simp [hs, Finset.mem_offDiag, cells]
  have hmap : ∀ e ∈ s, col e.1.1 e.2.2 ∈ otherUsed col pal j := by
    intro e he
    obtain ⟨h1, h2, h12⟩ := (mem_s e).1 he
    have haa' : e.1.1 ≠ e.2.1 := by
      intro h
      apply h12
      have := hrow _ _ _ (h1.trans (h ▸ h2.symm))
      exact Prod.ext h this
    simp only [otherUsed, Finset.mem_filter, Finset.mem_univ, true_and]
    exact ⟨⟨_, _, rfl⟩, pal_ne_of_cells col pal hrow hcol hint h1 h2 haa'⟩
  have hfib : ∀ l ∈ otherUsed col pal j, (s.filter fun e => col e.1.1 e.2.2 = l).card ≤ 2 := by
    intro l _
    by_cases hne : (s.filter fun e => col e.1.1 e.2.2 = l).Nonempty
    · obtain ⟨e0, he0⟩ := hne
      rw [Finset.mem_filter] at he0
      have hsub : (s.filter fun e => col e.1.1 e.2.2 = l) ⊆ {e0, (e0.2, e0.1)} := by
        intro e he
        rw [Finset.mem_filter] at he
        obtain ⟨h1, h2, h12⟩ := (mem_s e0).1 he0.1
        obtain ⟨h3, h4, h34⟩ := (mem_s e).1 he.1
        have key := fibre col pal hrow hcol hint (a := e0.1.1) (b := e0.1.2) (a' := e0.2.1)
          (b' := e0.2.2) (c := e.1.1) (d := e.1.2) (c' := e.2.1) (d' := e.2.2) h1 h2 h12 h3 h4 h34
          (he0.2.trans he.2.symm)
        simp only [Finset.mem_insert, Finset.mem_singleton]
        rcases key with h | h
        · left; simpa using h
        · right; simpa using h
      exact (Finset.card_le_card hsub).trans (Finset.card_insert_le _ _)
    · rw [Finset.not_nonempty_iff_eq_empty] at hne
      simp [hne]
  have hcard := Finset.card_le_mul_card_image_of_maps_to hmap 2 hfib
  rw [Finset.offDiag_card] at hcard
  rw [Nat.choose_two_right]
  apply Nat.div_le_of_le_mul
  rw [Nat.mul_sub_one]
  exact hcard

/-- Every colour class meets each row at most once, hence has at most `N` cells. -/
theorem cells_card_le (j : κ) : (cells col j).card ≤ N := by
  have := Finset.card_le_card_of_injOn (fun p : Fin N × Fin N => p.1) (s := cells col j)
    (t := Finset.univ)
    (fun _ _ => Finset.mem_univ _) (by
      intro p hp q hq h
      have hp' : col p.1 p.2 = j := by simpa [cells] using hp
      have hq' : col q.1 q.2 = j := by simpa [cells] using hq
      have h' : p.1 = q.1 := h
      exact Prod.ext h' (hrow _ _ _ (hp'.trans (h' ▸ hq'.symm))))
  simpa using this

/-- A colour class has at most `gstar N (#other palette)` cells. -/
theorem cells_card_le_gstar (j : κ) :
    (cells col j).card ≤ gstar N (otherUsed col pal j).card := by
  apply le_gstar (cells_card_le col pal hrow hcol hint j)
  rw [← Nat.choose_two_right]
  exact key_lemma col pal hrow hcol hint j

/-- **Counting inequality.** With `p = #(used colours with pal = true)` and
`q = #(used colours with pal = false)`: `N² ≤ q * gstar N p + p * gstar N q`. -/
theorem counting :
    N ^ 2 ≤ (usedPal col pal false).card * gstar N (usedPal col pal true).card
      + (usedPal col pal true).card * gstar N (usedPal col pal false).card := by
  classical
  set used : Finset κ := univ.filter fun c => ∃ a b, col a b = c with hused
  have hmaps : ∀ x ∈ (univ : Finset (Fin N × Fin N)), col x.1 x.2 ∈ used := by
    intro x _; simp [hused]
  have hsum := Finset.card_eq_sum_card_fiberwise (f := fun x : Fin N × Fin N => col x.1 x.2)
    (s := univ) (t := used) (fun x hx => hmaps x hx)
  have hN : (univ : Finset (Fin N × Fin N)).card = N ^ 2 := by simp [sq]
  rw [hN] at hsum
  rw [hsum, ← Finset.sum_filter_add_sum_filter_not used (fun c => pal c = true)]
  have hT : used.filter (fun c => pal c = true) = usedPal col pal true := by
    ext c; simp [hused, usedPal]
  have hS : used.filter (fun c => ¬ pal c = true) = usedPal col pal false := by
    ext c; simp [hused, usedPal]
  rw [hT, hS, add_comm]
  apply add_le_add
  · refine (Finset.sum_le_card_nsmul _ _ (gstar N (usedPal col pal true).card) ?_).trans
      (by simp)
    intro c hc
    have hpc : pal c = false := by simp [usedPal] at hc; exact hc.2
    have := cells_card_le_gstar col pal hrow hcol hint c
    rw [otherUsed_eq, hpc] at this
    exact this
  · refine (Finset.sum_le_card_nsmul _ _ (gstar N (usedPal col pal false).card) ?_).trans
      (by simp)
    intro c hc
    have hpc : pal c = true := by simp [usedPal] at hc; exact hc.2
    have := cells_card_le_gstar col pal hrow hcol hint c
    rw [otherUsed_eq, hpc] at this
    exact this

/-- **Corollary (trade-off, real form).** With `p = #T` (`pal = true`), `q = #S` (`pal = false`),
if `q ≤ p` then `N² ≤ 2 p (√(2q) + 1)`. -/
theorem tradeoff_real (hqp : (usedPal col pal false).card ≤ (usedPal col pal true).card) :
    (N : ℝ) ^ 2 ≤ 2 * (usedPal col pal true).card *
      (Real.sqrt (2 * (usedPal col pal false).card) + 1) := by
  set p := (usedPal col pal true).card
  set q := (usedPal col pal false).card
  have hc : ((N ^ 2 : ℕ) : ℝ) ≤ ((q * gstar N p + p * gstar N q : ℕ) : ℝ) := by
    exact_mod_cast counting col pal hrow hcol hint
  push_cast at hc
  have h1 := gstar_le_real_sqrt N p
  have h2 := gstar_le_real_sqrt N q
  have h3 := mul_sqrt_le (p := p) (q := q) hqp
  have hq0 : (0 : ℝ) ≤ q := by positivity
  have hp0 : (0 : ℝ) ≤ p := by positivity
  have hqp' : (q : ℝ) ≤ p := by exact_mod_cast hqp
  calc (N : ℝ) ^ 2 ≤ q * gstar N p + p * gstar N q := hc
    _ ≤ q * (Real.sqrt (2 * p) + 1) + p * (Real.sqrt (2 * q) + 1) := by
        gcongr
    _ ≤ _ := by nlinarith

end KeyLemma

end Intercalate

#print axioms Intercalate.key_lemma
#print axioms Intercalate.counting
#print axioms Intercalate.gstar_le_sqrt
#print axioms Intercalate.tradeoff_real
