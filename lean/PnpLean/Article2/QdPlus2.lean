import PnpLean.Article2.PerPoint

/-!
# The `Q_{d+2}` lemma (commutative class, one point per direction)

Source: `code/gen8/g8d/notes.md` §3.

Setting. For every vertex `z` of the cube `Fin (d+2) → Bool` a polynomial `P z ∈ ℚ[X]` of degree
`≤ d` with `P z (0) = 1`; pairwise distinct nonzero points `τ k`, one per direction `k`.
The edge `{z, z^k}` is *good* iff `(P z + P (z^k)) (τ k) = 0`, *bad* otherwise.
`qd2_lemma`: the cube has at least `4` bad edges.

* `two_splittings`: the polynomial core (the same `Ψ_Q` computed through two splittings).
* `four_le_card_bad`: facet-level abstract form (only uses the facet sums `Ψ a β`).
* `qd2_lemma`: the statement on the cube.
-/

open Polynomial

namespace QdPlus2

/-- Barrier B core: a polynomial of degree `≤ d`, nonzero at `0`, has at most `d` roots. -/
theorem card_le_of_roots (Ψ : ℚ[X]) (d : ℕ) (hdeg : Ψ.natDegree ≤ d) (h0 : Ψ.eval 0 ≠ 0)
    (Z : Finset ℚ) (hZ : ∀ t ∈ Z, Ψ.eval t = 0) : Z.card ≤ d := by
  by_contra h
  push Not at h
  have := eq_zero_of_natDegree_lt_card_of_eval_eq_zero' Ψ Z hZ (lt_of_le_of_lt hdeg h)
  exact h0 (by simp [this])

/-- **Polynomial core of the `Q_{d+2}` lemma.** Two splittings `Ψ₀ + Ψ₁ = Φ₀ + Φ₁` of the same
polynomial into pieces of degree `≤ d`, with `Ψ₀(0) = Φ₀(0)` and `Ψ₁(0) ≠ 0`; a common root set `W`
with `|W| ≥ d - 1`; `Ψ₀, Φ₀` also vanish at `τc`, `Ψ₁` at `τb`, `Φ₁` at `τa`, where
`τa ≠ τb` and `0, τa, τb, τc ∉ W`, `τc ≠ 0`. This is impossible. -/
theorem two_splittings (d : ℕ) (Ψ₀ Ψ₁ Φ₀ Φ₁ : ℚ[X])
    (hΨ₀ : Ψ₀.natDegree ≤ d) (hΨ₁ : Ψ₁.natDegree ≤ d) (hΦ₀ : Φ₀.natDegree ≤ d)
    (h00 : Ψ₀.eval 0 = Φ₀.eval 0) (h10 : Ψ₁.eval 0 ≠ 0)
    (hsum : Ψ₀ + Ψ₁ = Φ₀ + Φ₁)
    (W : Finset ℚ) (hW : d < W.card + 2) (τa τb τc : ℚ)
    (hab : τa ≠ τb) (haW : τa ∉ W) (hbW : τb ∉ W) (hcW : τc ∉ W) (h0W : (0 : ℚ) ∉ W)
    (hc0 : τc ≠ 0)
    (hΨ₀W : ∀ t ∈ W, Ψ₀.eval t = 0) (hΨ₀c : Ψ₀.eval τc = 0)
    (hΦ₀W : ∀ t ∈ W, Φ₀.eval t = 0) (hΦ₀c : Φ₀.eval τc = 0)
    (hΨ₁W : ∀ t ∈ W, Ψ₁.eval t = 0) (hΨ₁b : Ψ₁.eval τb = 0)
    (hΦ₁a : Φ₁.eval τa = 0) : False := by
  have h1 : Ψ₀ - Φ₀ = 0 := by
    apply eq_zero_of_natDegree_lt_card_of_eval_eq_zero' (Ψ₀ - Φ₀) (insert 0 (insert τc W))
    · intro t ht
      simp only [Finset.mem_insert] at ht
      rcases ht with rfl | rfl | ht
      · simp [h00]
      · simp [hΨ₀c, hΦ₀c]
      · simp [hΨ₀W t ht, hΦ₀W t ht]
    · have hc : (insert 0 (insert τc W)).card = W.card + 2 := by
        rw [Finset.card_insert_of_notMem, Finset.card_insert_of_notMem hcW]
        simp [Finset.mem_insert, h0W, hc0.symm]
      have := (natDegree_sub_le Ψ₀ Φ₀).trans (max_le hΨ₀ hΦ₀)
      omega
  have h2 : Ψ₁ = Φ₁ := by
    have h : Ψ₀ = Φ₀ := sub_eq_zero.1 h1
    rw [h] at hsum
    exact add_left_cancel hsum
  have h3 : Ψ₁ = 0 := by
    apply eq_zero_of_natDegree_lt_card_of_eval_eq_zero' Ψ₁ (insert τa (insert τb W))
    · intro t ht
      simp only [Finset.mem_insert] at ht
      rcases ht with rfl | rfl | ht
      · rw [h2]; exact hΦ₁a
      · exact hΨ₁b
      · exact hΨ₁W t ht
    · have hc : (insert τa (insert τb W)).card = W.card + 2 := by
        rw [Finset.card_insert_of_notMem, Finset.card_insert_of_notMem hbW]
        simp [Finset.mem_insert, haW, hab]
      omega
  exact h10 (by simp [h3])

/-- **Facet-level form of the `Q_{d+2}` lemma.** `S` is the set of `d+2` directions, `Ψ a β` the
sum of the vertex polynomials over the facet `{z_a = β}`. `B` is the set of bad edges; a bad
edge `e` has a direction `dir e` and lies, for every other direction `a`, in the facet
`(a, side e a)`. Hypothesis `hvan` is the only link: a facet polynomial vanishes at `τ k` when the
facet contains no bad `k`-edge. Then `|B| ≥ 4`. -/
theorem four_le_card_bad {S ε : Type*} [Fintype S] [DecidableEq S] [DecidableEq ε] (d : ℕ)
    (hS : Fintype.card S = d + 2)
    (τ : S → ℚ) (hτ : Function.Injective τ) (hτ0 : ∀ k, τ k ≠ 0)
    (Ψ : S → Bool → ℚ[X]) (Ψtot : ℚ[X]) (c : ℚ) (hc : c ≠ 0)
    (hdeg : ∀ a β, (Ψ a β).natDegree ≤ d) (h0 : ∀ a β, (Ψ a β).eval 0 = c)
    (hadd : ∀ a, Ψ a false + Ψ a true = Ψtot)
    (B : Finset ε) (dir : ε → S) (side : ε → S → Bool)
    (hvan : ∀ a β k, k ≠ a → (∀ e ∈ B, dir e = k → side e a ≠ β) → (Ψ a β).eval (τ k) = 0) :
    4 ≤ B.card := by
  by_contra hlt
  push Not at hlt
  -- every facet contains a bad edge (Barrier B)
  have L1 : ∀ a β, ∃ e ∈ B, dir e ≠ a ∧ side e a = β := by
    intro a β
    by_contra hne
    push Not at hne
    have hroots : ∀ t ∈ (Finset.univ.erase a).image τ, (Ψ a β).eval t = 0 := by
      intro t ht
      obtain ⟨k, hk, rfl⟩ := Finset.mem_image.1 ht
      have hka : k ≠ a := Finset.ne_of_mem_erase hk
      exact hvan a β k hka (fun e he hek => hne e he (hek ▸ hka))
    have := card_le_of_roots (Ψ a β) d (hdeg a β) (by rw [h0]; exact hc) _ hroots
    rw [Finset.card_image_of_injective _ hτ, Finset.card_erase_of_mem (Finset.mem_univ a),
      Finset.card_univ, hS] at this
    omega
  have hSne : Nonempty S := Fintype.card_pos_iff.1 (by omega)
  obtain ⟨a₀⟩ := hSne
  obtain ⟨e₁, he₁, -, -⟩ := L1 a₀ false
  obtain ⟨e₂, he₂, h2a, h2s⟩ := L1 (dir e₁) false
  obtain ⟨e₃, he₃, h3a, h3s⟩ := L1 (dir e₁) true
  have h12 : e₁ ≠ e₂ := fun h => h2a (h ▸ rfl)
  have h13 : e₁ ≠ e₃ := fun h => h3a (h ▸ rfl)
  have h23 : e₂ ≠ e₃ := by
    intro h
    rw [h, h3s] at h2s
    exact Bool.noConfusion h2s
  have hmem : ∀ e ∈ B, e = e₁ ∨ e = e₂ ∨ e = e₃ := by
    have hB : ({e₁, e₂, e₃} : Finset ε) = B := by
      apply Finset.eq_of_subset_of_card_le
      · intro x hx
        simp only [Finset.mem_insert, Finset.mem_singleton] at hx
        rcases hx with rfl | rfl | rfl <;> assumption
      · have : ({e₁, e₂, e₃} : Finset ε).card = 3 :=
          Finset.card_eq_three.2 ⟨e₁, e₂, e₃, h12, h13, h23, rfl⟩
        omega
    intro e he
    rw [← hB] at he
    simpa using he
  -- splitting along b = dir e₂: e₁ and e₃ lie in different b-facets, and dir e₃ ≠ b
  obtain ⟨f, hf, hfb, hfs⟩ := L1 (dir e₂) false
  obtain ⟨g, hg, hgb, hgs⟩ := L1 (dir e₂) true
  have hcb : dir e₃ ≠ dir e₂ ∧ side e₁ (dir e₂) ≠ side e₃ (dir e₂) := by
    rcases hmem f hf with h | h | h <;> rw [h] at hfb hfs
    · rcases hmem g hg with h' | h' | h' <;> rw [h'] at hgb hgs
      · rw [hfs] at hgs
        exact Bool.noConfusion hgs
      · exact absurd rfl hgb
      · exact ⟨hgb, by rw [hfs, hgs]; decide⟩
    · exact absurd rfl hfb
    · rcases hmem g hg with h' | h' | h' <;> rw [h'] at hgb hgs
      · exact ⟨hfb, by rw [hfs, hgs]; decide⟩
      · exact absurd rfl hgb
      · rw [hfs] at hgs
        exact Bool.noConfusion hgs
  obtain ⟨hcb1, hcb2⟩ := hcb
  -- vanishing criterion specialised to B = {e₁, e₂, e₃}
  have clean : ∀ a β k, k ≠ a → (dir e₁ = k → side e₁ a ≠ β) → (dir e₂ = k → side e₂ a ≠ β) →
      (dir e₃ = k → side e₃ a ≠ β) → (Ψ a β).eval (τ k) = 0 := by
    intro a β k hka c1 c2 c3
    apply hvan a β k hka
    intro e he
    rcases hmem e he with rfl | rfl | rfl <;> assumption
  set W : Finset ℚ := (((Finset.univ.erase (dir e₁)).erase (dir e₂)).erase (dir e₃)).image τ
    with hWdef
  have hWcard : W.card + 3 = d + 2 := by
    rw [hWdef, Finset.card_image_of_injective _ hτ, Finset.card_erase_of_mem,
      Finset.card_erase_of_mem, Finset.card_erase_of_mem (Finset.mem_univ _), Finset.card_univ, hS]
    · have : 3 ≤ Fintype.card S := by
        have h3 : ({dir e₁, dir e₂, dir e₃} : Finset S).card = 3 :=
          Finset.card_eq_three.2 ⟨_, _, _, h2a.symm, h3a.symm, hcb1.symm, rfl⟩
        rw [← h3]
        exact Finset.card_le_univ _
      omega
    · simp [h2a]
    · simp [h3a, hcb1]
  have hWmem : ∀ t ∈ W, ∃ k, k ≠ dir e₁ ∧ k ≠ dir e₂ ∧ k ≠ dir e₃ ∧ τ k = t := by
    intro t ht
    obtain ⟨k, hk, rfl⟩ := Finset.mem_image.1 ht
    simp only [Finset.mem_erase, Finset.mem_univ, and_true] at hk
    exact ⟨k, hk.2.2, hk.2.1, hk.1, rfl⟩
  have hnotW : ∀ k, (k = dir e₁ ∨ k = dir e₂ ∨ k = dir e₃) → τ k ∉ W := by
    intro k hk hkW
    obtain ⟨k', h1, h2, h3, h4⟩ := hWmem _ hkW
    have := hτ h4
    subst this
    rcases hk with h | h | h <;> contradiction
  have h0W : (0 : ℚ) ∉ W := by
    intro h
    obtain ⟨k, -, -, -, hk⟩ := hWmem _ h
    exact hτ0 k hk
  have vanW : ∀ a β, ∀ t ∈ W, (a = dir e₁ ∨ a = dir e₂ ∨ a = dir e₃) → (Ψ a β).eval t = 0 := by
    intro a β t ht ha
    obtain ⟨k, h1, h2, h3, rfl⟩ := hWmem t ht
    apply clean a β k
    · rintro rfl
      rcases ha with h | h | h <;> contradiction
    · exact fun h => absurd h.symm h1
    · exact fun h => absurd h.symm h2
    · exact fun h => absurd h.symm h3
  refine two_splittings d (Ψ (dir e₁) false) (Ψ (dir e₁) true)
    (Ψ (dir e₂) (side e₁ (dir e₂))) (Ψ (dir e₂) (side e₃ (dir e₂)))
    (hdeg _ _) (hdeg _ _) (hdeg _ _) (by rw [h0, h0]) (by rw [h0]; exact hc) ?_
    W (by omega) (τ (dir e₁)) (τ (dir e₂)) (τ (dir e₃))
    (fun h => h2a (hτ h).symm) (hnotW _ (Or.inl rfl)) (hnotW _ (Or.inr (Or.inl rfl)))
    (hnotW _ (Or.inr (Or.inr rfl))) h0W (hτ0 _)
    (fun t ht => vanW _ _ t ht (Or.inl rfl)) ?_
    (fun t ht => vanW _ _ t ht (Or.inr (Or.inl rfl))) ?_
    (fun t ht => vanW _ _ t ht (Or.inl rfl)) ?_ ?_
  · rw [hadd]
    cases h1 : side e₁ (dir e₂) <;> cases h3 : side e₃ (dir e₂)
    · exact absurd (h1.trans h3.symm) hcb2
    · exact (hadd _).symm
    · rw [add_comm]; exact (hadd _).symm
    · exact absurd (h1.trans h3.symm) hcb2
  · -- Ψ (a,false) at τ c
    apply clean _ _ _ h3a
    · exact fun h => absurd h.symm h3a
    · exact fun h => absurd h.symm hcb1
    · intro _; rw [h3s]; decide
  · -- Ψ (b, side e₁ b) at τ c
    apply clean _ _ _ hcb1
    · exact fun h => absurd h.symm h3a
    · exact fun h => absurd h.symm hcb1
    · exact fun _ => hcb2.symm
  · -- Ψ (a,true) at τ b
    apply clean _ _ _ h2a
    · exact fun h => absurd h.symm h2a
    · intro _; rw [h2s]; decide
    · exact fun h => absurd h hcb1
  · -- Ψ (b, side e₃ b) at τ a
    apply clean _ _ _ h2a.symm
    · exact fun _ => hcb2
    · exact fun h => absurd h h2a
    · exact fun h => absurd h h3a

section Cube

variable {n : ℕ}

/-- The neighbour `z^k` of `z` in direction `k`. -/
def flip (z : Fin n → Bool) (k : Fin n) : Fin n → Bool := Function.update z k (!z k)

theorem flip_same (z : Fin n → Bool) (k : Fin n) : flip z k k = !z k := by simp [flip]

theorem flip_ne (z : Fin n → Bool) {k a : Fin n} (h : a ≠ k) : flip z k a = z a := by
  simp [flip, Function.update_of_ne h]

theorem flip_flip (z : Fin n → Bool) (k : Fin n) : flip (flip z k) k = z := by
  funext a
  by_cases h : a = k
  · subst h; simp [flip_same]
  · simp [flip_ne _ h]

/-- Bad edges of the cube. The edge `{z, z^k}` is represented by `(z, k)` with `z k = false`;
it is bad iff `(P z + P (z^k)) (τ k) ≠ 0`. -/
noncomputable def badEdges (P : (Fin n → Bool) → ℚ[X]) (τ : Fin n → ℚ) :
    Finset ((Fin n → Bool) × Fin n) :=
  Finset.univ.filter fun e => e.1 e.2 = false ∧ eval (τ e.2) (P e.1 + P (flip e.1 e.2)) ≠ 0

/-- Sum of the vertex polynomials over the facet `{z | z a = β}`. -/
noncomputable def facetSum (P : (Fin n → Bool) → ℚ[X]) (a : Fin n) (β : Bool) : ℚ[X] :=
  ∑ z ∈ Finset.univ.filter (fun z : Fin n → Bool => z a = β), P z

theorem facetSum_add (P : (Fin n → Bool) → ℚ[X]) (a : Fin n) :
    facetSum P a false + facetSum P a true = ∑ z, P z := by
  unfold facetSum
  rw [← Finset.sum_filter_add_sum_filter_not Finset.univ (fun z : Fin n → Bool => z a = false)]
  congr 2
  ext z
  simp

/-- A facet sum vanishes at `τ k` if the facet contains no bad `k`-edge. -/
theorem facetSum_eval_eq_zero (P : (Fin n → Bool) → ℚ[X]) (τ : Fin n → ℚ) (a : Fin n) (β : Bool)
    (k : Fin n) (hka : k ≠ a)
    (hclean : ∀ e ∈ badEdges P τ, e.2 = k → e.1 a ≠ β) :
    (facetSum P a β).eval (τ k) = 0 := by
  have hsplit : facetSum P a β
      = ∑ z ∈ Finset.univ.filter (fun z : Fin n → Bool => z a = β ∧ z k = false),
          (P z + P (flip z k)) := by
    unfold facetSum
    rw [Finset.sum_add_distrib,
      ← Finset.sum_filter_add_sum_filter_not
        (Finset.univ.filter (fun z : Fin n → Bool => z a = β)) (fun z => z k = false),
      Finset.filter_filter, Finset.filter_filter]
    congr 1
    apply Finset.sum_nbij' (fun z => flip z k) (fun z => flip z k)
    · intro z hz
      simp only [Finset.mem_filter, Finset.mem_univ, true_and] at hz ⊢
      refine ⟨by rw [flip_ne _ hka.symm]; exact hz.1, ?_⟩
      rw [flip_same]
      simpa using hz.2
    · intro z hz
      simp only [Finset.mem_filter, Finset.mem_univ, true_and] at hz ⊢
      refine ⟨by rw [flip_ne _ hka.symm]; exact hz.1, ?_⟩
      rw [flip_same]
      simp [hz.2]
    · intro z _; exact flip_flip z k
    · intro z _; exact flip_flip z k
    · intro z _; rw [flip_flip]
  rw [hsplit, eval_finsetSum]
  apply Finset.sum_eq_zero
  intro z hz
  simp only [Finset.mem_filter, Finset.mem_univ, true_and] at hz
  by_contra hne
  have hmem : (z, k) ∈ badEdges P τ := by
    unfold badEdges
    rw [Finset.mem_filter]
    exact ⟨Finset.mem_univ _, hz.2, hne⟩
  exact hclean (z, k) hmem rfl hz.1

/-- **The `Q_{d+2}` lemma** (G8d §3). On the cube `Fin (d+2) → Bool`, with vertex polynomials
`P z` of degree `≤ d`, `P z (0) = 1`, and pairwise distinct nonzero points `τ k`, there are at
least `4` bad edges. -/
theorem qd2_lemma (d : ℕ) (P : (Fin (d + 2) → Bool) → ℚ[X]) (τ : Fin (d + 2) → ℚ)
    (hτ : Function.Injective τ) (hτ0 : ∀ k, τ k ≠ 0)
    (hdeg : ∀ z, (P z).natDegree ≤ d) (hP0 : ∀ z, (P z).eval 0 = 1) :
    4 ≤ (badEdges P τ).card := by
  have hcard : ∀ (a : Fin (d + 2)) (β : Bool),
      (Finset.univ.filter (fun z : Fin (d + 2) → Bool => z a = β)).card = 2 ^ (d + 1) := by
    intro a β
    have h := PerPoint.card_subcube (d + 2) {a} (fun _ => β)
    rw [Fintype.card_subtype] at h
    simpa using h
  refine four_le_card_bad d (by simp) τ hτ hτ0 (facetSum P) (∑ z, P z) (2 ^ (d + 1))
    (by positivity) ?_ ?_ (facetSum_add P) (badEdges P τ) (fun e => e.2) (fun e a => e.1 a) ?_
  · intro a β
    exact natDegree_sum_le_of_forall_le _ _ (fun z _ => hdeg z)
  · intro a β
    unfold facetSum
    rw [eval_finsetSum]
    simp [hP0, hcard a β]
  · intro a β k hka hclean
    exact facetSum_eval_eq_zero P τ a β k hka hclean

end Cube

end QdPlus2

#print axioms QdPlus2.two_splittings
#print axioms QdPlus2.four_le_card_bad
#print axioms QdPlus2.qd2_lemma
